# Create job list of all model iterations to be run
#
# @param models dataframe listing model types and specifications
# @param data dataframe containing data for models
# @param cv_cluster optional columns name to specify cv clusters
# @param fixed_effects_list optional list of fixed effects - corresponding to column names in data
# @return list of all model iterations
# traits <- c('N', 'C', 'SLA', 'DMC', 'HC', 'CL', 'LG')
# model_df <- data.frame(model_type = c('w', 'w', 'ne'),
#                        fixed_effects = c(TRUE, FALSE, FALSE),
#                        random_effects = c(FALSE, TRUE, FALSE))
# jobs_list <- create_jobs(models = model_df,
#                          data = df,
#                          cv_cluster = 'species_code',
#                          fixed_effects_list = traits)

create_jobs <- function (models, data, cv_cluster = NULL, fixed_effects_list = NULL) {

  if (is.null(fixed_effects_list) & nrow(models[isTRUE(models$fixed_effects), ]) > 0) {
    stop('fixed effects reported but no list provided')
  }

  if (!is.null(fixed_effects_list)) {

    # create list of traits as expressions with tilda and intercept terms
    formulas <- c('~ 1', sprintf("%s + %s", '~ 1', fixed_effects_list))
    n <- length(formulas)
    ne_formulas <- data.frame(model_type = rep('ne', n),
                              fixed_effects = rep('FE', n),
                              formula_k = formulas,
                              stringsAsFactors = FALSE)

    models <- merge(models, ne_formulas, by = c('model_type', 'fixed_effects'), all.x = TRUE)

    formula_grid <- expand.grid(formulas, formulas, stringsAsFactors = FALSE)
    g <- nrow(formula_grid)
    w_formulas <- data.frame(model_type = rep('w', g),
                             fixed_effects = rep('FE', g),
                             formula_alpha = formula_grid[,1],
                             formula_beta = formula_grid[,2],
                             stringsAsFactors = FALSE)

    models <- merge(models, w_formulas, by = c('model_type', 'fixed_effects'), all.x = TRUE)

  }

  if (is.null(fixed_effects_list)) {
    models$formula_k <- NA
    models$formula_alpha <- NA
    models$formula_beta <- NA
  }

  # give formulas to non fixed effects models
  models$formula_k[models$model_type == 'ne' & is.na(models$formula_k)] <- '~ 1'
  models$formula_alpha[models$model_type == 'w' & is.na(models$formula_alpha)] <- '~ 1'
  models$formula_beta[models$model_type == 'w' & is.na(models$formula_beta)] <- '~ 1'

  # number the models
  models$model <- seq(1:nrow(models))

  models <- models[, - which(names(models) %in% 'fixed_effects')]

  # change dataframe to list
  jobs <- split(models, seq(nrow(models)))

  # if cv cluster is provided, these jobs are cv. and therefore need to be expanded
  if (!is.null(cv_cluster)) {
    clusters <- unique(data[, cv_cluster])
    tmp <- lapply(jobs, expand_models, clusters)
    jobs <- unlist(tmp, recursive = FALSE)
  }

  jobs

}

expand_models <- function (input, clusters) {

  # repeat each model n = number of species times
  expand <- expand.grid(cv_cluster = clusters,
                        model = input$model)
  expand$model <- as.character(expand$model)

  # merge with formula dataframe so also have the parameter functions listed
  job <- merge(expand, input, by = 'model')

  list <- split(job, seq(nrow(job)))

  return(list)

}

