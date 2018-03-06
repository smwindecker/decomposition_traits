#' Create job list of all model iterations to be run
#'
#' @param model_structures dataframe listing model types and specifications
#' @param data dataframe containing data for models
#' @param folder folder path to save output
#' @param cv_cluster optional columns name to specify cv clusters
#' @param fixed_effects_list optional list of fixed effects - corresponding to column names in data
#' @return list of all model iterations
#' @importFrom utils write.csv
#' @examples
#' traits <- c('N', 'C', 'SLA', 'LDMC', 'HC', 'CL', 'LG')
#' model_df <- data.frame(model_type = c('w', 'w', 'ne'),
#'                        fixed_effects = c(TRUE, FALSE, FALSE),
#'                        random_effects = c(FALSE, TRUE, FALSE),
#'                        cv = c(TRUE, TRUE, TRUE))
#' jobs_list <- create_jobs(model_structures = model_df,
#'                          data = df,
#'                          folder = 'output/stan/',
#'                          cv_cluster = 'species_code',
#'                          fixed_effects_list = traits)
#'
#' @export

create_jobs <- function(model_structures, data, folder, cv_cluster = NULL, fixed_effects_list = NULL) {

  # modify structures df so in correct form
  model_structures$fixed_effects[model_structures$fixed_effects == TRUE] <- 'FE'
  model_structures$fixed_effects[model_structures$fixed_effects == FALSE] <- 'noFE'

  model_structures$random_effects[model_structures$random_effects == TRUE] <- 'RE'
  model_structures$random_effects[model_structures$random_effects == FALSE] <- 'noRE'

  model_structures$cv[model_structures$cv == TRUE] <- 'CV'
  model_structures$cv[model_structures$cv == FALSE] <- 'noCV'

  if (is.null(fixed_effects_list) & nrow(model_structures[model_structures$fixed_effects == 'FE', ]) > 0) {
    stop('fixed effects reported but no list provided')
  }

  models <- model_structures

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

  # write the formula dataframe
  utils::write.csv(models, paste0(folder, 'model_formulas.csv'))

  models[] <- lapply(models, as.character)

  # change dataframe to list
  jobs <- split(models, seq(nrow(models)))

  if (!is.null(cv_cluster)) {
    clusters <- unique(data[, cv_cluster])
  } else {
    clusters <- NULL
  }

  folder_path <- folder

  # expand models
  jobs <- lapply(jobs, expand_models, clusters, folder_path)
  jobs <- unlist(jobs, recursive = FALSE)

  # append correct data to each job in list
  data_appended <- lapply(jobs, append_data, data, cv_cluster)

  data_appended

}

