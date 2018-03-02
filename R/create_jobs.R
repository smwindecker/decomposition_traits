#' Create design matrix of unique model combinations
#'
#' @param trait_list a list of traits
#' @param model_structures dataframe listing model types and specifications
#' @param data dataframe containing data for models
#' @param cv_cluster optional specification for grouping of cv clusters
#' @param folder folder to save output
#' @return dataframe of all model structures
#' @examples
#' traits <- c('N', 'C', 'SLA', 'LDMC', 'HC', 'CL', 'LG')
#' model_df <- data.frame(model_type = c('w', 'w', 'ne'),
#'                        fixed_effects = c(TRUE, FALSE, FALSE),
#'                        random_effects = c(FALSE, TRUE, FALSE),
#'                        cv = c(TRUE, TRUE, TRUE))
#' jobs_list <- create_jobs(model_df, df,
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

  if (!is.null(fixed_effects_list)) {

    # create list of traits as expressions with tilda and intercept terms
    formulas <- c('~ 1', sprintf("%s + %s", '~ 1', fixed_effects_list))

    if (model_structures$model_type == 'ne') {
      formula_grid <- data.frame(formulas)
      colnames(formula_grid) <- 'formula_k'
    }

    if (model_structures$model_type == 'w') {
      # expand this formula list so that it includes all combinations of pairs
      formula_grid <- expand.grid(formulas, formulas, stringsAsFactors = FALSE)
      colnames(formula_grid) <- c('formula_alpha', 'formula_beta')
    }

    #specify that these correspond to those models with fixed effects
    formula_grid$fixed_effects <- 'FE'

    # merge, retaining rows that aren't fixed effects
    models <- merge(model_structures, formula_grid, all.x = TRUE)

  }

  if (is.null(fixed_effects_list)) {
    models <- model_structures
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
  write.csv(models, paste0(folder, 'model_formulas.csv'))

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
  lapply(jobs, append_data, data)

}





