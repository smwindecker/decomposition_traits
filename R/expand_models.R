#' Expand each cv model to have separate model run per cluster
#'
#' @param input job list item
#' @param clusters unique list of cv cluster id's
#' @param folder file path for output
#' @return list of jobs additional models added
#'
#' @export

expand_models <- function(input, clusters, folder) {

  if (!is.null(clusters)) {
    # repeat each model n = number of species times
    expand <- expand.grid(cv_cluster = clusters,
                          model = input$model)
    expand$model <- as.character(expand$model)

    # merge with formula dataframe so also have the parameter functions listed
    job <- merge(expand, input, by = 'model')

    job$cv_cluster <- as.character(job$cv_cluster)

    # create job id as mix of the model number and species code
    job$job_id <- paste0(job$model, '_', job$cv_cluster)
  }

  if (is.null(clusters)) {
    job <- input
    job$job_id <- paste0(job$model, '_noCV')
  }

  job$stan_group <- paste0(job$model_type, job$cv, job$random_effects)

  # add filename attribute
  job$filename <- NULL
  job$filename <- sprintf(paste0(folder, '%s.rds'), job$job_id)

  list <- split(job, seq(nrow(job)))

  return(list)

}



