#' Expand each cv model to have separate model run per cluster
#'
#' @param input job list item
#' @return list of jobs additional models added
#'
#' @export

expand_models <- function(input, cluster, folder) {

  if (!is.null(cluster)) {
    # repeat each model n = number of species times
    expand <- expand.grid(cv_cluster = cluster,
                          model = input$model)
    expand$model <- as.character(expand$model)

    # merge with formula dataframe so also have the parameter functions listed
    job <- merge(expand, input, by = 'model')

    job$cv_cluster <- as.character(job$cv_cluster)

    # create job id as mix of the model number and species code
    job$job_id <- paste0(job$model, '_', job$cv_cluster)
  }

  if (is.null(cluster)) {
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



