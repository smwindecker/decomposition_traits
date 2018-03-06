#' Wrapper for extract_output, runs all models, using the compiled files
#'
#' @param jobs_list job list
#' @param mass column name in data referring to logged mass remaining
#' @param initial_mass columne name in data referring to logged initial mass
#' @param time column name in data referring to time in years
#' @param group_id column name in data referring to random effects cluster
#' @return list with model type, iteration, neg log likelihood output, predicted and real for test data, and diagnostics
#' @importFrom foreach foreach %dopar%
#'
#' @export

run_models <- function(jobs_list, mass, initial_mass, time, group_id) {

  # run all jobs using compiled models as pre-fits
  output_list <- list()

  output_list <- foreach(i = 1:length(jobs_list)) %dopar% {
    output_i <- extract_output(jobs_list[[i]], mass, initial_mass, time, group_id)
    return(output_i)
  }

  return(output_list)
}
