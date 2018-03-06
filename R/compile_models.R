#' Select models to compile and compile them
#'
#' @param jobs_list job list
#' @param mass column name in data referring to logged mass remaining
#' @param initial_mass columne name in data referring to logged initial mass
#' @param time column name in data referring to time in years
#' @param group_id column name in data referring to random effects cluster
#' @return list with model type, iteration, neg log likelihood output, predicted and real for test data, and diagnostics
#' @import reshape2
#'
#' @export

compile_models <- function(jobs_list, mass, initial_mass, time, group_id) {

  # find the first instance of each model type, to run it with compile version
  all_stan_groups <- NULL

  get_stan_group <- function (x) {
    g <- x$stan_group
    c(g, all_stan_groups)
  }

  all <- unlist(lapply(jobs_list, get_stan_group))

  to_compile <- tapply(seq_along(all), all, identity)[unique(all)] %>%
    lapply(., head, 1) %>%
    unlist(.)

  foreach(i = to_compile) %dopar% {
    execute_stan(jobs_list[[i]],
                 mass = mass,
                 initial_mass = initial_mass,
                 time = time,
                 group_id = group_id,
                 compile_model = TRUE)
  }

}
