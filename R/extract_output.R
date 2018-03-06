#' Run single model, using the compiled files
#'
#' @param job job
#' @param mass column name in data referring to logged mass remaining
#' @param initial_mass columne name in data referring to logged initial mass
#' @param time column name in data referring to time in years
#' @param group_id column name in data referring to random effects cluster
#' @return list for single model with model type, iteration, neg log likelihood output, predicted and real for test data, and diagnostics
#' @importFrom reshape2 melt
#' @importFrom rstan summary
#'
#' @export

extract_output <- function(job, mass, initial_mass, time, group_id) {

  model <- job$model
  iter <- job$cv_cluster

  fit <- execute_stan(job,
                      mass = mass,
                      initial_mass = initial_mass,
                      time = time,
                      group_id = group_id,
                      compile_model = FALSE)

  # make 1 row dataframe of the neglogliklihood
  neg_loglik_df <- as.data.frame(fit, 'neg_loglik')
  neg_loglik <- data.frame(model = model,
                           iter = iter,
                           mean = mean(neg_loglik_df$neg_loglik, na.rm = TRUE),
                           stringsAsFactors = FALSE)

  # make dataframe of the test data predicted v. real
  mT_pred_wide <- as.data.frame(fit, 'mT_pred')
  mT_pred <- reshape2::melt(mT_pred_wide,
                            variable.name = 'data_point',
                            value.name = 'draw')
  pred_real <- data.frame(model = rep(model, nrow(mT_pred)),
                          iter = rep(iter, nrow(mT_pred)),
                          data_point = mT_pred$data_point,
                          draw = mT_pred$draw,
                          mT_real = rep(job$test[, mass], each = nrow(mT_pred) / nrow(job$test)),
                          stringsAsFactors = FALSE)

  # diagnostics
  fit_summary <- rstan::summary(fit)$summary
  abs_rhat <- max(abs(fit_summary[,'Rhat'] - 1))
  neff_min <- min(fit_summary[,'n_eff'])
  sampler_params <- get_sampler_params(fit, inc_warmup = FALSE)
  sum_div <- sum(sapply(sampler_params, function(x) sum(x[, "divergent__"])))
  max_treedepth <- sapply(sampler_params, function(x) max(x[, "treedepth__"]))
  diagnostics <- data.frame(model = model,
                            iter = iter,
                            abs_rhat = abs_rhat,
                            neff_min = neff_min,
                            sum_div = sum_div,
                            max_treedepth = max_treedepth,
                            stringsAsFactors = FALSE)

  list(model = model,
       iter = iter,
       neg_loglik = neg_loglik,
       pred_real = pred_real,
       diagnostics = diagnostics)

}
