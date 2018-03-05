#' Expand each cv model to have separate model run per cluster
#'
#' @param job job list item
#' @return list with model type, iteration, neg log likelihood output, and predicted and real for test data
#' @import reshape2
#'
#' @export

extract_output <- function(job) {

  model <- job$model
  iter <- job$cv_cluster

  fit <- run_stan_models(job,
                         mass = 'log_mean',
                         initial_mass = 'log_i',
                         time = 't',
                         group_id = 'group_id',
                         compile_model = FALSE)

  neg_loglik_df <- as.data.frame(fit, 'neg_loglik')
  neg_loglik_mean <- mean(neg_loglik_df$neg_loglik, na.rm = TRUE)
  quants <- quantile(neg_loglik_df$neg_loglik, c(0.025, 0.975))

  neg_loglik <- data.frame(model = model,
                           iter = iter,
                           mean = neg_loglik_mean,
                           lower = quants[[1]],
                           upper = quants[[2]])

  mT_pred_wide <- as.data.frame(fit, 'mT_pred')

  mT_pred <- melt(mT_pred_wide,
                  variable.name = 'data_point',
                  value.name = 'draw')
  mT_real <- rep(job$test$log_mean, each = 3000)
  pred_real <- cbind(mT_pred, mT_real)

  list(model = model,
       iter = iter,
       neg_loglik = neg_loglik,
       pred_real = pred_real)

}
