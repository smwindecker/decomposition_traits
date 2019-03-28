## Manuscript functions
# Deviance for model
dev <- function (model_file) {

  dev <- median(2*model_file$mod_specs$neg_loglik)
  dev

}

sigma_sp_beta <- function (mod) {

  median(rstan::extract(mod$fit, 'sigma_sp_beta')$sigma_sp_beta)

}
