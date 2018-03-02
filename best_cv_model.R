best_cv_model <- function (jobs_list, model) {
  best_items <- function(x) {
    identical(x, model)
  }
  best_rows <- NULL
  for (i in 1:length(jobs_list)) {
    if (best_items(jobs_list[[i]]$model) == TRUE) {
      best_rows <- c(best_rows, i)
    }
  }
  foreach(i = best_rows) %dopar% {
    rep_weibull(init_w_output, jobs_list[[i]], formula_grid)
  }
}


run_full_model <- function (data, model) {

  formula_alpha <- as.formula(formula_grid[model, 1])
  formula_beta <- as.formula(formula_grid[model, 2])

  data$species_level <- as.factor(as.numeric(as.factor(as.character(data$species_code))))

  if (formula_alpha == '~1' & formula_beta == '~1') {
    X_alpha <- as.matrix(model.matrix(formula_alpha, data)[1:29,])
    X_beta <- as.matrix(model.matrix(formula_beta, data)[1:29,])
  } else if (formula_alpha == '~1' & formula_beta != '~1') {
    X_alpha <- as.matrix(model.matrix(formula_alpha, data)[1:29,])
    X_beta <- unique(model.matrix(formula_beta, data))
  } else if (formula_alpha != '~1' & formula_beta == '~1') {
    X_alpha <- unique(model.matrix(formula_alpha, data))
    X_beta <- as.matrix(model.matrix(formula_beta, data)[1:29,])
  } else {
    X_alpha <- unique(model.matrix(formula_alpha, data))
    X_beta <- unique(model.matrix(formula_beta, data))
  }

  stan_data <- list(mT = data$mRem,
                    m0 = data$mInit,
                    time = data$t,
                    N = nrow(data),
                    sp = as.numeric(as.factor(data$species_level)),
                    J = nlevels(data$species_level),
                    X_alpha = X_alpha,
                    X_beta = X_beta,
                    P_alpha = ncol(X_alpha),
                    P_beta = ncol(X_beta))

  full_output <- stan(file = "R/Matrix_Weibull_NoCV.stan", data = stan_data,
                      iter = 2000, chains = 3, cores = 3)

  saveRDS(full_output, paste0('R/stanOUTPUT/full_model_', model, '.rds'))
  full_output

}


run_full_negexp <- function (data) {

  formula_k <- as.formula(~ 1)

  X <- as.matrix(model.matrix(formula_k, data)[1:29,])

  stan_data <- list(mT = data$mRem,
                    m0 = data$mInit,
                    time = data$t,
                    N = nrow(data),
                    sp = data$species_level,
                    J = nlevels(data$sp),
                    X = X,
                    P = ncol(X))

  full_output <- stan(file = "R/Matrix_NegExp_NoCV.stan", data = stan_data,
                      iter = 2000, chains = 3)

  saveRDS(full_output, 'R/stanOUTPUT/null_model_negexp.rds')
  full_output

}
