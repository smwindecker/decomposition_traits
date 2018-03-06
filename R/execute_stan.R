#' Run stan models from job list
#'
#' @param job job from job list
#' @param mass column name in data referring to logged mass remaining
#' @param initial_mass columne name in data referring to logged initial mass
#' @param time column name in data referring to time in years
#' @param group_id column name in data referring to random effects cluster
#' @param compile_model whether to fully compile model or not
#' @return stan output, and saved .rds file
#' @import rstan
#' @importFrom stats model.matrix
#' @export

execute_stan <- function (job, mass, initial_mass, time, group_id = NULL, compile_model = TRUE) {

  cv <- job$cv
  re <- job$random_effects
  model_type <- job$model_type

  # where 'train' refers either to training data in cv model or full dataset in non cv model
  train <- job$train
  test <- job$test

  if (!is.null(group_id)) {

    # set group levels
    train$group_level <- as.factor(as.numeric(as.factor(as.character(train[, group_id]))))

  }

  if (is.null(group_id)) {
    train$group_level <- 1
  }

  # number of levels
  n_train_levels <- length(unique(train$group_level))

  make_matrix <- function (parameter) {
    if (parameter == '~1') {
      matrix <- as.matrix(stats::model.matrix(parameter, train)[1:n_train_levels,])
    }
    if (parameter != '~1') {
      matrix <- unique(stats::model.matrix(parameter, train))
    }
    matrix
  }

  # model formulas for each possible model parameter
  if (!is.na(job$formula_k)) {
    formula_k <- as.formula(job$formula_k)
    X_k <- make_matrix(formula_k)
    X_k_test <- unique(model.matrix(formula_k, test))
  }
  if (!is.na(job$formula_alpha)) {
    formula_alpha <- as.formula(job$formula_alpha)
    X_alpha <- make_matrix(formula_alpha)
    X_alpha_test <- unique(model.matrix(formula_alpha, test))
  }
  if (!is.na(job$formula_beta)) {
    formula_beta <- as.formula(job$formula_beta)
    X_beta <- make_matrix(formula_beta)
    X_beta_test <- unique(model.matrix(formula_beta, test))
  }

  if (model_type == 'ne' & cv == 'CV') {

    stan_data <- list(mT = train[, mass],
                      mT_test = test[, mass],
                      m0 = train[, initial_mass],
                      m0_test = test[, initial_mass],
                      time = train[, time],
                      time_test = test[, time],
                      N = nrow(train),
                      N_test = nrow(test),
                      sp = as.numeric(as.factor(train$group_level)),
                      J = n_train_levels,
                      X_k = X_k,
                      X_k_test = X_k_test,
                      P = ncol(X_k))

  }

  if (model_type == 'ne' & cv == 'noCV') {

    stan_data <- list(mT = train[, mass],
                      m0 = train[, initial_mass],
                      time = train[, time],
                      N = nrow(train),
                      sp = as.numeric(as.factor(train$group_level)),
                      J = n_train_levels,
                      X_k = X_k,
                      P = ncol(X_k))

  }

  if (model_type == 'w' & cv == 'CV') {

    stan_data <- list(mT = train[, mass],
                      mT_test = test[, mass],
                      m0 = train[, initial_mass],
                      m0_test = test[, initial_mass],
                      time = train[, time],
                      time_test = test[, time],
                      N = nrow(train),
                      N_test = nrow(test),
                      sp = as.numeric(as.factor(train$group_level)),
                      J = n_train_levels,
                      X_alpha = X_alpha,
                      X_alpha_test = X_alpha_test,
                      X_beta = X_beta,
                      X_beta_test = X_beta_test,
                      P_alpha = ncol(X_alpha),
                      P_beta = ncol(X_beta))

  }

  if (model_type == 'w' & cv == 'noCV') {

    stan_data <- list(mT = train[, mass],
                      m0 = train[, initial_mass],
                      time = train[, time],
                      N = nrow(train),
                      sp = as.numeric(as.factor(train$group_level)),
                      J = n_train_levels,
                      X_alpha = X_alpha,
                      X_beta = X_beta,
                      P_alpha = ncol(X_alpha),
                      P_beta = ncol(X_beta))

  }

  if (compile_model == TRUE) {

    # run stan model
    cv_output <- rstan::stan(file = paste0('R/stan_', model_type, '_', cv, '_', re, '.stan'),
                             data = stan_data,
                             iter = 2000,
                             chains = 3,
                             control = list(adapt_delta = 0.99, max_treedepth = 15))

    saveRDS(cv_output, paste0('output/stan/', model_type, '_', cv, '_', re, '.rds'))

  }

  if (compile_model == FALSE) {

    cv_output <- rstan::stan(fit = readRDS(paste0('output/stan/', model_type, '_', cv, '_', re, '.rds')),
                             data = stan_data,
                             iter = 2000,
                             chains = 3)

    # save output as .rds
    saveRDS(cv_output, job$filename)

  }

  cv_output

}
