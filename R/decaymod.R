## Model analysis

# runs all models from jobs list
run_models <- function (jobs_list, data, initial_mass, removal_mass, time, group, n_cores) {

  registerDoMC(n_cores)

  # so my decaymod outputs the fit. i want to append it to the list.
  model_output <- foreach(i = 1:length(jobs_list)) %dopar% {
    output_i <- evaluate_decaymod(jobs_list[[i]], data, initial_mass, removal_mass, time, group)
    return(output_i)
  }

  return(model_output)
}

# call the function to fit the model, extracts the fit, diagnostics, and negative loglikelihood
evaluate_decaymod <- function (df, data, initial_mass, removal_mass, time, group) {

  if ('cv_cluster' %in% colnames(df)) {
    cross_validation <- TRUE
    group_id <- df$cv_cluster
  } else {
    cross_validation <- FALSE
    group_id <- NA
  }

  fit <- decaymod(data = data,
                  initial_mass = initial_mass,
                  removal_mass = removal_mass,
                  time = time,
                  group = group,
                  model_type = df$model_type,
                  random_effects = df$random_effects,
                  cross_validation = cross_validation,
                  group_id = group_id,
                  formula_k = df$formula_k,
                  formula_alpha = df$formula_alpha,
                  formula_beta = df$formula_beta)

  # make dataframe of the test data predicted v. real
  # mT_pred_wide <- as.data.frame(fit, 'mT_pred')
  # mT_pred <- reshape2::melt(mT_pred_wide,
  #                           variable.name = 'data_point',
  #                           value.name = 'draw')
  # pred <- data.frame(data_point = mT_pred$data_point,
  #                    draw = mT_pred$draw,
  #                    stringsAsFactors = FALSE)
  #
  # mT_real = rep(job$test[, mass], each = nrow(mT_pred) / nrow(job$test)),

  # diagnostics
  fit_summary <- rstan::summary(fit)$summary
  abs_rhat <- max(abs(fit_summary[,'Rhat'] - 1))
  neff_min <- min(fit_summary[,'n_eff'])
  sampler_params <- rstan::get_sampler_params(fit, inc_warmup = FALSE)
  sum_div <- sum(sapply(sampler_params, function(x) sum(x[, "divergent__"])))
  max_treedepth <- sapply(sampler_params, function(x) max(x[, "treedepth__"]))
  diagnostics <- data.frame(abs_rhat = abs_rhat,
                            neff_min = neff_min,
                            sum_div = sum_div,
                            max_treedepth = max_treedepth,
                            stringsAsFactors = FALSE)

  list <- unlist(df, recursive = FALSE)

  neg_loglik_df <- as.data.frame(fit, 'neg_loglik')
  neg_loglik <- mean(neg_loglik_df$neg_loglik, na.rm = TRUE)

  fit_list <- list(mod_specs = list,
                   fit = fit,
                   neg_loglik = neg_loglik,
                   diagnostics = diagnostics)
  fit_list

}

# conducts the model fit itself
decaymod <- function (data, initial_mass, removal_mass, time, group,
                      model_type, random_effects, cross_validation = FALSE, group_id = NA,
                      formula_k = NA, formula_alpha = NA, formula_beta = NA) {

  if (isTRUE(cross_validation)) {
    train <- data[data[, group] != group_id, ]
    test <- data[data[, group] == group_id, ]

    # assign level to each group in the training dataset (as now has one fewer item)
    train$group_level <- as.factor(as.numeric(as.factor(as.character(train[, group]))))

    # number of levels
    n_train_levels <- length(unique(train$group_level))

    if (model_type == 'ne') {
      X <- make_matrix(model_type = 'ne', X_type = 'train', n_train_levels. = n_train_levels,
                       train. = train, formula_k. = formula_k)
      X_test <- make_matrix(model_type = 'ne', X_type = 'test', n_train_levels. = n_train_levels,
                            test. = test, formula_k. = formula_k)

      # for model without random effects
      if (!isTRUE(random_effects)) {
        fit <- ne_CV_noRE_stan(mT = train[, removal_mass],
                               mT_test = test[, removal_mass],
                               m0 = train[, initial_mass],
                               m0_test = test[, initial_mass],
                               time = train[, time],
                               time_test = test[, time],
                               sp = as.numeric(as.factor(train$group_level)),
                               J = n_train_levels,
                               X = X,
                               X_test = X_test)
      }

      # for model with random effects
      if (isTRUE(random_effects)) {
        fit <- ne_CV_RE_stan(mT = train[, removal_mass],
                             mT_test = test[, removal_mass],
                             m0 = train[, initial_mass],
                             m0_test = test[, initial_mass],
                             time = train[, time],
                             time_test = test[, time],
                             sp = as.numeric(as.factor(train$group_level)),
                             J = n_train_levels,
                             X = X,
                             X_test = X_test)
      }
    }

    if (model_type == 'w') {
      X_alpha <- make_matrix(model_type = 'w', X_type = 'train', n_train_levels. = n_train_levels,
                             train. = train, formula_alpha. = formula_alpha)
      X_beta <- make_matrix(model_type = 'w', X_type = 'train', n_train_levels. = n_train_levels,
                            train. = train, formula_beta. = formula_beta)
      X_alpha_test <- make_matrix(model_type = 'w', X_type = 'test', n_train_levels. = n_train_levels,
                                  test. = test, formula_alpha. = formula_alpha)
      X_beta_test <- make_matrix(model_type = 'w', X_type = 'test', n_train_levels. = n_train_levels,
                                 test. = test, formula_beta. = formula_beta)

      # for weibull model without random effects
      if (!isTRUE(random_effects)) {
        fit <- w_CV_noRE_stan(mT = train[, removal_mass],
                              mT_test = test[, removal_mass],
                              m0 = train[, initial_mass],
                              m0_test = test[, initial_mass],
                              time = train[, time],
                              time_test = test[, time],
                              sp = as.numeric(as.factor(train$group_level)),
                              J = n_train_levels,
                              X_alpha = X_alpha,
                              X_alpha_test = X_alpha_test,
                              X_beta = X_beta,
                              X_beta_test = X_beta_test)
      }

      # for weibull model with random effects
      if (isTRUE(random_effects)) {
        fit <- w_CV_RE_stan(mT = train[, removal_mass],
                            mT_test = test[, removal_mass],
                            m0 = train[, initial_mass],
                            m0_test = test[, initial_mass],
                            time = train[, time],
                            time_test = test[, time],
                            sp = as.numeric(as.factor(train$group_level)),
                            J = n_train_levels,
                            X_alpha = X_alpha,
                            X_alpha_test = X_alpha_test,
                            X_beta = X_beta,
                            X_beta_test = X_beta_test)
      }
    }
  }

  if(!isTRUE(cross_validation)) {

    train <- data

    # assign level to each group in the training dataset (as now has one fewer item)
    train$group_level <- as.factor(as.numeric(as.factor(as.character(train[, group]))))

    # number of levels
    n_train_levels <- length(unique(train$group_level))

    if (model_type == 'ne') {
      X <- make_matrix(model_type = 'ne', X_type = 'train', n_train_levels. = n_train_levels,
                       train. = train, formula_k. = formula_k)

      if (isTRUE(random_effects)) {
        fit <- ne_noCV_RE_stan(mT = train[, removal_mass],
                               m0 = train[, initial_mass],
                               time = train[, time],
                               sp = as.numeric(as.factor(train$group_level)),
                               J = n_train_levels,
                               X = X)
      }

      if (!isTRUE(random_effects)) {
        fit <- ne_noCV_noRE_stan(mT = train[, removal_mass],
                                 m0 = train[, initial_mass],
                                 time = train[, time],
                                 sp = as.numeric(as.factor(train$group_level)),
                                 J = n_train_levels,
                                 X = X)
      }
    }

    if (model_type == 'w') {
      X_alpha <- make_matrix(model_type = 'w', X_type = 'train', n_train_levels. = n_train_levels,
                             train. = train, formula_alpha. = formula_alpha)
      X_beta <- make_matrix(model_type = 'w', X_type = 'train', n_train_levels. = n_train_levels,
                            train. = train, formula_beta. = formula_beta)

      if (isTRUE(random_effects)) {
        fit <- w_noCV_RE_stan(mT = train[, removal_mass],
                              m0 = train[, initial_mass],
                              time = train[, time],
                              sp = as.numeric(as.factor(train$group_level)),
                              J = n_train_levels,
                              X_alpha = X_alpha,
                              X_beta = X_beta)
      }

      if (!isTRUE(random_effects)) {
        fit <- w_noCV_noRE_stan(mT = train[, removal_mass],
                                m0 = train[, initial_mass],
                                time = train[, time],
                                sp = as.numeric(as.factor(train$group_level)),
                                J = n_train_levels,
                                X_alpha = X_alpha,
                                X_beta = X_beta)
      }
    }
  }

  fit
}

# creates model matrix
make_matrix <- function (model_type, X_type, n_train_levels. = n_train_levels,
                         train. = NULL, test. = NULL, formula_k. = NULL,
                         formula_alpha. = NULL, formula_beta. = NULL) {

  # negative exponential model matrices
  if (X_type == 'train') {
    if (model_type == 'ne') {
      if (formula_k. == '~1' | formula_k. == '~ 1') {
        mat <- as.matrix(model.matrix(as.formula(formula_k.), train.)[1:n_train_levels.,])
      }
      else {
        mat <- unique(model.matrix(as.formula(formula_k.), train.))
      }
    }

    if (model_type == 'w' & !is.null(formula_alpha.)) {

      # create alpha model matrix for training data
      if (formula_alpha. == '~1' | formula_alpha. == '~ 1') {
        mat <- as.matrix(model.matrix(as.formula(formula_alpha.), train.)[1:n_train_levels., ])
      }
      else {
        mat <- unique(model.matrix(as.formula(formula_alpha.), train.))
      }
    }

    if (model_type == 'w' & !is.null(formula_beta.)) {

      # create beta model matrix for training data
      if (formula_beta. == '~1' | formula_beta. == '~ 1') {
        mat <- as.matrix(model.matrix(as.formula(formula_beta.), train.)[1:n_train_levels., ])
      }
      else {
        mat <- unique(model.matrix(as.formula(formula_beta.), train.))
      }
    }
  }

  if (X_type == 'test') {
    if (model_type == 'ne') mat <- unique(model.matrix(as.formula(formula_k.), test.))
    if (model_type == 'w') {
      if (!is.null(formula_alpha.)) mat <- unique(model.matrix(as.formula(formula_alpha.), test.))
      if (!is.null(formula_beta.)) mat <- unique(model.matrix(as.formula(formula_beta.), test.))
    }
  }

  mat
}
