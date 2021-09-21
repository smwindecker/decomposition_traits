## need to model traits, no cv, no re.


create_jobs <- function (model_type,
                         data,
                         cv_cluster = NULL,
                         random_effects = TRUE,
                         fixed_effects = NULL,
                         trait_param = NULL) {

  if (is.null(fixed_effects)) {
    param_formula <- list('~ 1')
  }

  if (!is.null(fixed_effects)) {
    param_formula <- lapply(fixed_effects, mypaste)
  }

  # for rows that are FE and traits on one only
  if (!is.null(cv_cluster) & is.null(trait_param)) {

    clusters <- unique(data[, cv_cluster])
    tmp <- lapply(param_formula, expand_models, clusters)
    tmp_expanded <- unlist(tmp, recursive = FALSE)

    jobs <- lapply(tmp_expanded, model_details,
                   model_type = model_type,
                   random_effects = random_effects,
                   cv = TRUE)
  }

  # not FE traits on one only
  if (is.null(cv_cluster) & is.null(trait_param)) {

    jobs <- lapply(param_formula, model_details,
                   model_type = model_type,
                   random_effects = random_effects,
                   cv = FALSE)
  }

  # for rows that are FE traits on both
  if (!is.null(fixed_effects) & !is.null(trait_param)) {
    unlist_formulas <- unlist(param_formula)
    x <- gtools::combinations(n = length(unlist_formulas), r = 2, v = unlist_formulas, repeats.allowed = TRUE)
    ext_param_formulas <- tapply(x, rep(1:nrow(x), ncol(x)), function(i) i)
  }

  if (!is.null(cv_cluster) & !is.null(trait_param)) {

    clusters <- unique(data[, cv_cluster])
    tmp <- lapply(ext_param_formulas, expand_models, clusters, both = TRUE)
    tmp_expanded <- unlist(tmp, recursive = FALSE)

    jobs <- lapply(tmp_expanded, model_details,
                   model_type = model_type,
                   random_effects = random_effects,
                   cv = TRUE)
  }

  if (is.null(cv_cluster) & !is.null(trait_param)) {

    jobs <- lapply(ext_param_formulas, model_details,
                   model_type = model_type,
                   random_effects = random_effects,
                   cv = FALSE)
  }


  return(jobs)
}







sim_alpha <- rep(1, 10000)

t1 <- rep(seq(min(decay[,trt1]), max(decay[,trt1]), length.out = 100), 100)

if (!is.null(trt2)) {
  t2 <- rep(seq(min(decay[,trt2]), max(decay[,trt2]), length.out = 100), each = 100)
  sim_beta <- cbind(sim_alpha, t1, t2)
}

if (is.null(trt2)) {
  sim_beta <- cbind(sim_alpha, t1)
}

decay$group_level <- as.factor(as.numeric(as.factor(as.character(decay[, 'species_code']))))
J <- length(unique(decay$group_level))

X_alpha <- as.matrix(model.matrix(as.formula('~ 1'), decay)[1:J,])

if (is.null(trt2)) formula <- paste('~ 1 +', trt1)
if (!is.null(trt2)) formula <- paste('~ 1 +', trt1, '+', trt2)

X_beta <- unique(model.matrix(as.formula(as.character(formula)), decay))

X_alpha_sim <- as.matrix(model.matrix(as.formula('~ 1'), as.data.frame(sim_beta))[1:10000,])
#X_beta_sim <- unique(model.matrix(as.formula(as.character(beta_formula)), sim_traits))

fit <- w_sim_best_model_stan(mT = decay[, 'mRem'],
                             m0 = decay[, "mInit"],
                             time = decay[, 't'],
                             sp = as.numeric(as.factor(decay$group_level)),
                             J = J,
                             X_alpha = X_alpha,
                             X_beta = X_beta,
                             m0_sim = log(4100),
                             time_sim = 0.7,
                             J_sim = nrow(sim_beta),
                             X_alpha_sim = X_alpha_sim,
                             X_beta_sim = sim_beta,
                             sp_sim = 1:10000)

return(fit)
}


trait_nocvre_model <- function (traits) {

  sim_traits_model <- function (decay, trt1, trt2 = NULL) {

    sim_alpha <- rep(1, 10000)

    t1 <- rep(seq(min(decay[,trt1]), max(decay[,trt1]), length.out = 100), 100)

    if (!is.null(trt2)) {
      t2 <- rep(seq(min(decay[,trt2]), max(decay[,trt2]), length.out = 100), each = 100)
      sim_beta <- cbind(sim_alpha, t1, t2)
    }

    if (is.null(trt2)) {
      sim_beta <- cbind(sim_alpha, t1)
    }

    decay$group_level <- as.factor(as.numeric(as.factor(as.character(decay[, 'species_code']))))
    J <- length(unique(decay$group_level))

    X_alpha <- as.matrix(model.matrix(as.formula('~ 1'), decay)[1:J,])

    if (is.null(trt2)) formula <- paste('~ 1 +', trt1)
    if (!is.null(trt2)) formula <- paste('~ 1 +', trt1, '+', trt2)

    X_beta <- unique(model.matrix(as.formula(as.character(formula)), decay))

    X_alpha_sim <- as.matrix(model.matrix(as.formula('~ 1'), as.data.frame(sim_beta))[1:10000,])
    #X_beta_sim <- unique(model.matrix(as.formula(as.character(beta_formula)), sim_traits))

    fit <- w_sim_best_model_stan(mT = decay[, 'mRem'],
                                 m0 = decay[, "mInit"],
                                 time = decay[, 't'],
                                 sp = as.numeric(as.factor(decay$group_level)),
                                 J = J,
                                 X_alpha = X_alpha,
                                 X_beta = X_beta,
                                 m0_sim = log(4100),
                                 time_sim = 0.7,
                                 J_sim = nrow(sim_beta),
                                 X_alpha_sim = X_alpha_sim,
                                 X_beta_sim = sim_beta,
                                 sp_sim = 1:10000)

    return(fit)
  }
}

trait_cvnore_model <- function (traits) {

  comb_traits <- combn(traits, 2, simplify = FALSE)

  best.model_job <- create_jobs(model_type = 'w',
                                data = decay,
                                random_effects = TRUE,
                                cv_cluster = 'species_code',
                                fixed_effects = comb_traits)

  best.model <- run_models(best.model_job,
                           decay,
                           initial_mass = "mInit",
                           removal_mass = "mRem",
                           time = "t",
                           group = "species_code",
                           n_cores = 4,
                           trait_param = 'beta')
  best.model
}

sim_traits_model <- function (decay, trt1, trt2 = NULL) {

  sim_alpha <- rep(1, 10000)

  t1 <- rep(seq(min(decay[,trt1]), max(decay[,trt1]), length.out = 100), 100)

  if (!is.null(trt2)) {
    t2 <- rep(seq(min(decay[,trt2]), max(decay[,trt2]), length.out = 100), each = 100)
    sim_beta <- cbind(sim_alpha, t1, t2)
  }

  if (is.null(trt2)) {
    sim_beta <- cbind(sim_alpha, t1)
  }

  decay$group_level <- as.factor(as.numeric(as.factor(as.character(decay[, 'species_code']))))
  J <- length(unique(decay$group_level))

  X_alpha <- as.matrix(model.matrix(as.formula('~ 1'), decay)[1:J,])

  if (is.null(trt2)) formula <- paste('~ 1 +', trt1)
  if (!is.null(trt2)) formula <- paste('~ 1 +', trt1, '+', trt2)

  X_beta <- unique(model.matrix(as.formula(as.character(formula)), decay))

  X_alpha_sim <- as.matrix(model.matrix(as.formula('~ 1'), as.data.frame(sim_beta))[1:10000,])
  #X_beta_sim <- unique(model.matrix(as.formula(as.character(beta_formula)), sim_traits))

  fit <- w_sim_best_model_stan(mT = decay[, 'mRem'],
                               m0 = decay[, "mInit"],
                               time = decay[, 't'],
                               sp = as.numeric(as.factor(decay$group_level)),
                               J = J,
                               X_alpha = X_alpha,
                               X_beta = X_beta,
                               m0_sim = log(4100),
                               time_sim = 0.7,
                               J_sim = nrow(sim_beta),
                               X_alpha_sim = X_alpha_sim,
                               X_beta_sim = sim_beta,
                               sp_sim = 1:10000)

  return(fit)
}

sim_over_time <- function (decay, var_trait, beta_formula) {

  # test new sim stan
  if (var_trait == 'N'){
    N <- seq(exp(min(decay$N)), exp(max(decay$N)), length.out = 100)
    C <- rep(exp(ave(decay$C)), length.out = 100)
  }

  if (var_trait == 'C'){
    N <- rep(exp(ave(decay$N)), length.out = 100)
    C <- seq(exp(min(decay$C)), exp(max(decay$C)), length.out = 100)
  }

  sim_alpha <- rep(1, 100)
  sim_beta <- cbind(sim_alpha, log(N), log(C))

  decay$group_level <- as.factor(as.numeric(as.factor(as.character(decay[, 'species_code']))))
  J <- length(unique(decay$group_level))

  X_alpha <- as.matrix(model.matrix(as.formula('~ 1'), decay)[1:J,])
  X_beta <- unique(model.matrix(as.formula(as.character(beta_formula)), decay))

  X_alpha_sim <- as.matrix(model.matrix(as.formula('~ 1'), as.data.frame(sim_beta))[1:100,])
  #X_beta_sim <- unique(model.matrix(as.formula(as.character(beta_formula)), sim_traits))

  time_sim <- seq(0, 0.7, length.out = 100)

  fit <- w_sim_over_time_stan(mT = decay[, 'mRem'],
                              m0 = decay[, "mInit"],
                              time = decay[, 't'],
                              sp = as.numeric(as.factor(decay$group_level)),
                              J = J,
                              X_alpha = X_alpha,
                              X_beta = X_beta,
                              m0_sim = log(4100),
                              time_sim = time_sim,
                              J_sim = nrow(sim_beta),
                              X_alpha_sim = X_alpha_sim,
                              X_beta_sim = sim_beta,
                              sp_sim = 1:100)

  return(fit)
}
