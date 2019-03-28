
best_re <- function (trt1, trt2) {

  fe <- list(mod1 = c(trt1, trt2))

  best_re_job <- decaymod::create_jobs(model_type = 'w',
                                       data = decay,
                                       random_effects = TRUE,
                                       fixed_effects = fe[1])

  best_re_mod <- decaymod::run_models(best_re_job,
                                      decay,
                                      initial_mass = "mInit",
                                      removal_mass = "mRem",
                                      time = "t",
                                      group = "species_code",
                                      n_cores = 4,
                                      trait_param = 'beta',
                                      save_fit = TRUE)

}

# functions to run on finished models
residence_time <- function (alpha, beta) {

  mrt <- beta * gamma (1 + (1/alpha))
  names(mrt) <- 'res_times'
  mrt
}

half_life <- function (alpha, beta) {

  hl <- beta * (log(2)) ^ (1/alpha)
  names(hl) <- 'half_life'
  hl
}

prep_sim_df <- function (nocv_re_file) {

  fit <- nocv_re_file$fit

  mt_sim_vals <- sprintf('mT_sim[%s]', 1:5800)
  mt_sim_all <- rstan::summary(fit, pars = mt_sim_vals)$summary

  sim_df <- data.frame(time_sim = rep(seq(0, 0.7, length.out = 200), 29),
                       sp_sim = rep(1:29, each = 200),
                       sim_mean = exp(mt_sim_all[,'mean'])/4100,
                       sim_lower = exp(mt_sim_all[,"2.5%"])/4100,
                       sim_upper = exp(mt_sim_all[,"97.5%"])/4100)
  sim_df
}

# for each species extract all the

prep_sim_alpha_df <- function (nocv_re_file) {

  fit <- nocv_re_file$fit

  alpha_fit_vals <- sprintf('alpha_fit[%s]', 1:29)
  alpha_fit_all <- rstan::summary(fit, pars = alpha_fit_vals)$summary

  sim_df <- data.frame(time_sim = rep(seq(0, 0.7, length.out = 200), 29),
                       sp_sim = rep(1:29, each = 200),
                       sim_mean = exp(mt_sim_all[,'mean'])/4100,
                       sim_lower = exp(mt_sim_all[,"2.5%"])/4100,
                       sim_upper = exp(mt_sim_all[,"97.5%"])/4100)
  sim_df
}


# table of species with parameter estimates and calculations
param_output <- function (decay_data, ne_nocv_re_output, w_nocv_re_output) {

  k <- format_param(ne_nocv_re_output, 'k_fit')
  k$group_level <- as.factor(1:29)

  alphas <- format_param(w_nocv_re_output, 'alpha_fit')
  alphas$group_level <- as.factor(1:29)

  betas <- format_param(w_nocv_re_output, 'beta_fit')
  betas$group_level <- as.factor(1:29)

  params <- decay_data
  params$group_level <- as.factor(as.numeric(as.factor(as.character(params[, 'species_code']))))
  params1 <- merge(params, k, by = 'group_level', all = TRUE, suffixes = c("", ".y"))
  names(params1)[names(params1) == "num_mean"] <- "k_mean"
  names(params1)[names(params1) == "2.5%"] <- "k_lower"
  names(params1)[names(params1) == "97.5%"] <- "k_upper"
  names(params1)[names(params1) == "formatted"] <- "k_form"

  params2 <- merge(params1, alphas, by = 'group_level', all = TRUE, suffixes = c("", ".y"))
  names(params2)[names(params2) == "num_mean"] <- "a_mean"
  names(params2)[names(params2) == "2.5%"] <- "a_lower"
  names(params2)[names(params2) == "97.5%"] <- "a_upper"
  names(params2)[names(params2) == "formatted"] <- "a_form"

  params3 <- merge(params2, betas, by = 'group_level', all = TRUE, suffixes = c("", ".y"))
  names(params3)[names(params3) == "num_mean"] <- "b_mean"
  names(params3)[names(params3) == "2.5%"] <- "b_lower"
  names(params3)[names(params3) == "97.5%"] <- "b_upper"
  names(params3)[names(params3) == "formatted"] <- "b_form"

  params3$half_life <- round(half_life(params3$a_mean, params3$b_mean), 2)
  params3$res_time <- round(residence_time(params3$a_mean, params3$b_mean), 2)
  params4 <- params3[order(params3$b_mean),]

  return(params4)
}


# little function to format three params
format_param <- function (file, param) {

  fit <- file$fit
  values <- as.data.frame(rstan::summary(fit, param)$summary)[, c('mean', '2.5%', '97.5%')]
  values$num_mean <- values$mean
  values[, c('mean', '2.5%', '97.5%')] <- lapply(round(values[, c('mean', '2.5%', '97.5%')], 2), sprintf, fmt = "%.2f")
  values$formatted <- paste0(values$mean, '[',
                             values$`2.5%`, ', ',
                             values$`97.5%`, ']')

  values
}

# prepare posteriors
prep_post <- function (nocv_re, decay_data, parameter) {

  fit <- nocv_re$fit

  if (parameter == 'k') {
    means <- colMeans(rstan::extract(fit, 'k_fit')$k_fit)
    pars <- c('sigma_sp', sprintf('k_fit[%s]', 1:29))
  }

  if (parameter == 'alpha') {
    means <- colMeans(rstan::extract(fit, 'alpha_fit')$alpha_fit)
    pars <- c('sigma_sp_alpha', sprintf('alpha_fit[%s]', 1:29))
  }

  if (parameter == 'beta') {
    means <- colMeans(rstan::extract(fit, 'beta_fit')$beta_fit)
    pars <- c('sigma_sp_beta', sprintf('beta_fit[%s]', 1:29))
  }

  posterior <- as.array(fit)
  post <- posterior[,,pars]

  sp <- as.character(unique(decay_data[, 'sp_abrev']))
  dimnames(post)[3] <- list(parameters = c('sigma', sp))

  df_means <- data.frame(species = sp,
                         mean = as.numeric(means))

  ordered_sp <- df_means[order(df_means$mean), 1]

  library(memisc)
  post_ordered <- reorder(post, 3, ordered_sp)
  post_ordered
}


prep_post_gf <- function (nocv_re, decay_data, parameter) {

  fit <- nocv_re$fit

  if (parameter == 'k') {
    means <- colMeans(rstan::extract(fit, 'k_fit')$k_fit)
    pars <- sprintf('k_fit[%s]', 1:29)
  }

  if (parameter == 'alpha') {
    means <- colMeans(rstan::extract(fit, 'alpha_fit')$alpha_fit)
    pars <- sprintf('alpha_fit[%s]', 1:29)
  }

  if (parameter == 'beta') {
    means <- colMeans(rstan::extract(fit, 'beta_fit')$beta_fit)
    pars <- sprintf('beta_fit[%s]', 1:29)
  }

  posterior <- as.array(fit)
  post <- posterior[,,pars]

  sp <- unique(decay[, c('sp_abrev', 'gf')])
  dimnames(post)[3] <- list(parameters = sp[,'gf'])

  ### try and name by growth forms


  df_means <- data.frame(species = as.character(sp),
                         mean = as.numeric(means))

  ordered_sp <- df_means[order(df_means$mean), 1]

  library(memisc)
  post_ordered <- reorder(post, 3, ordered_sp)
  post_ordered
}

# prepare the long df of predicted values for cv no re models
pred_all <- function (model_type, cv_nore_file) {

  mod_specs_all <- dplyr::bind_rows(lapply(1:length(cv_nore_file), function(x) {
    return(cv_nore_file[[x]]$mod_specs)
  }))

  idx <- cv_nore_file[which(mod_specs_all$model_type == model_type)]

  pred_all <- dplyr::bind_rows(lapply(idx, pred_df))
  pred_all
}

# extract the pred values
pred_df <- function (list_item) {

  fit <- list_item$fit

  mt_pred_df <- as.data.frame(fit, pars = 'mT_pred')
  mt_pred_all <- reshape2::melt(mt_pred_df)

  mt_pred_all
}











values <- function (output_list, decay_data, type) {

  mod_specs_all <- dplyr::bind_rows(lapply(1:length(output_list), function(x) {
    return(output_list[[x]]$mod_specs)
  }))

  w_idx <- output_list[which(mod_specs_all$model_type == 'w')]
  w_fit <- w_idx[[1]]$fit
  w_values <- as.data.frame(rstan::summary(w_fit, type)$summary)[, c('mean', '2.5%', '97.5%')]
  colnames(w_values) <- c('w_mean', 'w_lower', 'w_upper')
  w_values$data_point <- rownames(w_values)

  n_idx <- output_list[which(mod_specs_all$model_type == 'ne')]
  n_fit <- n_idx[[1]]$fit
  n_values <- as.data.frame(rstan::summary(n_fit, type)$summary)[, c('mean', '2.5%', '97.5%')]
  colnames(n_values) <- c('n_mean', 'n_lower', 'n_upper')
  n_values$data_point <- rownames(n_values)

  fit_values <- plyr::join(w_values, n_values, by = 'data_point')
  fit_values[, - which (names(fit_values) %in% 'data_point')] <-
    exp(fit_values[, - which (names(fit_values) %in% 'data_point')])

  df <- cbind(decay_data, fit_values)

  df
}



param_estimates <- function (output_list, decay_data, summary = TRUE) {

  mod_specs_all <- dplyr::bind_rows(lapply(1:length(output_list), function(x) {
    return(output_list[[x]]$mod_specs)
  }))

  decay_data$group_level <- as.factor(as.numeric(as.factor(as.character(decay_data[, 'species_code']))))
  spp <- unique(decay_data[, c('group_level', 'species_code', 'species', 'sp_abrev', 'gf', 'gf_old', 'family')])

  # extract fits
  w_idx <- output_list[which(mod_specs_all$model_type == 'w')]
  w_fit <- w_idx[[1]]$fit
  n_idx <- output_list[which(mod_specs_all$model_type == 'ne')]
  n_fit <- n_idx[[1]]$fit

  k_int <- as.data.frame(rstan::summary(n_fit, 'b')$summary)[, 'mean']
  alpha_int <- as.data.frame(rstan::summary(w_fit, 'b_alpha')$summary)[, 'mean']
  beta_int <- as.data.frame(rstan::summary(w_fit, 'b_beta')$summary)[, 'mean']

  ## all estimates
  # extract w parameter estimates
  w_all <- as.data.frame(w_fit, c('alpha_fit', 'beta_fit'))
  w_long <- reshape2::melt(w_all)

  # n estimates
  n_all <- as.data.frame(n_fit, 'k_fit')
  n_long <- reshape2::melt(n_all)

  # bind estimates
  all_params <- rbind(w_long, n_long)
  all_params$param <- rep(c('alpha', 'beta', 'k'), each = 4000*29)
  all_params$group_level <- rep(rep(seq(1:29), each = 4000), 3)

  # all parameter estimates with species data
  df_all <- plyr::join(spp, all_params, by = 'group_level')
  df_all$k_int <- k_int
  df_all$alpha_int <- alpha_int
  df_all$beta_int <- beta_int

  ## summary
  # w summary values
  alpha_values <- as.data.frame(rstan::summary(w_fit, 'alpha_fit')$summary)[, c('mean', '2.5%', '97.5%')]
  colnames(alpha_values) <- c('alpha_mean', 'alpha_lower', 'alpha_upper')
  beta_values <- as.data.frame(rstan::summary(w_fit, 'beta_fit')$summary)[, c('mean', '2.5%', '97.5%')]
  colnames(beta_values) <- c('beta_mean', 'beta_lower', 'beta_upper')

  # k summary values
  k_values <- as.data.frame(rstan::summary(n_fit, 'k_fit')$summary)[, c('mean', '2.5%', '97.5%')]
  colnames(k_values) <- c('k_mean', 'k_lower', 'k_upper')

  # bind values
  param_values <- cbind(alpha_values, beta_values, k_values)
  param_values$group_level <- seq(1:29)

  # parameter summaries with species data
  df_sum <- plyr::join(spp, param_values, by = 'group_level')
  df_sum$k_int <- k_int
  df_sum$alpha_int <- alpha_int
  df_sum$beta_int <- beta_int

  df_sum$half_life <- half_life(df_sum$alpha_mean, df_sum$beta_mean)
  df_sum$res_time <- residence_time(df_sum$alpha_mean, df_sum$beta_mean)
  df_sum$prop <- df_sum$res_time / df_sum$half_life

  if (isTRUE(summary)) {
    return(df_sum)
  } else {
    return(df_all)
  }
}

old <- function (output_list, decay_data, summary = TRUE) {

  mod_specs_all <- dplyr::bind_rows(lapply(1:length(output_list), function(x) {
    return(output_list[[x]]$mod_specs)
  }))

  decay_data$group_level <- as.factor(as.numeric(as.factor(as.character(decay_data[, 'species_code']))))
  spp <- unique(decay_data[, c('group_level', 'species_code', 'sp_abrev', 'gf')])

  # extract fits
  w_idx <- output_list[which(mod_specs_all$model_type == 'w')]
  w_fit <- w_idx[[1]]$fit
  n_idx <- output_list[which(mod_specs_all$model_type == 'ne')]
  n_fit <- n_idx[[1]]$fit

  k_int <- as.data.frame(rstan::summary(n_fit, 'b')$summary)[, 'mean']
  alpha_int <- as.data.frame(rstan::summary(w_fit, 'b_alpha')$summary)[, 'mean']
  beta_int <- as.data.frame(rstan::summary(w_fit, 'b_beta')$summary)[, 'mean']

  ## all estimates
  # extract w parameter estimates
  w_all <- as.data.frame(w_fit, c('alpha_fit', 'beta_fit'))
  w_long <- reshape2::melt(w_all)

  # n estimates
  n_all <- as.data.frame(n_fit, 'k_fit')
  n_long <- reshape2::melt(n_all)

  # bind estimates
  all_params <- rbind(w_long, n_long)
  all_params$param <- rep(c('alpha', 'beta', 'k'), each = 4000*29)
  all_params$group_level <- rep(rep(seq(1:29), each = 4000), 3)

  # all parameter estimates with species data
  df_all <- plyr::join(spp, all_params, by = 'group_level')
  df_all$k_int <- k_int
  df_all$alpha_int <- alpha_int
  df_all$beta_int <- beta_int

  ## summary
  # w summary values
  alpha_values <- as.data.frame(rstan::summary(w_fit, 'alpha_fit')$summary)[, c('mean', '2.5%', '97.5%')]
  colnames(alpha_values) <- c('alpha_mean', 'alpha_lower', 'alpha_upper')
  beta_values <- as.data.frame(rstan::summary(w_fit, 'beta_fit')$summary)[, c('mean', '2.5%', '97.5%')]
  colnames(beta_values) <- c('beta_mean', 'beta_lower', 'beta_upper')

  # k summary values
  k_values <- as.data.frame(rstan::summary(n_fit, 'k_fit')$summary)[, c('mean', '2.5%', '97.5%')]
  colnames(k_values) <- c('k_mean', 'k_lower', 'k_upper')

  # bind values
  param_values <- cbind(alpha_values, beta_values, k_values)
  param_values$group_level <- seq(1:29)

  # parameter summaries with species data
  df_sum <- plyr::join(spp, param_values, by = 'group_level')
  df_sum$k_int <- k_int
  df_sum$alpha_int <- alpha_int
  df_sum$beta_int <- beta_int

  df_sum$half_life <- half_life(df_sum$alpha_mean, df_sum$beta_mean)
  df_sum$res_time <- residence_time(df_sum$alpha_mean, df_sum$beta_mean)
  df_sum$prop <- df_sum$res_time / df_sum$half_life

  if (isTRUE(summary)) {
    return(df_sum)
  } else {
    return(df_all)
  }
}





