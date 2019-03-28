
library(decaymod)
library(mixchar)
library(magrittr)
library(rstan)
library(memisc)

path <- 'chapters/Chapter5_decomp/'
R.utils::sourceDirectory(paste0(path, 'R/'))

# ==================
# prepare decomposition data

# species <- load_species("chapters/Chapter3_traits/data-raw/species.csv")
#
# combined_traits <- traits_combine(species_data = species,
#                                   leco_data = "chapters/Chapter3_traits/data-raw/leco.csv",
#                                   trait_data = "chapters/Chapter3_traits/data-raw/traits.csv",
#                                   tga_data_folder = "chapters/Chapter3_traits/data-raw/tga/",
#                                   ranseed = 1)
#
# mean_traits <- traits_mean_only(combined_traits)

decay <- prepare_decay_data(initial_weight = paste0(path, "data-raw/pre_weight.csv"),
                            removal_weight = paste0(path, "data-raw/decomp_weight.csv"),
                            trait_data = mean_traits)

decay <- read.csv(paste0(path, 'decay_data.csv'))

# ================
# null models
ne_cv_nore <- readRDS(paste0(path, 'model_output/ne.cv.nore.RDS'))
w_cv_nore <- readRDS(paste0(path, 'model_output/w.cv.nore.RDS'))

# species estimates in null models
ne_nocv_re <- readRDS(paste0(path, 'model_output/ne.nocv.re.RDS'))
w_nocv_re <- readRDS(paste0(path, 'model_output/w.nocv.re.RDS'))

best_re_mod <- best_re('N', 'C')

ne_sim_df <- prep_sim_df(ne_nocv_re)
w_sim_df <- prep_sim_df(w_nocv_re)

# decay_and_params <- param_output(decay, ne_nocv_re, w_nocv_re)

alphas <- format_param(w_nocv_re, 'alpha_fit')
betas <- format_param(w_nocv_re, 'beta_fit')

a_post <- prep_post(w_nocv_re, decay, 'alpha')
b_post <- prep_post(w_nocv_re, decay, 'beta')
k_post <- prep_post(ne_nocv_re, decay, 'k')

all_trait_models_mean <- get_deviance(path)
all_trait_models_median <- get_deviance(path, fun = 'median')

png(paste0(path, 'figs/mod_four.png'), width = 1000, height = 980)
mod_plot_four(decay, w_sim_df, ne_sim_df, alphas, betas,
              species_names = c('C', 'BB', 'CC', 'L'))
dev.off()

png(paste0(path, 'figs/mods_ar.png'), width = 1500, height = 980)
mod_plot_ar(decay, w_sim_df, ne_sim_df, subfig = 'a')
dev.off()

png(paste0(path, 'figs/mods_at.png'), width = 1500, height = 1940)
mod_plot_at(decay, w_sim_df, ne_sim_df, subfig = 'b')
dev.off()

png(paste0(path, 'figs/mods_tda.png'), width = 1500, height = 1460)
mod_plot_tda(decay, w_sim_df, ne_sim_df, subfig = 'c')
dev.off()

png(paste0(path, 'figs/mods_tdr.png'), width = 1500, height = 980)
mod_plot_tdr(decay, w_sim_df, ne_sim_df, subfig = 'd')
dev.off()

png(paste0(path, 'figs/trait_plots.png'))
traits_nocvnore_model(decay)
dev.off()

png('figs/kb_plot.png', width = 1000, height = 980)
kb_param_plot(k_post, b_post)
dev.off()

png(paste0(path, 'figs/sigma_sp.png'), width = 1000, height = 800)
sigma_sp_plot(w_nocv_re, best_re_mod)
dev.off()

png(paste0(path, 'figs/a_plot.png'), width = 1000, height = 980)
a_param_plot(a_post)
dev.off()

png(paste0(path, 'figs/ab_cor_plot.png'), width = 1000, height = 980)
kb_cor(a_post, b_post, 'alpha')
dev.off()

png('figs/kb_cor_plot.png', width = 1000, height = 980)
kb_cor(k_post, b_post, 'kappa')
dev.off()

# png('figs/best_mod_trait_effects.png', width = 1000, height = 980)
# sim_best_traits_plot(decay, sim_best_model)
# dev.off()

# ------ tables -------
param_table(decay, ne_nocv_re, w_nocv_re, paste0(path, 'figs/species_parameters.tex'))

# trait models on w
deviance_trait_models(all_trait_models_mean, paste0(path, 'figs/deviance_mean.tex'))
deviance_trait_models(all_trait_models_median, paste0(path, 'figs/deviance_median.tex'))
deviance_trait_models(all_trait_models_median,
                      paste0(path, 'figs/short_deviance_median.tex'),
                      short = TRUE)



knitr::knit('manuscript.Rnw')
tinytex::pdflatex("manuscript.tex", clean = TRUE)



# IN PROGRESS ================
alpha <- rstan::extract(w.nocv.re$fit, pars = 'alpha_fit')$alpha_fit
beta <- rstan::extract(w.nocv.re$fit, pars = 'beta_fit')$beta_fit

hl <- sapply(1:ncol(alpha), function(i) half_life(alpha[,i], beta[,i]))
mrt <- sapply(1:ncol(alpha), function(i) residence_time(alpha[,i], beta[,i]))



traits <- c('N', 'C')

all_fe <- c(traits, combn(traits, 2, simplify = FALSE))

best.traits.nocv.re_jobs <- create_jobs(model_type = 'w',
                                        data = decay,
                                        random_effects = TRUE,
                                        fixed_effects = all_fe)
best.traits.nocv.re <- run_models(best.traits.nocv.re_jobs[3],
                                  decay,
                                  initial_mass = "mInit",
                                  removal_mass = "mRem",
                                  time = "t",
                                  group = "species_code",
                                  n_cores = 4,
                                  trait_param = 'beta',
                                  save_fit = TRUE)
saveRDS(best.traits.nocv.re, 'model_output/best.traits.nocv.re.RDS')
fit_ss <- rstan::extract(best.traits.nocv.re$fit, permuted = TRUE)
mt_sim <- fit_ss$mT_sim

sim_df <- data.frame(m0_sim = rep(log(4100), 5800),
                     time_sim = rep(seq(0, 0.7, length.out = 200), 29),
                     sp_sim = rep(1:29, each = 200))

plot(sim_df$time_sim, exp(colMean(mt_sim)), pch = 20, cex = .2)
points(decay$t, decay$mass_rem, pch = 20, cex = .4, col = 'red')

sim_NC <- sim_traits_model(decay, 'N', 'C')
fit_nc <- rstan::extract(sim_NC, permuted = TRUE)

time_sim <- seq(0, 0.7, length.out = 200)
mt_nc <- 4100*exp(b)

lines(sim_NC)


# best trait model on ne
ne.traits.cv.nore <- readRDS('model_output/ne.traits.cv.nore.RDS')

# simulate traits from best model - without re, without cv
sim_NC <- sim_traits_model(decay, 'N', 'C')
sim_NLG <- sim_traits_model(decay, 'N', 'LG')
sim_DMC <- sim_traits_model(decay, 'DMC')

sim_NDMC <- sim_traits_model(decay, 'N', 'DMC')
sim_NHC <- sim_traits_model(decay, 'N', 'HC')
sim_C <- sim_traits_model(decay, beta_formula = '~ 1 + C')

sim_plot(decay, sim_NC, 'sim_decay', 'N', 'C')
sim_plot(decay, sim_NC, 'hl', 'N', 'C')
sim_plot(decay, sim_NC, 'mrt', 'N', 'C')

sim_plot(decay, sim_NLG, 'sim_decay', 'N', 'LG')
sim_plot(decay, sim_NLG, 'hl', 'N', 'LG')
sim_plot(decay, sim_NLG, 'mrt', 'N', 'LG')

sim_plot(decay, sim_NDMC, 'sim_decay', 'N', 'DMC')
sim_plot(decay, sim_NDMC, 'hl', 'N', 'DMC')
sim_plot(decay, sim_NDMC, 'mrt', 'N', 'DMC')


plot(sim_NC, pars = c('b_beta[1]', 'b_beta[2]', 'b_beta[3]', 'sigma_obs'))
plot(sim_DMC, pars = c('b_beta[1]', 'b_beta[2]', 'sigma_obs'))

# what do these plots tell me?
# since traits are all linear on beta, linear effect on hl and mrt. in same direction.
# c traits pos. affect beta. high c traits low rate. beta pos. affects hl and mrt.
# low rate is high time at half mass and high mean reisdence time.
# similar parameter effects.
par(mfrow = c(1, 3),
    mar = c(8, 6, 4, 6))
my_plot <- function () {

  sim_plot(decay, sim_NC, 'sim_decay', 'N', 'C')
  sim_plot(decay, sim_NLG, 'sim_decay', 'N', 'LG')
  sim_plot(decay, sim_NDMC, 'sim_decay', 'N', 'DMC')
}
png('test.png', width = 1000, height = 300)
sim_plot(decay, sim_NC, 'sim_decay', 'N', 'C')
dev.off()


# best model with w, traits, and re
w.traits.cv.re <- readRDS('model_output/w.traits.cv.re.RDS')


subset <- readRDS('model_output/beta.traits.cv.nore400.RDS')$mod_specs
subset$param_formula <- as.factor(subset$param_formula)
my_formula <- levels(subset$param_formula)[9]
w.traits.cv.nore <- subset[subset$param_formula == my_formula, ]

# this looks bad. :(
deviance_reduction(ne.cv.nore, w.cv.nore, w.traits.cv.nore, w.traits.cv.re)
#


































plot_df <- data.frame(model = c('n_null', 'w_null', 'best'),
                      neg_loglik = c(11.12, 10.05, 3.26))

levels(plot_df) <- c('n_null', 'w_null', 'best')
plot(plot_df$model, plot_df$neg_loglik, pch = 20, ylim = c(0, 12))

# figure out how to extract best model.
# run it without cv
# run traits on ne
# run two traits on each param weibull






## for others it misses a significant early phase dip.
## beta next to alpha will show that the two aren't perfectly corr.




diagnostics_table(nocv.re, length_data = nrow(decay),
                  'figs/diagnostics.csv')

params_all <- param_estimates(nocv.re_list, decay, summary = FALSE)








fit_df <- values(nocv.re_list, decay, 'mT_fit')


sigma_obs <- as.matrix(nocv.re_list[[1]]$fit, 'sigma_obs')

test <- tapply(params_all, params_all$group_level, mT_sim, sigma_obs)
subset <- params_all[params_all$group_level == '1',]

test <- mT_sim(subset, sigma_obs)

mT_sim <- function (parameter_draws, sigma_obs) {
  alpha_fit <- parameter_draws$value[parameter_draws$param == 'alpha']
  beta_fit <- parameter_draws$value[parameter_draws$param == 'beta']
  sp <- as.numeric(parameter_draws$group_level)
  mT_sim <- weibull_sim_rng(beta_fit, alpha_fit, sp, sigma_obs)
}



































# simulate data
set.seed(123)

i <- rnorm(609, 4111, 156)
log_i <- log(i)
t <- rep(seq(0, 0.7, .035), 29)
a <- rnorm(609, 0.6, .1)
b <- rnorm(609, 1.2, .2)
k <- rnorm(609, 1.3, .4)
group <- rep(1:29, each = 21)

# ne
log_mean_k <- log_i - (k*t)
plot(t, log_mean_k, ylim= c(0, max(log_mean_k)))
mean_k <- i*exp(-k*t)
plot(t, mean_k, ylim = c(0, max(mean_k)))

sim <- data.frame(t = t, log_i = log_i, log_mean = log_mean_k, group_id = group)
model_df <- data.frame(model_type = c('ne'),
                       fixed_effects = c(FALSE),
                       random_effects = c(FALSE),
                       cv = c(TRUE))
jobs_list <- create_jobs(model_df, sim,
                         folder = 'output/stan/',
                         cv_cluster = 'group_id')

# w
log_mean_w <- log_i - ((t/b)^a)
plot(t, log_mean_w, ylim= c(0, max(log_mean_w)))
mean_w <- exp(log_mean_w)
plot(t, mean_w, ylim = c(0, max(mean_w)))

sim <- data.frame(t = t, log_i = log_i, log_mean = log_mean_w, group_id = group)



