
library(parallel)
library(plyr)
library(dplyr)
library(foreach)
library(doMC)
library(rstan)
library(decay)

# using my own data
n <- read.table('munge/inundated_data.txt', header = TRUE)
n$species_code <- as.character(n$species_code)

traits <- c('N', 'C', 'SLA', 'LDMC', 'HC', 'CL', 'LG')

model_df <- data.frame(model_type = c('w', 'w', 'ne'),
                       fixed_effects = c(TRUE, FALSE, FALSE),
                       random_effects = c(FALSE, TRUE, FALSE),
                       cv = c(TRUE, TRUE, TRUE))



# simulate data to check each model
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
registerDoMC(3)
foreach(i = 1:3) %dopar% { #length(jobs_list)
  run_stan_models(jobs_list[[i]],
                  mass = 'log_mean',
                  initial_mass = 'log_i',
                  time = 't',
                  group_id = 'group_id',
                  compile_model = TRUE)
}

# w
log_mean_w <- log_i - ((t/b)^a)
plot(t, log_mean_w, ylim= c(0, max(log_mean_w)))
mean_w <- exp(log_mean_w)
plot(t, mean_w, ylim = c(0, max(mean_w)))

sim <- data.frame(t = t, log_i = log_i, log_mean = log_mean_w, group_id = group)

model_df <- data.frame(model_type = c('w', 'ne'),
                       fixed_effects = c(FALSE, FALSE),
                       random_effects = c(FALSE, FALSE),
                       cv = c(TRUE, TRUE))
jobs_list <- create_jobs(model_df, sim,
                         folder = 'output/stan/',
                         cv_cluster = 'group_id')
registerDoMC(3)

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
  run_stan_models(jobs_list[[i]],
                  mass = 'log_mean',
                  initial_mass = 'log_i',
                  time = 't',
                  group_id = 'group_id',
                  compile_model = TRUE)
}

# run all jobs using compiled models as pre-fits
output_list <- list()

output_list <- foreach(i = 1:2) %dopar% { #length(jobs_list)
  output_i <- extract_output(jobs_list[[i]], 'log_mean', 'log_i', 't', 'group_id')
  return(output_i)
}




# negative log likelihood results
neg_ll <- dplyr::bind_rows(lapply(1:length(output_list), function(x) {
  return(output_list[[x]]$neg_loglik)
}))

summary_neg_ll <- neg_ll %>%
  group_by(model) %>%
  summarise(mn = mean(mean),
            lwr = quantile(mean, 0.025),
            upr = quantile(mean, 0.975))

# predictions v. real for all folds of model
pred <- dplyr::bind_rows(lapply(1:length(output_list), function(x) {
  return(output_list[[x]]$pred_real)
}))

mean_pred <- pred %>%
  group_by(model, iter, data_point, mT_real) %>%
  summarise(mean = mean(draw)) %>%
  as.data.frame()

for (i in unique(mean_pred$model)) {
  i_pred <- pred[pred$model == i, ]
  i_mean_pred <- mean_pred[mean_pred$model == i, ]
  png(paste0('output/stan/eval_plots/model_', i, '.png'))
  boxplot(draw ~ as.numeric(mT_real), i_pred,
          ylab = 'Posterior predicted distrbitions',
          xlab = 'Real test data')#,
  #axes = FALSE)
  fit <- lm(i_mean_pred$mean ~ as.numeric(i_mean_pred$mT_real))
  legend('topleft', bty = 'n', legend = paste('R2 is', format(summary(fit)$adj.r.squared, digits = 4)))
  dev.off()

}

# diagnostics
diag <- dplyr::bind_rows(lapply(1:length(output_list), function(x) {
  return(output_list[[x]]$diagnostics)
}))

not_converged <- diag[(diag$abs_rhat > 1.1 |
                         diag$neff_min/nrow(n) < 0.001 |
                         diag$sum_div > 0 |
                         diag$max_tree_c1 > 10 |
                         diag$max_tree_c2 > 10 |
                         diag$max_tree_c3 > 10) ,
                      ]

# what's next?









jobs_list <- create_jobs(model_df, n,
                         folder = 'output/stan/',
                         cv_cluster = 'species_code',
                         fixed_effects_list = traits)




all_stan_groups <- NULL

get_stan_group <- function (x) {
  g <- x$stan_group
  c(g, all_stan_groups)
}

all <- unlist(lapply(jobs_list, get_stan_group))

to_compile <- tapply(seq_along(all), all, identity)[unique(all)] %>%
  lapply(., head, 1) %>%
  unlist(.)

# set number of cores
registerDoMC(10)

# run jobs that are first of their kind, fully compiling models
foreach(i = to_compile) %dopar% {
  run_stan_models(jobs_list[[i]], 'mRem', 'mInit', 't', 'species_code', compile_model = TRUE)
}

# run all jobs using compiled models as pre-fits
foreach(i = length(jobs_list)) %dopar% {
  run_stan_models(jobs_list[[i]], 'mRem', 'mInit', 't', 'species_code', compile_model = FALSE)
}
```

