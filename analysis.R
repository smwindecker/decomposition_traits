
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

jobs_list <- create_jobs(model_df, n,
                         folder = 'output/stan/',
                         cv_cluster = 'species_code')

registerDoMC(20)

compile_models(jobs_list, 'mRem', 'mInit', 't', 'species_code')

output_list <- run_models(jobs_list, 'mRem', 'mInit', 't', 'species_code')

outcome <- evaluate_models(output_list,
                           predR2_path = 'output/stan/eval_plots/')


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

model_df <- data.frame(model_type = c('w', 'ne'),
                       fixed_effects = c(FALSE, FALSE),
                       random_effects = c(FALSE, FALSE),
                       cv = c(TRUE, TRUE))
jobs_list <- create_jobs(model_df, sim,
                         folder = 'output/stan/',
                         cv_cluster = 'group_id')
registerDoMC(3)

compile_models(jobs_list, 'log_mean', 'log_i', 't', 'group_id')

output_list <- run_models(jobs_list, 'log_mean', 'log_i', 't', 'group_id')

outcome <- evaluate_models(output_list,
                           predR2_path = 'output/stan/eval_plots/')

