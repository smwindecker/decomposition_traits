extract_diagnostics <- function (data, formula_grid) {

  n_models <- 29*nrow(formula_grid)
  diagnostics <- matrix(NA, n_models, 8)
  diagnostics[, 1] <- rep(1:nrow(formula_grid), each = 29)
  diagnostics[, 2] <- rep(unique(data[,'species_code']), times = nrow(formula_grid))

  for (i in 1:nrow(formula_grid)) {
    for (j in unique(data[,'species_code'])) {

      f <- paste0('R/stanOUTPUT/', i, '_', j, '.rds')
      if(file.exists(f)) {

        fit <- readRDS(f)
        fit_summary <- summary(fit)$summary
        abs_rhat <- max(abs(fit_summary[,'Rhat'] - 1))
        neff_min <- min(fit_summary[,'n_eff'])
        sampler_params <- get_sampler_params(fit, inc_warmup = FALSE)
        sum_div <- sum(sapply(sampler_params, function(x) sum(x[, "divergent__"])))
        max_treedepth <- sapply(sampler_params, function(x) max(x[, "treedepth__"]))
        diagnostics[(diagnostics[, 1] == i & diagnostics[, 2] == j), ] <-
          c(i, j, abs_rhat, neff_min, sum_div, max_treedepth)
      }
    }
  }

  colnames(diagnostics) <- c('model', 'test_sp', 'abs_rhat', 'neff_min', 'sum_div',
                             'max_tree_c1', 'max_tree_c2', 'max_tree_c3')

  diag <- merge(formula_grid, diagnostics)

  write.csv(diag, 'output/cv_models_diagnostics.csv')

  cols <- c('abs_rhat', 'neff_min', 'sum_div', 'max_tree_c1', 'max_tree_c2', 'max_tree_c3')
  diag[,cols] <- apply(diag[,cols], 2, function(x) as.numeric(as.character(x)))

  not_converged <- diag[(diag$abs_rhat > 1.1 |
                           diag$neff_min/nrow(n) < 0.001 |
                           diag$sum_div > 0 |
                           diag$max_tree_c1 > 10 |
                           diag$max_tree_c2 > 10 |
                           diag$max_tree_c3 > 10) ,
                        ]
  not_converged

}
