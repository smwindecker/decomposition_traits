extract_loglik <- function (data, formula_grid) {

  nmod <- 66
  output <- data.frame(model = 1:nmod, lower = NA, mean = NA, upper = NA)
  for (i in 1:nmod) {
    input <- data.frame(species_code = unique(data[,'species_code']), nll_ave = NA)
    for (j in unique(data[,'species_code'])) {
      f <- paste0('R/stanOUTPUT/', i, '_', j, '.rds')
      if(file.exists(f)) {
        fit <- readRDS(f)

        neg_loglik <- unlist(extract(fit, 'neg_loglik'))
        input$nll_ave[input$species_code == j] <- mean(neg_loglik, na.rm = TRUE)
      }
    }
    output$mean[output$model == i] <- mean(input$nll_ave)
    quants <- quantile(input$nll_ave, c(0.025, 0.975))
    output$lower[output$model == i] <- quants[[1]]
    output$upper[output$model == i] <- quants[[2]]
  }
  ll <- merge(formula_grid[1:nmod, ], output)
  ll$deviance <- 2 * ll$mean
  #write.csv(ll, 'output/rank.csv')
  ll
}
extract_loglik <- function (formula_grid) {

  neg_loglik_matrix <- matrix(NA, 36, 29)

  for (i in 1:36) {
    for (j in 1:29) {
      f <- paste0('R/stanOUTPUT/', i, '_', j, '.rds')
      fit <- readRDS(f)
      neg_loglik <- unlist(extract(fit, 'neg_loglik'))
      neg_loglik_matrix[i, j] <- mean(neg_loglik, na.rm = TRUE)
      neg_loglik_ave <- rowMeans(neg_loglik_matrix, na.rm = TRUE)
    }
  }

  ll <- cbind(formula_grid, neg_loglik_ave)
  ll$deviance <- 2 * ll$neg_loglik_ave
  ll
}

extract_loglik <- function (fit, i, j) {
  log_lik <- as.data.frame(extract(fit, 'neg_loglik')[[1]])
  log_lik$sum <- rowSums(log_lik, na.rm = TRUE)
  log_lik_matrix[i, j] <- mean(log_lik$sum)
  return(log_lik_matrix)
}
