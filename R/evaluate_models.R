#' Evaulate the contents of the model outputs
#'
#' @param output_list list of output of run_models
#' @param predR2_path file path for where to save pred v. real test data plot
#' @return list with deviance per model type and unconverged models
#' @importFrom dplyr bind_rows
#' @importFrom stats coef lm quantile
#' @importFrom graphics boxplot legend
#' @importFrom grDevices png dev.off
#'
#' @export

evaluate_models <- function(output_list, predR2_path) {

  # negative log likelihood results
  neg_ll <- dplyr::bind_rows(lapply(1:length(output_list), function(x) {
    return(output_list[[x]]$neg_loglik)
  }))

  deviance <- neg_ll %>%
    group_by(model) %>%
    summarise(mn = mean(2*mean),
              lwr = stats::quantile(2*mean, 0.025),
              upr = stats::quantile(2*mean, 0.975))

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
    grDevices::png(paste0(predR2_path, 'model_', i, '.png'))
    graphics::boxplot(draw ~ as.numeric(mT_real), i_pred,
                      ylab = 'Posterior predicted distrbitions',
                      xlab = 'Real test data')#,
    #axes = FALSE)
    fit <- stats::lm(i_mean_pred$mean ~ as.numeric(i_mean_pred$mT_real))
    R2 <- paste('R2 is', format(summary(fit)$adj.r.squared, digits = 4))
    int <- paste('Intercept is', format(stats::coef(fit)["(Intercept)"], digits = 4))
    graphics::legend('topleft', bty = 'n', legend = c(R2,
                                                      int))
    grDevices::dev.off()
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

  outcome_list <- list(deviance = deviance,
                       not_converged = not_converged)

  return(outcome_list)

}
