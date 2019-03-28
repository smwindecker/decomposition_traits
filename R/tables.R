## Tables

# prepare table of deviance results
deviance_table <- function (output_list, fun = 'mean', type = NULL) {

  # negative log likelihood results
  neg_ll <- dplyr::bind_rows(lapply(1:length(output_list), function(x) {
    return(output_list[[x]]$mod_specs)
  }))

  if (!is.null(type)) {

    if (type == 'alpha') {
      names(neg_ll)[names(neg_ll) == 'param_formula'] <- 'alpha_formula'
      neg_ll$beta_formula <- '~ 1'
    }

    if (type == 'beta') {
      names(neg_ll)[names(neg_ll) == 'param_formula'] <- 'beta_formula'
      neg_ll$alpha_formula <- '~ 1'
    }
  }

  if (fun == 'mean') {
    deviance <- neg_ll %>%
      dplyr::group_by(model_type, alpha_formula, beta_formula) %>%
      dplyr::summarise(mn = mean(2*neg_loglik),
                       lwr = stats::quantile(2*neg_loglik, 0.025),
                       upr = stats::quantile(2*neg_loglik, 0.975)) %>%
      as.data.frame() %>%
      dplyr::arrange(., mn)
  }

  if (fun == 'median') {
    deviance <- neg_ll %>%
      dplyr::group_by(model_type, alpha_formula, beta_formula) %>%
      dplyr::summarise(mn = median(2*neg_loglik),
                       lwr = stats::quantile(2*neg_loglik, 0.025),
                       upr = stats::quantile(2*neg_loglik, 0.975)) %>%
      as.data.frame() %>%
      dplyr::arrange(., mn)
  }

  deviance
}

param_table <- function (decay_data, ne_nocv_re_output, w_nocv_re_output, output_file) {

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
  species_params <- params3[order(params3$b_mean),]

  # modify growth form labels
  species_params$gf_old <- as.character(species_params$gf_old)
  species_params$gf_old[species_params$gf_old == 'G'] <- 'Graminoid'
  species_params$gf_old[species_params$gf_old == 'F'] <- 'Forb'
  species_params$gf_old[species_params$species == 'Marsilea_drumondii'] <- 'Fern'
  species_params$gf_old[species_params$gf_old == 'NV'] <- 'Non-vascular'
  species_params$gf_old[species_params$gf_old == 'S'] <- 'Shrub'
  species_params$gf_old[species_params$gf_old == 'T'] <- 'Tree'

  # add italics latex code
  species_params$species <- gsub('_', '\\ ', species_params$species)
  species_params$species <- paste0('\\textit{', species_params$species, '}')

  # exclude_rows <- c('group_level', 'species_code', 'sp_abrev')
  # out_table <- species_params[, -which(names(species_params) %in% exclude_rows)]
  out_table <- unique(species_params[c('family', 'species', 'gf', 'gf_old', 'k_form',
                                       'a_form', 'b_form', 'half_life', 'res_time')])
  out_table <- out_table[with(out_table, order(family, species)),]

  # produce xtable
  table <- xtable::xtable(out_table)

  # save as .tex file
  print(table,
        include.rownames = FALSE,
        include.colnames = FALSE,
        only.contents = TRUE,
        comment = FALSE,
        sanitize.text.function = identity,
        hline.after = NULL,
        file = output_file)
}

x_table <- function (models, output_file) {

  rank_tex <- xtable::xtable(models)

  print(rank_tex,
        include.rownames = FALSE,
        include.colnames = FALSE,
        only.contents = TRUE,
        comment = FALSE,
        hline.after = NULL,
        file = output_file)
}

# extract diagnostics
diagnostics_table <- function (output_list, length_data, output_file) {

  diag <- dplyr::bind_rows(lapply(1:length(output_list), function(x) {
    return(output_list[[x]]$diagnostics)
  }))

  not_converged <- diag[(diag$abs_rhat > 1.1 |
                           diag$neff_min/length_data < 0.001 |
                           diag$sum_div > 0 |
                           diag$max_tree_c1 > 10 |
                           diag$max_tree_c2 > 10 |
                           diag$max_tree_c3 > 10), ]

  write.csv(not_converged, output_file)
}

deviance_trait_models <- function (deviance_df, output_file, short = FALSE) {

  deviance_df <- deviance_df[order(deviance_df$mn),]

  deviance_df$deviance <- paste0(signif(deviance_df$mn, 3), '[',
                                 signif(deviance_df$lwr, 3), ', ',
                                 signif(deviance_df$upr, 3), ']')

  deviance_df <- deviance_df[, c('alpha_formula',
                                 'beta_formula',
                                 'deviance')]

  # create xtable
  dev_table <- xtable::xtable(deviance_df)

  if (isTRUE(short)) {
    dev_table <- xtable::xtable(deviance_df[1:5,])
  }

  print(dev_table,
        include.rownames = FALSE,
        include.colnames = FALSE,
        only.contents = TRUE,
        comment = FALSE,
        sanitize.text.function = identity,
        hline.after = NULL,
        file = output_file)
}
