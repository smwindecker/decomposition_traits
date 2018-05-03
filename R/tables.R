## Tables

# Produce pca loadings table
pca_loadings <- function (prin, output_file) {

  # isolate first two axes
  loadings <- as.data.frame(prin$loadings[,1:2])

  # expand trait names
  loadings$trait[row.names(loadings) == 'SLA'] <- 'Specific litter area'
  loadings$trait[row.names(loadings) == 'DMC'] <- 'Litter dry matter content'
  loadings$trait[row.names(loadings) == 'N'] <- 'Litter nitrogen'
  loadings$trait[row.names(loadings) == 'C'] <- 'Litter carbon'
  loadings$trait[row.names(loadings) == 'HC'] <- 'Litter hemicelluloses'
  loadings$trait[row.names(loadings) == 'CL'] <- 'Litter cellulose'
  loadings$trait[row.names(loadings) == 'LG'] <- 'Litter lignin'

  loadings <- loadings[,c('trait', 'Comp.1', 'Comp.2')]

  # create xtable
  pca_loadings <- xtable::xtable(loadings)

  print(pca_loadings,
        include.rownames = FALSE,
        include.colnames = FALSE,
        only.contents = TRUE,
        comment = FALSE,
        hline.after = NULL,
        file = output_file)
}

# Produce output trait table
traits_table <- function (traits_df, output_file) {

  # modify growth form labels
  traits_df$gf <- as.character(traits_df$gf)
  traits_df$gf[traits_df$gf == 'G'] <- 'graminoid'
  traits_df$gf[traits_df$gf == 'F'] <- 'forb'
  traits_df$gf[traits_df$gf == 'NV'] <- 'nonvascular'
  traits_df$gf[traits_df$gf == 'S'] <- 'shrub'
  traits_df$gf[traits_df$gf == 'T'] <- 'tree'

  # round trait values
  traits_df$SLA <- sprintf("%.2f", round(traits_df$SLA, 2))
  traits_df$DMC <- sprintf("%.0f", round(traits_df$DMC, 0))
  traits_df[,c('N', 'C')] <- lapply(round(traits_df[,c('N', 'C')], 1), sprintf, fmt = "%.1f")
  traits_df[,c('HC_1', 'HC_2', 'CL', 'LG')] <- round(traits_df[,c('HC_1', 'HC_2', 'CL', 'LG')], 1)

  # for each biomass trait, paste mean and ci of weights together
  for (i in c('HC_1', 'HC_2', 'CL', 'LG')) {

    for (j in unique(as.character(traits_df$species))) {

      if (!is.na(traits_df[(traits_df$wt_type == 'mean' & traits_df$species == j), i])) {

        traits_df[(traits_df$wt_type == 'mean' & traits_df$species == j), i] <- paste0(
          traits_df[(traits_df$wt_type == 'mean' & traits_df$species == j), i], ' [',
          traits_df[(traits_df$wt_type == '2.5%' & traits_df$species == j), i], ', ',
          traits_df[(traits_df$wt_type == '97.5%' & traits_df$species == j), i], ']')

      }
    }
  }

  # subset to only include modified mean rows
  trt_table <- traits_df[traits_df$wt_type == 'mean', ]
  exclude_rows <- c('species_code', 'sp_abrev', 'wt_type')
  trt_table <- trt_table[, -which(names(trt_table) %in% exclude_rows)]

  trt_table <- trt_table[c('species', 'family', 'gf', 'SLA', 'DMC', 'N', 'C', 'HC_1', 'HC_2', 'CL', 'LG')]

  # add italics specification for latex
  trt_table$species <- paste0('\\textit{', trt_table$species, '}')

  # create xtable
  trt_table <- xtable::xtable(trt_table)

  print(trt_table,
        include.rownames = FALSE,
        include.colnames = FALSE,
        only.contents = TRUE,
        comment = FALSE,
        sanitize.text.function = identity,
        hline.after = NULL,
        file = output_file)
}

# prepare table of deviance results
deviance_table <- function (output_list, output_file) {

  # negative log likelihood results
  neg_ll <- dplyr::bind_rows(lapply(1:length(output_list), function(x) {
    return(output_list[[x]]$neg_loglik)
  }))

  deviance <- neg_ll %>%
    dplyr::group_by(model) %>%
    dplyr::summarise(mn = mean(2*mean),
                     lwr = stats::quantile(2*mean, 0.025),
                     upr = stats::quantile(2*mean, 0.975)) %>%
    as.data.frame() %>%
    xtable::xtable(.)

  print(deviance,
        include.rownames = FALSE,
        include.colnames = FALSE,
        only.contents = TRUE,
        comment = FALSE,
        hline.after = NULL,
        file = output_file)
}
