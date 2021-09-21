# decay_data$group_level <- as.factor(as.numeric(as.factor(as.character(decay_data[, 'species_code']))))
# group_level <- unique(decay_data$group_level[decay_data$sp_abrev == sp_abrev])

mod_plot_four <- function (species_data, w_sim_df, n_sim_df, alphas, betas) {
  
  layout(matrix(c(1,2,3,4,5,5), nrow = 3, ncol = 2, byrow = TRUE), heights = c(0.8, 0.8, 0.2))
  par(oma = c(3, 6, 0, 2), mar = c(3, 3, 2, 0))
  
  for (i in 1:29) {
    sp_abrev[i] <- unique(species_data$sp_abrev[species_data$species_code == species_names[i]])
    b <- round(betas$num_mean[idx], 2)
  }
  
  species_data$group_level <- as.factor(as.numeric(as.factor(as.character(species_data[, 'species_code']))))
  
  
  
  label <- function (x, data = species_data, alpha_vals = alphas, beta_vals = betas) {
    # alphas <- list$species_alphas
    # betas <- list$species_betas
    idx <- as.numeric(unique(data$group_level[data$species_code == species_names[x]]))
    a <- round(alpha_vals$num_mean[idx], 2)
    b <- round(beta_vals$num_mean[idx], 2)
    text <- bquote(alpha ~ '=' ~ .(a) ~ ',' ~ beta ~ '=' ~ .(b))
    text
  }
  
  mod_plot(sp_abrev_1, species_data, w_sim_df, n_sim_df)
  axis(side = 2, at = c(0, 0.5, 1), cex.axis = 2.8,
       labels = c(0, 0.5, 1))
  legend_subfig_cex('a', cex = 2.8)
  text(x = 0.7, y = 0.9,
       adj = c(1, NA),
       labels = label(1),
       cex = 2.6)
  
  mod_plot(sp_abrev_2, species_data, w_sim_df, n_sim_df)
  legend_subfig_cex('b', cex = 2.8)
  text(x = 0.7, y = 0.9,
       adj = c(1, NA),
       labels = label(2),
       cex = 2.6)
  
  mod_plot(sp_abrev_3, species_data, w_sim_df, n_sim_df)
  axis(side = 1, at = c(0, 0.35, 0.7), cex.axis = 2.8, labels = c(0, 0.35, 0.7))
  axis(side = 2, at = c(0, 0.5, 1), cex.axis = 2.8,
       labels = c(0, 0.5, 1))
  legend_subfig_cex('c', cex = 2.8)
  text(x = 0.7, y = 0.9,
       adj = c(1, NA),
       labels = label(3),
       cex = 2.6)
  
  mod_plot(sp_abrev_4, species_data, w_sim_df, n_sim_df)
  axis(side = 1, at = c(0, 0.35, 0.7), cex.axis = 2.8, labels = c(0, 0.35, 0.7))
  legend_subfig_cex('d', cex = 2.8)
  text(x = 0.7, y = 0.9,
       adj = c(1, NA),
       labels = label(4),
       cex = 2.6)
  
  mtext(text = 'Time (years)',
        side = 1,
        line = 0,
        outer = TRUE,
        cex = 2.5)
  mtext(text = 'Proportion of initial mass remaining (mg)',
        side = 2,
        line = 2,
        outer = TRUE,
        cex = 2.5,
        adj = 0.65)
  
  legend_decay_curves_horizontal(cex = 2.2)
  
}

# All species decay curves
mod_plot_all <- function (species_data, w_sim_df, n_sim_df, alphas, betas) {
  
  layout(matrix(c(1,2,3,4,5,6,7,8,9,10,
                  11,12,13,14,15,16,17,18,19,20,
                  21,22,23,24,25,26,27,28,29,30), 
                nrow = 10, ncol = 3, byrow = TRUE), 
         heights = c(0.8, 0.8, 0.8, 0.8, 0.8, 
                     0.8, 0.8, 0.8, 0.8, 0.8))
  par(oma = c(3, 6, 0, 2), mar = c(3, 3, 2, 0))
  
  

  label <- function (x, data = species_data, alpha_vals = alphas, beta_vals = betas) {
    idx <- as.numeric(unique(data$group_level[data$species_code == species_names[x]]))
    a <- round(alpha_vals$num_mean[idx], 2)
    b <- round(beta_vals$num_mean[idx], 2)
    text <- bquote(alpha ~ '=' ~ .(a) ~ ',' ~ beta ~ '=' ~ .(b))
    text
  }
  
  sorted_species <- species_data[order(species_data$species),]
  
  
  
}