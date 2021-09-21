## Figures

# Produce parameter simulation of Fraser-Suzuki function
simulate_fraser_suzuki <- function () {

  x <- seq(200, 700)
  height_1 <- mixchar::fs_function(x, 0.004, -0.25, 400, 60)
  height_2 <- mixchar::fs_function(x, 0.006, -0.25, 400, 60)
  height_3 <- mixchar::fs_function(x, 0.008, -0.25, 400, 60)
  h4 <- mixchar::fs_function(x, 0.010, -0.25, 400, 60)

  skew_1 <- mixchar::fs_function(x, 0.010, -0.55, 400, 60)
  skew_2 <- mixchar::fs_function(x, 0.010, -0.25, 400, 60)
  skew_3 <- mixchar::fs_function(x, 0.010, 0.25, 400, 60)
  s4 <- mixchar::fs_function(x, 0.010, 0.55, 400, 60)

  position_1 <- mixchar::fs_function(x, 0.010, -0.25, 350, 60)
  position_2 <- mixchar::fs_function(x, 0.010, -0.25, 400, 60)
  position_3 <- mixchar::fs_function(x, 0.010, -0.25, 450, 60)
  p4 <- mixchar::fs_function(x, 0.010, -0.25, 500, 60)

  width_1 <- mixchar::fs_function(x, 0.010, -0.25, 400, 30)
  width_2 <- mixchar::fs_function(x, 0.010, -0.25, 400, 60)
  width_3 <- mixchar::fs_function(x, 0.010, -0.25, 400, 90)
  w4 <- mixchar::fs_function(x, 0.010, -0.25, 400, 120)

  par(oma = c(5, 3, 0, 2), mar = c(1, 3, 2, 0), mfrow = c(2, 2))

  plot(x, h4, type = 'l', lty = 4, xaxt = 'n', yaxt = 'n', cex = 1.6,
       xlab = '',
       ylab = '')
  axis(side = 2, at = c(0.00, 0.005, 0.009), cex.axis = 1.8,
       labels = c(0.001, 0.005, 0.009))
  lines(x, height_2, lty = 2)
  lines(x, height_3, lty = 3)
  lines(x, height_1, lty = 1)
  legend('topleft', legend = '(a)', bty = 'n', cex = 2)
  legend('topright', legend = c(expression(paste('h = 0.004 C'^'-1')),
                                expression(paste('h = 0.006 C'^'-1')),
                                expression(paste('h = 0.008 C'^'-1')),
                                expression(paste('h = 0.010 C'^'-1')),
                                's = -0.25',
                                'p = 400 C',
                                'w = 60 C'),
         bty = 'n', cex = 1.2,
         lty = c(1, 2, 3, 4, NA, NA, NA)
  )

  plot(x, skew_1, type = 'l', lty = 1, xaxt = 'n', yaxt = 'n', cex = 1.6,
       xlab = '',
       ylab = '')
  lines(x, skew_2, lty = 2)
  lines(x, skew_3, lty = 3)
  lines(x, s4, lty = 4)
  legend('topleft', legend = '(b)', bty = 'n', cex = 2)
  legend('topright', legend = c(expression(paste('h = 0.010 C'^'-1')),
                                's = -0.5',
                                's = -0.25',
                                's = 0.25',
                                's = 0.55',
                                'p = 400 C',
                                'w = 60 C'),
         bty = 'n', cex = 1.2,
         lty = c(NA, 1, 2, 3, 4, NA, NA)
  )

  plot(x, position_1, type = 'l', lty = 1, xaxt = 'n', yaxt = 'n', cex = 1.6,
       xlab = '',
       ylab = '')
  axis(side = 1, at = c(200, 400, 600), cex.axis = 1.8,
       labels = c(200, 400, 600))
  axis(side = 2, at = c(0.00, 0.005, 0.009), cex.axis = 1.8,
       labels = c(0.001, 0.005, 0.009))
  lines(x, position_2, lty = 2)
  lines(x, position_3, lty = 3)
  lines(x, p4, lty = 4)
  legend('topleft', legend = '(c)', bty = 'n', cex = 2)
  legend('topright', legend = c(expression(paste('h = 0.010 C'^'-1')),
                                's = -0.25',
                                'p = 350 C',
                                'p = 400 C',
                                'p = 450 C',
                                'p = 500 C',
                                'w = 60 C'),
         bty = 'n', cex = 1.2,
         lty = c(NA, NA, 1, 2, 3, 4, NA)
  )

  plot(x, width_1, type = 'l', lty = 1, xaxt = 'n', yaxt = 'n', cex = 1.6,
       xlab = '',
       ylab = '')
  axis(side = 1, at = c(200, 400, 600), cex.axis = 1.8,
       labels = c(200, 400, 600))
  lines(x, width_2, lty = 2)
  lines(x, width_3, lty = 3)
  lines(x, w4, lty = 4)
  legend('topleft', legend = '(d)', bty = 'n', cex = 2)
  legend('topright', legend = c(expression(paste('h = 0.010 C'^'-1')),
                                's = -0.25',
                                'p = 400 C',
                                'w = 30 C',
                                'w = 60 C',
                                'w = 90 C',
                                'w = 120'),
         bty = 'n', cex = 1.2,
         lty = c(NA, NA, NA, 1, 2, 3, 4))

  mtext(text = 'Temperature (C)',
        side = 1,
        line = 2.1,
        outer = TRUE,
        cex = 1.9)
  mtext(text = expression(paste('Rate of mass loss (C'^'-1', ')')),
        side = 2,
        line = 0.2,
        outer = TRUE,
        cex = 1.9)
}

# Individual curves for TGA theory explanation figure
tga_theory_plots <- function (tga_data) {

  # read raw TGA
  tmp <- process_raw_tga(tga_data)

  # plot TG curve
  layout(matrix(c(1,2,3), nrow = 1, ncol = 3, byrow = TRUE), heights = c(0.8, 0.2))
  par(oma = c(8, 5, 0, 2), mar = c(3, 6, 3, 3))

  plot(tmp$data$temp_C, tmp$data$mass_T, yaxs = 'i', ylim = c(0, 22), xlim = c(33, 800),
       xaxs = 'i', ylab = 'Mass (mg)', xlab = '', xaxt = 'n', yaxt = 'n',
       pch = 20, cex = 0.3, cex.lab = 3.3)
  axis(side = 1, at = c(200, 400, 600), cex.axis = 2.8,
       labels = c(200, 400, 600),
       padj = 1)
  axis(side = 2, at = c(0, 10, 20), cex.axis = 2.8,
       labels = c(0, 10, 20))
  legend('topleft',
         legend = '(a)',
         bty = 'n',
         cex = 3)

  arrows(x0 = 266, y0 = 19, x1 = 266, y1 = 17, lwd = 2, length = 0.1)
  arrows(x0 = 317, y0 = 15, x1 = 317, y1 = 13, lwd = 2, length = 0.1)
  arrows(x0 = 340, y0 = 11.5, x1 = 340, y1 = 9.5, lwd = 2, length = 0.1)

  # plot DTG curve
  plot(tmp$data$temp_C, tmp$data$deriv, yaxs = 'i', ylim = c(0, 0.01),
       xaxs = 'i', ylab = expression(paste('Rate of mass loss (C'^'-1', ')')), xlab = '', xaxt = 'n', yaxt = 'n',
       pch = 20, cex.lab = 3.3, cex = 0.9)
  axis(side = 1, at = c(200, 400, 600), cex.axis = 2.8,
       labels = c(200, 400, 600),
       padj = 1)
  axis(side = 2, at = c(0, 0.004, 0.008), cex.axis = 2.8,
       labels = c(0, 0.004, 0.008))
  legend('topleft',
         legend = '(b)',
         bty = 'n',
         cex = 3)

  segments(x0 = 40, y0 = 0.0018, x1 = 40, y1 = 0.002, lwd = 3.5, col = 'darkgrey')
  segments(x0 = 40, y0 = 0.002, x1 = 120, y1 = 0.002, lwd = 3.5, col = 'darkgrey')
  segments(x0 = 120, y0 = 0.0018, x1 = 120, y1 = 0.002, lwd = 3.5, col = 'darkgrey')
  text(x = ((120-40)/2+40), y = 0.0024, '1', cex = 3.3, col = 'darkgrey')

  segments(x0 = 120, y0 = 0.0078, x1 = 120, y1 = 0.008, lwd = 3.5, col = 'darkgrey')
  segments(x0 = 120, y0 = 0.008, x1 = 650, y1 = 0.008, lwd = 3.5, col = 'darkgrey')
  segments(x0 = 650, y0 = 0.0078, x1 = 650, y1 = 0.008, lwd = 3.5, col = 'darkgrey')
  text(x = ((650-120)/2+120), y = 0.0084, '2', cex = 3.3, col = 'darkgrey')

  segments(x0 = 650, y0 = 0.0008, x1 = 650, y1 = 0.001, lwd = 3.5, col = 'darkgrey')
  segments(x0 = 650, y0 = 0.001, x1 = 790, y1 = 0.001, lwd = 3.5, col = 'darkgrey')
  segments(x0 = 790, y0 = 0.0008, x1 = 790, y1 = 0.001, lwd = 3.5, col = 'darkgrey')
  text(x = ((790-650)/2+650), y = 0.0014, '3', cex = 3.3, col = 'darkgrey')

  # deconvolve data
  output <- mixchar::deconvolve(tmp, upper_temp = 650, n_peaks = NULL)
  temp <- seq(output$temp_bounds[1], output$temp_bounds[2], length.out = nrow(output$data))
  fit <- output$model_fit
  params <- as.data.frame(summary(fit)$coefficients[,1])

  # plot mixture model outcome on DTG data
  plot(output$data$temp_C, output$data$deriv, yaxs = 'i', ylim = c(0, 0.01),
       ylab = expression(paste('Rate of mass loss (C'^'-1', ')')), xlab = '',
       xaxt = 'n', yaxt = 'n', pch = 20, cex = 0.9, cex.lab = 3.3)
  axis(side = 1, at = c(200, 400, 600), cex.axis = 3,
       labels = c(200, 400, 600),
       padj = 1)
  axis(side = 2, at = c(0, 0.004, 0.008), cex.axis = 3,
       labels = c(0, 0.004, 0.008))
  arrows(x0 = 266, y0 = 0.0062, x1 = 266, y1 = 0.0055, lwd = 3, length = 0.1)
  arrows(x0 = 317, y0 = 0.0087, x1 = 317, y1 = 0.008, lwd = 3, length = 0.1)
  arrows(x0 = 365, y0 = 0.0022, x1 = 365, y1 = 0.0015, lwd = 3, length = 0.1)

  y1 <- mixchar::fs_mixture(temp = temp,
                            height_1 = params['height_1',], skew_1 = params['skew_1',],
                            position_1 = params['position_1',], width_1 = params['width_1',],
                            height_2 = params['height_2',], skew_2 = params['skew_2',],
                            position_2 = params['position_2',], width_2 = params['width_2',],
                            height_3 = params['height_3',], skew_3 = params['skew_3',],
                            position_3 = params['position_3',], width_3 = params['width_3',])

  y2 <- mixchar::fs_function(temp = temp,
                             height = params['height_1',], skew = params['skew_1',],
                             position = params['position_1',], width = params['width_1',])

  y3 <- mixchar::fs_function(temp = temp,
                             height = params['height_2',], skew = params['skew_2',],
                             position = params['position_2',], width = params['width_2',])

  y4 <- mixchar::fs_function(temp = temp,
                             height = params['height_3',], skew = params['skew_3',],
                             position = params['position_3',], width = params['width_3',])

  lines(temp, y1, lty = 1, lwd = 2)
  lines(temp, y2, lty = 3, lwd = 3.5, col = '#440154FF')
  lines(temp, y3, lty = 4, lwd = 3.5, col = '#B8DE29FF')
  lines(temp, y4, lty = 5, lwd = 3.5, col = '#3CBB75FF')

  legend('topright',
         legend = c('DTG data', 'DTG modelled', 'Hemicelluloses', 'Cellulose', 'Lignin'),
         ncol = 1,
         cex = 2.8,
         bty = 'n',
         lty = c(NA, 1, 3, 4, 5),
         pch = c(20, NA, NA, NA, NA),
         col = c('black', 'black', '#440154FF', '#B8DE29FF', '#3CBB75FF'),
         lwd = 2)
  legend('topleft',
         legend = '(c)',
         bty = 'n',
         cex = 3)

  mtext(text = 'Temperature (C)',
        side = 1,
        line = 5,
        outer = TRUE,
        cex = 3)

}


# Produce boxplot
box_plot <- function (df) {

  # create functions to specify how to round
  mfloor <- function (x, base) {
    base*floor(x/base)
  }
  mround <- function (x, base) {
    base*round(x/base)
  }
  mceiling <- function (x, base) {
    base*ceiling(x/base)
  }

  par(oma = c(2, 3, 0, 2), mar = c(4, 6, 1, 1), mfrow = c(3, 3))

  low <- mfloor(min(df$LAM), .05)
  high <- mceiling(max(df$LAM), .05)
  mid <- mround((low + high)/2, .01)
  plot(df$gf, df$LAM, ylab = expression(paste('Litter area per mass (m'^'2', '/g)')),
       xlab = '', yaxt = 'n',
       ylim = c(0.99*low, (high + 0.15*(high-low))),
       cex.axis = 1.8, cex.lab = 2.2)
  points(df$gf, df$LAM)
  axis(side = 2, at = c(low, mid, high), cex.axis = 2,
       labels = sprintf("%.2f", c(low, mid, high)))
  legend('topleft', '(a)', bty = 'n', cex = 2)

  low <- mfloor(min(df$DMC), 5)
  high <- mceiling(max(df$DMC), 5)
  mid <- mround((low + high)/2, 1)
  plot(df$gf, df$DMC, ylab = 'Litter dry matter content (mg/g)',
       xlab = '', yaxt = 'n',
       ylim = c(0.99*low, (high + 0.15*(high-low))),
       cex.axis = 1.8, cex.lab = 2.2)
  points(df$gf, df$DMC)
  axis(side = 2, at = c(low, mid, high), cex.axis = 2,
       labels = sprintf("%.0f", c(low, mid, high)))
  legend('topleft', '(b)', bty = 'n', cex = 2)

  low <- mfloor(min(df$N), .05)
  high <- mceiling(max(df$N), .05)
  mid <- mround((low + high)/2, .01)
  plot(df$gf, df$N, ylab = 'Litter nitrogen content (wt%)',
       yaxt = 'n',
       ylim = c(0.99*low, (high + 0.15*(high-low))),
       cex.axis = 1.8, cex.lab = 2.2)
  points(df$gf, df$N)
  axis(side = 2, at = c(low, mid, high), cex.axis = 2,
       labels = sprintf("%.2f", c(low, mid, high)))
  legend('topleft', '(c)', bty = 'n', cex = 2)

  low <- mfloor(min(df$C), 1)
  high <- mceiling(max(df$C), 1)
  mid <- mround((low + high)/2, 1)
  plot(df$gf, df$C, ylab = 'Litter carbon content (wt%)',
       xlab = '', yaxt = 'n',
       ylim = c(0.99*low, (high + 0.15*(high-low))),
       cex.axis = 1.8, cex.lab = 2.2)
  points(df$gf, df$C)
  axis(side = 2, at = c(low, mid, high), cex.axis = 2,
       labels = sprintf("%.0f", c(low, mid, high)))
  legend('topleft', '(d)', bty = 'n', cex = 2)

  low <- mfloor(min(df$HC), .05)
  high <- mceiling(max(df$HC), .05)
  mid <- mround((low + high)/2, .01)
  plot(df$gf, df$HC, ylab = 'Litter hemicelluloses (wt%)',
       yaxt = 'n',
       ylim = c(0.99*low, (high + 0.15*(high-low))),
       cex.axis = 1.8, cex.lab = 2.2)
  points(df$gf, df$HC)
  axis(side = 2, at = c(low, mid, high), cex.axis = 2,
       labels = sprintf("%.2f", c(low, mid, high)))
  legend('topleft', '(e)', bty = 'n', cex = 2)

  low <- mfloor(min(df$CL), .05)
  high <- mceiling(max(df$CL), .05)
  mid <- mround((low + high)/2, .01)
  plot(df$gf, df$CL, ylab = 'Litter cellulose (wt%)',
       yaxt = 'n',
       ylim = c(0.99*low, (high + 0.15*(high-low))),
       cex.axis = 1.8, cex.lab = 2.2)
  points(df$gf, df$CL)
  axis(side = 2, at = c(low, mid, high), cex.axis = 2,
       labels = sprintf("%.2f", c(low, mid, high)))
  legend('topleft', '(f)', bty = 'n', cex = 2)

  low <- mfloor(min(df$LG), .05)
  high <- mceiling(max(df$LG), .05)
  mid <- mround((low + high)/2, .01)
  plot(df$gf, df$LG, ylab = 'Litter lignin (wt%)',
       yaxt = 'n',
       ylim = c(0.99*low, (high + 0.15*(high-low))),
       cex.axis = 1.8, cex.lab = 2.2)
  points(df$gf, df$LG)
  axis(side = 2, at = c(low, mid, high), cex.axis = 2,
       labels = sprintf("%.2f", c(low, mid, high)))
  legend('topleft', '(g)', bty = 'n', cex = 2)

}

# Produce pair plot of traits
pair_plot <- function (df) {

  panel.cor <- function (x, y, digits = 2, prefix = "", cex.cor = 1.8, ...)
  {
    usr <- par("usr"); on.exit(par(usr))
    par(usr = c(0, 1, 0, 1))
    r <- stats::cor(x, y)
    txt <- format(c(r, 0.123456789), digits = digits)[1]
    txt <- paste0(prefix, txt)
    if(missing(cex.cor)) cex.cor <- 0.7/strwidth(txt)
    text(0.5, 0.5, txt, cex = cex.cor * abs(r))

    p <- stats::cor.test(x, y)$p.value
    if (p < 0.05) sym <- 8
    if (p < 0.01) sym <- c(8,8)
    if (p <0.001) sym <- c(8,8,8)
    if (p < 0.05) legend('topright', legend = '', pch = sym, bty = 'n')
  }

  # Customize upper panel
  upper.panel<-function(x, y){
    points(x, y, xlab = '', ylab = '', cex = 2.2)
    mylm <- lm(y ~ x)
    abline(mylm, col = 'red', cex = 2.2)
    newx <- seq(min(x), max(x), length.out = 500)
    prd <- predict(mylm, newdata = data.frame(x = newx), interval = c('confidence'),
                   level = 0.90, type = 'response')
    lines(newx, prd[, 2], col = 'black', lty = 2, cex = 2.2)
    lines(newx, prd[, 3], col = 'black', lty = 2, cex = 2.2)

  }

  # Create the plot
  par(cex.axis = 2.7)
  pairs(df,
        lower.panel = panel.cor,
        upper.panel = upper.panel,
        cex.labels = 5)

}

# Produce phylo plot
phylo_plot <- function (phylo, tips) {

  tips[] <- scale(tips)

  n <- length(tips)
  PRGn <- c('#762a83', '#af8dc3', '#e7d4e8', '#f7f7f7',
            '#d9f0d3', '#7fbf7b', '#1b7837')

  viri <- c('#440154FF', '#453781FF', '#33638DFF', '#238A8DFF', '#20A387FF', '#55C667FF', '#DCE319FF')

  # change colour scheme
  colors <- grDevices::colorRampPalette(viri)(n)

  phytools::phylo.heatmap(phylo, tips, fsize = c(1.5, 1.5, 1), colors = colors)
}

# Call TGA plot for AR
tga_plot_ar <- function (species_deconvolved_list, species_data, gf) {

  sorted_species <- species_data[order(species_data$species),]
  arp_species <- as.character(unique(sorted_species$species_code[sorted_species$gf == gf]))

  layout(matrix(c(1,2,3,4,5,0,6,6,6), nrow = 3, ncol = 3, byrow = TRUE), heights = c(0.8, 0.8, 0.2))
  par(oma = c(3, 8, 0, 2), mar = c(3, 3, 2, 0))

  tga_plot(arp_species[1], species_deconvolved_list, species_data)
  axis(side = 2, at = c(0.001, 0.005, 0.009), cex.axis = 3.2,
       labels = c(sprintf("%.3f", c(0.001, 0.005, 0.009))))
  tga_plot(arp_species[2], species_deconvolved_list, species_data)
  tga_plot(arp_species[3], species_deconvolved_list, species_data)

  tga_plot(arp_species[4], species_deconvolved_list, species_data)
  axis(side = 2, at = c(0.001, 0.005, 0.009), cex.axis = 3.2,
       labels = c(sprintf("%.3f", c(0.001, 0.005, 0.009))))
  axis(side = 1, at = c(150, 400, 650), cex.axis = 3.2, labels = c(150, 400, 650), padj = 1)
  tga_plot(arp_species[5], species_deconvolved_list, species_data)
  axis(side = 1, at = c(150, 400, 650), cex.axis = 3.2, labels = c(150, 400, 650), padj = 1)

  mtext(text = 'Temperature (C)',
        side = 1,
        line = 1,
        outer = TRUE,
        cex = 3.5)
  mtext(text = expression(paste('Rate of mass loss (-dm/dT) (C'^'-1', ')')),
        side = 2,
        line = 2,
        outer = TRUE,
        cex = 3.5,
        adj = 0.55)

  legend_four_curves_horizontal()

}

# Call TGA plot for AT
tga_plot_at <- function (species_deconvolved_list, species_data, gf) {

  sorted_species <- species_data[order(species_data$species),]
  ate_species <- as.character(unique(sorted_species$species_code[sorted_species$gf == gf]))

  layout(matrix(c(1,2,3,4,5,6,7,8,9,10,0,0,11,11,11), nrow = 5, ncol = 3, byrow = TRUE), heights = c(0.8, 0.8, 0.8, 0.8, 0.2))
  par(oma = c(3, 8, 0, 2), mar = c(3, 3, 2, 0))

  tga_plot(ate_species[1], species_deconvolved_list, species_data)
  axis(side = 2, at = c(0.001, 0.005, 0.009), cex.axis = 3.2,
       labels = c(sprintf("%.3f", c(0.001, 0.005, 0.009))))
  tga_plot(ate_species[2], species_deconvolved_list, species_data)
  tga_plot(ate_species[3], species_deconvolved_list, species_data)

  tga_plot(ate_species[4], species_deconvolved_list, species_data)
  axis(side = 2, at = c(0.001, 0.005, 0.009), cex.axis = 3.2,
       labels = c(sprintf("%.3f", c(0.001, 0.005, 0.009))))
  tga_plot(ate_species[5], species_deconvolved_list, species_data)
  tga_plot(ate_species[6], species_deconvolved_list, species_data)

  tga_plot(ate_species[7], species_deconvolved_list, species_data)
  axis(side = 2, at = c(0.001, 0.005, 0.009), cex.axis = 3.2,
       labels = c(sprintf("%.3f", c(0.001, 0.005, 0.009))))
  tga_plot(ate_species[8], species_deconvolved_list, species_data)
  tga_plot(ate_species[9], species_deconvolved_list, species_data)

  tga_plot(ate_species[10], species_deconvolved_list, species_data)
  axis(side = 1, at = c(150, 400, 650), cex.axis = 3.2, labels = c(150, 400, 650), padj = 1)
  axis(side = 2, at = c(0.001, 0.005, 0.009), cex.axis = 3.2,
       labels = c(sprintf("%.3f", c(0.001, 0.005, 0.009))))

  mtext(text = 'Temperature (C)',
        side = 1,
        line = 1,
        outer = TRUE,
        cex = 3.5)
  mtext(text = expression(paste('Rate of mass loss (-dm/dT) (C'^'-1', ')')),
        side = 2,
        line = 2,
        outer = TRUE,
        cex = 3.5,
        adj = 0.55)

  legend_four_curves_horizontal()

}

# Call TGA plot for Tda
tga_plot_tda <- function (species_deconvolved_list, species_data, gf) {

  sorted_species <- species_data[order(species_data$species),]
  tda_species <- as.character(unique(sorted_species$species_code[sorted_species$gf == gf]))

  layout(matrix(c(1,2,3,4,5,6,7,8,9,10,10,10), nrow = 4, ncol = 3, byrow = TRUE), heights = c(0.8, 0.8, 0.8, 0.2))
  par(oma = c(3, 8, 0, 2), mar = c(3, 3, 2, 0))

  tga_plot(tda_species[1], species_deconvolved_list, species_data)
  axis(side = 2, at = c(0.001, 0.005, 0.009), cex.axis = 3.2,
       labels = c(sprintf("%.3f", c(0.001, 0.005, 0.009))))
  tga_plot(tda_species[2], species_deconvolved_list, species_data)
  tga_plot(tda_species[3], species_deconvolved_list, species_data)

  tga_plot(tda_species[4], species_deconvolved_list, species_data)
  axis(side = 2, at = c(0.001, 0.005, 0.009), cex.axis = 3.2,
       labels = c(sprintf("%.3f", c(0.001, 0.005, 0.009))))
  tga_plot(tda_species[5], species_deconvolved_list, species_data)
  tga_plot(tda_species[6], species_deconvolved_list, species_data)

  tga_plot(tda_species[7], species_deconvolved_list, species_data)
  axis(side = 1, at = c(150, 400, 650), cex.axis = 3.2, labels = c(150, 400, 650), padj = 1)
  axis(side = 2, at = c(0.001, 0.005, 0.009), cex.axis = 3.2,
       labels = c(sprintf("%.3f", c(0.001, 0.005, 0.009))))
  tga_plot(tda_species[8], species_deconvolved_list, species_data)
  axis(side = 1, at = c(150, 400, 650), cex.axis = 3.2, labels = c(150, 400, 650), padj = 1)
  tga_plot(tda_species[9], species_deconvolved_list, species_data)
  axis(side = 1, at = c(150, 400, 650), cex.axis = 3.2, labels = c(150, 400, 650), padj = 1)

  mtext(text = 'Temperature (C)',
        side = 1,
        line = 1,
        outer = TRUE,
        cex = 3.5)
  mtext(text = expression(paste('Rate of mass loss (-dm/dT) (C'^'-1', ')')),
        side = 2,
        line = 2,
        outer = TRUE,
        cex = 3.5,
        adj = 0.55)

  legend_four_curves_horizontal()

}

# Call TGA plot for TDr
tga_plot_tdr <- function (species_deconvolved_list, species_data, gf) {

  sorted_species <- species_data[order(species_data$species),]
  tdr_species <- as.character(unique(sorted_species$species_code[sorted_species$gf == gf]))

  layout(matrix(c(1,2,3,4,5,0,6,6,6), nrow = 3, ncol = 3, byrow = TRUE), heights = c(0.8, 0.8, 0.2))
  par(oma = c(3, 8, 0, 2), mar = c(3, 3, 2, 0))

  tga_plot(tdr_species[1], species_deconvolved_list, species_data)
  axis(side = 2, at = c(0.001, 0.005, 0.009), cex.axis = 3.2,
       labels = c(sprintf("%.3f", c(0.001, 0.005, 0.009))))
  tga_plot(tdr_species[2], species_deconvolved_list, species_data)
  tga_plot(tdr_species[3], species_deconvolved_list, species_data)

  tga_plot(tdr_species[4], species_deconvolved_list, species_data)
  axis(side = 2, at = c(0.001, 0.005, 0.009), cex.axis = 3.2,
       labels = c(sprintf("%.3f", c(0.001, 0.005, 0.009))))
  axis(side = 1, at = c(150, 400, 650), cex.axis = 3.2, labels = c(150, 400, 650), padj = 1)
  tga_plot(tdr_species[5], species_deconvolved_list, species_data)
  axis(side = 1, at = c(150, 400, 650), cex.axis = 3.2, labels = c(150, 400, 650), padj = 1)

  mtext(text = 'Temperature (C)',
        side = 1,
        line = 1,
        outer = TRUE,
        cex = 3.5)
  mtext(text = expression(paste('Rate of mass loss (-dm/dT) (C'^'-1', ')')),
        side = 2,
        line = 2,
        outer = TRUE,
        cex = 3.5,
        adj = 0.55)

  legend_four_curves_horizontal()

}

# Call TGA plot for three emblem species
tga_plot_three <- function (species_deconvoluted_list, species_data, species_names) {

  layout(matrix(c(1,2,3,4,4,4), nrow = 2, ncol = 3, byrow = TRUE), heights = c(0.8, 0.2))
  par(oma = c(3, 6, 0, 2), mar = c(3, 3, 2, 0))

  tga_plot(species_names[1], species_deconvoluted_list, species_data, sp_legend = FALSE)
  axis(side = 1, at = c(150, 400, 650), cex.axis = 2.5, labels = c(150, 400, 650), padj = 1)
  axis(side = 2, at = c(0.001, 0.005, 0.009), cex.axis = 2.5,
       labels = c(sprintf("%.3f", c(0.001, 0.005, 0.009))))
  legend_subfig('a', cex = 2.6)

  tga_plot(species_names[2], species_deconvoluted_list, species_data, sp_legend = FALSE)
  axis(side = 1, at = c(150, 400, 650), cex.axis = 2.5, labels = c(150, 400, 650), padj = 1)
  legend_subfig('b', cex = 2.6)

  tga_plot(species_names[3], species_deconvoluted_list, species_data, sp_legend = FALSE)
  axis(side = 1, at = c(150, 400, 650), cex.axis = 2.5, labels = c(150, 400, 650), padj = 1)
  legend_subfig('c', cex = 2.6)

  mtext(text = 'Temperature (C)',
        side = 1,
        line = 0,
        outer = TRUE,
        cex = 2.5)
  mtext(text = expression(paste('Rate of mass loss (-dm/dT) (C'^'-1', ')')),
        side = 2,
        line = 2,
        outer = TRUE,
        cex = 2.2,
        adj = 0.75)

  legend_three_curves_horizontal()

}

# Plot single species' TGA data
tga_plot <- function (species_code, species_deconvolved_list, species_data, sp_legend = TRUE) {

  x <- species_code

  output <- species_deconvolved_list[[x]]
  spname <- species_data$sp_abrev[species_data$species_code == x][1]

  # extract parameters from mixture model fit
  fit <- output$model_fit
  params <- as.data.frame(summary(fit)$coefficients[,1])

  # temperature bounds for plots
  temp <- seq(output$temp_bounds[1], output$temp_bounds[2], length.out = nrow(output$data))

  # isolate data
  data <- output$data

  # plot
  plot(data$temp_C, data$deriv, yaxs = 'i', ylim = c(0, 0.01),
       ylab = '', xlab = '', xaxt = 'n', yaxt = 'n', pch = 20, cex = 0.3)

  if (output$n_peaks == 4) {

    y1 <- mixchar::fs_mixture(temp = temp,
                              height_1 = params['height_1',], skew_1 = params['skew_1',],
                              position_1 = params['position_1',], width_1 = params['width_1',],
                              height_2 = params['height_2',], skew_2 = params['skew_2',],
                              position_2 = params['position_2',], width_2 = params['width_2',],
                              height_3 = params['height_3',], skew_3 = params['skew_3',],
                              position_3 = params['position_3',], width_3 = params['width_3',],
                              height_0 = params['height_0',], skew_0 = params['skew_0',],
                              position_0 = params['position_0',], width_0 = params['width_0',])

    y5 <- mixchar::fs_function(temp = temp,
                               height = params['height_0',], skew = params['skew_0',],
                               position = params['position_0',], width = params['width_0',])

    lines(temp, y5, lty = 6, lwd = 2.5, col = '#33638DFF')

  }

  if (output$n_peaks == 3) {

    y1 <- mixchar::fs_mixture(temp = temp,
                              height_1 = params['height_1',], skew_1 = params['skew_1',],
                              position_1 = params['position_1',], width_1 = params['width_1',],
                              height_2 = params['height_2',], skew_2 = params['skew_2',],
                              position_2 = params['position_2',], width_2 = params['width_2',],
                              height_3 = params['height_3',], skew_3 = params['skew_3',],
                              position_3 = params['position_3',], width_3 = params['width_3',])
  }

  y2 <- mixchar::fs_function(temp = temp,
                             height = params['height_1',], skew = params['skew_1',],
                             position = params['position_1',], width = params['width_1',])

  y3 <- mixchar::fs_function(temp = temp,
                             height = params['height_2',], skew = params['skew_2',],
                             position = params['position_2',], width = params['width_2',])

  y4 <- mixchar::fs_function(temp = temp,
                             height = params['height_3',], skew = params['skew_3',],
                             position = params['position_3',], width = params['width_3',])

  lines(temp, y1, lty = 1, lwd = 2)
  lines(temp, y2, lty = 3, lwd = 3.5, col = '#440154FF')
  lines(temp, y3, lty = 4, lwd = 3.5, col = '#B8DE29FF')
  lines(temp, y4, lty = 5, lwd = 3.5, col = '#3CBB75FF')

  if (isTRUE(sp_legend)) {legend_species_tga(spname)}
}

# Produce pca plot and loadings table
pca <- function (prin, df, species_data) {

  # first two axes' scores
  pc12 <- prin$scores[, 1:2]
  df_pc12 <- data.frame(pc12)
  df_pc12$sp_abrev <- rownames(df_pc12)

  # label with abreviations
  pc12_labeled <- merge(df_pc12, species_data[,c('sp_abrev', 'sp_pca_label', 'gf', 'gf_old')])
  rownames(pc12_labeled) <- pc12_labeled[,'sp_pca_label']

  # get length of axes
  fit <- vegan::envfit(pc12, na.omit(df))

  vars <- prin$sdev^2
  prop_vars <- vars/sum(vars)

  par(oma = c(2, 2, 0, 2))
  plot(pc12_labeled[, c('Comp.1', 'Comp.2')], ylab = '', xlab = '', xaxt = 'n', yaxt = 'n',
       ylim = c(-4, 4), xlim = c(-4.2, 5),
       cex.axis = 1, cex = 1.8, pch = c(2, 20, 8, 0)[as.numeric(pc12_labeled$gf)])
  plot(fit, cex = 2, col = 1, labels = list(vectors = c('LAM', 'DMC', 'N', 'C', 'HC', 'CL', 'LG')))

  mtext(text = paste0('Axis 1 (', (100*round(prop_vars[[1]], 2)), '%)'),
        side = 1,
        line = 0,
        outer = TRUE,
        cex = 2)
  mtext(text = paste0('Axis 2 (', (100*round(prop_vars[[2]], 2)), '%)'),
        side = 2,
        line = 0,
        outer = TRUE,
        cex = 2)

  axis(side = 1, at = c(-4, -2, 0, 2, 4), cex.axis = 1.8, labels = c(-4, -2, 0, 2, 4))
  axis(side = 2, at = c(-4, -2, 0, 2, 4), cex.axis = 1.8, labels = c(-4, -2, 0, 2, 4))
  legend(-4, 4,
         c('Amph. fluct-responders', 'Amph. fluct-tolerators', 'Terr. damp', 'Terr. dry'),
         bty = 'n',
         pch = c(2, 20, 8, 0),
         cex = 1.8)

  reg_labels <- pc12_labeled[!rownames(pc12_labeled) %in% c('Sph',
                                                            'P.pro',
                                                            'A.den',
                                                            'M.cri',
                                                            'L.aus',
                                                            'J.ama',
                                                            'P.dis'), ]

  text(x = reg_labels[, 'Comp.1'], y = reg_labels[, 'Comp.2'],
       labels = row.names(reg_labels), vfont = c('sans serif', 'bold italic'),
       cex = 1.5, pos = 4, col = 'black')

  up_labels <- pc12_labeled[c('A.den', 'M.cri', 'L.aus', 'J.ama', 'P.dis'), ]

  text(x = up_labels[, 'Comp.1'], y = up_labels[, 'Comp.2']+0.16,
       labels = row.names(up_labels), vfont = c('sans serif', 'bold italic'),
       cex = 1.5, pos = 4, col = 'black')

  down_labels <- pc12_labeled[c('Sph', 'P.pro'), ]

  text(x = down_labels[, 'Comp.1'], y = down_labels[, 'Comp.2']-0.16,
       labels = row.names(down_labels), vfont = c('sans serif', 'bold italic'),
       cex = 1.5, pos = 4, col = 'black')
}

## Custom legends

# Species name legend
legend_species_tga <- function (spname) {
  legend(650, .01,
         xjust = 1,
         legend = spname,
         text.font = 3,
         cex = 4.2,
         bty = 'n')
}

# Four curves legend
legend_four_curves <- function () {
  legend(120, 0.008,
         xjust = 0,
         legend = c('data', 'total DTG', 'Hemicelluloses-1', 'Hemicelluloses-2', 'Cellulose', 'Lignin'),
         ncol = 2,
         cex = 1.8,
         bty = 'n',
         lty = c(NA, 1, 6, 3, 4, 5),
         pch = c(20, NA, NA, NA, NA, NA),
         col = c('black', 'black', '#33638DFF', '#440154FF', '#B8DE29FF', '#3CBB75FF'),
         lwd = 2)
}

# Subfigure legend
legend_subfig <- function (subfig, cex = 2.2) {
  legend('topleft', paste0('(', subfig, ')'), bty = 'n', cex = cex)
}

# Three four horizontal legend
legend_four_curves_horizontal <- function () {

  plot(1, type = 'n', axes = FALSE, xlab = '', ylab = '')

  legend(x = "top", inset = 0,
         legend = c('DTG data', 'DTG modelled', 'Hemicelluloses-1', 'Hemicelluloses-2', 'Cellulose', 'Lignin'),
         horiz = TRUE,
         cex = 2.8,
         bty = 'n',
         lty = c(NA, 1, 6, 3, 4, 5),
         pch = c(20, NA, NA, NA, NA, NA),
         col = c('black', 'black', '#33638DFF', '#440154FF', '#B8DE29FF', '#3CBB75FF'),
         lwd = 2)
}

# Three three horizontal legend
legend_three_curves_horizontal <- function () {

  plot(1, type = 'n', axes = FALSE, xlab = '', ylab = '')

  legend(x = "top", inset = 0,
         legend = c('DTG data', 'DTG modelled', '', 'Hemicelluloses', 'Cellulose', 'Lignin'),
         horiz = TRUE,
         cex = 2.8,
         bty = 'n',
         lty = c(NA, 1, NA, 3, 4, 5),
         pch = c(20, NA, NA, NA, NA, NA),
         col = c('black', 'black', 'black', '#440154FF', '#B8DE29FF', '#3CBB75FF'),
         lwd = 2)
}

# Deconvolve raw materials plot
tga_raw_plots <- function (item_1, item_2) {

  layout(matrix(c(1,2,3,3), nrow = 2, ncol = 2, byrow = TRUE), heights = c(0.8, 0.2))
  par(oma = c(4, 4, 0, 2), mar = c(3, 6, 2, 0))

  plot(item_1$temp, item_1$obs,
       xlab = '',
       ylab = '',
       yaxs = 'i',
       yaxt = 'n',
       xaxt = 'n',
       ylim = c(0, 0.012),
       pch = 20,
       cex = 0.3)
  axis(side = 1, at = c(200, 400, 600), cex.axis = 2.2,
       labels = c(200, 400, 600), padj = 1)
  axis(side = 2, at = c(0.002, 0.006, 0.010), cex.axis = 2.2,
       labels = c(0.002, 0.006, 0.010))

  y1 <- mixchar::fs_function(item_1$temp,
                             item_1$height,
                             item_1$skew,
                             item_1$position,
                             item_1$width)

  lines(item_1$temp, y1, lty = 1, lwd = 3)

  legend_subfig('a', cex = 2.7)
  legend('topright',
         legend = c(paste('h =', round(item_1$height, digits = 4)),
                    paste('s =', round(item_1$skew, digits = 3)),
                    paste('p =', round(item_1$position, digits = 0)),
                    paste('w =', round(item_1$width, digits = 0))),
         bty = 'n',
         ncol = 1,
         cex = 2.2)

  plot(item_2$temp, item_2$obs,
       xlab = '',
       ylab = '',
       yaxs = 'i',
       yaxt = 'n',
       xaxt = 'n',
       ylim = c(0, 0.012),
       pch = 20,
       cex = 0.3)
  axis(side = 1, at = c(200, 400, 600), cex.axis = 2.2,
       labels = c(200, 400, 600), padj = 1)

  y2 <- mixchar::fs_function(item_2$temp,
                             item_2$height,
                             item_2$skew,
                             item_2$position,
                             item_2$width)

  lines(item_2$temp, y2, lty = 1, lwd = 3)

  legend_subfig('b', cex = 2.7)
  legend('topright',
         legend = c(paste('h =', round(item_2$height, digits = 4)),
                    paste('s =', round(item_2$skew, digits = 3)),
                    paste('p =', round(item_2$position, digits = 0)),
                    paste('w =', round(item_2$width, digits = 0))),
         bty = 'n',
         ncol = 1,
         cex = 2.2)

  mtext(text = 'Temperature (C)',
        side = 1,
        line = 0,
        outer = TRUE,
        cex = 2.7)
  mtext(text = expression(paste('Rate of mass loss (C'^'-1', ')')),
        side = 2,
        line = 0,
        outer = TRUE,
        cex = 2.7,
        adj = 0.75)

  # empty plot to get the legend on the bottom
  plot(1, type = 'n', axes = FALSE, xlab = '', ylab = '')
  legend(x = "top", inset = 0,
         legend = c('DTG data', 'DTG modelled'),
         horiz = TRUE,
         cex = 2.5,
         bty = 'n',
         lty = c(NA, 1),
         pch = c(20, NA),
         lwd = 2)
}

## Figures

# Produce parameter simulation of Fraser-Suzuki function
simulate_weibull <- function () {

  # simulated model results with 'real' mean parameter values
  x <- seq(0, 0.7, length.out = 1000)

  a.low <- 0.3
  a.high <- 1
  b.low <- 0.3
  b.high <- 2

  layout(matrix(c(1,2,3,4), nrow = 2, ncol = 2, byrow = TRUE), heights = c(0.8, 0.8))
  par(oma = c(5, 5, 0, 2), mar = c(3, 3, 2, 0))

  # first plot is high alpha, low beta
  simulate_plot(a.high, b.low, 'a')
  # axis(side = 2, at = c(0, 0.5, 1), cex.axis = 2.8,
  #      labels = c(0, 0.5, 1))

  # second plot is high alpha, high beta
  simulate_plot(a.high, b.high, 'b')

  # third plot is low alpha, low beta
  simulate_plot(a.low, b.low, 'c')
  # axis(side = 2, at = c(0, 0.5, 1), cex.axis = 2.8,
  #      labels = c(0, 0.5, 1))
  # axis(side = 1, at = c(0, 0.35, 0.7), cex.axis = 2.8, labels = c(0, 0.35, 0.7))

  # fourth plot is low alpha, high beta
  simulate_plot(a.low, b.high, 'd')
  # axis(side = 1, at = c(0, 0.35, 0.7), cex.axis = 2.8, labels = c(0, 0.35, 0.7))

  mtext(text = 'Beta',
        side = 1,
        line = 3,
        outer = TRUE,
        cex = 2.8)
  mtext(text = 'Alpha',
        side = 2,
        line = 0.9,
        outer = TRUE,
        cex = 2.8)
}

# Produce single curve parameter simulation of Fraser-Suzuki function
simulate_plot <- function (alpha, beta, subfig, top = TRUE) {

  # simulated model results with 'real' mean parameter values
  x <- seq(0, 0.7, length.out = 1000)

  y <- weibull(1, x, alpha, beta)
  hl <- round(half_life(alpha, beta), 2)

  text_alpha <- bquote(alpha ~ '=' ~ .(alpha))
  text_beta <- bquote(beta ~ '=' ~ .(beta))
  text_hl <- bquote(hl ~ '=' ~ .(hl))

  plot(x, y, type = 'l', lty = 2, ylim = c(0, 1.05), xlim = c(-.02, 0.8),
       xaxt = 'n', yaxt = 'n', xlab = 'Time', ylab = 'Mass', cex.lab = 2, cex = 3)
  lines(x, y, cex = 3)
  #abline(v = hl, lty = 4, lwd = 3)
  legend_subfig_cex(subfig)

  if (isTRUE(top)) {
    x <- 0.65
    y1 <- 1
  }
  if (!isTRUE(top)) {
    x <- 0.03
    y1 <- 0.3
  }

  # text(x = x,
  #      y = y1,
  #      labels = text_alpha,
  #      cex = 2.5,
  #      adj = 0)
  # text(x = x,
  #      y = y1 - 0.1,
  #      labels = text_beta,
  #      cex = 2.5,
  #      adj = 0)

  # text(x = x,
  #      y = y1,
  #      labels = text_hl,
  #      cex = 2.2,
  #      adj = 0)

  # legend(x = x - 0.07,
  #        y = y1 + 0.02,
  #        legend = NA,
  #        lty = 4,
  #        lwd = 3,
  #        bty = 'n')
}

# Weibull function
weibull <- function (m0, t, a, b) {
  m0 * exp(-(t/b)^a)
}

# Individual model plot
mod_plot <- function (sp_abrev, decay_data, w_sim_df, n_sim_df) {

  decay_data$group_level <- as.factor(as.numeric(as.factor(as.character(decay_data[, 'species_code']))))
  group_level <- unique(decay_data$group_level[decay_data$sp_abrev == sp_abrev])

  # real data for the species
  time <- decay_data$t[decay_data$group_level == group_level]
  mr <- decay_data$mass_rem[decay_data$group_level == group_level]
  mi <- decay_data$mass_init[decay_data$group_level == group_level]

  finish <- as.numeric(5800/29*as.numeric(group_level))
  start <- as.numeric(finish - 199)
  time_sim <- w_sim_df$time_sim[1:200]

  plot(time, mr/mi, pch = 20, ylim = c(0, 1.08), xlim = c(0, 0.72), xaxt = 'n', yaxt = 'n', cex = 2)
  lines(time_sim, w_sim_df$sim_mean[start:finish], col = 'red', lwd = 2.5)
  lines(time_sim, w_sim_df$sim_upper[start:finish], col = 'red', lty = 2, lwd = 2.5)
  lines(time_sim, w_sim_df$sim_lower[start:finish], col = 'red', lty = 2, lwd = 2.5)
  lines(time_sim, n_sim_df$sim_mean[start:finish], lwd = 2.5)
  lines(time_sim, n_sim_df$sim_upper[start:finish], lty = 2, lwd = 2.5)
  lines(time_sim, n_sim_df$sim_lower[start:finish], lty = 2, lwd = 2.5)

  legend_species(sp_abrev)

}

# Model plot for AR
mod_plot_ar <- function (species_data, w_sim_df, n_sim_df, subfig) {

  sorted_species <- species_data[order(species_data$species),]
  ar_species <- as.character(unique(sorted_species$sp_abrev[sorted_species$gf == 'AR']))

  layout(matrix(c(1,2,3,4,5,0,6,6,6), nrow = 3, ncol = 3, byrow = TRUE), heights = c(0.8, 0.8, 0.2))
  par(oma = c(3, 8, 0, 2), mar = c(3, 3, 2, 0))

  mod_plot(ar_species[1], species_data, w_sim_df, n_sim_df)
  axis(side = 2, at = c(0, 0.5, 1), cex.axis = 2.8,
       labels = c(0, 0.5, 1))
  legend_subfig_cex(subfig, cex = 3)
  mod_plot(ar_species[2], species_data, w_sim_df, n_sim_df)
  mod_plot(ar_species[3], species_data, w_sim_df, n_sim_df)

  mod_plot(ar_species[4], species_data, w_sim_df, n_sim_df)
  axis(side = 2, at = c(0, 0.5, 1), cex.axis = 2.8,
       labels = c(0, 0.5, 1))
  axis(side = 1, at = c(0, 0.35, 0.7), cex.axis = 2.8, labels = c(0, 0.35, 0.7))
  mod_plot(ar_species[5], species_data, w_sim_df, n_sim_df)
  axis(side = 1, at = c(0, 0.35, 0.7), cex.axis = 2.8, labels = c(0, 0.35, 0.7))

  mtext(text = 'Time (years)',
        side = 1,
        line = 0,
        outer = TRUE,
        cex = 2.8)
  mtext(text = 'Proportion of initial mass remaining (mg)',
        side = 2,
        line = 2,
        outer = TRUE,
        cex = 2.8,
        adj = 0.55)

  legend_decay_curves_horizontal()

}

# Call mod plot for AT
mod_plot_at <- function (species_data, w_sim_df, n_sim_df, subfig) {

  sorted_species <- species_data[order(species_data$species),]
  at_species <- as.character(unique(sorted_species$sp_abrev[sorted_species$gf == 'AT']))

  layout(matrix(c(1,2,3,4,5,6,7,8,9,10,0,0,11,11,11), nrow = 5, ncol = 3, byrow = TRUE), heights = c(0.8, 0.8, 0.8, 0.8, 0.2))
  par(oma = c(3, 8, 0, 2), mar = c(3, 3, 2, 0))

  mod_plot(at_species[1], species_data, w_sim_df, n_sim_df)
  axis(side = 2, at = c(0, 0.5, 1), cex.axis = 2.8,
       labels = c(0, 0.5, 1))
  legend_subfig_cex(subfig, cex = 2.8)
  mod_plot(at_species[2], species_data, w_sim_df, n_sim_df)
  mod_plot(at_species[3], species_data, w_sim_df, n_sim_df)

  mod_plot(at_species[4], species_data, w_sim_df, n_sim_df)
  axis(side = 2, at = c(0, 0.5, 1), cex.axis = 2.8,
       labels = c(0, 0.5, 1))
  mod_plot(at_species[5], species_data, w_sim_df, n_sim_df)
  mod_plot(at_species[6], species_data, w_sim_df, n_sim_df)

  mod_plot(at_species[7], species_data, w_sim_df, n_sim_df)
  axis(side = 2, at = c(0, 0.5, 1), cex.axis = 2.8,
       labels = c(0, 0.5, 1))
  mod_plot(at_species[8], species_data, w_sim_df, n_sim_df)
  mod_plot(at_species[9], species_data, w_sim_df, n_sim_df)

  mod_plot(at_species[10], species_data, w_sim_df, n_sim_df)
  axis(side = 1, at = c(0, 0.35, 0.7), cex.axis = 2.8, labels = c(0, 0.35, 0.7))
  axis(side = 2, at = c(0, 0.5, 1), cex.axis = 2.8,
       labels = c(0, 0.5, 1))

  mtext(text = 'Time (years)',
        side = 1,
        line = 0,
        outer = TRUE,
        cex = 2.8)
  mtext(text = 'Proportion of initial mass remaining (mg)',
        side = 2,
        line = 2,
        outer = TRUE,
        cex = 2.8,
        adj = 0.55)

  legend_decay_curves_horizontal()

}

# Call model plot for Tda
mod_plot_tda <- function (species_data, w_sim_df, n_sim_df, subfig) {

  sorted_species <- species_data[order(species_data$species),]
  tda_species <- as.character(unique(sorted_species$sp_abrev[sorted_species$gf == 'Tda']))

  layout(matrix(c(1,2,3,4,5,6,7,8,0,9,9,9), nrow = 4, ncol = 3, byrow = TRUE), heights = c(0.8, 0.8, 0.8, 0.2))
  par(oma = c(3, 8, 0, 2), mar = c(3, 3, 2, 0))

  mod_plot(tda_species[1], species_data, w_sim_df, n_sim_df)
  axis(side = 2, at = c(0, 0.5, 1), cex.axis = 2.8,
       labels = c(0, 0.5, 1))
  legend_subfig_cex(subfig, cex = 2.8)
  mod_plot(tda_species[2], species_data, w_sim_df, n_sim_df)
  mod_plot(tda_species[3], species_data, w_sim_df, n_sim_df)

  mod_plot(tda_species[4], species_data, w_sim_df, n_sim_df)
  axis(side = 2, at = c(0, 0.5, 1), cex.axis = 2.8,
       labels = c(0, 0.5, 1))
  mod_plot(tda_species[5], species_data, w_sim_df, n_sim_df)
  mod_plot(tda_species[6], species_data, w_sim_df, n_sim_df)

  mod_plot(tda_species[7], species_data, w_sim_df, n_sim_df)
  axis(side = 1, at = c(0, 0.35, 0.7), cex.axis = 2.8, labels = c(0, 0.35, 0.7))
  axis(side = 2, at = c(0, 0.5, 1), cex.axis = 2.8,
       labels = c(0, 0.5, 1))
  mod_plot(tda_species[8], species_data, w_sim_df, n_sim_df)
  axis(side = 1, at = c(0, 0.35, 0.7), cex.axis = 2.8, labels = c(0, 0.35, 0.7))

  mtext(text = 'Time (years)',
        side = 1,
        line = 0,
        outer = TRUE,
        cex = 2.8)
  mtext(text = 'Proportion of initial mass remaining (mg)',
        side = 2,
        line = 2,
        outer = TRUE,
        cex = 2.8,
        adj = 0.55)

  legend_decay_curves_horizontal()

}

# Call mod plot for TDr
mod_plot_tdr <- function (species_data, w_sim_df, n_sim_df, subfig) {

  sorted_species <- species_data[order(species_data$species),]
  tdr_species <- as.character(unique(sorted_species$sp_abrev[sorted_species$gf == 'Tdr']))

  layout(matrix(c(1,2,3,4,5,0,6,6,6), nrow = 3, ncol = 3, byrow = TRUE), heights = c(0.8, 0.8, 0.2))
  par(oma = c(3, 8, 0, 2), mar = c(3, 3, 2, 0))

  mod_plot(tdr_species[1], species_data, w_sim_df, n_sim_df)
  axis(side = 2, at = c(0, 0.5, 1), cex.axis = 2.8,
       labels = c(0, 0.5, 1))
  legend_subfig_cex(subfig, cex = 2.8)
  mod_plot(tdr_species[2], species_data, w_sim_df, n_sim_df)
  mod_plot(tdr_species[3], species_data, w_sim_df, n_sim_df)

  mod_plot(tdr_species[4], species_data, w_sim_df, n_sim_df)
  axis(side = 2, at = c(0, 0.5, 1), cex.axis = 2.8,
       labels = c(0, 0.5, 1))
  axis(side = 1, at = c(0, 0.35, 0.7), cex.axis = 2.8, labels = c(0, 0.35, 0.7))
  mod_plot(tdr_species[5], species_data, w_sim_df, n_sim_df)
  axis(side = 1, at = c(0, 0.35, 0.7), cex.axis = 2.8, labels = c(0, 0.35, 0.7))

  mtext(text = 'Time (years)',
        side = 1,
        line = 0,
        outer = TRUE,
        cex = 2.8)
  mtext(text = 'Proportion of initial mass remaining (mg)',
        side = 2,
        line = 2,
        outer = TRUE,
        cex = 2.8,
        adj = 0.55)

  legend_decay_curves_horizontal()

}

# Call mod plot for three emblem species
mod_plot_four <- function (species_data, w_sim_df, n_sim_df, alphas, betas, species_names) {

  layout(matrix(c(1,2,3,4,5,5), nrow = 3, ncol = 2, byrow = TRUE), heights = c(0.8, 0.8, 0.2))
  par(oma = c(3, 6, 0, 2), mar = c(3, 3, 2, 0))

  sp_abrev_1 <- unique(species_data$sp_abrev[species_data$species_code == species_names[1]])
  sp_abrev_2 <- unique(species_data$sp_abrev[species_data$species_code == species_names[2]])
  sp_abrev_3 <- unique(species_data$sp_abrev[species_data$species_code == species_names[3]])
  sp_abrev_4 <- unique(species_data$sp_abrev[species_data$species_code == species_names[4]])

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
       labels = c(0, 0.5, 1), las = 2)
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
  axis(side = 1, at = c(0, 0.35, 0.7), cex.axis = 2.8,
       labels = c(0, 0.35, 0.7), padj = 0.5)
  axis(side = 2, at = c(0, 0.5, 1), cex.axis = 2.8,
       labels = c(0, 0.5, 1), las = 2)
  legend_subfig_cex('c', cex = 2.8)
  text(x = 0.7, y = 0.9,
       adj = c(1, NA),
       labels = label(3),
       cex = 2.6)

  mod_plot(sp_abrev_4, species_data, w_sim_df, n_sim_df)
  axis(side = 1, at = c(0, 0.35, 0.7), cex.axis = 2.8,
       labels = c(0, 0.35, 0.7), padj = 0.5)
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
        line = 3,
        outer = TRUE,
        cex = 2.5,
        adj = 0.65)

  legend_decay_curves_horizontal(cex = 2.2)

}

# plot posteriors
post_plot <- function (posterior, y_text = NULL) {

  bayesplot::color_scheme_set('darkgray')
  plot <- bayesplot::mcmc_intervals(posterior) +
    ggplot2::theme_bw() +
    ggplot2::theme(panel.grid.major = ggplot2::element_blank(),
                   panel.grid.minor = ggplot2::element_blank(),
                   axis.text.x = ggplot2::element_text(size = 18))

  if (!is.null(y_text)) {
    plot <- plot + ggplot2::theme(axis.text.y = ggplot2::element_text(face = 'italic',
                                                                      size = 18))
  } else {
    plot <- plot + ggplot2::theme(axis.text.y = ggplot2::element_blank())
  }
  plot
}

# weibull posteriors two side by side
a_param_plot <- function (a_post) {

  a_plot <- post_plot(a_post, y_text = 'yes') +
    ggplot2::annotate('text', x = 1.3, y = 27.5, label = "alpha", parse = TRUE, size = 16) +
    ggplot2::geom_vline(xintercept = 1, linetype = 'dashed', colour = 'red', size = 1.5)
  a_plot

}

# correlation plot underlying function
cor_plot <- function (x, y, ...) {

  plot(x, y, cex = 3, cex.axis = 2, ...)
  mylm <- lm(y ~ x)
  abline(mylm, col = 'red', lwd = 3)
  newx <- seq(min(x), max(x), length.out = 500)
  prd <- predict(mylm, newdata = data.frame(x = newx), interval = c('confidence'),
                 level = 0.90, type = 'response')
  lines(newx, prd[, 2], col = 'black', lty = 2, lwd = 3)
  lines(newx, prd[, 3], col = 'black', lty = 2, lwd = 3)
  r <- stats::cor(x, y)
  p <- stats::cor.test(x, y)$p.value
  if (p < 0.05) sym <- 8
  if (p < 0.01) sym <- c(8,8)
  if (p < 0.001) sym <- c(8,8,8)
  if (p > 0.05) sym <- 20
  legend('topright',
         legend = format(c(r, 0.123456789),
                         digits = 2)[1],
         pch = sym, bty = 'n',
         cex = 3)

}

# alpha and beta correlation plot
kb_cor <- function (y_vals, b_post, param) {

  b_ordered_sp <- dimnames(b_post)[3]$parameters
  y_reordered <- reorder(y_vals, 3, b_ordered_sp)

  beta <- apply(b_post, 3, median)
  y <- apply(y_reordered, 3, median)

  if (param == 'kappa') {
    y_lab <- expression(paste('log(', kappa, ')'))
  }
  if (param == 'alpha') {
    y_lab <- expression(paste('log(', alpha, ')'))
  }

  par(mar = c(5, 7, 2, 2))
  cor_plot(log(beta), log(y), cex.lab = 3,
           ylab = y_lab,
           xlab = expression(paste('log(', beta, ')')))

}

# neg exp posteriors
kb_param_plot <- function (k_post, b_post) {

  b_plot <- post_plot(b_post, y_text = 'yes') +
    ggplot2::annotate('text', x = 5.6, y = 27.5, label = "beta", parse = TRUE, size = 16) +
    ggplot2::geom_vline(xintercept = 0, colour = 'white', size = 0.5)

  b_ordered_sp <- dimnames(b_post)[3]$parameters
  k_reordered <- reorder(k_post, 3, b_ordered_sp)
  k_plot <- post_plot(k_reordered) +
    ggplot2::annotate('text', x = 5.6, y = 27.5, label = "kappa", parse = TRUE, size = 16) +
    ggplot2::geom_vline(xintercept = 0, colour = 'white', size = 0.5)

  egg::ggarrange(b_plot, k_plot, nrow = 1)

}

sim_plot <- function (decay, fit, param, trt1, trt2) {

  mt_sim_vals <- sprintf('mT_sim[%s]', 1:10000)
  mt_sim_all <- rstan::summary(fit, pars = mt_sim_vals)$summary[,'mean']

  alpha_sim_vals <- sprintf('alpha_sim[%s]', 1:10000)
  alpha_sim <- unique(rstan::summary(fit, pars = alpha_sim_vals)$summary[,'mean'])

  t1 <- rep(seq(min(decay[,trt1]), max(decay[,trt1]), length.out = 100), 100)
  t2 <- rep(seq(min(decay[,trt2]), max(decay[,trt2]), length.out = 100), each = 100)

  int <- mean(rstan::extract(fit, pars = 'b_beta[1]')$`b_beta[1]`)
  t1_beta <- mean(rstan::extract(fit, pars = 'b_beta[2]')$`b_beta[2]`)
  t2_beta <- mean(rstan::extract(fit, pars = 'b_beta[3]')$`b_beta[3]`)

  beta_sim <- exp(int + t1_beta*t1 + t2_beta*t2)

  sim_df <- data.frame(sim_decay = exp(mt_sim_all)/4100,
                       hl = half_life(alpha_sim, beta_sim),
                       mrt = residence_time(alpha_sim, beta_sim),
                       t1 = t1,
                       t2 = t2)

  test <- reshape2::dcast(sim_df, t2 ~ t1, value.var = param)
  rownames(test) <- test[,1]
  test <- as.matrix(test[,2:101])

  palette <- colorRampPalette(c("black", 'lightgrey'))

  if (param == 'sim_decay') titletext <- 'Proportion\nmass remaining'
  if (param == 'hl') titletext <- 'Time to\nhalf mass'
  if (param == 'mrt') titletext <- 'Mean\nresidence time'

  filled.contour(x = sort(as.numeric(as.character(unique(sim_df$t2)))),
                 y = sort(as.numeric(as.character(unique(sim_df$t1)))),
                 z = test,
                 color.palette = palette,
                 key.title = title(titletext, adj = 0, cex = 1.5),
                 key.axes = axis(4, cex.axis = 1.5),
                 plot.axes = {
                   axis(1, cex.axis = 2)
                   axis(2, cex.axis = 2)
                 })
  points(x = unique(decay[,trt2]), y = unique(decay[,trt1]),
         pch = 20, col = 'white', cex = 2.5)

  full_name <- function (trt) {

    if (trt == 'HC') full_name <- 'Hemicelluloses'
    if (trt == 'CL') full_name <- 'Cellulose'
    if (trt == 'LG') full_name <- 'Lignin'
    if (trt == 'N') full_name <- 'Nitrogen'
    if (trt == 'C') full_name <- 'Carbon'
    if (trt == 'LAM') full_name <- 'Litter area mass'
    if (trt == 'DMC') full_name <- 'Dry matter content'

    return(full_name)

  }

  mtext(text = full_name(trt2),
        side = 1,
        line = 4,
        cex = 2.5)
  mtext(text = full_name(trt1),
        side = 2,
        line = 3,
        cex = 2.5)
}




deviance_reduction <- function (ne.null, w.null, ne.trait, w.trait, w.trait.re) {

  dev_1 <- 2*ne.null$mod_specs$neg_loglik
  dev_2 <- 2*w.null$mod_specs$neg_loglik
  dev_3 <- 2*ne.trait$mod_specs$neg_loglik
  dev_4 <- 2*w.trait$neg_loglik
  dev_5 <- 2*w.trait.re$mod_specs$neg_loglik

  boxplot(dev_1, dev_2, dev_3, dev_4, dev_5,
          labels = c('negative exponential, no traits, no re',
                     'weibull, no traits, no re',
                     'negative exponential, traits, no re',
                     'weibull, traits, no re',
                     'weibull, traits, re'), ylim = c(-30, 30))

}














# pred v. real
pred_real_plot <- function (pred_values_w, pred_values_ne, decay_data) {

  pred_values <- cbind(pred_values_w, pred_values_ne[,2])
  colnames(pred_values) <- c('pred_value', 'w_pred', 'n_pred')
  pred_values$real <- rep(decay_data$mRem, each = 4000)

  ### not sure how this works.
  boxplot(pred_values$w_pred ~ as.numeric(pred_values$real), horizontal = TRUE)

}


sina_plot <- function (decay_data, a_post) {


  # plot posteriors
  bayesplot::color_scheme_set('darkgray')
  plot <- bayesplot::mcmc_intervals(a_post) +
    ggplot2::theme_bw() +
    ggplot2::theme(panel.grid.major = ggplot2::element_blank(),
                   panel.grid.minor = ggplot2::element_blank(),
                   axis.text.x = ggplot2::element_text(size = 18))

  if (!is.null(y_text)) {
    plot <- plot + ggplot2::theme(axis.text.y = ggplot2::element_text(face = 'italic',
                                                                      size = 18))
  } else {
    plot <- plot + ggplot2::theme(axis.text.y = ggplot2::element_blank())
  }
  plot


  a_plot <- post_plot(a_post, y_text = 'yes') +
    ggplot2::annotate('text', x = 1.3, y = 27.5, label = "alpha", parse = TRUE, size = 16) +
    ggplot2::geom_vline(xintercept = 1, linetype = 'dashed', colour = 'red', size = 1.5)
  a_plot

}












# rank parameters
param_rank <- function (estimates, parameter) {

  subset <- estimates[estimates$param == parameter, ]

  if (parameter == 'beta') {
    ordered_sp <- with(subset, reorder(sp_abrev, value, median, na.rm = TRUE))
  }
  if (parameter == 'alpha' | parameter == 'k') {
    ordered_sp <- with(subset, reorder(sp_abrev, -value, median, na.rm = TRUE))
  }

  column <- paste0(parameter, '_int')
  v_value <- subset[1, column]

  if (v_value < min(subset$value)) {
    min <- v_value
  }
  if (min(subset$value) < v_value) {
    min <- min(subset$value)
  }

  par(font.axis = 3, mar = c(3, 7, 2, 2))
  boxplot(subset$value ~ ordered_sp, horizontal = TRUE, las = 1, xaxt = 'n',
          ylab = '', outline = FALSE)
  abline(v = v_value, lty = 2, col = 'red', cex = 2)
  text(x = 0.5, y = 2, labels = bquote(beta ~ '=' ~ .(v_value)), col = 'red', cex = 2, pos = 4)
  axis(side = 1, at = c(0:6), labels = c(0:6), font.axis = 1, cex = 2)

}

sigma_sp_plot <- function (null_mod, trt_mod) {

  par(mar = c(5, 5, 2, 2))
  null_fit <- as.data.frame(null_mod$fit)
  trt_fit <- as.data.frame(trt_mod$fit)
  plot(density(null_fit$sigma_sp_alpha),
       lwd = 3,
       xlim = c(0, 1.7), ylim = c(0, 8.7),
       xaxt = 'n', yaxt = 'n',
       xaxs = 'i', yaxs = 'i',
       bty = 'L', main = '',
       xlab = 'Species-level sigma',
       ylab = 'Density', cex.lab = 2.5)
  lines(density(trt_fit$sigma_sp_alpha), col = 'red', lwd = 3)
  lines(density(null_fit$sigma_sp_beta), lty = 4, lwd = 3)
  lines(density(trt_fit$sigma_sp_beta), lty = 4, col = 'red', lwd = 3)

  axis(side = 1, at = c(0, 0.5, 1, 1.5), labels = c(0, 0.5, 1, 1.5), cex.axis = 2)
  axis(side = 2, las = 2, at = c(0, 4, 8), labels = c(0, 4, 8), cex.axis = 2)

  legend('topright',
         legend = c('Null Weibull', 'Trait Weibull'),
         text.col = c('black', 'red', 'red'),
         bty = 'n',
         cex = 2)

  segments(x0 = 0.15, y0 = 7.8, x1 = 0.15, y1 = 8, lwd = 5, col = 'darkgrey')
  segments(x0 = 0.15, y0 = 8, x1 = 0.4, y1 = 8, lwd = 5, col = 'darkgrey')
  segments(x0 = 0.4, y0 = 7.8, x1 = 0.4, y1 = 8, lwd = 5, col = 'darkgrey')
  text(x = ((0.4-.15)/2+.15), y = 8.4, bquote(alpha), cex = 3, col = 'darkgrey')

  segments(x0 = 0.5, y0 = 3.8, x1 = 0.5, y1 = 4, lwd = 5, col = 'darkgrey')
  segments(x0 = 0.5, y0 = 4, x1 = 1.1, y1 = 4, lwd = 5, col = 'darkgrey')
  segments(x0 = 1.1, y0 = 3.8, x1 = 1.1, y1 = 4, lwd = 5, col = 'darkgrey')
  text(x = ((1.1-0.5)/2+0.5), y = 4.5, bquote(beta), cex = 3, col = 'darkgrey')

}

# Negative exponential function
negexp <- function (m0, t, k) {
  m0 * exp(-k*t)
}

## trait plots
# Produce boxplot
box_plot <- function (df) {

  # create functions to specify how to round
  mfloor <- function (x, base) {
    base*floor(x/base)
  }
  mround <- function (x, base) {
    base*round(x/base)
  }
  mceiling <- function (x, base) {
    base*ceiling(x/base)
  }

  par(oma = c(2, 3, 0, 2), mar = c(4, 6, 1, 1), mfrow = c(3, 3))

  low <- mfloor(min(df$LAM), .05)
  high <- mceiling(max(df$LAM), .05)
  mid <- mround((low + high)/2, .01)
  plot(df$gf, df$LAM, ylab = expression(paste('Litter area per mass (m'^'2', '/g)')),
       xlab = '', yaxt = 'n',
       ylim = c(0.99*low, (high + 0.15*(high-low))),
       cex.axis = 1.8, cex.lab = 2.2)
  points(df$gf, df$LAM)
  axis(side = 2, at = c(low, mid, high), cex.axis = 2,
       labels = sprintf("%.2f", c(low, mid, high)))
  legend('topleft', '(a)', bty = 'n', cex = 2)

  low <- mfloor(min(df$DMC), 5)
  high <- mceiling(max(df$DMC), 5)
  mid <- mround((low + high)/2, 1)
  plot(df$gf, df$DMC, ylab = 'Litter dry matter content (mg/g)',
       xlab = '', yaxt = 'n',
       ylim = c(0.99*low, (high + 0.15*(high-low))),
       cex.axis = 1.8, cex.lab = 2.2)
  points(df$gf, df$DMC)
  axis(side = 2, at = c(low, mid, high), cex.axis = 2,
       labels = sprintf("%.0f", c(low, mid, high)))
  legend('topleft', '(b)', bty = 'n', cex = 2)

  low <- mfloor(min(df$N), .05)
  high <- mceiling(max(df$N), .05)
  mid <- mround((low + high)/2, .01)
  plot(df$gf, df$N, ylab = 'Litter nitrogen content (wt%)',
       yaxt = 'n',
       ylim = c(0.99*low, (high + 0.15*(high-low))),
       cex.axis = 1.8, cex.lab = 2.2)
  points(df$gf, df$N)
  axis(side = 2, at = c(low, mid, high), cex.axis = 2,
       labels = sprintf("%.2f", c(low, mid, high)))
  legend('topleft', '(c)', bty = 'n', cex = 2)

  low <- mfloor(min(df$C), 1)
  high <- mceiling(max(df$C), 1)
  mid <- mround((low + high)/2, 1)
  plot(df$gf, df$C, ylab = 'Litter carbon content (wt%)',
       xlab = '', yaxt = 'n',
       ylim = c(0.99*low, (high + 0.15*(high-low))),
       cex.axis = 1.8, cex.lab = 2.2)
  points(df$gf, df$C)
  axis(side = 2, at = c(low, mid, high), cex.axis = 2,
       labels = sprintf("%.0f", c(low, mid, high)))
  legend('topleft', '(d)', bty = 'n', cex = 2)

  low <- mfloor(min(df$HC), .05)
  high <- mceiling(max(df$HC), .05)
  mid <- mround((low + high)/2, .01)
  plot(df$gf, df$HC, ylab = 'Litter hemicelluloses (wt%)',
       yaxt = 'n',
       ylim = c(0.99*low, (high + 0.15*(high-low))),
       cex.axis = 1.8, cex.lab = 2.2)
  points(df$gf, df$HC)
  axis(side = 2, at = c(low, mid, high), cex.axis = 2,
       labels = sprintf("%.2f", c(low, mid, high)))
  legend('topleft', '(e)', bty = 'n', cex = 2)

  low <- mfloor(min(df$CL), .05)
  high <- mceiling(max(df$CL), .05)
  mid <- mround((low + high)/2, .01)
  plot(df$gf, df$CL, ylab = 'Litter cellulose (wt%)',
       yaxt = 'n',
       ylim = c(0.99*low, (high + 0.15*(high-low))),
       cex.axis = 1.8, cex.lab = 2.2)
  points(df$gf, df$CL)
  axis(side = 2, at = c(low, mid, high), cex.axis = 2,
       labels = sprintf("%.2f", c(low, mid, high)))
  legend('topleft', '(f)', bty = 'n', cex = 2)

  low <- mfloor(min(df$LG), .05)
  high <- mceiling(max(df$LG), .05)
  mid <- mround((low + high)/2, .01)
  plot(df$gf, df$LG, ylab = 'Litter lignin (wt%)',
       yaxt = 'n',
       ylim = c(0.99*low, (high + 0.15*(high-low))),
       cex.axis = 1.8, cex.lab = 2.2)
  points(df$gf, df$LG)
  axis(side = 2, at = c(low, mid, high), cex.axis = 2,
       labels = sprintf("%.2f", c(low, mid, high)))
  legend('topleft', '(g)', bty = 'n', cex = 2)

}

## Custom legends

# Species name legend
legend_species <- function (spname) {

  legend(0.7, 1.1,
         xjust = 1,
         legend = spname,
         text.font = 3,
         cex = 3,
         bty = 'n')

}

# Three curves legend
legend_three_curves <- function () {

  legend(120, 0.008,
         xjust = 0,
         legend = c('data', 'total DTG', 'HC-1', 'HC-2', 'CL', 'LG'),
         ncol = 2,
         cex = 1.4,
         bty = 'n',
         lty = c(NA, 1, 6, 3, 4, 5),
         pch = c(20, NA, NA, NA, NA, NA),
         col = c('black', 'black', 'orange', 'red', 'green3', 'blue'),
         lwd = 2)

}

# Subfigure legend
legend_subfig_cex <- function (subfig, cex = 2.2) {

  legend('topleft', paste0('(', subfig, ')'), bty = 'n', cex = cex)

}

# Three curves horizontal legend
legend_three_curves_horizontal <- function () {

  plot(1, type = 'n', axes = FALSE, xlab = '', ylab = '')

  legend(x = "top", inset = 0,
         legend = c('data', 'total DTG', 'HC-1', 'HC-2', 'CL', 'LG'),
         horiz = TRUE,
         cex = 2.4,
         bty = 'n',
         lty = c(NA, 1, 6, 3, 4, 5),
         pch = c(20, NA, NA, NA, NA, NA),
         col = c('black', 'black', 'orange', 'red', 'green3', 'blue'),
         lwd = 2)
}

# legend for decay plots
legend_decay_curves_horizontal <- function (cex = 2.6) {

  plot(1, type = 'n', axes = FALSE, xlab = '', ylab = '')

  legend(x = "top", inset = 0,
         legend = c('data', 'Weibull', 'Negative exponential'),
         horiz = TRUE,
         cex = cex,
         bty = 'n',
         lty = c(NA, 1, 1),
         pch = c(20, NA, NA),
         col = c('black', 'red', 'black'),
         lwd = 3)
}

# traits_nocvnore_model <- function (df) {
#
#   mod1 <- c('N', 'C')
#   mod2 <- c('N', 'HC')
#   mod3 <- c('N', 'LG')
#   mod4 <- c('C')
#   mod5 <- c('N', 'DMC')
#
#   all_fe <- list(mod1, mod2, mod3, mod4, mod5)
#
#   jobs <- decaymod::create_jobs(model_type = 'w',
#                                 data = df,
#                                 random_effects = FALSE,
#                                 fixed_effects = all_fe)
#
#   mod1 <- decaymod::run_models(jobs[1],
#                                df,
#                                initial_mass = "mInit",
#                                removal_mass = "mRem",
#                                time = "t",
#                                group = "species_code",
#                                n_cores = 4,
#                                trait_param = 'beta',
#                                save_fit = TRUE)
#   mod2 <- decaymod::run_models(jobs[2],
#                                df,
#                                initial_mass = "mInit",
#                                removal_mass = "mRem",
#                                time = "t",
#                                group = "species_code",
#                                n_cores = 4,
#                                trait_param = 'beta',
#                                save_fit = TRUE)
#   mod3 <- decaymod::run_models(jobs[3],
#                                df,
#                                initial_mass = "mInit",
#                                removal_mass = "mRem",
#                                time = "t",
#                                group = "species_code",
#                                n_cores = 4,
#                                trait_param = 'beta',
#                                save_fit = TRUE)
#   mod4 <- decaymod::run_models(jobs[4],
#                                df,
#                                initial_mass = "mInit",
#                                removal_mass = "mRem",
#                                time = "t",
#                                group = "species_code",
#                                n_cores = 4,
#                                trait_param = 'beta',
#                                save_fit = TRUE)
#   mod5 <- decaymod::run_models(jobs[5],
#                                df,
#                                initial_mass = "mInit",
#                                removal_mass = "mRem",
#                                time = "t",
#                                group = "species_code",
#                                n_cores = 4,
#                                trait_param = 'beta',
#                                save_fit = TRUE)
#
#   names(mod1$fit)[1:4] <- c('alpha_intercept', 'beta_intercept', 'N', 'C')
#   names(mod2$fit)[1:4] <- c('alpha_intercept', 'beta_intercept', 'N', 'HC')
#   names(mod3$fit)[1:4] <- c('alpha_intercept', 'beta_intercept', 'N', 'LG')
#   names(mod4$fit)[1:3] <- c('alpha_intercept', 'beta_intercept', 'C')
#   names(mod5$fit)[1:4] <- c('alpha_intercept', 'beta_intercept', 'N', 'DMC')
#
#   par(mfrow = c(2, 3))
#   plot(mod1$fit, pars = c('alpha_intercept', 'beta_intercept', 'N', 'C'), cex = 2)
#   plot(mod2$fit, pars = c('alpha_intercept', 'beta_intercept', 'N', 'HC'), cex = 2)
#   plot(mod3$fit, pars = c('alpha_intercept', 'beta_intercept', 'N', 'LG'), cex = 2)
#   plot(mod4$fit, pars = c('alpha_intercept', 'beta_intercept', 'C'), cex = 2)
#   plot(mod5$fit, pars = c('alpha_intercept', 'beta_intercept', 'N', 'DMC'), cex = 2)
#
# }


### pre v real plot
# pred_real <- function (output_list) {
#
#   # extract appropriate values and plot
#   # but what if I want to plot that for every model?
#   # well that won't be a step in the remake. it can just be a step in my exploration.
#   # predictions v. real for all folds of model
#   pred <- dplyr::bind_rows(lapply(1:length(output_list), function(x) {
#     return(output_list[[x]]$pred)
#   }))
#
#   ##### now need to add real data
#
#   mean_pred <- pred %>%
#     dplyr::group_by(model, iter, data_point, mT_real) %>%
#     dplyr::summarise(mean = mean(draw)) %>%
#     as.data.frame()
#
#   for (i in unique(mean_pred$model)) {
#     i_pred <- pred[pred$model == i, ]
#     i_mean_pred <- mean_pred[mean_pred$model == i, ]
#     grDevices::png(paste0(predR2_path, 'model_', i, '.png'))
#     graphics::boxplot(draw ~ as.numeric(mT_real), i_pred,
#                       ylab = 'Posterior predicted distrbitions',
#                       xlab = 'Real test data')#,
#     #axes = FALSE)
#     fit <- stats::lm(i_mean_pred$mean ~ as.numeric(i_mean_pred$mT_real))
#     R2 <- paste('R2 is', format(summary(fit)$adj.r.squared, digits = 4))
#     int <- paste('Intercept is', format(stats::coef(fit)["(Intercept)"], digits = 4))
#     graphics::legend('topleft', bty = 'n', legend = c(R2,
#                                                       int))
#     grDevices::dev.off()
#   }
#
# }
