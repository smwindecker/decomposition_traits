## Figures

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

  low <- mfloor(min(df$SLA), .05)
  high <- mceiling(max(df$SLA), .05)
  mid <- mround((low + high)/2, .01)
  plot(df$gf, df$SLA, ylab = expression(paste('Specific litter area (m'^'2', '/g)')),
       xlab = '', yaxt = 'n',
       ylim = c(0.99*low, (high + 0.15*(high-low))),
       cex.axis = 1.8, cex.lab = 2.2)
  points(df$gf, df$SLA)
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

# Produce parameter simulation of Fraser-Suzuki function
simulate_weibull_4 <- function () {

  # simulated model results with 'real' mean parameter values
  x <- seq(0, 0.7, length.out = 1000)

  a.low <- 0.3
  a.high <- 1
  b.low <- 0.3
  b.high <- 2.5

  y.ah.bl <- weibull(4000, x, a.high, b.low)
  y.ah.bh <- weibull(4000, x, a.high, b.high)
  y.al.bl <- weibull(4000, x, a.low, b.low)
  y.al.bh <- weibull(4000, x, a.low, b.high)

  # simulated data
  m0.range <- rnorm(90, 4000, 20)
  x.range <- rep(seq(0, 0.7, length.out = 30), each = 3)

  a.low.range <- rnorm(90, 0.3, 0.05)
  a.high.range <- rnorm(90, 1, 0.1)
  b.low.range <- rnorm(90, 0.3, 0.05)
  b.high.range <- rnorm(90, 2.2, 0.3)

  pred.ah.bl <- weibull(m0.range, x.range, a.high.range, b.low.range)
  pred.ah.bh <- weibull(m0.range, x.range, a.high.range, b.high.range)
  pred.al.bl <- weibull(m0.range, x.range, a.low.range, b.low.range)
  pred.al.bh <- weibull(m0.range, x.range, a.low.range, b.high.range)

  layout(matrix(c(1,2,3,4,5,5), nrow = 3, ncol = 2, byrow = TRUE), heights = c(0.8, 0.8, 0.2))
  par(oma = c(5, 5, 0, 2), mar = c(1, 3, 2, 0))

  # first plot is high alpha, low beta
  plot(x.range, pred.al.bh, ylim = c(0, 4200), xaxt = 'n', yaxt = 'n', xlab = '', ylab = '', pch = 20)
  lines(x, y.al.bh, lty = 1)
  legend('topleft', legend = '(a)', bty = 'n', cex = 1.6)
  legend('topright',
         legend = c(expression(paste(alpha, ' = 0.3')),
                    expression(paste(beta, ' = 1.7'))),
         bty = 'n',
         cex = 1.5)
  axis(side = 2, at = c(0, 1000, 2000, 3000, 4000), cex.axis = 1.2,
       labels = c(0, 1000, 2000, 3000, 4000))

  # second plot is high alpha, high beta
  plot(x.range, pred.ah.bh, ylim = c(0, 4200), xaxt = 'n', yaxt = 'n', xlab = '', ylab = '', pch = 20)
  lines(x, y.ah.bh, lty = 1)
  legend('topleft', legend = '(b)', bty = 'n', cex = 1.6)
  legend('topright',
         legend = c(expression(paste(alpha, ' = 1')),
                    expression(paste(beta, ' = 1.7'))),
         bty = 'n',
         cex = 1.5)

  # third plot is low alpha, low beta
  plot(x.range, pred.al.bl, ylim = c(0, 4200), xaxt = 'n', yaxt = 'n', xlab = '', ylab = '', pch = 20)
  lines(x, y.al.bl, lty = 1)
  legend('topleft', legend = '(c)', bty = 'n', cex = 1.6)
  legend('topright',
         legend = c(expression(paste(alpha, ' = 0.3')),
                    expression(paste(beta, ' = 0.3'))),
         bty = 'n',
         cex = 1.5)
  axis(side = 2, at = c(0, 1000, 2000, 3000, 4000), cex.axis = 1.2,
       labels = c(0, 1000, 2000, 3000, 4000))
  axis(side = 1, at = c(0, 0.25, 0.5, 0.75, 1), cex.axis = 1.2,
       labels = c(0, 0.25, 0.5, 0.75, 1))

  # fourth plot is low alpha, high beta
  plot(x.range, pred.ah.bl, ylim = c(0, 4200), xaxt = 'n', yaxt = 'n', xlab = '', ylab = '', pch = 20)
  lines(x, y.ah.bl, lty = 1)
  legend('topleft', legend = '(d)', bty = 'n', cex = 1.6)
  legend('topright',
         legend = c(expression(paste(alpha, ' = 1')),
                    expression(paste(beta, ' = 0.3'))),
         bty = 'n',
         cex = 1.5)
  axis(side = 1, at = c(0, 0.25, 0.5, 0.75, 1), cex.axis = 1.2,
       labels = c(0, 0.25, 0.5, 0.75, 1))

  mtext(text = 'Time (years)',
        side = 1,
        line = 2.1,
        outer = TRUE,
        cex = 1.8)
  mtext(text = 'Mass Remaining (mg)',
        side = 2,
        line = 0.9,
        outer = TRUE,
        cex = 1.8)

  plot(1, type = 'n', axes = FALSE, xlab = '', ylab = '')

  legend(x = "top", inset = 0,
         legend = c('Simulated data', 'Weibull model'),
         horiz = TRUE,
         cex = 2.4,
         bty = 'n',
         lty = c(NA, 1),
         pch = c(20, NA),
         lwd = 2)

}

# Weibull function
weibull <- function (m0, t, a, b) {
  m0 * exp(-(t/b)^a)
}

# what about one plot. neg exp medium level rate. then show four other lines -
# alpha low beta high, alpha high beta high, alpha low beta low, alpha high beta low to demo? maybe.
#
# or alpha = 1 with high and low beta. then low alapha with both. on same plot.
# can label the other neg exp.
#
# we predict that alpha is lower than 1 for our species and therefore that weibull
# will be better suited. and that . yessss


# 2 plot : Produce parameter simulation of Fraser-Suzuki function
simulate_weibull <- function (constant_param, constant_value, low, high) {

  # simulated model results with 'real' mean parameter values
  x <- seq(0, 0.7, length.out = 1000)

  if (constant_param == 'beta'){
    y.low <- weibull(4000, x, low, constant_value)
    y.high <- weibull(4000, x, high, constant_value)
    subfig <- 'a'
    text_high <- bquote(alpha ~ '=' ~ .(high))
    text_low <- bquote(alpha ~ '=' ~ .(low))
  }

  if (constant_param == 'alpha') {
    y.low <- weibull(4000, x, constant_value, low)
    y.high <- weibull(4000, x, constant_value, high)
    subfig <- 'b'
    text_high <- bquote(beta ~ '=' ~ .(high))
    text_low <- bquote(beta ~ '=' ~ .(low))
  }

  par(oma = c(2, 2, 0, 0), mar = c(2, 2, 2, 2))
  plot(x, y.high, type = 'l', lty = 2, ylim = c(0, 4350),
       xlim = c(0, 0.9), xaxt = 'n', yaxt = 'n', xlab = '', ylab = '')
  lines(x, y.low, lty = 4)
  legend('topleft', legend = paste0('(', subfig, ')'), bty = 'n', cex = 1.6)
  text(x = max(x)+.12,
       y = y.high[1000],
       labels = text_high,
       cex = 1.6)
  text(x = max(x)+.12,
       y = y.low[1000],
       labels = text_low,
       cex = 1.6)

  mtext(text = 'Mass Remaining (mg)',
        side = 2,
        line = 0,
        outer = TRUE,
        cex = 1.6)

  mtext(text = 'Time (years)',
        side = 1,
        line = 0,
        outer = TRUE,
        cex = 1.6)
}


# 2 plot : Produce parameter simulation of Fraser-Suzuki function
simulate_weibull_new <- function (alpha.low, alpha.high, beta.low, beta.high) {

  internal.plot <- function (a, subfig, beta.low, beta.high) {

    x <- seq(0, 0.7, length.out = 1000)

    y.high <- weibull(4000, x, 1, beta.high)
    y.high.label <- bquote(beta ~ '=' ~ .(beta.high))
    y.low <- weibull(4000, x, 1, beta.low)
    y.low.label <- bquote(beta ~ '=' ~ .(beta.low))

    a1 <- weibull(4000, x, a, beta.high)
    a1.label <- bquote(alpha ~ '=' ~ .(a))
    a2 <- weibull(4000, x, a, beta.low)
    a2.label <- bquote(alpha ~ '=' ~ .(a))

    plot(x, y.high, type = 'l', lty = 1, ylim = c(0, 4350),
         xlim = c(0, 0.9), xaxt = 'n', yaxt = 'n', xlab = '', ylab = '')
    lines(x, y.low, lty = 1)
    lines(x, a1, lty = 2)
    lines(x, a2, lty = 4)
    legend('topleft', legend = paste0('(', subfig, ')'), bty = 'n', cex = 1.6)
    text(x = max(x)+.12,
         y = y.high[1000],
         labels = y.high.label,
         cex = 1.6)
    text(x = max(x)+.12,
         y = y.low[1000],
         labels = y.low.label,
         cex = 1.6)
    text(x = max(x)+.12,
         y = a1[1000],
         labels = a1.label,
         cex = 1.6)
    text(x = max(x)+.12,
         y = a1[1000]-100,
         labels = y.high.label,
         cex = 1.6)
    text(x = max(x)+.12,
         y = a2[1000],
         labels = a2.label,
         cex = 1.6)
    text(x = max(x)+.12,
         y = a2[1000]-100,
         labels = y.low.label,
         cex = 1.6)
  }

  par(oma = c(2, 2, 0, 0), mar = c(2, 2, 2, 2), mfrow = c(1,2))
  internal.plot(alpha.low, 'a', beta.low, beta.high)
  internal.plot(alpha.high, 'b', beta.low, beta.high)

  mtext(text = 'Mass Remaining (mg)',
        side = 2,
        line = 0,
        outer = TRUE,
        cex = 1.6)

  mtext(text = 'Time (years)',
        side = 1,
        line = 0,
        outer = TRUE,
        cex = 1.6)
}

# Produce pca plot
pca <- function (prin, df, species_data) {

  # first two axes' scores
  pc12 <- prin$scores[, 1:2]
  df_pc12 <- data.frame(pc12)
  df_pc12$sp_abrev <- rownames(df_pc12)

  # label with abreviations
  pc12_labeled <- merge(df_pc12, species_data[,c('sp_abrev', 'sp_pca_label', 'gf')])
  rownames(pc12_labeled) <- pc12_labeled[,'sp_pca_label']

  # get length of axes
  fit <- vegan::envfit(pc12, na.omit(df))

  vars <- prin$sdev^2
  prop_vars <- vars/sum(vars)

  par(oma = c(2, 2, 0, 2))
  plot(pc12_labeled[, c('Comp.1', 'Comp.2')], ylab = '', xlab = '', xaxt = 'n', yaxt = 'n',
       ylim = c(-4, 4), xlim = c(-4.2, 5),
       cex.axis = 1, cex = 1.8, pch = c(2, 20, 3, 8, 0)[as.numeric(pc12_labeled$gf)])
  plot(fit, cex = 2, col = 1, labels = list(vectors = c('SLA', 'DMC', 'N', 'C', 'HC', 'CL', 'LG')))

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
  legend('topright',
         c('Forb', 'Graminoid', 'Non-vascular', 'Shrub', 'Tree'),
         bty = 'n',
         pch = c(2, 20, 3, 8, 0),
         cex = 1.8)

  reg_labels <- pc12_labeled[!rownames(pc12_labeled) %in% c('Sph',
                                                            'P.pro',
                                                            'A.den',
                                                            'M.cri',
                                                            'L.aus',
                                                            'J.ama',
                                                            'P.dis'), ]

  text(x = reg_labels[, 'Comp.1'], y = reg_labels[, 'Comp.2'],
       labels = row.names(reg_labels), vfont = c('sans serif', 'italic'),
       cex = 1.5, pos = 4)

  up_labels <- pc12_labeled[c('A.den', 'M.cri', 'L.aus', 'J.ama', 'P.dis'), ]

  text(x = up_labels[, 'Comp.1'], y = up_labels[, 'Comp.2']+0.16,
       labels = row.names(up_labels), vfont = c('sans serif', 'italic'),
       cex = 1.5, pos = 4)

  down_labels <- pc12_labeled[c('Sph', 'P.pro'), ]

  text(x = down_labels[, 'Comp.1'], y = down_labels[, 'Comp.2']-0.16,
       labels = row.names(down_labels), vfont = c('sans serif', 'italic'),
       cex = 1.5, pos = 4)

}

## Custom legends

# Species name legend
legend_species <- function (spname) {

  legend(650, .01,
         xjust = 1,
         legend = spname,
         text.font = 3,
         cex = 3.1,
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
legend_subfig <- function (subfig, cex = 2.2) {

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
