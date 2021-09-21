get_deviance <- function (path, fun = 'mean') {

  # trait models
  alpha.traits.cv.nore200 <- readRDS(paste0(path, 'model_output/alpha.traits.cv.nore200.RDS'))
  alpha.traits.cv.nore400 <- readRDS(paste0(path, 'model_output/alpha.traits.cv.nore400.RDS'))
  alpha.traits.cv.nore600 <- readRDS(paste0(path, 'model_output/alpha.traits.cv.nore600.RDS'))
  alpha.traits.cv.nore812 <- readRDS(paste0(path, 'model_output/alpha.traits.cv.nore812.RDS'))
  beta.traits.cv.nore400 <- readRDS(paste0(path, 'model_output/beta.traits.cv.nore400.RDS'))
  beta.traits.cv.nore812 <- readRDS(paste0(path, 'model_output/beta.traits.cv.nore812.RDS'))
  both.traits.cv.nore400 <- readRDS(paste0(path, 'model_output/both.traits.cv.nore400.RDS'))
  both.traits.cv.nore812 <- readRDS(paste0(path, 'model_output/both.traits.cv.nore812.RDS'))

  alpha.cv <- list(alpha.traits.cv.nore200,
                   alpha.traits.cv.nore400,
                   alpha.traits.cv.nore600,
                   alpha.traits.cv.nore812)
  alpha.deviance <- deviance_table(alpha.cv, fun, 'alpha')

  beta.cv <- list(beta.traits.cv.nore400,
                  beta.traits.cv.nore812)
  beta.deviance <- deviance_table(beta.cv, fun, 'beta')

  both.cv <- list(both.traits.cv.nore400,
                  both.traits.cv.nore812)
  both.deviance <- deviance_table(both.cv, fun)

  out <- rbind(alpha.deviance, beta.deviance, both.deviance)
  return(out)
}

