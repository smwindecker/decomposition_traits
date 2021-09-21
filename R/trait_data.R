## Load data

# Prepare C and N trait data
load_leco_traits <- function (leco_data) {

  # read leco data
  leco <- read.csv(leco_data, header = F, skip = 12)

  # exclude malformed sample and unused species. select columns for C and N
  leco <- leco[!(leco[, 1] == 'donotuse' | leco[, 1] == 'AA'), c(1,6,7)]
  colnames(leco) <- c('species_code', 'C', 'N')
  leco

}

# Prepare LES traits
load_les_traits <- function (trait_data, species_data) {

  # read raw trait data
  traits <- read.csv(trait_data, header = T)

  # merge with the species data file
  trait <- merge(traits[!traits$species_code == 'AA',],
                 species_data[!species_data$species_code == 'AA',
                              c('species_code', 'sp_abrev', 'species', 'family', 'plant_part', 'gf', 'gf_old')],
                 by = c('species_code', 'plant_part'))

  # calculate LAM (m2/g) and DMC (mg/g) for each of ten samples
  trait$longLAM <- (trait$area/100)/trait$dry_weight
  trait$longDMC <- (trait$dry_weight*1000)/trait$wet_weight

  # calculate mean LAM and DMC of ten samples
  trait_1 <- trait %>%
    dplyr::group_by(species_code) %>%
    dplyr::summarize(
      LAM = mean(longLAM),
      DMC = mean(longDMC)
    )

  # merge with species data
  trt <- merge(trait_1, unique(trait[,c('species_code', 'sp_abrev', 'species', 'family', 'gf', 'gf_old')]))
  trt

}

# Prepare species data
load_species <- function (species_data) {

  species <- read.csv(species_data, header = T) %>%
    dplyr::filter(species_code != 'AA')

  species

}

# Wrapper for all species for TGA functions
tga_wrapper <- function (species_data, function_name, ...) {

  input <- as.character(species_data$species_code)
  sapply(input, FUN = function_name, ..., simplify = FALSE, USE.NAMES = TRUE)

}

# Function to deconvolve species' TGA data
tga_deconvolve <- function (species_code, data_folder) {

  # extract species code
  x <- species_code

  three <- c('KK', 'M', 'X')

  if (x %in% three) {
    n_peaks <- 3
  } else if (x == 'CC') {
    n_peaks <- 4
  } else {
    n_peaks = NULL
  }

  # read raw TGA
  file <- paste0(data_folder, x, '_TGA.csv')
  tmp <- process_raw_tga(file)

  # deconvolve TGA data
  output <- mixchar::deconvolve(tmp, upper_temp = 650, n_peaks = n_peaks)

  # extract weight estimates
  weights <- output$weights

  # if 3-curves, add HC_1 and change name of HC to HC_2
  if (ncol(weights) == 4) {
    weights$HC_1 <- NA
    colnames(weights)[1] <- 'HC_2'
  }

  # set species of weight outputs
  weights$species_code <- x
  weights <- weights[, c('species_code', 'HC_1', 'HC_2', 'CL', 'LG', 'value_type')]

  # set species of parameter outputs and make data wide
  params <- as.data.frame(summary(output$model_fit)$coefficients[,1])
  params$species_code <- x
  params$parameter <- row.names(params)
  colnames(params) <- c('coefficient', 'species_code', 'parameter')
  params <- tidyr::spread(params, parameter, coefficient)

  return(x = list(data = output$data,
                  bounds = output$bounds,
                  minpack.lm = output$model_fit,
                  n_peaks = output$n_peaks,
                  weights = weights,
                  params = params))

}

# Prepare weights
load_tga_traits <- function (species_deconvolved_list) {

  # combine weight estimates from all deconvolved species data
  weight_estimates <- dplyr::bind_rows(lapply(1:length(species_deconvolved_list), function(x) {
    return(species_deconvolved_list[[x]]$weights)
  }))

  weight_estimates

}

# Load and process raw TGA data
process_raw_tga <- function (raw_file) {

  df <- read.csv(raw_file, header = FALSE, skip = 29)
  names(df) <- c('temp', 'time', 'mass_loss')
  init_mass <- read.csv(raw_file, nrows = 1, header = FALSE, skip = 17)[1,2]
  tmp <- mixchar::process(df,
                          init_mass = init_mass,
                          temp = 'temp',
                          mass_loss = 'mass_loss')

  tmp
}

## Analyse data
# Merge trait datasets
traits_combine <- function (species_data, leco_data, trait_data, tga_data_folder, ranseed) {

    # here first load the appropriate ones from load_data.R
    species_df <- species_data
    carbon_nitrogen <- load_leco_traits(leco_data = leco_data)
    economic_traits <- load_les_traits(trait_data = trait_data,
                                       species_data = species_df)
    tga_output <- tga_wrapper(species_data = species_df,
                              function_name = tga_deconvolve,
                              data_folder = tga_data_folder)
    biomass_traits <- load_tga_traits(tga_output)

    t <- dplyr::full_join(economic_traits, carbon_nitrogen, by = 'species_code') %>%
      dplyr::full_join(., biomass_traits, by = 'species_code') %>%
      .[, c(1,4:8,2,3,10,9,11:15)]

  t
}

# Select mean of traits
traits_mean_only <- function (traits) {

  t_mean <- traits[traits$value_type == 'mean',]
  t_mean$HC <- t_mean$HC_1 + t_mean$HC_2
  t_mean$HC[is.na(t_mean$HC)] <- t_mean$HC_2[is.na(t_mean$HC)]

  t_mean
}

# Produce matrix of logged traits
traits_log <- function (t_mean) {

  cov <- t_mean %>%
    magrittr::set_rownames(t_mean[, 'sp_abrev']) %>%
    dplyr::select(LAM, DMC, N, C, HC, CL, LG)
  cov[] <- log(cov[])

  cov
}

# Conduct PCA
pca_data <- function (df) {

  # calculate loadings
  prin <- stats::princomp((stats::na.omit(df)), cor = TRUE, scores = TRUE)

  prin
}

