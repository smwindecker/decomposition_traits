## Prepare decomposition data

prepare_decay_data <- function (initial_weight, removal_weight, trait_data) {

  # weights of the litter before decomposition
  leaves <- read.csv(initial_weight, header = T)

  # weights at removal of baggies (mid decomposition)
  weights <- read.csv(removal_weight, header = T)

  # merge the initial and removal weights
  decomp <- merge(leaves, weights)

  # add zero weight amounts
  bzero <- unique(decomp[,c('label_code', 'rep', 'species_code', 'treatment')])
  bzero$baggie <- 0
  bzero$pre_weight <- 6
  bzero$removal_date <- NA
  bzero$weeks <- 0
  bzero$days <- 0
  bzero$removal_weight <- 6
  bzero <- bzero[,c(1,2,5,3,4,6:10)]
  decomp_0 <- rbind(decomp, bzero)

  # remove rows where removal weight is NA - a +G bag was filled with foreign plant material
  decomp_1 <- decomp_0[!(is.na(decomp_0$removal_weight) | decomp_0$species_code == 'J'),]

  # remove average weight of baggie and tag and convert to mg
  decomp_1$mass_init <- (decomp_1$pre_weight - 2.13) * 1000
  decomp_1$mass_rem <- (decomp_1$removal_weight - 2.13) * 1000

  # remove datapoint where the baggie was calculated as having negative mass
  decomp_1 <- decomp_1[!decomp_1$mass_rem < 0 ,]

  # make treatment a factor
  decomp_1$fTreatment <- factor(decomp_1$treatment, levels = c('variable', 'inundated'))

  # necessary to keep the replicates organised
  decomp_1$fTub <- as.factor(paste(decomp_1$species_code, decomp_1$rep))

  # subset
  decomp_2 <- decomp_1[,c(4,9,11:14)]

  # tidy mislabelled species code
  decomp_2$species_code <- gsub('AA', 'J', decomp_2$species_code)

  # merge with trait data
  m_1 <- merge(decomp_2, trait_data, by = 'species_code')

  # add underscore between species' genus and species
  m_1$species <- gsub(' ', '_', m_1$species, fixed = TRUE)

  # take the log of the masses
  m_1$mRem <- log(m_1$mass_rem)
  m_1$mInit <- log(m_1$mass_init)

  # put time in years
  m_1$t <- m_1$days/365

  # scale and log traits
  trait_list <- c('SLA', 'DMC', 'N', 'C', 'HC', 'CL', 'LG')
  for (i in unique(trait_list)) {
    m_1[, i] <- scale(log(m_1[, i]))
  }

  # reorder by species
  m_2 <- m_1[order(m_1$species_code),]

  # subset inundated data set
  n <- m_2[m_2$fTreatment == 'inundated', ]
  n$species_code <- as.character(n$species_code)

  n
  # variable dataset
  #v_sp <- unique(m_1$species[m_2$fTreatment == 'variable'])
  #v <- m_1[m_1$species %in% v_sp, ]
  #write.table(m_2, 'munge/all_data.txt')

}








