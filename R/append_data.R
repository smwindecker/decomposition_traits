#' Append data subsets to job list
#'
#' @param input job list
#' @param data data you will subset
#' @param fold_id list item that corresponds to what to divide the test and training data by
#' @return list of jobs with subsetted data appended
#'
#' @export

append_data <- function(input, data) {

  list <- as.list(unlist(input[1,]))

  if (list$cv == 'CV') {

    fold <- list$cv_cluster

    # specify training data and append to list item called 'train'
    list[['train']] <- data[data$species_code != fold, ]

    # specify test data and append to list item called 'test'
    list[['test']] <- data[data$species_code == fold, ]

  }

  if (list$cv == 'noCV') {
    list[['train']] <- data
  }

  return(list)

}
