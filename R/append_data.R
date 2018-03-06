#' Append data subsets to job list
#'
#' @param input job list item
#' @param data dataframe containing data for models
#' @param cv_cluster optional columns name to specify cv clusters
#' @return list of jobs with subsetted data appended
#'
#' @export

append_data <- function(input, data, cv_cluster) {

  list <- as.list(unlist(input[1,]))

  if (list$cv == 'CV') {

    fold <- list$cv_cluster

    # specify training data and append to list item called 'train'
    list[['train']] <- data[data[, cv_cluster] != fold, ]

    # specify test data and append to list item called 'test'
    list[['test']] <- data[data[, cv_cluster] == fold, ]

  }

  if (list$cv == 'noCV') {
    list[['train']] <- data
  }

  return(list)

}
