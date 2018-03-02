context('create_design_matrix')

test_that('formula list is correct', {
  trait_list <- c('N', 'C', 'SLA')
  formula_output <- create_design_matrix(trait_list, extra_null = c('ne', 're'))
  expect_match(unique(formula_output$formula_alpha), c('~ 1',
                                 '~ 1 + N',
                                 '~ 1 + C',
                                 '~ 1 + SLA'))


  n_params <- 1 + length(trait_list)
  n_rows <- n_params*n_params + 2
  expect_equal(dim(formula_output), c(n_rows, 4))

})
