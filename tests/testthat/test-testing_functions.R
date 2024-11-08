# Test if check_dataset function works
test_that("column names are correct and body temp values are biologically sound", {
  df <- tibble(id = c(1,1,1), 
               body_tempx = c(10, 20, 25), 
               torpor_onset = c(30, 30, 30))
  expect_error(check_dataset(df), "column body_temp is missing")
})


# Test if add_temp_lag_cols function works
test_that("adding the three columns body_temp_diff, body_temp_diff_2 and below_threshold works", {
  df <- tibble(id = c(1,1,1), 
               body_temp = c(10, 20, 25), 
               torpor_onset = c(30, 30, 30))
  expect_equal(add_temp_lag_cols(df), tibble(id = c(1,1,1), 
                                            body_temp = c(10, 20, 25), 
                                            torpor_onset = c(30, 30, 30),
                                            body_temp_diff = c(NA, 10, 5),
                                            body_temp_diff_2 = c(NA, NA, 15),
                                            below_threshold = c(TRUE, TRUE, TRUE)) |> group_by(id) )
})


# Test if myrleid function works with numerical values
test_that("myrleid function working with numerical", {
  x = c(1,1,2,2,2)
  expect_equal(myrleid(x), c(1,1,2,2,2))
})


# Test if myrleid function works with TRUE/FALSE
test_that("myrleid function working with TRUE/FALSE", {
  x = c(TRUE,FALSE,TRUE,TRUE)
  expect_equal(myrleid(x), c(1,2,3,3))
})


# Test if add_section_numbers function works, adding 
test_that("adding the three columns body_temp_diff, body_temp_diff_2 and below_threshold works", {
  df <- tibble(below_threshold = c(TRUE,FALSE,TRUE,TRUE))
  expect_equal(add_section_numbers(df), tibble(below_threshold = c(TRUE,FALSE,TRUE,TRUE), 
                                               num = c(1,2,3,3),
                                               count_length = c(1,1,2,2),
                                               torpor = c("torpor", "not torpor", "torpor", "torpor")))
})

