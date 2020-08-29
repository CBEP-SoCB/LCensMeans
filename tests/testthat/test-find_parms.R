library(dplyr)

# Epsilon provides a factor that defines if two numeric values should be
# considered equivalent.  THis is for deviations numerical representations.
# Some functions are stochastic, and so do not provide identical results every
# time.  Results vary from slightly from run to run. It is also useful
# if comparing to  prior value that is represented less accurately in printed
# output than internally.


# Could shift to using the "tolerance" argument of expect_equal (passed through
# to compare). instead of using expect_lt and epsilon.
epsilon <- 10^5


context('find_parms')
test_that("find_parms provides correct output", {
  the_dat <- mussels  %>% filter(PAH == "Fluoranthene")
  res <- find_parms(the_dat$Concentration,
                    the_dat$Flag)
  expect_lt(abs(res[[1]] - 3.699337 ), epsilon)
  expect_lt(abs(res[[2]] - 1.311283 ), epsilon)
})

test_that("find_parms raises an error if parms are different length", {
  conc <- rep(10,20)
  ff   <- c(rep(TRUE,10), rep(FALSE, 20))
  expect_error(find_parms(conc, ff))
})

test_that("find_parms raises an error if no data", {
  conc <- rep(NA_real_,20)
  ff   <- c(rep(TRUE,7), rep(FALSE, 17))
  expect_error(find_parms(conc, ff))

  conc <- rep(1,20)
  ff   <- rep(NA,20)
  expect_error(find_parms(conc, ff))

  conc <- c(rep(NA_real_, 5), rep(1,15))
  ff <-  c(rep(TRUE, 3), rep(FALSE,2), rep(NA, 15))
  expect_error(find_parms(conc, ff))
})
