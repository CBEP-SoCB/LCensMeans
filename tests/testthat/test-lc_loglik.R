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

test_that("lc_loglik returns correct values", {
  conc <- c(rep(1,10), 1:20)
  ff <- c(rep(TRUE,10), rep(FALSE, 20))
  expect_lt(abs(lc_loglik(c(2,2),  conc, ff) - (-94.59033)), epsilon)
  expect_lt(abs(lc_loglik(c(0,2),  conc, ff) - (-94.27940)), epsilon)
  expect_lt(abs(lc_loglik(c(2,3),  conc, ff) - (-97.16272)), epsilon)
  expect_lt(abs(lc_loglik(c(2,10), conc, ff) - (-115.4876)), epsilon)
})

test_that("lc_loglik accepts missing values in Concentration", {
  the_dat <- mussels  %>% filter(PAH == "Benzo(b+k)fluoranthene")
  expect_error(lc_loglik(c(2,2),
                         the_dat$Concentration,
                         the_dat$Flag), NA)
  expect_equal(lc_loglik(c(2,2),
                         the_dat$Concentration,
                         the_dat$Flag),  -63.81571)
})

test_that("lc_loglik returns NA if no data", {
  conc <- rep(NA_real_,30)
  ff <- c(rep(TRUE,10), rep(FALSE, 20))
  expect_equal(lc_loglik(c(2,2), conc, ff), NA_real_)
})

test_that("lc_loglik produces error when parameters are inappropriate", {
  conc <- rep(10,20)
  ff <- c(rep(TRUE,10), rep(FALSE, 20))
  expect_error(lc_loglik(c(2,2), conc, ff))
  conc <- rep(10,30)
  expect_error(lc_loglik(c(2),   conc, ff))
  expect_error(lc_loglik(c(2,2,2),   conc, ff))
})
