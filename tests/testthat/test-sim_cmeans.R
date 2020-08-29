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

# The current implementation uses Vectorize() to generate a vectorized
# version of .sim_single() , so is is complicated to add conditions.

test_that('length of output from sim_cmeans matches length of input', {
  the_dat <- mussels %>% filter(PAH == "Fluoranthene")
  parms <- find_parms(the_dat$Concentration, the_dat$Flag)
  res <-  sim_cmeans(parms[[1]], parms[[2]], the_dat$Concentration)
  expect_equal(length(res), length(the_dat$Concentration))
})

test_that('length of output from sim_cmeans with missing values matches length of input', {
  the_dat <- mussels  %>% filter(PAH == "Benzo(b+k)fluoranthene")
  parms <- find_parms(the_dat$Concentration, the_dat$Flag)
  res <-  sim_cmeans(parms[[1]], parms[[2]], the_dat$Concentration)
  expect_equal(length(res), length(the_dat$Concentration))
})

test_that('sim_cmeans returns values as expected',{
  # This is a little tricky because we need to see that multiple
  # values are close to what we have calculated previously.
  # This test will have to change for more efficient code.
  the_dat <- mussels %>% filter(PAH == "Fluoranthene")
  parms <- find_parms(the_dat$Concentration, the_dat$Flag)
  original <- c(49.120045, 49.548623, 47.000921, 48.453291, 24.672966,
                23.519998, 22.007736, 23.558695, 24.523340, 16.695698,
                5.891711, 8.098373, 24.438973, 24.854111, 24.061155,
                20.711412, 61.711595, 50.748058, 59.135875, 43.065625,
                6.307027, 5.854042, 6.191214, 5.743561, 5.789952,
                9.925236, 5.455754, 9.186314, 10.282295, 9.192600)
  res <-  sim_cmeans(parms[[1]], parms[[2]], the_dat$Concentration)
  expect_equal(original, res, tolerance = 1)
})


test_that('Results of sim_cmeans are less than or equal to cutoff',{
  the_dat <- mussels %>% filter(PAH == "Fluoranthene")
  parms <- find_parms(the_dat$Concentration, the_dat$Flag)
  res <-  sim_cmeans(parms[[1]], parms[[2]], the_dat$Concentration)
  expect_true(all(res <= the_dat$Concentration))
})

test_that('sim_cmeans returns NAs when fed NAs',{
  conc <- 1:20
  conc[c(5,10,15)] <- NA_real_
  ff   <- c(rep(TRUE,3), rep(FALSE,17))
  parms <- find_parms(conc, ff)
  res <-  sim_cmeans(parms[[1]], parms[[2]], conc)
  expect_equal(sum(is.na(res)), 3)
  expect_equal(sum(is.na(res[c(5,10,15)])), 3)
})

test_that('sim_cmeans returns NAs when fed some negative concentrations',{
  conc <- 1:20
  conc[c(5,10,15)] <- - conc[c(5,10,15)]
  ff   <- c(rep(TRUE,3), rep(FALSE,17))
  parms <- find_parms(conc, ff)
  res <-  sim_cmeans(parms[[1]], parms[[2]], conc)
  expect_equal(sum(is.na(res)), 3)
  expect_equal(sum(is.na(res[c(5,10,15)])), 3)
})

test_that('sim_cmeans throws an error when expected',{})

