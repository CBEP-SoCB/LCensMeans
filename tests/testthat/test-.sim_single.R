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


context('.sim_single')
test_that('.sim_single returns NA for cutoff <= 0', {
  expect_equal(.sim_single(lmu=0,
                           lsigma=1,
                           cutoff = 0), NA_real_)
  expect_equal(.sim_single(lmu=0,
                           lsigma=1,
                           cutoff = -1), NA_real_)
})

test_that('.sim_single returns correct values',{
  expect_equal(.sim_single(lmu=0,
                           lsigma=0,
                           cutoff = 2), 1)
  expect_lt(.sim_single(lmu=0,
                        lsigma=1,
                        cutoff = 1), 1)
  expect_lt(.sim_single(lmu=0,
                        lsigma=1,
                        cutoff = 0.1), 0.1)
})

test_that('.sim_single returns NA when it encounters difficulties', {
  expect_equal(.sim_single(lmu=0, lsigma=1, cutoff = 0.01), NA_real_)
})

