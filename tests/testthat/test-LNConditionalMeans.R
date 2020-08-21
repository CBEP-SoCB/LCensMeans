test_that("multiplication works", {
  expect_equal(2 * 2, 4)
})


test_that("left_censored_loglik returns expected values", {
  tt = c(rep(1,10), 1:20)
  ff = c(rep(TRUE,10), rep(FALSE, 20))
  epsilon <- 10^5
  expect_true(left_censored_loglik(c(2,2), tt,ff) - (-94.59033) < epsilon)
  expect_true(left_censored_loglik(c(0,2), tt,ff) - (-94.27940) < epsilon)
  expect_true(left_censored_loglik(c(2,3), tt,ff) - (-97.16272) < epsilon)
  expect_true(left_censored_loglik(c(2,10), tt,ff) - ( -115.4876) < epsilon)
})
