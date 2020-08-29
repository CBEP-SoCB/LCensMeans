library(dplyr)

test_that('length of output from sub_cmeans matches length of input', {
  the_dat <- mussels %>% filter(PAH == "Fluoranthene")
  res <-  sub_cmeans(the_dat$Concentration, the_dat$Flag)
  expect_equal(length(res), length(the_dat$Concentration))
})

test_that('length of output from sub_cmeans with missing values matches length of input', {
  the_dat <- mussels  %>% filter(PAH == "Benzo(b+k)fluoranthene")
  res <-  sub_cmeans(the_dat$Concentration, the_dat$Flag)
  expect_equal(length(res), length(the_dat$Concentration))
})

test_that('sub_cmeans returns values as expected',{
  # This is a little tricky because we need to see that multiple
  # values are close to what we have calculated previously.
  # This test will have to change for more efficient code.
  the_dat <- mussels %>% filter(PAH == "Fluoranthene")
  original = c(200, 193, 178, 189, 59,
               58, 57, 58, 61, 38,
               5.74, 15, 62, 65, 57,
               48, 345, 228, 294, 148,
               11, 5.82, 11, 5.91, 5.67,
               19, 9.5, 17, 19.4, 16.7)
  res <-  sub_cmeans(the_dat$Concentration, the_dat$Flag)
  expect_equal(original, res, tolerance = 1)
})

test_that('Results of sub_cmeans are less than or equal to cutoff',{
  the_dat <- mussels %>% filter(PAH == "Fluoranthene")
  res <-  sub_cmeans(the_dat$Concentration, the_dat$Flag)
  expect_true(all(res <= the_dat$Concentration))
})
