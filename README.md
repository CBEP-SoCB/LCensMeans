See the R Notebook "Conditional Means of Censored Distributions"
for an explanation of the following functions.

# Introduction
This tiny package exposes a couple of functions for estimating the
conditional mean of censored values for left censored environmental data.

It is still common practice in envrionemtnal analysis to replace
"non-detects" in laboratory data with essentially arbitrary values,
usually  zero, one half of the detection limit, or the detection limit. None of
these approaches rests on any statistical foundation, and they ignore the
information available  from knowing that a sample concentration was
below a certain value.  The frequency of non-detects itself contains information
relevant to estimating teh value of censored data.

We can do better.

The goal of this package is to offer a simple interface to replace non-detects
with a maximum likelihood-derived estimator of the mean of unobserved
(censored) observations.

It estimates those means by first calculating maximum likelihood estimates of
parameters of a log normal distribution that could have produced the observed
censored data, and then drawing a sample from that distribution, and calculating
the mean of values that would have been censored.

The package can handle different levels of censoring for different observations.

# Caveats
This package assumes data are distributed according to a log normal
distribution. While that is a good first guess for many environmetnal data sets,
envrionmetnal data are often substantially more right skewed than a lognormal
distribution implies.  The user is cautioned to evaluate whether making the
assumption of a log lormal distribution is reasonable for tehir data or not.


# Explanation
The idea is to fit a maximum likelihood model to the data
assuming we are looking at a censored lognormal distribution.

With a lognormal density in hand, we can estimate a conditional
mean of "unobserved" observations below a detection or reporting limit
by sampling from the truncated distribution and colculating a mean.

The first function ("left_censored_loglik") provides a likelihood function
for a left-censored lognormal distribution. It is used with maximum
likelihood estimation to estimate parameters for the underlying
lognormal distribution.  The parameters consist of a vector of
concentrations and a vector of flags indicating which
observations are censored, and therefore represent the applicable
detection limit.

The second function estimates the conditional mean for a truncated
lognormal distribution, given parameters of the lognormal distribution
and a cutoff value (here, the detection limit). The third function
vectorizes the second, allowing us to pass a vector of detection
limits, rather than calculating this for each sample.

The last function manages the process of combining the other functions to generate a "corrected" vector of estimated concentrations, where
censored values are replaced by conditional means.
