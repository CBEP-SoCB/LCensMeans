# LCensMeans

This tiny package exposes a couple of functions for estimating the
conditional mean of censored values for left censored environmental data.

# Installation
## Install `devtools`
If you have not already installed the `devtools`package, 
you will need to do so.

On the R Command Line, enter

`install.packages("devtools")`

Alternatively, you can use the GUI menus in RStudio:

`Tools -> Install Packages... `

Then select the `devtools` package from  the enormous list 
of packages available on CRAN, and the package will be 
installed.

## Install Package

To install a package of interest from GitHub, you need 
to specify both the "Author" and the "Package".  If you 
have a URL to the GitHub Repo, the form of that URL is 
`https://github.com/<Author>/<Package>`, so you have all
the information you need.

```	
library(devtools)
install_github("CBEP-SoCB/LCensMeans")
```

Or, alternatively, if you want to avoid polluting your 
search path with an unnecessary environment associated
with `devtools`, you can accomplish the same thing with
the following.

```
devtools::install_github("CBEP-SoCB/LCensMeans")
```

# Package Contents
It is still common practice in environmental analysis to replace "non-detects"
in laboratory data with essentially arbitrary values, usually  zero, one half of
the detection limit, or the detection limit. None of these approaches rests on
any statistical foundation, and they ignore the information available  from
knowing that a sample concentration was below a certain value.  The frequency of
non-detects itself contains information relevant to estimating teh value of
censored data.

We can do better.

The goal of this package is to offer a simple interface to replace non-detects
with a maximum likelihood-derived estimator of the mean of unobserved
(censored) observations.

It estimates those means by first calculating maximum likelihood estimates of
parameters of a log normal distribution that could have produced the observed
censored data, then drawing a sample from that distribution, and calculating
the mean of values that would have been censored.

The package can handle different levels of censoring for different observations.

# Caveats
This package assumes data are distributed according to a lognormal
distribution. While that is a good first guess for many environmental data sets,
environmental data are often substantially more right skewed than a lognormal
distribution implies.  The user is cautioned to evaluate whether making the
assumption of a lognormal distribution is reasonable for their data or not.


# Explanation
The idea is to fit a maximum likelihood model to the data
assuming we are looking at a censored lognormal distribution.

With a lognormal density in hand, we can estimate a conditional
mean of "unobserved" observations below detection or reporting limits
by sampling from the truncated distribution and calculating a mean.

The first function (`lc_loglik()`) provides a likelihood function for a
left-censored lognormal distribution. It is used with maximum likelihood
estimation in the function `find_parms()` to estimate parameters for the
underlying lognormal distribution.  

The next function (`sim_cmeans()`) estimates conditional means of values from a
specified lognormal distribution (with specific parameters) that fall below
(variable) cutoff values. Estimated means of censored values are determined via 
simulation.

The last function (`sub_cmeans()`) manages the process of combining the other
functions to generate a "corrected" vector of estimated concentrations, where
censored values are replaced by the estimated conditional means of censored 
values.

See the R Notebook ["Conditional Means of Censored Distributions"](https://github.com/ccb60/PortlandHarborToxics/blob/master/Analysis/Conditional_means_of_censored_distributions.Rmd)
for additional consideration of these functions.

# Data Format
The Functions in this package expect data to be presented in a two-column
format.  The first column consists of a vector of concentrations. The
Second contains a vector of flags that indicate whether the observations was
censored or not.  Where the censoring flag is `TRUE`, the value in the 
concentration vector represents the detection or reporting limit.  Where the
flag is `FALSE`, the value in the concentration vector is an observed 
concentration.
