#' Calculate the log likelihood of left censored lognormal data.
#'
#' Calculate the log likelihood of left censored lognormal data.  The
#' log likelihood has two parts -- a term for observed values over the
#' detection limit, and term for censored values (where all we know is that
#' the value falls below the detection limit).
#'
#' @param params A two item list or vector containing the mean and sd of the
#'   distribution on the log scale (i.e., the parameters of the related
#'   normal distribution).
#'
#' @param cc A vector of data values, including both observed values, where they
#'   exist, or the detection limits, where data was censored.
#'
#' @param flg A vector of TRUE or FALSE values, of the same length as cc, that
#'    indicates which values are detection limits (TRUE) and which are measured
#'    values (FALSE)
#'
#' @return The calculated log likelihood under a censored lognormal distribution with
#'    the specific parameters.
#'
#' @examples
#' df <- data.frame(sim = sort(stats::rlnorm(25,2,3)),
#'                   cens=c(rep.int(TRUE,5), rep.int(FALSE,20)))
#' df$sim[1:4] <- df$sim[5]
#' left_censored_loglik(c(2,5), df$sim, df$cens)
#' left_censored_loglik(c(2,3), df$sim, df$cens)
#' left_censored_loglik(c(1,3), df$sim, df$cens)
#' @export
left_censored_loglik <- function(params, cc, flg) {
    stopifnot(length(params)==2, length(cc)==length(flg))
    lmu    <- params[[1]]
    lsigma <- params[[2]]
    if (lsigma < 0) return(NA)
    ll <- sum(dplyr::if_else(flg,
                      stats::plnorm(cc, lmu, lsigma, log.p = TRUE),  # Density below DL
                      stats::dlnorm(cc, lmu, lsigma, log = TRUE)))   # And normal obs.
    return(ll)
}

#' Calculate maximum likelihood parameter estimates from censored data.
#'
#' Calculate maximum likelihood parameter estimates for underlying (uncensored)
#' distribution based on left censored data.  This is a thin wrapper around
#' maxLik::maxLik() that strips any missing values and sets it up with
#' left_censored_loglik() and a simple default starting point for optimization.
#' Unlike maxLIk(), the return value is ONLY the parameter values, not the whole
#' maxLik object.
#'
#' @references
#'  Henningsen, Arne and Toomet, Ott (2011). maxLik: A package for maximum
#'  likelihood estimation in R. Computational Statistics 26(3), 443-458. DOI
#'  10.1007/s00180-010-0217-1.
#'
#' @param cc A vector of data values, including both observed values, where they
#'   exist, or the detection limits, where data was censored.
#'
#' @param flg A vector of TRUE or FALSE values, of the same length as cc, that
#'   indicates which values are detection limits (TRUE) and which are measured
#'   values (FALSE).
#'
#' @param start  A list or vector containing parameters of the underlying
#'   uncensored probability distribution (e.g., mean and SD of the related
#'   normal distribution for the default lognormal distribution).  This is used
#'   to initialize the numerical search for maximum likelihood estimates.  A
#'   good starting value may help speed the convergence. If convergence is not
#'   achieved, consider providing a better starting point for the optimization.
#'
#' @return Vector of (named) Maximum likelihood estimators of parameters to the
#'   underlying (uncensored) probability distribution.  (By default, a lognormal
#'   distribution). The calculated log likelihood under a censored lognormal
#'   distribution with the specific parameters.
#'
#' @examples
#'  df <- data.frame(sim = sort(stats::rlnorm(25,2,3)),
#'                   cens=c(rep.int(TRUE,5), rep.int(FALSE,20)))
#'  df$sim[1:4] <- df$sim[5]
#'  find_parms(df$sim, df$cens, c(2,5))
#'  find_parms(df$sim, df$cens, c(lmu=2,lsigma=5))
#'  find_parms(df$sim, df$cens, c(1,1))
#' @export
find_parms <- function(cc, flg, start = c(1,1))  {
  stopifnot(length(cc) == length(flg))
  cc2  <- cc [! (is.na(cc) | is.na(flg))]
  flg2 <- flg[! (is.na(cc) | is.na(flg))]
  res <- maxLik::maxLik(left_censored_loglik, start = start,
                        cc = cc2,
                        flg = flg2)
  parms <- res$estimate
  if (length(names(start))>0) {
    names(parms) <- names(start)
  }
  else {
    names(parms) <- paste0('p', 1:length(parms))
  }
  return(parms)
}


.simulate_mean_censored_single <- function(lmu, lsigma, cutoff,  sz = 1000) {
  estsamplesz <- sz + sz / stats::plnorm(cutoff, lmu, lsigma)  # Guess how many
  rawsample <- stats::rlnorm(estsamplesz, lmu, lsigma)         # Calculate extras
  sample <- rawsample[rawsample < cutoff]               # Throw out if too large
  sample <- sample[1:sz]                                # Take the sample
  return(mean(sample))                                  # Calculate the mean
}

#' Estimate the mean of unobserved left censored values
#'
#' \code{simulate_mean_censored} estimates conditional means of
#' values from specified lognormal distribution (with specific parameters) that
#' fall below (variable) cutoff values.  Estimated means are determined via
#' simulation.
#'
#' This functions in inefficient in the current context, as it simulates means
#' for all values, not just censored values.  That will make it slow for large
#' data sets or high values of sz.
#'
#' @param lmu The mean (on the log scale) of the lognormal distribution.
#'
#' @param lsigma The standard deviation of the lognormal distribution.
#'
#' @param cutoff Vector of detection limits (or observed values).
#'
#' @param sz The target size of the (simulated) sample for calculating
#'  conditional means (default = 1000). The current approach simulates a large
#'  draw from the underlying uncensored lognormal distribution, and retains
#'  only values below the detection limits.  The function estimates the size of
#'  the oversample needed, but because this is based on probabalistic reasoning,
#'  the sample is not guaranteed to always be exactly the size requested.
#'
#' @returns a vector of estimated conditional means.  Note that this function
#' provides conditional means for all observations, not only the censored ones.
#' This is wasteful, and potentially confusing to users.  See examples.
#'
#' @examples
#' df <- data.frame(sim = sort(stats::rlnorm(25,2,3)),
#'                 cens=c(rep.int(TRUE,5), rep.int(FALSE,20)))
#' df$sim[1:4] <- df$sim[5]
#' est <- simulate_mean_censored(lmu = 2, lsigma = 3, cutoff = df$sim)
#' library(ggplot2)
#' ggplot(df, aes(x = 1:25)) +
#' geom_line(aes(y = sim, color=cens)) +
#'   geom_point(aes(y = est)) +
#'   scale_y_log10() +
#'   scale_color_discrete(name = 'Censored') +
#'   theme_minimal() +
#'   xlab('Rank Order') +
#'   ylab('Raw Data (Line) and Conditional Means (points)')
#' rm(est)
#' @export
simulate_mean_censored <-
  Vectorize(.simulate_mean_censored_single, "cutoff")

#' Replace censored values with estimated conditional means
#'
#' \code{sub_conditional_means} replaces left censored values with estimated
#' conditional means.  The means are conditioned upon the fact that the
#' observation was censored.  All one knows about the value of a censored
#' observation is that the true value lies below a detection limit or other
#' threshold.  The estimated conditional means use information on the number
#' of observations below the detection limit and the distribution of values
#' above the detection limit to estimate censored values.
#'
#' An assumption of the method, however, is that all observations come from a
#' single underlying distribution.
#'
#' Thus if the goal is to compare concentrations of contaminants from different
#' populations, the correction should be applied separately to each population
#' before conducting additional analyses.
#'
#' These procedures may not be well suited for use where a covariate may alter
#' conditional means. For example, where rainfall or river discharge has a large
#' effect on concentrations of pollutants of interest, the assumption that all
#' observations are drawn from a single lognormal distribution may be untenable.
#' In practice, however, if censored observations are infrequent, the effect on
#' further analyses is likely to be small, and these methods may still be
#' preferable to making arbitrary choices about what value to use to replace
#' non-detects.
#'
#' The method simulates conditional means by drawing from a best fit lognormal
#' distribution, selected based on maximum likelihood.  Because the analysis is
#' based on simulation, results will not be identical for subsequent runs.
#'
#' @param cc A vector of data values, including both observed values, where they
#'   exist, or applicable detection limits, where data was censored.
#'
#' @param flg A vector of TRUE or FALSE values, of the same length as cc, that
#'    indicates which values are detection limits (TRUE) and which are measured
#'    values (FALSE).
#'
#' @param sz The target size of the (simulated) sample for calculating
#'  conditional means (default = 1000).
#'
#' @param start  A list or vector containing parameters of the underlying
#'   uncensored probability distribution (e.g., mean and SD of the related normal
#'   distribution for the default lognormal distribution).  This is used to
#'   initialize the numerical search for maximum likelihood estimates.  A good
#'   starting value may help speed the convergence. If convergence is not
#'   achieved, consider providing a better starting point for the optimization.
#'
#' @returns a vector containing original uncensored values, and estimates of the
#' conditional means (expected value) of censored observations. Detection limits
#' for  censored observations can differ.
#'
#' @examples
#' df = data.frame(sim = sort(stats::rlnorm(25,2,3)),
#'                 cens=c(rep.int(TRUE,5), rep.int(FALSE,20)))
#' df$sim[1:4] <- df$sim[5]
#' vals <- sub_conditional_means(cc = df$sim, flg =df$cens)
#' library(ggplot2)
#' ggplot(df, aes(x = 1:25)) +
#' geom_line(aes(y = sim, color=cens)) +
#'   geom_point(aes(y = vals)) +
#'   scale_y_log10() +
#'   scale_color_discrete(name = 'Censored') +
#'   theme_minimal() +
#'   xlab('Rank Order') +
#'   ylab('Raw Data (Line) and Data with Substitutions (points)')
#' rm(vals)
#' @export
sub_conditional_means <- function(cc, flg, sz = 1001, start= c(1,1)) {
  # Calculate maximum likelihood parameter estimates
  es <- find_parms(cc = cc,
                    flg = flg,
                    start = start)
  lmu <- es[1]
  lsigma <- es[2]
  # Calculate replacement values for censored observations
  res <- dplyr::if_else(flg,
                simulate_mean_censored(lmu, lsigma, cc, sz = sz),
                cc)
  return(res)
}
