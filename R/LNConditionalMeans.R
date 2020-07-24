#' Calculate the log likelihood of right censored lognormal data.
#'
#' Calculate the log likelihood of right censored lognormal data.  The
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
#' @return
#'  The calculated log likelihood under a censored lognormal distribution with
#'  the specific parameters.
#'
#' @examples
#' df = data.frame(sim = sort(rlnorm(25,2,3)),
#'                   cens=c(rep.int(TRUE,5), rep.int(FALSE,20)))
#' df$sim[1:4] <- df$sim[5]
#' left_censored_loglik(c(2,5), df$sim, df$cens)
#' left_censored_loglik(c(2,3), df$sim, df$cens)
#' left_censored_loglik(c(1,3), df$sim, df$cens)
left_censored_loglik <- function(params, cc, flg) {
    stopifnot(length(params)==2, length(cc)==length(flg))
    lmu    <- params[[1]]
    lsigma <- params[[2]]
    if (lsigma < 0) return(NA)
    ll <- sum(dplyr::if_else(flg,
                      plnorm(cc, lmu, lsigma, log.p = TRUE),  # Density below DL
                      dlnorm(cc, lmu, lsigma, log = TRUE)))   # And normal obs.
    return(ll)
}

.simulate_mean_censored_single <- function(lmu, lsigma, cutoff,  sz = 1000) {
  estsamplesz <- sz + sz / plnorm(cutoff, lmu, lsigma)  # Guess how many
  rawsample <- rlnorm(estsamplesz, lmu, lsigma)         # Calculate extras
  sample <- rawsample[rawsample < cutoff]               # Throw out if too large
  sample <- sample[1:sz]                                # Take the sample
  return(mean(sample))                                  # Calculate the mean
}

simulate_mean_censored <-
  Vectorize(.simulate_mean_censored_single, "cutoff")

sub_conditional_means <- function(cc, flg) {
  cat(flg)
  # Calculate maximum likelihood parameter estimates
  r <- maxLik::maxLik(left_censored_loglik, start = c(-2, 3), cc = cc, flg = flg)
  es <- r$estimate
  lmu <- es[1]
  lsigma <- es[2]
  # Calculate replacement values for censored observations
  cat(flg)
  res <- ifelse(flg,
                simulate_mean_censored(lmu, lsigma, cc),
                cc)
  return(res)
}
