#' Concentrations of PAH contaminants in Casco Bay blue mussels.
#'
#' A dataset containing concentrations of PAH contaminants from
#' blue mussels (\emph{Mytilus edulis}) from Casco Bay, Maine, USA.
#' All concentrations are in ng/g, or (approximately) ppb.
#'
#'  Data are derived from NOAA Gulfwatch data collected from
#'  1993 through 2000.
#'
#' Data exhibits many of the challenges of working with non-detects,
#' in a small, compact form.
#' \itemize{
#' \item{Some compounds were observed all or almost all the time.}
#' \item{Some compounds were not sampled every year, generating
#' missing values}
#' \item{Some compounds were detected rarely, or even never detected.}
#' }
#'
#' Note that format of data follows a specific convention for handlign non-detects.
#' One data column contains a column of concentrations (including
#' detection/reporting limits for "Non-detects" and a column if flags
#' (logical TRUE/FALSE values) indicating which values represent
#' non-detects.
#'
#' @format A data frame with 780 rows and 5 variables:
#' \describe{
#'   \item{Site}{A four-letter site code.  Casco Bay sites include
#'   "MEPH" for Portland Harbor, "MEPR" for the Presumpscot River
#'   estuary, "MYRY" for the Royal River Estuary, and "MEBC" for
#'   Broad Cove.}
#'   \item{Year}{Year the sample was collected, as an integer}
#'   \item{PAH}{Chemical name, as a string.  Note that many of
#'   the chemical names can not be used as syntactic names in R.}
#'   \item{Flag}{Logical (TRUE/FALSE) value indicating whether the
#'   specific contaminant was detected or not. A TRUE value means
#'   the associated compound was not detected.  The Concentration
#'   value then represents the associated detection limit.}
#'   \item{Concentration}{Concentration expressed in ng/g.  Value
#'   either represents the observed value (if Flag == FALSE) or a
#'   related detection or reporting limit (if Flag == TRUE).}
#' }
#' @source \url{https://gulfofmaine.org/public/gulfwatch-contaminants-monitoring/findings/}
"mussels"
