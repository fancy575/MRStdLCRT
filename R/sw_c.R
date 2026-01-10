#' Example of stepped wedge CRT dataset for continuous outcome
#'
#' A toy dataset with cluster, period, and individual records for illustrating
#' estimands in stepped wedge CRT with a continuous outcome.
#'
#' @format A data frame with columns:
#' \describe{
#'   \item{cluster}{Cluster identifier (integer).}
#'   \item{period}{Period index (integer).}
#'   \item{id}{Individual identifier within cluster-period (integer).}
#'   \item{trt}{Treatment indicator (0/1).}
#'   \item{x1}{Auxiliary covariate (0/1).}
#'   \item{x2}{Auxiliary covariate (numeric).}
#'   \item{y}{Outcome (numeric, continuous).}
#' }
#'
#' @usage data(sw_c)
#' @keywords datasets
#' @examples
#' data(sw_c)
#' head(sw_c)
"sw_c"
