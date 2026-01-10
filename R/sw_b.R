#' Example stepped wedge CRT dataset with binary outcome
#'
#' A toy dataset with cluster, period, and individual records for illustrating
#' estimands in stepped wedge CRT with a binary outcome.
#'
#' @format A data frame with columns:
#' \describe{
#'   \item{cluster}{Cluster identifier (integer).}
#'   \item{period}{Period index (integer).}
#'   \item{id}{Individual identifier within cluster-period (integer).}
#'   \item{trt}{Treatment indicator (0/1).}
#'   \item{x1}{Auxiliary covariate (0/1).}
#'   \item{x2}{Auxiliary covariate (numeric).}
#'   \item{y}{Outcome (0/1, binary).}
#' }
#'
#' @usage data(sw_b)
#' @keywords datasets
#' @examples
#' data(sw_b)
#' head(sw_b)
"sw_b"
