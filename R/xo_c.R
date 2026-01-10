#' Example crossover cluster-randomized trial dataset with continuous outcome
#'
#' A small simulated 2×2 crossover trial dataset with a continuous outcome.
#'
#' @format A tibble/data.frame with one row per subject and the following columns:
#' \describe{
#'   \item{h}{Integer cluster ID (hospital)}
#'   \item{p}{Integer period (1 or 2)}
#'   \item{k}{Integer subject index within cluster-period}
#'   \item{trt}{Treatment indicator (0 = control, 1 = treatment)}
#'   \item{x_c01, x_c02}{Continuous covariates}
#'   \item{x_b01}{Binary covariate (0/1)}
#'   \item{x_cat1_2, x_cat1_3}{Dummy variables for a 3-level categorical covariate (level 1 is reference)}
#'   \item{y_cont}{Observed continuous outcome}
#' }
#'
#' @examples
#' data(xo_c)
#' str(xo_c)
#' head(xo_c)
"xo_c"
