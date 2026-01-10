#' @importFrom rlang .data
NULL


utils::globalVariables(c(
  ".data",
  "Nij","Zij","Ybar_ij","nZ","minZ","maxZ","is_mixed",
  "N_j","N_i","m0","m1","mz","pot_z","omega"
))

# =============================================================================
# Small utilities
# =============================================================================

mrs_inv_logit <- function(eta) 1/(1 + exp(-eta))
mrs_logit     <- function(p) log(p/(1 - p))
mrs_sdiv      <- function(a,b) ifelse(b == 0, NA_real_, a/b)

mrs_jk_cov <- function(mat_Ix4) {
  I <- nrow(mat_Ix4)
  theta_bar <- colMeans(mat_Ix4)
  D <- t(t(mat_Ix4) - theta_bar)
  D <- t(D)
  ((I - 1) / I) * (D %*% t(D))
}

# Remove random effects terms for GEE formulas
mrs_strip_re_from_formula <- function(fml) {
  f_txt <- paste(deparse(fml), collapse = " ")
  f_txt <- gsub("\\([^\\)]*\\|[^\\)]*\\)", "", f_txt) # remove ( ... | ... )
  f_txt <- gsub("\\+\\s*\\+", "+", f_txt)
  f_txt <- gsub("\\s*\\+\\s*$", "", f_txt)
  f_txt <- gsub("\\s{2,}", " ", f_txt)
  f_txt <- trimws(f_txt)
  stats::as.formula(f_txt, env = parent.frame())
}

# =============================================================================
# Cluster-period table + period mixture rule
# =============================================================================

mrs_make_cp_table <- function(data, cluster_id, period, trt, y_name) {
  cl_sym  <- rlang::sym(cluster_id)
  per_sym <- rlang::sym(period)
  trt_sym <- rlang::sym(trt)
  y_sym   <- rlang::sym(y_name)

  cp0 <- data |>
    dplyr::group_by(!!cl_sym, !!per_sym) |>
    dplyr::summarise(
      Nij     = dplyr::n(),
      Zij     = dplyr::first(!!trt_sym),
      nZ      = dplyr::n_distinct(!!trt_sym),
      Ybar_ij = mean(!!y_sym),
      .groups = "drop"
    )

  if (any(cp0$nZ > 1, na.rm = TRUE)) {
    bad <- cp0 |> dplyr::filter(.data$nZ > 1)
    stop(
      "Within at least one cluster-period cell, 'trt' is not constant.\n",
      "Please ensure treatment is constant within each cluster-period.\n",
      "Example problematic cells:\n",
      paste(utils::capture.output(utils::head(bad, 10)), collapse = "\n"),
      call. = FALSE
    )
  }

  cp0 |> dplyr::select(-.data$nZ)
}

mrs_period_mixture_table <- function(cp_tab, period, z_col = "Zij") {
  per_sym <- rlang::sym(period)
  z_sym   <- rlang::sym(z_col)

  cp_tab |>
    dplyr::group_by(!!per_sym) |>
    dplyr::summarise(
      minZ = min(!!z_sym, na.rm = TRUE),
      maxZ = max(!!z_sym, na.rm = TRUE),
      .groups = "drop"
    ) |>
    dplyr::mutate(is_mixed = (.data$minZ == 0 & .data$maxZ == 1)) |>
    dplyr::arrange(!!per_sym)
}

mrs_kept_periods_from_mixture <- function(per_mix_tbl, period) {
  per_sym <- rlang::sym(period)
  per_mix_tbl |>
    dplyr::filter(.data$is_mixed) |>
    dplyr::pull(!!per_sym)
}

mrs_detect_sw_edge_pattern <- function(per_mix_tbl, period) {
  if (is.null(per_mix_tbl) || !nrow(per_mix_tbl)) {
    return(list(is_sw_edge = NA, first_period = NA, last_period = NA))
  }
  per_mix_tbl <- per_mix_tbl[order(per_mix_tbl[[period]]), , drop = FALSE]
  first_row <- per_mix_tbl[1, , drop = FALSE]
  last_row  <- per_mix_tbl[nrow(per_mix_tbl), , drop = FALSE]

  first_all_control <- isTRUE(first_row$minZ == 0 && first_row$maxZ == 0)
  last_all_treated  <- isTRUE(last_row$minZ  == 1 && last_row$maxZ  == 1)

  list(
    is_sw_edge   = first_all_control && last_all_treated,
    first_period = first_row[[period]],
    last_period  = last_row[[period]]
  )
}

# =============================================================================
# Working model fit + prediction (rank-deficiency safe)
# =============================================================================

mrs_fit_working_model <- function(data_fit, formula, cluster_id, method, family, corstr) {
  method <- match.arg(method, c("gee","lmer","glmer"))
  family <- match.arg(family, c("gaussian","binomial"))

  if (method == "lmer"  && family != "gaussian") stop("lmer requires family='gaussian'.")
  if (method == "glmer" && family != "binomial") stop("glmer requires family='binomial'.")

  if (method == "gee") {
    fam <- if (family == "gaussian") stats::gaussian() else stats::binomial(link = "logit")
    f_use <- mrs_strip_re_from_formula(formula)

    fit <- NULL
    utils::capture.output({
      fit <- suppressWarnings(gee::gee(
        formula   = f_use,
        id        = data_fit[[cluster_id]],
        data      = data_fit,
        family    = fam,
        corstr    = corstr,
        na.action = stats::na.omit
      ))
    }, type = "output")

    # IMPORTANT: store the evaluated formula (avoid symbol-only call storage)
    fit$formula <- f_use
    fit$call$formula <- f_use

    # Optional but helpful for factors in newdata
    mf_train <- stats::model.frame(f_use, data = data_fit, na.action = stats::na.omit)
    fit$xlevels <- attr(mf_train, "xlevels")

    return(fit)
  }


  if (method == "lmer")  return(lme4::lmer(formula, data = data_fit))
  if (method == "glmer") return(lme4::glmer(formula, data = data_fit,
                                            family = stats::binomial(link = "logit")))
  stop("Unknown method.")
}

mrs_re_var_total <- function(fit) {
  vc <- lme4::VarCorr(fit)
  out <- 0
  for (g in names(vc)) {
    M <- as.matrix(vc[[g]])
    out <- out + sum(diag(M))
  }
  out
}

# Predict mean response under newdata (fixed-effects only, aligned columns)
mrs_predict_mean <- function(fit, newdata, method, family) {
  method <- match.arg(method, c("gee","lmer","glmer"))
  family <- match.arg(family, c("gaussian","binomial"))

  if (inherits(fit, "merMod")) {
    # fixed-effect formula
    f_fix <- reformulas::nobars(stats::formula(fit))
    TT <- stats::terms(f_fix, data = newdata)

    # factor levels as in training
    xlev <- NULL
    if (!is.null(fit@frame) && !is.null(fit@frame$xlevels)) xlev <- fit@frame$xlevels

    mf <- stats::model.frame(TT, data = newdata, xlev = xlev, na.action = stats::na.pass)

    # carry contrasts (helps consistency)
    contr <- NULL
    Xtrain <- tryCatch(lme4::getME(fit, "X"), error = function(e) NULL)
    if (!is.null(Xtrain)) contr <- attr(Xtrain, "contrasts")

    X <- stats::model.matrix(TT, mf, contrasts.arg = contr)

    b <- lme4::fixef(fit)
    keep <- intersect(colnames(X), names(b))
    if (!length(keep)) stop("No overlapping fixed-effect columns between newdata and fitted model.")

    eta <- as.numeric(X[, keep, drop = FALSE] %*% b[keep])

    if (inherits(fit, "lmerMod")) return(eta)

    # glmer binomial: approximate marginal mean using logistic-normal scale adjustment
    sig2 <- mrs_re_var_total(fit)
    sc   <- sqrt(1 + (3/pi^2) * sig2)
    return(mrs_inv_logit(eta / sc))
  }

  if (method == "gee") {
    fml <- fit$formula
    if (is.null(fml)) fml <- fit$call$formula

    # If still a symbol/name, we cannot safely recover it later; require stored formula
    if (is.symbol(fml) || is.name(fml)) {
      stop("GEE fit contains only a symbolic formula reference. Store the evaluated formula in fit$formula when fitting.")
    }

    fml <- stats::as.formula(fml)

    TT  <- stats::delete.response(stats::terms(fml, data = newdata))
    mf  <- stats::model.frame(TT, data = newdata, xlev = fit$xlevels, na.action = stats::na.pass)
    X   <- stats::model.matrix(TT, mf)

    b <- stats::coef(fit)
    keep <- intersect(colnames(X), names(b))
    if (!length(keep)) stop("No overlapping columns between model matrix and GEE coefficients.")

    eta <- as.numeric(X[, keep, drop = FALSE] %*% b[keep])
    if (family == "gaussian") return(eta)
    return(mrs_inv_logit(eta))
  }




  stop("Unknown fit class/method.")
}

# =============================================================================
# Estimand scales / contrasts
# =============================================================================

mrs_contrast <- function(mu1, mu0, family, scale = c("RD","RR","OR")) {
  family <- match.arg(family, c("gaussian","binomial"))
  if (family == "gaussian") return(mu1 - mu0)

  scale <- match.arg(scale)
  eps <- 1e-12
  mu1c <- pmin(pmax(mu1, eps), 1 - eps)
  mu0c <- pmin(pmax(mu0, eps), 1 - eps)

  if (scale == "RD") return(mu1 - mu0)
  if (scale == "RR") return(log(mu1c/mu0c))                  # log RR
  if (scale == "OR") return(mrs_logit(mu1c) - mrs_logit(mu0c))# log OR
  stop("Unknown scale.")
}

# =============================================================================
# Period-specific unadjusted & augmented means, then aggregate
# =============================================================================

mrs_mu_period_tables <- function(cp_k, period, w_col, omega_col) {
  per_sym <- rlang::sym(period)
  w_sym   <- rlang::sym(w_col)
  om_sym  <- rlang::sym(omega_col)

  # unadjusted mu_{j}(z)
  mu_unadj_jz <- cp_k |>
    dplyr::group_by(!!per_sym, .data$Zij) |>
    dplyr::summarise(
      num = sum(!!w_sym * .data$Ybar_ij),
      den = sum(!!w_sym),
      mu  = mrs_sdiv(.data$num, .data$den),
      .groups = "drop"
    ) |>
    dplyr::rename(Z = .data$Zij) |>
    dplyr::select(!!per_sym, .data$Z, .data$mu)

  # augmented mu_{j}(z): term2 + term1
  mu_aug_jz <- cp_k |>
    tidyr::crossing(pot_z = c(0L, 1L)) |>
    dplyr::mutate(mz = ifelse(.data$pot_z == 0L, .data$m0, .data$m1)) |>
    dplyr::group_by(!!per_sym, .data$pot_z) |>
    dplyr::summarise(
      sum_w_mz = sum(!!w_sym * .data$mz),
      num_t1   = sum(!!w_sym * (.data$Zij == .data$pot_z) * (.data$Ybar_ij - .data$mz)),
      den_t1   = sum(!!w_sym * (.data$Zij == .data$pot_z)),
      omega    = dplyr::first(!!om_sym),
      .groups  = "drop"
    ) |>
    dplyr::mutate(
      term1  = mrs_sdiv(.data$num_t1, .data$den_t1),
      term2  = mrs_sdiv(.data$sum_w_mz, .data$omega),
      mu_aug = .data$term1 + .data$term2
    ) |>
    dplyr::rename(Z = .data$pot_z) |>
    dplyr::select(!!per_sym, .data$Z, mu = .data$mu_aug)

  list(mu_unadj_jz = mu_unadj_jz, mu_aug_jz = mu_aug_jz)
}

mrs_aggregate_mu <- function(mu_jz, omega_j, period, omega_col) {
  per_sym <- rlang::sym(period)
  om_sym  <- rlang::sym(omega_col)

  mu_jz |>
    dplyr::left_join(
      omega_j |> dplyr::select(!!per_sym, !!om_sym),
      by = stats::setNames(period, period)
    ) |>
    dplyr::rename(omega = !!om_sym) |>
    dplyr::group_by(.data$Z) |>
    dplyr::summarise(
      mu = sum(.data$omega * .data$mu, na.rm = TRUE) / sum(.data$omega, na.rm = TRUE),
      .groups = "drop"
    ) |>
    dplyr::arrange(.data$Z)
}

mrs_mu_pair <- function(mu_tbl) {
  c(mu0 = mu_tbl$mu[mu_tbl$Z == 0], mu1 = mu_tbl$mu[mu_tbl$Z == 1])
}

# =============================================================================
# Core computation on full data fit; aggregate on mixed periods only
# =============================================================================

mrs_fit_core <- function(
    data, formula,
    cluster_id = "cluster",
    period     = "period",
    trt        = "trt",
    method     = c("gee","lmer","glmer"),
    family     = c("gaussian","binomial"),
    corstr     = "independence",
    scale      = c("RD","RR","OR")
) {
  method <- match.arg(method)
  family <- match.arg(family)

  if (!all(c(cluster_id, period, trt) %in% names(data))) {
    stop("data must contain columns: cluster_id, period, trt.")
  }

  # coerce trt to 0/1 integer
  data[[trt]] <- as.integer(data[[trt]])
  if (!all(stats::na.omit(unique(data[[trt]])) %in% c(0L, 1L))) {
    stop("Treatment column must be coded as 0/1 (or coercible to 0/1).")
  }

  y_name <- all.vars(formula)[1]
  if (!length(y_name) || !(y_name %in% names(data))) stop("Outcome variable not found in data.")

  # complete cases for variables in formula + IDs
  need_vars <- unique(c(all.vars(formula), cluster_id, period, trt))
  need_vars <- intersect(need_vars, names(data))
  data_fit <- data[stats::complete.cases(data[, need_vars, drop = FALSE]), , drop = FALSE]
  if (!nrow(data_fit)) stop("No complete cases after removing missing values.")

  # CP table on full data (for mixture + SW edge detection)
  cp_full <- mrs_make_cp_table(data_fit, cluster_id, period, trt, y_name)
  per_mix_full <- mrs_period_mixture_table(cp_full, period, z_col = "Zij")
  kept_periods <- mrs_kept_periods_from_mixture(per_mix_full, period)
  if (!length(kept_periods)) stop("No periods have mixed treatment assignment (both 0 and 1). Cannot proceed.")

  sw_edge <- mrs_detect_sw_edge_pattern(per_mix_full, period)

  # Fit working model on FULL data
  fit <- mrs_fit_working_model(data_fit, formula, cluster_id, method, family, corstr)

  # Restrict to kept periods for aggregation
  data_k <- data_fit[data_fit[[period]] %in% kept_periods, , drop = FALSE]
  cp_k   <- mrs_make_cp_table(data_k, cluster_id, period, trt, y_name)

  cl_sym  <- rlang::sym(cluster_id)
  per_sym <- rlang::sym(period)

  # period totals N_j on kept periods
  N_j_tab <- cp_k |>
    dplyr::group_by(!!per_sym) |>
    dplyr::summarise(N_j = sum(.data$Nij), .groups = "drop")

  # cluster totals N_i on kept periods (for h-cATE)
  N_i_tab <- cp_k |>
    dplyr::group_by(!!cl_sym) |>
    dplyr::summarise(N_i = sum(.data$Nij), .groups = "drop")

  cp_k <- cp_k |>
    dplyr::left_join(N_j_tab, by = stats::setNames(period, period)) |>
    dplyr::left_join(N_i_tab, by = stats::setNames(cluster_id, cluster_id))

  # Predictions on kept-period individual data (fit from full data)
  nd0 <- data_k; nd0[[trt]] <- 0L
  nd1 <- data_k; nd1[[trt]] <- 1L

  mu0_i <- mrs_predict_mean(fit, nd0, method = method, family = family)
  mu1_i <- mrs_predict_mean(fit, nd1, method = method, family = family)

  pred_i <- data_k
  pred_i$.mu0 <- mu0_i
  pred_i$.mu1 <- mu1_i

  cp_pred <- pred_i |>
    dplyr::group_by(!!cl_sym, !!per_sym) |>
    dplyr::summarise(
      m0 = mean(.data$.mu0),
      m1 = mean(.data$.mu1),
      .groups = "drop"
    )

  cp_k <- cp_k |>
    dplyr::left_join(cp_pred, by = stats::setNames(c(cluster_id, period), c(cluster_id, period)))

  # weights (kept periods)
  cp_k <- cp_k |>
    dplyr::mutate(
      w_hi = .data$Nij,
      w_hc = mrs_sdiv(.data$Nij, .data$N_i),
      w_vi = mrs_sdiv(.data$Nij, .data$N_j),
      w_vc = 1
    )

  # period weights omega_j for each estimand scheme
  omega_j <- cp_k |>
    dplyr::group_by(!!per_sym) |>
    dplyr::summarise(
      omega_hi = sum(.data$w_hi),
      omega_hc = sum(.data$w_hc),
      omega_vi = 1,
      omega_vc = 1,
      .groups  = "drop"
    )

  # attach omega columns to cp_k for period-specific augmented formula
  cp_k <- cp_k |>
    dplyr::left_join(omega_j, by = stats::setNames(period, period))

  # period tables
  tabs_hi <- mrs_mu_period_tables(cp_k, period, "w_hi", "omega_hi")
  tabs_hc <- mrs_mu_period_tables(cp_k, period, "w_hc", "omega_hc")
  tabs_vi <- mrs_mu_period_tables(cp_k, period, "w_vi", "omega_vi")
  tabs_vc <- mrs_mu_period_tables(cp_k, period, "w_vc", "omega_vc")

  # aggregate mu(z)
  mu_hi_unadj <- mrs_aggregate_mu(tabs_hi$mu_unadj_jz, omega_j, period, "omega_hi")
  mu_hc_unadj <- mrs_aggregate_mu(tabs_hc$mu_unadj_jz, omega_j, period, "omega_hc")
  mu_vi_unadj <- mrs_aggregate_mu(tabs_vi$mu_unadj_jz, omega_j, period, "omega_vi")
  mu_vc_unadj <- mrs_aggregate_mu(tabs_vc$mu_unadj_jz, omega_j, period, "omega_vc")

  mu_hi_aug <- mrs_aggregate_mu(tabs_hi$mu_aug_jz, omega_j, period, "omega_hi")
  mu_hc_aug <- mrs_aggregate_mu(tabs_hc$mu_aug_jz, omega_j, period, "omega_hc")
  mu_vi_aug <- mrs_aggregate_mu(tabs_vi$mu_aug_jz, omega_j, period, "omega_vi")
  mu_vc_aug <- mrs_aggregate_mu(tabs_vc$mu_aug_jz, omega_j, period, "omega_vc")

  mu_unadj <- rbind(
    "h-iATE" = mrs_mu_pair(mu_hi_unadj),
    "h-cATE" = mrs_mu_pair(mu_hc_unadj),
    "v-iATE" = mrs_mu_pair(mu_vi_unadj),
    "v-cATE" = mrs_mu_pair(mu_vc_unadj)
  )
  mu_aug <- rbind(
    "h-iATE" = mrs_mu_pair(mu_hi_aug),
    "h-cATE" = mrs_mu_pair(mu_hc_aug),
    "v-iATE" = mrs_mu_pair(mu_vi_aug),
    "v-cATE" = mrs_mu_pair(mu_vc_aug)
  )

  scale_use <- if (family == "gaussian") "RD" else match.arg(scale)

  ate_unadj <- c(
    "h-iATE" = mrs_contrast(mu_unadj["h-iATE","mu1"], mu_unadj["h-iATE","mu0"], family, scale_use),
    "h-cATE" = mrs_contrast(mu_unadj["h-cATE","mu1"], mu_unadj["h-cATE","mu0"], family, scale_use),
    "v-iATE" = mrs_contrast(mu_unadj["v-iATE","mu1"], mu_unadj["v-iATE","mu0"], family, scale_use),
    "v-cATE" = mrs_contrast(mu_unadj["v-cATE","mu1"], mu_unadj["v-cATE","mu0"], family, scale_use)
  )
  ate_aug <- c(
    "h-iATE" = mrs_contrast(mu_aug["h-iATE","mu1"], mu_aug["h-iATE","mu0"], family, scale_use),
    "h-cATE" = mrs_contrast(mu_aug["h-cATE","mu1"], mu_aug["h-cATE","mu0"], family, scale_use),
    "v-iATE" = mrs_contrast(mu_aug["v-iATE","mu1"], mu_aug["v-iATE","mu0"], family, scale_use),
    "v-cATE" = mrs_contrast(mu_aug["v-cATE","mu1"], mu_aug["v-cATE","mu0"], family, scale_use)
  )

  # counts on kept periods (for summary)
  counts_period <- cp_k |>
    dplyr::group_by(!!per_sym) |>
    dplyr::summarise(
      n_obs      = sum(.data$Nij),
      n_clusters = dplyr::n_distinct(.data[[cluster_id]]),
      n_cp_cells = dplyr::n(),
      .groups = "drop"
    ) |>
    dplyr::arrange(!!per_sym)

  counts_cluster <- cp_k |>
    dplyr::group_by(!!cl_sym) |>
    dplyr::summarise(
      n_obs     = sum(.data$Nij),
      n_periods = dplyr::n_distinct(.data[[period]]),
      .groups = "drop"
    ) |>
    dplyr::arrange(!!cl_sym)

  list(
    fit = fit,
    kept_periods = kept_periods,
    period_mix_table = per_mix_full,
    sw_edge = sw_edge,
    omega_j = omega_j,
    mu_unadj = mu_unadj,
    mu_aug   = mu_aug,
    ate_unadj = ate_unadj,
    ate_aug   = ate_aug,
    counts = list(period = counts_period, cluster = counts_cluster),
    method = method,
    family = family,
    corstr = corstr,
    scale = scale_use
  )
}

# =============================================================================
# Public API
# =============================================================================

#' Fit model-robust standardization for longitudinal CRTs
#'
#' @param data data.frame with outcome, treatment, period, cluster, covariates.
#' @param formula model formula; may include interactions and random effects.
#' @param cluster_id cluster id column name.
#' @param period period column name.
#' @param trt treatment column name (0/1).
#' @param method "gee","lmer","glmer".
#' @param family "gaussian","binomial".
#' @param corstr gee correlation.
#' @param scale For binomial: "RD","RR","OR" (RR/OR are on log scale).
#'
#' @return Object of class "mrs".
#' @export
mrstdlcrt_fit <- function(
    data, formula,
    cluster_id = "cluster",
    period     = "period",
    trt        = "trt",
    method     = c("gee","lmer","glmer"),
    family     = c("gaussian","binomial"),
    corstr     = "independence",
    scale      = c("RD","RR","OR")
) {
  method <- match.arg(method)
  family <- match.arg(family)

  core <- mrs_fit_core(
    data = data, formula = formula,
    cluster_id = cluster_id, period = period, trt = trt,
    method = method, family = family, corstr = corstr, scale = scale
  )

  cl_sym <- rlang::sym(cluster_id)
  cl_ids <- unique(dplyr::pull(data, !!cl_sym))
  I <- length(cl_ids)
  if (I < 2) stop("Need at least 2 clusters for delete-1 jackknife.")

  reps_unadj <- matrix(NA_real_, nrow = I, ncol = 4, dimnames = list(NULL, names(core$ate_unadj)))
  reps_aug   <- reps_unadj

  for (g in seq_len(I)) {
    d_sub <- dplyr::filter(data, (!!cl_sym) != cl_ids[g])
    res_g <- mrs_fit_core(
      data = d_sub, formula = formula,
      cluster_id = cluster_id, period = period, trt = trt,
      method = method, family = family, corstr = corstr, scale = scale
    )
    reps_unadj[g, ] <- res_g$ate_unadj
    reps_aug[g, ]   <- res_g$ate_aug
  }

  jk_cov_unadj <- mrs_jk_cov(reps_unadj)
  jk_cov_aug   <- mrs_jk_cov(reps_aug)

  se_unadj <- sqrt(diag(jk_cov_unadj))
  se_aug   <- sqrt(diag(jk_cov_aug))

  estimates <- dplyr::tibble(
    method_type = c("unadjusted","adjusted"),
    `h-iATE` = c(core$ate_unadj["h-iATE"], core$ate_aug["h-iATE"]),
    `h-cATE` = c(core$ate_unadj["h-cATE"], core$ate_aug["h-cATE"]),
    `v-iATE` = c(core$ate_unadj["v-iATE"], core$ate_aug["v-iATE"]),
    `v-cATE` = c(core$ate_unadj["v-cATE"], core$ate_aug["v-cATE"])
  )

  jk_se <- dplyr::tibble(
    method_type  = c("unadjusted","adjusted"),
    `SE(h-iATE)` = c(se_unadj[1], se_aug[1]),
    `SE(h-cATE)` = c(se_unadj[2], se_aug[2]),
    `SE(v-iATE)` = c(se_unadj[3], se_aug[3]),
    `SE(v-cATE)` = c(se_unadj[4], se_aug[4])
  )

  meta <- list(
    call    = match.call(),
    formula = formula,
    method  = core$method,
    family  = core$family,
    scale   = core$scale,
    corstr  = if (method == "gee") corstr else NA_character_,
    n_clusters = I,
    n_periods_total = length(unique(data[[period]])),
    n_periods_kept  = length(core$kept_periods),
    kept_periods    = core$kept_periods,
    n_subjects = nrow(data),
    cluster_id = cluster_id,
    period     = period,
    trt        = trt,
    sw_edge_pattern = isTRUE(core$sw_edge$is_sw_edge),
    sw_first_period = core$sw_edge$first_period,
    sw_last_period  = core$sw_edge$last_period
  )

  res <- list(
    estimates    = estimates,
    jk_se        = jk_se,
    jk_cov_unadj = jk_cov_unadj,
    jk_cov_aug   = jk_cov_aug,
    reps         = core,
    meta         = meta
  )
  class(res) <- "mrs"
  res
}

# =============================================================================
# S3: print / summary / plot
# =============================================================================

#' Print method for mrs objects
#'
#' @param x An object of class \code{"mrs"}.
#' @param ... Unused.
#' @return \code{x} invisibly.
#' @export
print.mrs <- function(x, ...) {
  stopifnot(inherits(x, "mrs"))
  cat("\nMRS-LCRT fit\n")
  cat(rep("=", 72), "\n", sep = "")
  cat(sprintf("Method: %s   Family: %s   Scale: %s\n",
              x$meta$method, x$meta$family, x$meta$scale))
  if (!is.na(x$meta$corstr)) cat(sprintf("GEE corstr: %s\n", x$meta$corstr))
  cat(sprintf("Clusters: %d   Periods (total): %d   Periods (kept): %d\n",
              x$meta$n_clusters, x$meta$n_periods_total, x$meta$n_periods_kept))
  cat(sprintf("Kept periods: %s\n", paste(x$meta$kept_periods, collapse = ", ")))

  if (!is.null(x$meta$sw_edge_pattern) && !is.na(x$meta$sw_edge_pattern)) {
    if (isTRUE(x$meta$sw_edge_pattern)) {
      cat(sprintf("Stepped-wedge edge pattern: YES (first period %s all control; last period %s all treatment)\n",
                  as.character(x$meta$sw_first_period),
                  as.character(x$meta$sw_last_period)))
    } else {
      cat("Stepped-wedge edge pattern: NO (first period not all control and/or last period not all treatment)\n")
    }
  }

  cat(sprintf("N (rows): %d\n\n", x$meta$n_subjects))
  cat("Point estimates:\n")
  print(x$estimates)
  invisible(x)
}

#' Summarize an mrs fit
#'
#' @param object An object of class \code{"mrs"}.
#' @param level Confidence level.
#' @param estimand Optional subset of estimands.
#' @param digits Digits to print.
#' @param show_counts Print counts tables.
#' @param ... Unused.
#' @return Invisibly returns a list of printed tables and metadata.
#' @export
summary.mrs <- function(object,
                        level = 0.95,
                        estimand = NULL,
                        digits = 6,
                        show_counts = TRUE,
                        ...) {
  stopifnot(inherits(object, "mrs"))

  est  <- object$estimates
  se   <- object$jk_se
  meta <- object$meta
  reps <- object$reps

  valid <- c("h-iATE","h-cATE","v-iATE","v-cATE")
  sel <- if (is.null(estimand)) valid else intersect(estimand, valid)
  if (!length(sel)) sel <- valid

  df <- max(1L, meta$n_clusters - 1L)
  alpha <- 1 - level
  crit <- stats::qt(1 - alpha/2, df = df)

  cat("\nMRS-LCRT Summary\n")
  cat(rep("=", 72), "\n", sep = "")
  cat(sprintf("Method: %s   Family: %s   Scale: %s\n", meta$method, meta$family, meta$scale))
  if (!is.na(meta$corstr)) cat(sprintf("GEE corstr: %s\n", meta$corstr))
  cat(sprintf("Clusters: %d   Periods kept: %d (of %d total)\n",
              meta$n_clusters, meta$n_periods_kept, meta$n_periods_total))
  cat(sprintf("Kept periods: %s\n", paste(meta$kept_periods, collapse = ", ")))
  cat(sprintf("CIs: %.1f%% (t, df = %d)\n", 100*level, df))

  if (!is.null(meta$sw_edge_pattern) && !is.na(meta$sw_edge_pattern)) {
    if (isTRUE(meta$sw_edge_pattern)) {
      cat(sprintf("Stepped-wedge edge pattern: YES (first period %s all control; last period %s all treatment)\n",
                  as.character(meta$sw_first_period),
                  as.character(meta$sw_last_period)))
    } else {
      cat("Stepped-wedge edge pattern: NO (first period not all control and/or last period not all treatment)\n")
    }
  }
  cat("\n")

  if (!is.null(reps$period_mix_table) && nrow(reps$period_mix_table)) {
    cat("Period mixture table (min/max of Z_ij by period)\n")
    cat(rep("-", 72), "\n", sep = "")
    print(reps$period_mix_table)
    cat("\n")
  }

  if (isTRUE(show_counts) && !is.null(reps$counts)) {
    if (!is.null(reps$counts$period) && nrow(reps$counts$period)) {
      cat("Counts used in aggregation (kept periods only)\n")
      cat(rep("-", 72), "\n", sep = "")
      print(reps$counts$period)
      cat("\n")
    }
    if (!is.null(reps$counts$cluster) && nrow(reps$counts$cluster)) {
      cat("Per-cluster counts (kept periods only)\n")
      cat(rep("-", 72), "\n", sep = "")
      print(reps$counts$cluster)
      cat("\n")
    }
  }

  se_col <- function(nm) paste0("SE(", nm, ")")
  ratio_scale <- (meta$family == "binomial" && meta$scale %in% c("RR","OR"))

  printed <- list()

  for (nm in sel) {
    vals <- est[[nm]]
    ses  <- se[[se_col(nm)]]
    lcl  <- vals - crit * ses
    ucl  <- vals + crit * ses

    tab <- data.frame(
      Estimate = vals,
      SE       = ses,
      LCL      = lcl,
      UCL      = ucl,
      row.names = est$method_type,
      check.names = FALSE
    )

    cat(nm, "\n", sep = "")
    print(round(tab, digits))
    cat("\n")
    printed[[nm]] <- tab

    if (isTRUE(ratio_scale)) {
      tab2 <- data.frame(
        Ratio     = exp(vals),
        Ratio_LCL = exp(lcl),
        Ratio_UCL = exp(ucl),
        row.names = est$method_type,
        check.names = FALSE
      )
      cat(if (meta$scale == "RR") "Exponentiated (RR)\n" else "Exponentiated (OR)\n")
      print(round(tab2, digits))
      cat("\n")
      printed[[paste0(nm, "_ratio")]] <- tab2
    }
  }

  invisible(list(
    estimates = est,
    jk_se = se,
    shown = sel,
    level = level,
    df = df,
    crit = crit,
    meta = meta,
    tables = printed
  ))
}

#' Plot method for mrs objects
#'
#' @param x An object of class \code{"mrs"}.
#' @param level Confidence level.
#' @param estimand Subset of estimands to plot.
#' @param point_size Point size.
#' @param ... Unused.
#' @return ggplot object invisibly.
#' @export
plot.mrs <- function(x, level = 0.95, estimand = NULL, point_size = 2.8, ...) {
  if (!requireNamespace("ggplot2", quietly = TRUE))
    stop("Package 'ggplot2' is required for plotting.", call. = FALSE)
  if (!requireNamespace("dplyr", quietly = TRUE) || !requireNamespace("tidyr", quietly = TRUE))
    stop("Packages 'dplyr' and 'tidyr' are required for plotting.", call. = FALSE)

  est <- x$estimates
  se  <- x$jk_se

  valid <- c("h-iATE","h-cATE","v-iATE","v-cATE")
  sel <- if (is.null(estimand)) valid else intersect(estimand, valid)
  if (!length(sel)) sel <- valid

  df <- max(1L, x$meta$n_clusters - 1L)
  alpha <- 1 - level
  crit  <- stats::qt(1 - alpha/2, df = df)

  est_long <- tidyr::pivot_longer(est, cols = tidyselect::all_of(valid),
                                  names_to = "Estimand", values_to = "Estimate")
  se_long  <- tidyr::pivot_longer(se, cols = tidyselect::starts_with("SE("),
                                  names_to = "Estimand", values_to = "SE")
  se_long$Estimand <- sub("^SE\\((.*)\\)$", "\\1", se_long$Estimand)

  plotdf <- dplyr::left_join(est_long, se_long, by = c("method_type","Estimand")) |>
    dplyr::filter(.data$Estimand %in% sel) |>
    dplyr::mutate(
      Type = dplyr::case_when(
        .data$method_type == "unadjusted" ~ "Unadjusted",
        .data$method_type == "adjusted"   ~ "Adjusted",
        TRUE                              ~ .data$method_type
      ),
      LCL = .data$Estimate - crit * .data$SE,
      UCL = .data$Estimate + crit * .data$SE,
      Type = factor(.data$Type, levels = c("Unadjusted","Adjusted"))
    )

  ylab <- if (x$meta$family == "gaussian") {
    "Effect (mean difference)"
  } else {
    if (x$meta$scale == "RD") "Effect (risk difference)"
    else if (x$meta$scale == "RR") "Effect (log risk ratio)"
    else "Effect (log odds ratio)"
  }

  p <- ggplot2::ggplot(plotdf, ggplot2::aes(x = .data$Type, y = .data$Estimate)) +
    ggplot2::geom_hline(yintercept = 0, linewidth = 0.4, linetype = 2, alpha = 0.6) +
    ggplot2::geom_pointrange(ggplot2::aes(ymin = .data$LCL, ymax = .data$UCL), linewidth = 0.5) +
    ggplot2::geom_point(size = point_size) +
    ggplot2::facet_wrap(~ Estimand, scales = "free_y", nrow = 1) +
    ggplot2::labs(
      x = NULL, y = ylab,
      title = "Model-robust standardization (Unadjusted vs Adjusted)",
      subtitle = sprintf("Confidence level: %.0f%%; df = %d; kept periods: %d/%d",
                         100*level, df, x$meta$n_periods_kept, x$meta$n_periods_total)
    ) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      strip.text = ggplot2::element_text(face = "bold"),
      plot.title = ggplot2::element_text(face = "bold")
    )

  print(p)
  invisible(p)
}
