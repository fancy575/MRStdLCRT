#' Model-Robust Standardization for Crossover, and Stepped-wedge cluster-randomized trial
#'
#' @description
#' Core fitting utilities to compute horizontal individual ATE (h-iATE),
#' horizontal cluster ATE (h-cATE), vertical individual ATE (v-iATE),
#' and vertical cluster-period ATE (v-cATE) in multi-period crossover CRTs.
#' Provides unadjusted & augmented (model-robust) estimators and delete-1-cluster
#' jackknife SEs. If \code{SW=TRUE}, the first and last periods are excluded from
#' model fitting, unadjusted means, augmentation, and treatment effect coefficients. The design
#' matrix columns for period (\code{pK}) and period-by-treatment (\code{pzK}) keep
#' the original period numbering and simply omit the excluded periods.
#'
#' @importFrom dplyr %>% group_by summarise mutate left_join arrange n distinct pull
#' @importFrom dplyr ungroup select across bind_cols bind_rows count filter first everything rename
#' @importFrom tidyr pivot_longer pivot_wider crossing
#' @importFrom tidyselect all_of any_of starts_with
#' @importFrom rlang sym
#' @importFrom lme4 lmer glmer fixef VarCorr
#' @importFrom gee gee
#' @importFrom stats model.frame model.matrix terms delete.response as.formula
#' @importFrom stats gaussian binomial qt
#' @importFrom utils head
#' @keywords internal
#' @name utilities

utils::globalVariables(c(
  ".", "Yij", "Z", "Nij", "N_i", "N_j",
  "wij_hi","wij_hc","wij_vi","wij_vc",
  "wj_hi","wj_hc","wj_vi","wj_vc",
  "w_j","mz","pot_z","num_term1","den_term1","sum_mz",
  "scheme","w_scheme_cp",".per_idx"
))

.inv_logit <- function(eta) 1/(1+exp(-eta))
.sdiv      <- function(a,b) ifelse(b==0, NA_real_, a/b)

# --------- Build pK and pzK (all levels, original numbering) ---------
add_p_and_pz <- function(df, period = "period", trt = "trt") {
  per_levels <- sort(unique(df[[period]]))
  per_fac    <- factor(df[[period]], levels = per_levels)
  J <- nlevels(per_fac)
  # p1..pJ
  P <- model.matrix(~ 0 + per_fac)
  colnames(P) <- paste0("p", seq_len(J))
  # pz1..pzJ
  tr <- as.integer(df[[trt]] == 1L)
  PZ <- sapply(seq_len(J), function(k) as.integer(per_fac == levels(per_fac)[k]) * tr)
  colnames(PZ) <- paste0("pz", seq_len(J))
  bind_cols(df, as.data.frame(P), as.data.frame(PZ))
}

# --------- Ensure pzK matches current trt and pK in a data.frame ----------
.sync_pz <- function(df, trt = "trt") {
  p_cols  <- grep("^p\\d+$",  names(df), value = TRUE)
  pz_cols <- grep("^pz\\d+$", names(df), value = TRUE)
  if (!length(p_cols) || !length(pz_cols)) return(df)
  k_p  <- sub("^p",  "", p_cols)
  k_pz <- sub("^pz", "", pz_cols)
  common <- intersect(k_p, k_pz)
  if (!length(common)) return(df)
  for (k in common) {
    pk  <- paste0("p",  k)
    pzk <- paste0("pz", k)
    df[[pzk]] <- as.integer(df[[pk]] == 1L) * as.integer(df[[trt]] == 1L)
  }
  df
}

.print_nice_table <- function(x, title) {
  cat(title, "\n")
  cat(rep("-", nchar(title)), "\n", sep = "")
  print(utils::head(x, n = max(6, nrow(x))))
  if (nrow(x) > 50) cat(sprintf("\n[%d rows total]\n", nrow(x)))
  cat("\n")
}

# --------- Rewrite a user formula to use pK / pzK explicitly ----------
rewrite_formula_to_p <- function(formula, data, period = "period", trt = "trt") {
  txt <- paste(deparse(formula), collapse = "")
  # Determine J
  per_levels <- sort(unique(data[[period]])); J <- length(per_levels)
  p_terms  <- paste0("p",  seq_len(J), collapse = " + ")
  pz_terms <- paste0("pz", seq_len(J), collapse = " + ")

  # Replace factor(period)*trt  ->  0 + p1+...+pJ + pz1+...+pzJ
  txt <- gsub(paste0("factor\\(",period,"\\)\\s*\\*\\s*", trt),
              paste0("0 + ", p_terms, " + ", pz_terms), txt)
  txt <- gsub(paste0(trt,"\\s*\\*\\s*factor\\(",period,"\\)"),
              paste0("0 + ", p_terms, " + ", pz_terms), txt)

  # Replace factor(period):trt  ->  0 + p1+...+pJ + pz1+...+pzJ
  txt <- gsub(paste0("factor\\(",period,"\\)\\s*:\\s*", trt),
              paste0("0 + ", p_terms, " + ", pz_terms), txt)
  txt <- gsub(paste0(trt,"\\s*:\\s*factor\\(",period,"\\)"),
              paste0("0 + ", p_terms, " + ", pz_terms), txt)

  # Replace lone factor(period) with 0 + p1+...+pJ
  txt <- gsub(paste0("factor\\(",period,"\\)"),
              paste0("0 + ", p_terms), txt)

  # Remove any accidental ' + 0 +' duplication
  txt <- gsub("\\+\\s*0\\s*\\+", "+", txt)
  as.formula(txt, env = environment(formula))
}

# Remove (1|...) terms for GEE (if user passed them)
strip_re_from_formula <- function(fml) {
  f_txt <- paste(deparse(fml), collapse = " ")
  f_txt <- gsub("\\(\\s*1\\s*\\|[^\\)]*\\)", "", f_txt) # remove (1|group)
  f_txt <- gsub("\\+\\s*\\+", "+", f_txt)
  f_txt <- gsub("\\s*\\+\\s*$", "", f_txt)
  f_txt <- gsub("\\s{2,}", " ", f_txt)
  f_txt <- trimws(f_txt)
  stats::as.formula(f_txt, env = parent.frame())
}

# Build model matrix columns aligned with a fitted object
mm_fixed <- function(fit, newdata, method) {
  if (method == "gee") {
    TT <- delete.response(terms(formula(fit), data = newdata))
    mf <- model.frame(TT, data = newdata)
    X  <- model.matrix(TT, mf)
    keep <- intersect(colnames(X), names(coef(fit)))
    X[, keep, drop = FALSE]
  } else {
    TT <- delete.response(terms(fit))
    mf <- model.frame(TT, data = newdata)
    X  <- model.matrix(TT, mf)
    keep <- intersect(colnames(X), names(lme4::fixef(fit)))
    X[, keep, drop = FALSE]
  }
}

# Predict m0/m1 at CP means (gaussian) or individual-level (binomial)
pred_cp_means <- function(fit, method, family, cp_base, trt, period) {
  rebuild_pz <- function(df, z) {
    df2 <- df
    df2[[trt]] <- z
    df2 <- .sync_pz(df2, trt = trt)
    df2
  }
  new0 <- rebuild_pz(cp_base, 0L)
  new1 <- rebuild_pz(cp_base, 1L)

  X0 <- mm_fixed(fit, new0, method); X1 <- mm_fixed(fit, new1, method)

  if (method == "gee") {
    b <- stats::coef(fit)
    eta0 <- drop(X0 %*% b[colnames(X0)]); eta1 <- drop(X1 %*% b[colnames(X1)])
    if (family=="gaussian") { m0 <- eta0; m1 <- eta1 } else { m0 <- .inv_logit(eta0); m1 <- .inv_logit(eta1) }
  } else if (inherits(fit,"lmerMod")) {
    b <- lme4::fixef(fit)
    m0 <- drop(X0 %*% b[colnames(X0)]); m1 <- drop(X1 %*% b[colnames(X1)])
  } else {
    b <- lme4::fixef(fit)
    eta0 <- drop(X0 %*% b[colnames(X0)]); eta1 <- drop(X1 %*% b[colnames(X1)])
    vc    <- lme4::VarCorr(fit)
    sigma <- sum(sapply(vc, function(x) attr(x, "stddev")^2))
    scale <- sqrt((sigma + pi^2/3) / (pi^2/3))
    m0 <- .inv_logit(eta0/scale); m1 <- .inv_logit(eta1/scale)
  }
  list(m0 = m0, m1 = m1)
}

# ---------- Core: compute ATEs with SW only in aggregation ----------
MRS_ATEs <- function(
    data, formula,
    cluster_id = "cluster",
    period     = "period",
    trt        = "trt",
    method     = c("gee","lmer","glmer"),
    family     = c("gaussian","binomial"),
    corstr     = "independence",
    SW         = FALSE
){
  method <- match.arg(method); family <- match.arg(family)
  if (method == "lmer"  && family != "gaussian") stop("lmer requires gaussian (use glmer/gee for binomial).")
  if (method == "glmer" && family != "binomial") stop("glmer requires binomial.")

  # Detect "treatment-only" working model (no covariates, no period terms)
  rhs_terms <- attr(terms(formula), "term.labels")
  rhs_terms <- gsub("`", "", rhs_terms)
  is_trt_only <- (length(rhs_terms) == 1L && identical(rhs_terms, trt))

  # Add p/pz to FULL data (safe even if not used)
  df <- add_p_and_pz(data, period, trt)
  per_all <- sort(unique(df[[period]]))
  if (length(per_all) < 2) stop("Need at least 2 periods.")

  # keep for aggregation (SW affects which periods are averaged)
  keep_levels <- if (SW && length(per_all) >= 3) per_all[2:(length(per_all)-1)] else per_all

  # Rewrite formula to p/pz form (full set); model pruning happens below
  f_fit <- rewrite_formula_to_p(formula, df, period, trt)

  # SW-aware pruning of edge-period regressors from the WORKING MODEL
  if (SW && length(per_all) >= 3) {
    J <- length(per_all)
    drop_terms <- c(paste0("p",  c(1, J)), paste0("pz", c(1, J)))
    drop_terms <- intersect(drop_terms, all.vars(delete.response(terms(f_fit))))
    if (length(drop_terms)) {
      f_fit <- stats::update(f_fit, paste(". ~ . -", paste(drop_terms, collapse = " - ")))
    }
  }

  # Prepare sizes & weights on FULL data (but N_i may be SW-aware)
  cl_sym <- rlang::sym(cluster_id); per_sym <- rlang::sym(period)
  cp_sizes <- df %>% dplyr::count(!!cl_sym, !!per_sym, name = "Nij")

  # Counts used in aggregation (respect SW)
  counts_period <- cp_sizes %>%
    dplyr::filter(.data[[period]] %in% keep_levels) %>%
    dplyr::group_by(!!per_sym) %>%
    dplyr::summarise(
      n_obs       = sum(Nij),
      n_clusters  = dplyr::n_distinct(.data[[cluster_id]]),
      .groups     = "drop"
    ) %>%
    dplyr::arrange(!!per_sym)

  counts_cluster <- cp_sizes %>%
    dplyr::filter(.data[[period]] %in% keep_levels) %>%
    dplyr::group_by(!!cl_sym) %>%
    dplyr::summarise(
      n_obs      = sum(Nij),
      n_periods  = dplyr::n_distinct(.data[[period]]),
      .groups    = "drop"
    ) %>%
    dplyr::arrange(!!cl_sym)

  counts_cp <- cp_sizes %>%
    dplyr::filter(.data[[period]] %in% keep_levels) %>%
    tidyr::pivot_wider(
      names_from  = tidyselect::all_of(period),
      values_from = Nij,
      values_fill = 0
    ) %>%
    dplyr::arrange(!!cl_sym)

  per_levels <- sort(unique(cp_sizes[[period]]))
  per_map    <- setNames(seq_along(per_levels), per_levels)

  # SW-aware N_i for h-cATE weights (sum CP sizes over kept periods only)
  N_i <- cp_sizes %>%
    { if (SW && length(per_all) >= 3) dplyr::filter(., .data[[period]] %in% keep_levels) else . } %>%
    dplyr::group_by(!!cl_sym) %>%
    dplyr::summarise(N_i = sum(Nij), .groups="drop")

  # Per-period totals (used by v-schemes)
  N_j <- cp_sizes %>% dplyr::group_by(!!per_sym) %>%
    dplyr::summarise(N_j = sum(Nij), .groups="drop")

  df <- df %>%
    dplyr::left_join(cp_sizes, by = setNames(c(cluster_id, period), c(cluster_id, period))) %>%
    dplyr::left_join(N_i,       by = setNames(cluster_id, cluster_id)) %>%
    dplyr::left_join(N_j,       by = setNames(period, period)) %>%
    dplyr::mutate(
      wijk_hi = 1,        # h-iATE
      wijk_hc = 1/N_i,    # h-cATE (SW-aware)
      wijk_vi = 1/N_j,    # v-iATE
      wijk_vc = 1/Nij     # v-cATE (CP)
    )

  # Unadjusted means (aggregation possibly SW-filtered)
  y_name <- all.vars(update(f_fit, . ~ 0))[1]

  cp_unadj <- df %>%
    dplyr::group_by(!!cl_sym, !!per_sym) %>%
    dplyr::summarise(
      Yij    = mean(.data[[y_name]]),
      Z      = dplyr::first(.data[[trt]]),
      wij_hi = sum(wijk_hi), wij_hc = sum(wijk_hc),
      wij_vi = sum(wijk_vi), wij_vc = sum(wijk_vc),
      .groups = "drop"
    ) %>%
    dplyr::mutate(.per_idx = per_map[as.character(.data[[period]])])

  p_w <- cp_unadj %>%
    dplyr::group_by(!!per_sym) %>%
    dplyr::summarise(
      wj_hi = sum(wij_hi), wj_hc = sum(wij_hc),
      wj_vi = sum(wij_vi), wj_vc = sum(wij_vc),
      .groups="drop"
    ) %>%
    dplyr::mutate(.per_idx = per_map[as.character(.data[[period]])])

  cp_z <- cp_unadj %>%
    dplyr::group_by(!!per_sym, Z) %>%
    dplyr::summarise(
      num_hi = sum(wij_hi * Yij), den_hi = sum(wij_hi),
      num_hc = sum(wij_hc * Yij), den_hc = sum(wij_hc),
      num_vi = sum(wij_vi * Yij), den_vi = sum(wij_vi),
      num_vc = sum(wij_vc * Yij), den_vc = sum(wij_vc),
      .groups="drop"
    ) %>%
    dplyr::left_join(p_w, by = setNames(period, period)) %>%
    dplyr::mutate(
      Ybar_j_hi = .sdiv(num_hi, den_hi),
      Ybar_j_hc = .sdiv(num_hc, den_hc),
      Ybar_j_vi = .sdiv(num_vi, den_vi),
      Ybar_j_vc = .sdiv(num_vc, den_vc)
    )

  mu_unadj <- cp_z %>% dplyr::filter(.data[[period]] %in% keep_levels) %>%
    dplyr::group_by(Z) %>%
    dplyr::summarise(
      mu_hi_unadj = sum(wj_hi * Ybar_j_hi)/sum(wj_hi),
      mu_hc_unadj = sum(wj_hc * Ybar_j_hc)/sum(wj_hc),
      mu_vi_unadj = sum(wj_vi * Ybar_j_vi)/sum(wj_vi),
      mu_vc_unadj = sum(wj_vc * Ybar_j_vc)/sum(wj_vc),
      .groups="drop"
    ) %>% dplyr::arrange(Z)

  if (family=="gaussian") {
    ATEs_unadj <- mu_unadj %>% dplyr::summarise(
      `h-iATE` = diff(mu_hi_unadj),
      `h-cATE` = diff(mu_hc_unadj),
      `v-iATE` = diff(mu_vi_unadj),
      `v-cATE` = diff(mu_vc_unadj)
    )
  } else {
    logit <- function(p) log(p/(1-p))
    ATEs_unadj <- mu_unadj %>% dplyr::summarise(
      `h-iATE` = logit(mu_hi_unadj[2]) - logit(mu_hi_unadj[1]),
      `h-cATE` = logit(mu_hc_unadj[2]) - logit(mu_hc_unadj[1]),
      `v-iATE` = logit(mu_vi_unadj[2]) - logit(mu_vi_unadj[1]),
      `v-cATE` = logit(mu_vc_unadj[2]) - logit(mu_vc_unadj[1])
    )
  }

  # Fit model EVEN IF treatment-only, so Coef/SE are available
  fit <- switch(method,
                "gee"  = {
                  fam <- if (family=="gaussian") gaussian() else binomial(link="logit")
                  f_use <- strip_re_from_formula(f_fit)
                  fit_gee <- NULL
                  capture.output({
                    fit_gee <- gee::gee(formula = f_use, id = df[[cluster_id]], data = df,
                                        family = fam, corstr = corstr, na.action = na.omit)
                  }, type = "output")
                  fit_gee
                },
                "lmer"  = lme4::lmer(f_fit, data = df),
                "glmer" = lme4::glmer(f_fit, data = df, family = binomial(link="logit"))
  )

  # Generate predictions for augmentation
  if (family == "binomial") {
    nd0 <- df; nd0[[trt]] <- 0L; nd0 <- .sync_pz(nd0, trt = trt)
    nd1 <- df; nd1[[trt]] <- 1L; nd1 <- .sync_pz(nd1, trt = trt)
    X0 <- mm_fixed(fit, nd0, method); X1 <- mm_fixed(fit, nd1, method)
    if (method == "gee") {
      b <- coef(fit)
      eta0 <- drop(X0 %*% b[colnames(X0)]); eta1 <- drop(X1 %*% b[colnames(X1)])
      mu0_s <- .inv_logit(eta0); mu1_s <- .inv_logit(eta1)
    } else {
      b <- lme4::fixef(fit)
      eta0 <- drop(X0 %*% b[colnames(X0)]); eta1 <- drop(X1 %*% b[colnames(X1)])
      vc    <- lme4::VarCorr(fit)
      sigma <- sum(sapply(vc, function(x) attr(x, "stddev")^2))
      scale <- sqrt((sigma + pi^2/3) / (pi^2/3))
      mu0_s <- .inv_logit(eta0/scale); mu1_s <- .inv_logit(eta1/scale)
    }

    data_sum <- df %>%
      dplyr::mutate(mu0_s = mu0_s, mu1_s = mu1_s) %>%
      dplyr::group_by(!!cl_sym, !!per_sym) %>%
      dplyr::summarise(mu0 = mean(mu0_s), mu1 = mean(mu1_s), .groups="drop")

    cp_pred <- df %>%
      dplyr::distinct(!!cl_sym, !!per_sym) %>%
      dplyr::left_join(data_sum, by = setNames(c(cluster_id, period), c(cluster_id, period))) %>%
      dplyr::mutate(m0 = mu0, m1 = mu1)

  } else {
    cp_proto <- df %>%
      dplyr::arrange(!!cl_sym, !!per_sym) %>% dplyr::group_by(!!cl_sym, !!per_sym) %>% dplyr::slice(1) %>% dplyr::ungroup()

    all_used <- all.vars(delete.response(terms(fit)))
    extras   <- grep("^p\\d+$|^pz\\d+$", names(df), value = TRUE)
    never_avg <- unique(c(y_name, trt, cluster_id, period, "Nij","N_i","N_j"))
    cand <- setdiff(unique(c(all_used, extras)), never_avg)
    cand <- intersect(cand, names(df))
    num_cols <- cand[sapply(df[cand], is.numeric)]

    cp_means <- if (length(num_cols)) {
      df %>% dplyr::group_by(!!cl_sym, !!per_sym) %>%
        dplyr::summarise(dplyr::across(tidyselect::all_of(num_cols), mean), .groups="drop")
    } else {
      df %>% dplyr::distinct(!!cl_sym, !!per_sym)
    }

    cp_base <- cp_proto %>%
      dplyr::select(-tidyselect::any_of(num_cols)) %>%
      dplyr::left_join(cp_means, by = setNames(c(cluster_id, period), c(cluster_id, period)))

    preds  <- pred_cp_means(fit, method, family, cp_base, trt, period)
    cp_pred <- cp_base %>% dplyr::mutate(m0 = preds$m0, m1 = preds$m1) %>%
      dplyr::filter(.data[[period]] %in% keep_levels)
  }

  # Attach weights & observed CP means for kept periods
  cp_unadj_min <- cp_unadj %>%
    dplyr::select(dplyr::all_of(c(cluster_id, period)), .per_idx, wij_hi, wij_hc, wij_vi, wij_vc, Yij, Z)

  cp_pred <- cp_pred %>%
    dplyr::left_join(cp_unadj_min, by = setNames(c(cluster_id, period), c(cluster_id, period)))

  # Long format for augmentation
  cp_long <- cp_pred %>%
    tidyr::crossing(pot_z = c(0,1)) %>%
    dplyr::mutate(mz = ifelse(pot_z==0, m0, m1)) %>%
    tidyr::pivot_longer(
      cols = c(wij_hi, wij_hc, wij_vi, wij_vc),
      names_to = "w_scheme_cp", values_to = "w_ij"
    ) %>%
    dplyr::mutate(
      scheme = dplyr::case_when(
        w_scheme_cp=="wij_hi" ~ "h-iATE",
        w_scheme_cp=="wij_hc" ~ "h-cATE",
        w_scheme_cp=="wij_vi" ~ "v-iATE",
        w_scheme_cp=="wij_vc" ~ "v-cATE"
      ),
      w_j = dplyr::case_when(
        scheme=="h-iATE" ~ p_w$wj_hi[match(.per_idx, p_w$.per_idx)],
        scheme=="h-cATE" ~ p_w$wj_hc[match(.per_idx, p_w$.per_idx)],
        scheme=="v-iATE" ~ p_w$wj_vi[match(.per_idx, p_w$.per_idx)],
        scheme=="v-cATE" ~ p_w$wj_vc[match(.per_idx, p_w$.per_idx)]
      )
    )

  aug_data <- cp_long %>%
    dplyr::group_by(!!per_sym, .per_idx, scheme, pot_z) %>%
    dplyr::summarise(
      num_term1 = sum(w_ij * (Z==pot_z) * (Yij - mz)),
      den_term1 = sum(w_ij * (Z==pot_z)),
      sum_mz    = sum(w_ij * mz),
      w_j       = dplyr::first(w_j),
      .groups   = "drop"
    ) %>%
    dplyr::mutate(
      term1 = .sdiv(num_term1, den_term1),
      term2 = .sdiv(sum_mz, w_j),
      mu_aug = term1 + term2
    )

  final_agg <- aug_data %>% dplyr::filter(.data[[period]] %in% keep_levels) %>%
    dplyr::group_by(scheme, pot_z) %>%
    dplyr::summarise(
      numerator   = sum(w_j * mu_aug),
      denominator = sum(w_j),
      mu_aug_overall = numerator/denominator,
      .groups="drop"
    )

  mu_aug_table <- final_agg %>%
    dplyr::select(pot_z, scheme, mu_aug_overall) %>%
    tidyr::pivot_wider(names_from = scheme, values_from = mu_aug_overall) %>%
    dplyr::arrange(pot_z) %>% dplyr::rename(Z = pot_z)

  if (family=="gaussian") {
    ATEs_aug <- tibble::tibble(
      `h-iATE` = mu_aug_table$`h-iATE`[mu_aug_table$Z==1] - mu_aug_table$`h-iATE`[mu_aug_table$Z==0],
      `h-cATE` = mu_aug_table$`h-cATE`[mu_aug_table$Z==1] - mu_aug_table$`h-cATE`[mu_aug_table$Z==0],
      `v-iATE` = mu_aug_table$`v-iATE`[mu_aug_table$Z==1] - mu_aug_table$`v-iATE`[mu_aug_table$Z==0],
      `v-cATE` = mu_aug_table$`v-cATE`[mu_aug_table$Z==1] - mu_aug_table$`v-cATE`[mu_aug_table$Z==0]
    )
  } else {
    logit <- function(p) log(p/(1-p))
    ATEs_aug <- tibble::tibble(
      `h-iATE` = logit(mu_aug_table$`h-iATE`[mu_aug_table$Z==1]) - logit(mu_aug_table$`h-iATE`[mu_aug_table$Z==0]),
      `h-cATE` = logit(mu_aug_table$`h-cATE`[mu_aug_table$Z==1]) - logit(mu_aug_table$`h-cATE`[mu_aug_table$Z==0]),
      `v-iATE` = logit(mu_aug_table$`v-iATE`[mu_aug_table$Z==1]) - logit(mu_aug_table$`v-iATE`[mu_aug_table$Z==0]),
      `v-cATE` = logit(mu_aug_table$`v-cATE`[mu_aug_table$Z==1]) - logit(mu_aug_table$`v-cATE`[mu_aug_table$Z==0])
    )
  }

  # TE coefficients (weighted across kept periods)
  coef_vec <- if (method == "gee") coef(fit) else lme4::fixef(fit)
  vc_full  <- if (method == "gee") fit[["robust.variance"]] else as.matrix(vcov(fit))
  cnames   <- names(coef_vec)

  te_names_ps <- grep("^pz\\d+$", cnames, value = TRUE)
  if (!length(te_names_ps)) {
    te_names_ps <- grep(paste0("^factor\\(",period,"\\)\\d+:", trt), cnames, value = TRUE)
  }

  w_per <- p_w %>% dplyr::select(.per_idx, wj_hi, wj_hc, wj_vi, wj_vc)

  idx_from_name <- function(nm) as.integer(sub("^pz","", nm))
  if (length(te_names_ps)) {
    k_idx <- vapply(te_names_ps, idx_from_name, integer(1))
    weighted_TE <- function(wcol) {
      ww <- w_per[[wcol]]
      pos <- match(k_idx, w_per$.per_idx); ok <- !is.na(pos)
      if (!any(ok)) return(c(estimate=NA_real_, variance=NA_real_))
      wvec <- ww[pos[ok]]; wvec <- wvec / sum(wvec)
      keeps <- te_names_ps[ok]
      est <- sum(wvec * coef_vec[keeps])
      V   <- vc_full[keeps, keeps, drop=FALSE]
      var <- as.numeric(t(wvec) %*% V %*% wvec)
      c(estimate=est, variance=var)
    }
    te_hi <- weighted_TE("wj_hi")
    te_hc <- weighted_TE("wj_hc")
    te_vi <- weighted_TE("wj_vi")
    te_vc <- weighted_TE("wj_vc")
    te_table <- tibble::tibble(
      scheme = c("h-iATE","h-cATE","v-iATE","v-cATE"),
      te_estimate = c(te_hi["estimate"], te_hc["estimate"], te_vi["estimate"], te_vc["estimate"]),
      te_variance = c(te_hi["variance"], te_hc["variance"], te_vi["variance"], te_vc["variance"])
    )
  } else {
    trt_nm <- intersect(cnames, c(trt, paste0(trt,"1"), paste0("`",trt,"`"), "Z", "Zij"))
    if (!length(trt_nm)) {
      te_est <- rep(NA_real_, 4); te_var <- rep(NA_real_, 4)
    } else {
      te <- unname(coef_vec[trt_nm[1]]); va <- unname(vc_full[trt_nm[1], trt_nm[1]])
      te_est <- rep(te, 4); te_var <- rep(va, 4)
    }
    te_table <- tibble::tibble(
      scheme = c("h-iATE","h-cATE","v-iATE","v-cATE"),
      te_estimate = te_est, te_variance = te_var
    )
  }

  list(
    fit            = fit,
    mu_unadj       = mu_unadj,
    ATEs_unadj     = ATEs_unadj,
    mu_aug_table   = mu_aug_table,
    ATEs_aug       = ATEs_aug,
    te_table       = te_table,
    p_weights      = p_w,
    SW             = SW,
    kept_periods   = keep_levels,
    counts         = list(
      period  = counts_period,
      cluster = counts_cluster,
      cp      = counts_cp
    )
  )
}

.extract_ates_vec <- function(res, which = c("aug","unadj")) {
  which <- match.arg(which)
  tab <- if (which == "aug") res$ATEs_aug else res$ATEs_unadj
  v <- unlist(tab[1, ], use.names = TRUE)
  v[c("h-iATE","h-cATE","v-iATE","v-cATE")]
}

.jk_var_scalar <- function(theta_reps) {
  I <- length(theta_reps); m <- mean(theta_reps)
  ((I - 1) / I) * sum((theta_reps - m)^2)
}

#' Fit model-robust standardization for crossover CRTs
#'
#' Fits the model-robust standardization (MRS) estimators for 2x2 or multi-period
#' crossover cluster randomized trials (CRTs), producing unadjusted and augmented
#' estimates for four estimands: \emph{h-iATE}, \emph{h-cATE}, \emph{v-iATE}, and \emph{v-cATE}.
#' In the augmented path, a working model is fit using \code{lmer}, \code{glmer}, or \code{gee},
#' and predictions are standardized to obtain model-robust estimates. Uncertainty
#' is quantified using delete-1-cluster jackknife.
#'
#' @param data A \code{data.frame} containing the outcome, treatment, period, cluster,
#'   and covariates referenced by \code{formula}.
#' @param formula Analysis formula. May include \code{factor(period)} and, for mixed
#'   models, random intercepts such as \code{(1 | h)} or \code{(1 | h:p)}.
#'   For \code{gee}, any random-effects terms are internally removed.
#' @param cluster_id Character string naming the cluster ID column (default \code{"cluster"}).
#' @param period Character string naming the period column (default \code{"period"}).
#' @param trt Character string naming the treatment indicator column (default \code{"trt"}; 0/1 expected).
#' @param method One of \code{"gee"}, \code{"lmer"}, or \code{"glmer"} specifying the working model.
#' @param family Outcome family: \code{"gaussian"} or \code{"binomial"}. For \code{lmer},
#'   the family must be \code{"gaussian"}; for \code{glmer}, it must be \code{"binomial"}.
#' @param corstr Working correlation for \code{gee} (default \code{"independence"}).
#'   Ignored for \code{lmer}/\code{glmer}.
#' @param SW Logical; if \code{TRUE}, excludes the first and last periods from
#'   the \emph{aggregation} step (i.e., period weights and averaging). Model fitting
#'   and design matrices always use all observed periods.
#'
#' @details
#' \strong{Estimands.} Horizontal estimands treat clusters as the primary unit:
#' \code{"h-iATE"} (hospital-weighted individual ATE) and \code{"h-cATE"} (hospital-weighted
#' cluster-period ATE). Vertical estimands weight over periods: \code{"v-iATE"} and \code{"v-cATE"}.
#'
#' \strong{Binary outcomes.} For \code{family = "binomial"}, results are reported on the
#' log-odds-ratio scale. In the \code{glmer} path, a latent-scale adjustment is applied
#' using the fitted random-effect variance.
#'
#' \strong{GEE.} Any random-effects terms in \code{formula} are dropped automatically.
#'
#' @return An object of class \code{"mrs"} with components:
#' \describe{
#'   \item{estimates}{2-row data frame (Unadjusted, Adjusted) with columns
#'     \code{h-iATE}, \code{h-cATE}, \code{v-iATE}, \code{v-cATE}.}
#'   \item{jk_se}{2-row data frame of jackknife SEs for each estimand.}
#'   \item{jk_cov_unadj}{4x4 jackknife covariance matrix for unadjusted estimands.}
#'   \item{jk_cov_aug}{4x4 jackknife covariance matrix for augmented estimands.}
#'   \item{reps}{List of intermediate objects (means, augmented tables, etc.).}
#'   \item{meta}{List with call, method/family, correlation structure, sample sizes, and variable names.}
#'   \item{te_table}{Optional table of period-weighted treatment-effect coefficients and variances.}
#'   \item{p_weights}{Period weights used in aggregation.}
#' }
#'
#' @seealso \link{summary.mrs}, \link{plot.mrs}
#'
#' @examples
#' \dontrun{
#' ## Continuous outcome (lmer)
#' data(ex_c, package = "MRStdLCRT")
#' fml_c <- y~ 0 + factor(period)*trt +
#'         x1 + x2 + (1 | cluster) + (1 | cluster:period)
#' fit_c <- mrstdlcrt_fit(
#'   data = ex_c, formula = fml_c,
#'   cluster_id = "cluster", period = "period", trt = "trt",
#'   method = "lmer", family = "gaussian", SW = FALSE
#' )
#' summary(fit_c, estimand = c("h-iATE","h-cATE"))
#'
#' ## Binary outcome (gee)
#' data(ex_b, package = "MRStdLCRT")
#' fml_b <- y~ 0 + factor(period)*trt +
#'           x1 + x2 + (1 | cluster) + (1 | cluster:period)
#' fit_b <- mrstdlcrt_fit(
#'   data = ex_b, formula = fml_b,
#'   cluster_id = "cluster", period = "period", trt = "trt",
#'   method = "gee", family = "binomial",
#'   corstr = "exchangeable", SW = FALSE
#' )
#' summary(fit_b, estimand = c("h-iATE","h-cATE"))
#' }
#' @export
mrstdlcrt_fit <- function(
    data, formula,
    cluster_id = "cluster",
    period     = "period",
    trt        = "trt",
    method  = c("gee","lmer","glmer"),
    family  = c("gaussian","binomial"),
    corstr  = "independence",
    SW      = TRUE
){
  method <- match.arg(method); family <- match.arg(family)

  full <- MRS_ATEs(
    data, formula, cluster_id, period, trt, method, family, corstr, SW = SW
  )
  full_aug   <- .extract_ates_vec(full, "aug")    # length 4: h-iATE, h-cATE, v-iATE, v-cATE
  full_unadj <- .extract_ates_vec(full, "unadj")  # length 4

  cl_sym <- rlang::sym(cluster_id)
  cl_ids <- unique(dplyr::pull(data, !!cl_sym)); I <- length(cl_ids)

  # delete-1 cluster replicates (I x 4) for adjusted and unadjusted
  reps_aug   <- matrix(NA_real_, nrow = I, ncol = 4,
                       dimnames = list(NULL, c("h-iATE","h-cATE","v-iATE","v-cATE")))
  reps_unadj <- reps_aug

  for (g in seq_len(I)) {
    d_sub <- dplyr::filter(data, (!!cl_sym) != cl_ids[g])
    res_g <- MRS_ATEs(
      data = d_sub, formula = formula,
      cluster_id = cluster_id, period = period, trt = trt,
      method = method, family = family, corstr = corstr, SW = SW
    )
    reps_aug[g, ]   <- .extract_ates_vec(res_g, "aug")
    reps_unadj[g, ] <- .extract_ates_vec(res_g, "unadj")
  }

  # JK covariance matrices a la var_est()
  .jk_cov <- function(mat_Ix4) {
    theta_bar <- colMeans(mat_Ix4)  # 1x4
    D <- t(t(mat_Ix4) - theta_bar)  # I x 4 deviations
    D <- t(D)                       # 4 x I
    ((I - 1) / I) * (D %*% t(D))   # 4 x 4 covariance
  }

  jk_cov_unadj <- .jk_cov(reps_unadj)
  jk_cov_aug   <- .jk_cov(reps_aug)

  # SEs are sqrt of diagonal
  se_unadj <- sqrt(diag(jk_cov_unadj))
  se_aug   <- sqrt(diag(jk_cov_aug))

  # point estimates table
  estimates <- dplyr::tibble(
    method_type = c("unadjusted","adjusted"),
    `h-iATE`    = c(full_unadj["h-iATE"], full_aug["h-iATE"]),
    `h-cATE`    = c(full_unadj["h-cATE"], full_aug["h-cATE"]),
    `v-iATE`    = c(full_unadj["v-iATE"], full_aug["v-iATE"]),
    `v-cATE`    = c(full_unadj["v-cATE"], full_aug["v-cATE"])
  )

  # SE table (diagonals only)
  jk_se <- dplyr::tibble(
    method_type     = c("unadjusted","adjusted"),
    `SE(h-iATE)`    = c(se_unadj[1], se_aug[1]),
    `SE(h-cATE)`    = c(se_unadj[2], se_aug[2]),
    `SE(v-iATE)`    = c(se_unadj[3], se_aug[3]),
    `SE(v-cATE)`    = c(se_unadj[4], se_aug[4])
  )

  meta <- list(
    call    = match.call(),
    formula = formula,
    method  = method,
    family  = family,
    corstr  = if (method=="gee") corstr else NA_character_,
    n_clusters = I,
    n_periods  = length(unique(data[[period]])),
    n_subjects = nrow(data),
    cluster_id = cluster_id,
    period     = period,
    trt        = trt,
    SW         = SW
  )

  res <- list(
    estimates    = estimates,
    jk_se        = jk_se,
    jk_cov_unadj = jk_cov_unadj,
    jk_cov_aug   = jk_cov_aug,
    reps         = full,
    meta         = meta,
    te_table     = full$te_table,
    p_weights    = full$p_weights,
    counts       = full$counts
  )
  class(res) <- "mrs"
  res
}

#' Summarize an mrs fit
#'
#' Produces a compact table of point estimates, jackknife standard errors, and
#' t-based confidence intervals for the crossover CRT estimands in \code{object}.
#' Optionally filter to a subset of estimands. If available, also prints the
#' corresponding model coefficient ("Coef") and its variance for each estimand.
#'
#' @method summary mrs
#' @export
#'
#' @param object An object of class \code{"mrs"} returned by \link{mrstdlcrt_fit}.
#' @param level Confidence level for the CI (default \code{0.95}).
#' @param estimand Optional character vector specifying which estimands to show.
#'   Must be a subset of \code{c("h-iATE","h-cATE","v-iATE","v-cATE")}.
#'   If \code{NULL} (default), all four are shown.
#' @param ... Unused.
summary.mrs <- function(object, level = 0.95, estimand = NULL, ...) {
  stopifnot(inherits(object, "mrs"))
  est  <- object$estimates
  se   <- object$jk_se
  meta <- object$meta
  te   <- object$te_table

  valid_names <- c("h-iATE","h-cATE","v-iATE","v-cATE")
  sel <- if (is.null(estimand)) valid_names else intersect(estimand, valid_names)
  if (!length(sel)) sel <- valid_names

  df <- max(1, meta$n_clusters - 1L)
  alpha <- 1 - level
  crit  <- stats::qt(1 - alpha/2, df = df)

  add_ci <- function(val, se) {
    lo <- val - crit * se
    hi <- val + crit * se
    cbind(Estimate = val, SE = se, LCL = lo, UCL = hi)
  }

  cat("\nMRS-XO Summary",
      if (isTRUE(meta$SW)) " (SW: first & last periods dropped in aggregation only)" else " (all periods)",
      "\n", sep = "")
  cat(rep("=", 72), "\n", sep = "")
  cat(sprintf("Method: %s   (family: %s)\n", meta$method, meta$family))
  if (!is.na(meta$corstr)) cat(sprintf("GEE corstr: %s\n", meta$corstr))
  cat(sprintf("Clusters: %d   Periods: %d   N: %d\n",
              meta$n_clusters, meta$n_periods, meta$n_subjects))
  cat(sprintf("Cluster id: %s   Period: %s   Trt: %s\n",
              meta$cluster_id, meta$period, meta$trt))
  cat("Formula used for fitting: "); print(meta$formula)
  cat(sprintf("CIs: %.1f%% (t, df = %d)\n\n", 100*level, df))

  # Counts used in aggregation (respecting SW)
  if (!is.null(object$counts)) {
    kp <- if (!is.null(object$reps$kept_periods)) object$reps$kept_periods else NA
    sw_note <- if (isTRUE(meta$SW)) {
      sprintf(" (SW=TRUE; kept periods: %s)", paste(kp, collapse = ", "))
    } else {
      " (SW=FALSE; all periods kept)"
    }

    cat("\nCounts used to construct estimators", sw_note, "\n", sep = "")
    cat(rep("=", 72), "\n", sep = "")

    if (!is.null(object$counts$period) && nrow(object$counts$period)) {
      per_tab <- object$counts$period
      colnames(per_tab)[colnames(per_tab) == meta$period] <- "period"
      .print_nice_table(per_tab, "Per-period counts (n_obs, n_clusters)")
    }

    if (!is.null(object$counts$cluster) && nrow(object$counts$cluster)) {
      cl_tab <- object$counts$cluster
      colnames(cl_tab)[colnames(cl_tab) == meta$cluster_id] <- "cluster"
      .print_nice_table(cl_tab, "Per-cluster counts (n_obs, n_periods)")
    }

    if (!is.null(object$counts$cp) && nrow(object$counts$cp)) {
      .print_nice_table(object$counts$cp, "Cluster x Period cell sizes (N_ij)")
    }
  }

  # Map estimand -> column names in 'est' and 'se'
  col_est <- list(
    "h-iATE" = "h-iATE",
    "h-cATE" = "h-cATE",
    "v-iATE" = "v-iATE",
    "v-cATE" = "v-cATE"
  )
  col_se <- list(
    "h-iATE" = "SE(h-iATE)",
    "h-cATE" = "SE(h-cATE)",
    "v-iATE" = "SE(v-iATE)",
    "v-cATE" = "SE(v-cATE)"
  )

  print_block <- function(short_name) {
    ce <- col_est[[short_name]]
    cs <- col_se[[short_name]]
    tab <- as.data.frame(add_ci(est[[ce]], se[[cs]]))
    rownames(tab) <- est$method_type
    cat(short_name, "\n", sep = "")
    print(round(tab, 6)); cat("\n")
    if (!is.null(te)) {
      row <- te[match(short_name, te$scheme), , drop = FALSE]
      if (nrow(row)) {
        te_se <- sqrt(row$te_variance)
        cat("Coef (model):\n")
        print(round(data.frame(TE = row$te_estimate, SE = te_se, Var = row$te_variance), 6))
        cat("\n")
      }
    }
  }

  for (nm in sel) print_block(nm)

  invisible(list(
    estimates = est,
    jk_se     = se,
    te_table  = te,
    p_weights = object$p_weights,
    meta      = meta,
    level     = level,
    shown     = sel
  ))
}

#' Plot method for mrs objects
#'
#' Draws one panel per estimand (or selected subset), with three points:
#' Unadjusted, Adjusted, and the model Coef; all with CIs.
#'
#' @param x        An object of class "mrs" (from mrstdlcrt_fit()).
#' @param level    Confidence level (default 0.95) for the error bars.
#' @param estimand Character vector of estimands to show. Any of
#'                 c("h-iATE","h-cATE","v-iATE","v-cATE"). If NULL, show all.
#' @param point_size Point size for the markers (default 2.8).
#' @param ci_width   Width of the CI end caps (default 0.2).
#' @param ...      Ignored (for S3 compatibility).
#' @return A ggplot object (invisibly) and draws the plot.
#' @export
plot.mrs <- function(x,
                     level = 0.95,
                     estimand = NULL,
                     point_size = 2.8,
                     ci_width = 0.2,
                     ...) {
  if (!requireNamespace("ggplot2", quietly = TRUE))
    stop("Package 'ggplot2' is required for plotting.", call. = FALSE)
  if (!requireNamespace("dplyr", quietly = TRUE) ||
      !requireNamespace("tidyr", quietly = TRUE))
    stop("Packages 'dplyr' and 'tidyr' are required for plotting.", call. = FALSE)

  est <- x$estimates
  se  <- x$jk_se
  te  <- x$te_table

  valid <- c("h-iATE","h-cATE","v-iATE","v-cATE")
  sel <- if (is.null(estimand)) valid else intersect(estimand, valid)
  if (!length(sel)) sel <- valid

  df <- max(1, x$meta$n_clusters - 1L)
  alpha <- 1 - level
  crit  <- stats::qt(1 - alpha/2, df = df)

  est_long <- tidyr::pivot_longer(
    est,
    cols = tidyselect::all_of(valid),
    names_to = "Estimand",
    values_to = "Estimate"
  )

  se_long <- tidyr::pivot_longer(
    se,
    cols = tidyselect::starts_with("SE("),
    names_to = "Estimand",
    values_to = "SE"
  ) |> dplyr::mutate(Estimand = sub("^SE\\((.*)\\)$", "\\1", Estimand))

  base_long <- est_long |>
    dplyr::left_join(se_long, by = c("method_type","Estimand")) |>
    dplyr::mutate(Type = dplyr::case_when(
      method_type == "unadjusted" ~ "Unadjusted",
      method_type == "adjusted"   ~ "Adjusted",
      TRUE                        ~ method_type
    )) |>
    dplyr::select(Estimand, Type, Estimate, SE)

  coef_long <- NULL
  if (!is.null(te) && nrow(te)) {
    coef_long <- te |>
      dplyr::transmute(
        Estimand = scheme,
        Type     = "Coef",
        Estimate = te_estimate,
        SE       = sqrt(te_variance)
      )
  }

  plotdf <- dplyr::bind_rows(base_long, coef_long) |>
    dplyr::filter(Estimand %in% sel) |>
    dplyr::mutate(
      Type = factor(Type, levels = c("Unadjusted","Adjusted","Coef")),
      LCL = Estimate - crit * SE,
      UCL = Estimate + crit * SE
    )

  plotdf$LCL[!is.finite(plotdf$LCL)] <- NA_real_
  plotdf$UCL[!is.finite(plotdf$UCL)] <- NA_real_

  p <- ggplot2::ggplot(plotdf, ggplot2::aes(x = Type, y = Estimate)) +
    ggplot2::geom_hline(yintercept = 0, linewidth = 0.4, linetype = 2, alpha = 0.6) +
    ggplot2::geom_pointrange(
      ggplot2::aes(ymin = LCL, ymax = UCL),
      position = ggplot2::position_dodge(width = 0.0),
      size = point_size / 4,
      linewidth = 0.5
    ) +
    ggplot2::geom_point(size = point_size) +
    ggplot2::facet_wrap(~ Estimand, scales = "free_y", nrow = 1) +
    ggplot2::labs(
      x = NULL, y = "Estimate (with CI)",
      title = "Model-standardization: Unadjusted, Adjusted, and Model Coefficient",
      subtitle = sprintf("Confidence level: %.0f%%; df = %d", 100*level, df)
    ) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_text(margin = ggplot2::margin(t = 2)),
      strip.text = ggplot2::element_text(face = "bold"),
      plot.title = ggplot2::element_text(face = "bold")
    )

  print(p)
  invisible(p)
}
