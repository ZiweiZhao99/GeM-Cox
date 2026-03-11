###############################################################
## GeM-Cox v3: Key revisions to GeMCox_improved_v2.R
##
## This file sources v2 and adds/overrides:
##   1. preprocess_features() — log-transform + standardize
##   2. shared_breslow() — shared-baseline Breslow estimator
##   3. cox_poisson_weighted() — Poisson regression M-step for shared baseline
##   4. cv_select_K_gamma_v3() — with built-in 1-SE rule
##   5. gemcox_pipeline() — end-to-end wrapper
##
## Usage: source("GeMCox_improved_v3.R")
##        (automatically sources GeMCox_improved_v2.R)
###############################################################

source("GeMCox_improved_v2.R")   # base implementation

###############################################################
## 1. Preprocessing utilities
###############################################################

#' Log-transform skewed features
#' @param X  numeric matrix
#' @param skew_threshold  features with skewness > this get log(1+x)
#' @return list(X_transformed, log_cols = indices of transformed columns)
preprocess_log_transform <- function(X, skew_threshold = 1.0) {
  X <- as.matrix(X)
  skewness <- function(x) {
    x <- x[is.finite(x)]
    n <- length(x)
    if (n < 3) return(0)
    m <- mean(x); s <- sd(x)
    if (s < 1e-12) return(0)
    (sum((x - m)^3) / n) / (s^3)
  }
  
  skew_vals <- apply(X, 2, skewness)
  log_cols <- which(skew_vals > skew_threshold & apply(X, 2, function(x) all(x[is.finite(x)] >= 0)))
  
  X_out <- X
  if (length(log_cols) > 0) {
    X_out[, log_cols] <- log1p(pmax(X[, log_cols], 0))
  }
  
  list(X = X_out, log_cols = log_cols, skew_before = skew_vals)
}

#' Full preprocessing: log-transform + impute (median) + standardize
#' @param X_gmm  GMM feature matrix (n x q)
#' @param X_cox  Cox feature matrix (n x p), may include X_gmm columns + covariates
#' @param log_transform  whether to apply log-transform
#' @param skew_threshold  skewness threshold for log-transform
#' @return list with transformed matrices and scalers for test-time application
preprocess_train <- function(X_gmm, X_cox, log_transform = TRUE, skew_threshold = 1.0) {
  X_gmm <- as.matrix(X_gmm)
  X_cox <- as.matrix(X_cox)
  
  # Log-transform
  gmm_log <- if (log_transform) preprocess_log_transform(X_gmm, skew_threshold)
             else list(X = X_gmm, log_cols = integer(0))
  cox_log <- if (log_transform) preprocess_log_transform(X_cox, skew_threshold)
             else list(X = X_cox, log_cols = integer(0))
  
  # Median imputation
  impute_median <- function(X) {
    for (j in seq_len(ncol(X))) {
      na_idx <- which(is.na(X[, j]))
      if (length(na_idx) > 0) {
        X[na_idx, j] <- median(X[, j], na.rm = TRUE)
      }
    }
    X
  }
  X_gmm_imp <- impute_median(gmm_log$X)
  X_cox_imp <- impute_median(cox_log$X)
  
  # Standardize
  gmm_scaler <- make_x_scaler(X_gmm_imp)
  cox_scaler <- make_x_scaler(X_cox_imp)
  
  list(
    X_gmm = apply_x_scaler(X_gmm_imp, gmm_scaler),
    X_cox = X_cox_imp,  # Cox standardization handled inside glmnet
    gmm_scaler = gmm_scaler,
    cox_scaler = cox_scaler,
    gmm_log_cols = gmm_log$log_cols,
    cox_log_cols = cox_log$log_cols
  )
}

#' Apply training preprocessing to test data
preprocess_test <- function(X_gmm_new, X_cox_new, train_prep) {
  X_gmm_new <- as.matrix(X_gmm_new)
  X_cox_new <- as.matrix(X_cox_new)
  
  # Log-transform same columns
  if (length(train_prep$gmm_log_cols) > 0)
    X_gmm_new[, train_prep$gmm_log_cols] <- log1p(pmax(X_gmm_new[, train_prep$gmm_log_cols], 0))
  if (length(train_prep$cox_log_cols) > 0)
    X_cox_new[, train_prep$cox_log_cols] <- log1p(pmax(X_cox_new[, train_prep$cox_log_cols], 0))
  
  # Median impute with training medians (approximate — proper MI would be better)
  for (j in seq_len(ncol(X_gmm_new))) {
    na_idx <- which(is.na(X_gmm_new[, j]))
    if (length(na_idx) > 0) X_gmm_new[na_idx, j] <- train_prep$gmm_scaler$center[j]
  }
  for (j in seq_len(ncol(X_cox_new))) {
    na_idx <- which(is.na(X_cox_new[, j]))
    if (length(na_idx) > 0) X_cox_new[na_idx, j] <- train_prep$cox_scaler$center[j]
  }
  
  list(
    X_gmm = apply_x_scaler(X_gmm_new, train_prep$gmm_scaler),
    X_cox = X_cox_new
  )
}


###############################################################
## 2. Shared-baseline Breslow estimator
###############################################################

#' Compute the shared Breslow baseline hazard across all clusters
#' @param time     observed times
#' @param status   event indicators
#' @param eta_list list of length C, each an n-vector of linear predictors
#' @param tau      n x C matrix of cluster weights
#' @return list(time, jump, cumhaz) — same format as breslow_weighted
shared_breslow <- function(time, status, eta_list, tau,
                           denom_floor = 1e-8, max_jump = 10, max_cumhaz = 500) {
  n <- length(time)
  C <- ncol(tau)
  status <- as.integer(status)
  
  uniqT <- sort(unique(time[status == 1]))
  if (length(uniqT) == 0) {
    return(list(time = numeric(0), jump = numeric(0), cumhaz = numeric(0)))
  }
  
  # Precompute exp(eta) for each cluster
  exp_eta_list <- lapply(eta_list, function(eta) exp(clamp(eta, -30, 30)))
  
  jump <- numeric(length(uniqT))
  for (m in seq_along(uniqT)) {
    t_m <- uniqT[m]
    # Numerator: unweighted event count (since sum_c tau_ic = 1)
    num <- sum(status[time == t_m])
    # Denominator: sum over clusters of weighted risk contributions
    at_risk <- (time >= t_m)
    denom <- 0
    for (c in seq_len(C)) {
      denom <- denom + sum(tau[at_risk, c] * exp_eta_list[[c]][at_risk])
    }
    denom <- max(denom, denom_floor)
    jump[m] <- min(num / denom, max_jump)
  }
  
  cumhaz <- pmin(cumsum(jump), max_cumhaz)
  list(time = uniqT, jump = jump, cumhaz = cumhaz)
}


###############################################################
## 3. Poisson regression M-step for shared baseline
###############################################################

#' Fit weighted Poisson regression for cluster c with shared baseline
#' @param time     observed times
#' @param status   event indicators
#' @param X        feature matrix (n x p)
#' @param w        cluster weights tau[, c]
#' @param baseline shared Breslow baseline (from shared_breslow)
#' @param lambda   penalty parameter
#' @param alpha    elastic net mixing (0.5 recommended)
#' @return list with beta_s, scaler, etc.
cox_poisson_weighted <- function(time, status, X, w, baseline,
                                 lambda = 0.1, alpha = 0.5,
                                 min_eff_events = 0.5,
                                 fallback_fit = NULL) {
  X <- as.matrix(X)
  status <- as.integer(status)
  w <- as.numeric(w)
  w[!is.finite(w) | w < 0] <- 0
  
  scaler <- make_x_scaler(X)
  Xs <- apply_x_scaler(X, scaler)
  
  n_eff_event <- sum(w * status)
  if (!is.finite(n_eff_event) || n_eff_event < min_eff_events) {
    if (!is.null(fallback_fit)) {
      out <- fallback_fit; out$ok <- FALSE; out$reason <- "fallback_low_events"
      return(out)
    }
    return(list(beta_s = rep(0, ncol(X)), beta_orig = rep(0, ncol(X)),
                scaler = scaler, baseline = baseline,
                lambda = NA, ok = FALSE, reason = "low_events",
                eff_events = n_eff_event))
  }
  
  # Compute offset = log(Lambda_0(t_i))
  Lambda_at_t <- get_cumhaz_at(baseline, time)
  Lambda_at_t <- pmax(Lambda_at_t, 1e-10)
  offset <- log(Lambda_at_t)
  
  fit <- tryCatch({
    gfit <- glmnet(
      x = Xs,
      y = status,   # response is delta_i
      family = "poisson",
      weights = w,
      offset = offset,
      alpha = alpha,
      lambda = lambda,
      standardize = FALSE,
      maxit = 1e5
    )
    beta_s <- as.numeric(coef(gfit, s = lambda))[-1]  # drop intercept
    list(beta_s = beta_s, lambda = lambda, ok = TRUE, reason = "ok")
  }, error = function(e) NULL)
  
  if (is.null(fit)) {
    if (!is.null(fallback_fit)) {
      out <- fallback_fit; out$ok <- FALSE; out$reason <- "fallback_poisson_failed"
      return(out)
    }
    fit <- list(beta_s = rep(0, ncol(X)), lambda = lambda,
                ok = FALSE, reason = "poisson_failed")
  }
  
  beta_s <- fit$beta_s
  beta_s[!is.finite(beta_s)] <- 0
  
  list(
    beta_s = beta_s,
    beta_orig = beta_s / scaler$scale,
    scaler = scaler,
    baseline = baseline,  # shared baseline stored in each coxfit for interface compatibility
    lambda = fit$lambda,
    ok = isTRUE(fit$ok),
    reason = fit$reason,
    eff_events = n_eff_event
  )
}


###############################################################
## 4. CV with built-in 1-SE rule
###############################################################

#' Cross-validate (K, gamma) with 1-SE rule
#' @param ... same arguments as cv_select_K_gamma
#' @param use_1se_rule  if TRUE, apply 1-SE rule for parsimony
#' @return same as cv_select_K_gamma but $best uses 1-SE rule
cv_select_K_gamma_v3 <- function(..., use_1se_rule = TRUE) {
  cv_res <- cv_select_K_gamma(...)
  
  if (!use_1se_rule || is.null(cv_res$summary)) return(cv_res)
  
  summary <- cv_res$summary
  ok_rows <- which(is.finite(summary$mean_score))
  if (length(ok_rows) == 0) return(cv_res)
  
  best_idx <- ok_rows[which.max(summary$mean_score[ok_rows])]
  best_score <- summary$mean_score[best_idx]
  best_se <- summary$sd_score[best_idx]
  if (!is.finite(best_se)) best_se <- 0
  
  threshold <- best_score - best_se
  eligible <- summary[ok_rows, , drop = FALSE]
  eligible <- eligible[eligible$mean_score >= threshold, , drop = FALSE]
  eligible <- eligible[order(eligible$K, eligible$gamma), , drop = FALSE]
  
  if (nrow(eligible) > 0) {
    cv_res$best <- eligible[1, , drop = FALSE]
    cv_res$best_method <- "1se_rule"
  }
  
  cv_res
}


###############################################################
## 5. End-to-end pipeline wrapper
###############################################################

#' Full GeM-Cox pipeline: preprocess -> CV select -> fit -> predict
#'
#' @param X_gmm        GMM features (n x q), may contain NAs
#' @param X_cox        Cox features (n x p), may contain NAs
#' @param time         observed times
#' @param status       event indicators
#' @param K_grid       candidate K values (default 1:3)
#' @param gamma_grid   candidate gamma values
#' @param alpha        elastic net mixing (default 0.5)
#' @param log_transform  apply log(1+x) to skewed features
#' @param use_1se_rule use parsimony rule for K selection
#' @param nfolds       CV folds (default 5)
#' @param n_starts     EM restarts for CV (default 5)
#' @param n_starts_final EM restarts for final fit (default 10)
#' @param max_iter     max EM iterations
#' @param verbose      print progress
#'
#' @return list with fit, cv_result, preprocessing info
gemcox_pipeline <- function(
    X_gmm, X_cox, time, status,
    K_grid       = 1:3,
    gamma_grid   = c(0, 0.25, 0.5, 1.0, 2.0),
    alpha        = 0.5,
    log_transform = TRUE,
    use_1se_rule = TRUE,
    nfolds       = 5,
    n_starts     = 5,
    n_starts_final = 10,
    max_iter     = 50,
    seed         = 1,
    verbose      = TRUE
) {
  # Preprocess
  if (verbose) cat("Preprocessing...\n")
  prep <- preprocess_train(X_gmm, X_cox,
                           log_transform = log_transform)
  
  # CV model selection
  if (verbose) cat("Cross-validating (K, gamma)...\n")
  cv_res <- cv_select_K_gamma_v3(
    X_gmm  = prep$X_gmm,
    X_cox  = prep$X_cox,
    time   = time,
    status = status,
    K_grid = K_grid,
    gamma_grid = gamma_grid,
    nfolds = nfolds,
    alpha  = alpha,
    max_iter = max_iter,
    n_starts = n_starts,
    seed   = seed,
    verbose = verbose,
    use_1se_rule = use_1se_rule
  )
  
  if (is.null(cv_res$best)) stop("CV failed to select any valid (K, gamma).")
  
  K_best     <- cv_res$best$K[1]
  gamma_best <- cv_res$best$gamma[1]
  if (verbose) cat(sprintf("Selected K=%d, gamma=%.2f\n", K_best, gamma_best))
  
  # Final fit
  if (verbose) cat("Fitting final model...\n")
  lambda_final <- select_lambda(time, status, prep$X_cox, alpha = alpha)
  
  final_fit <- gemcox_full_multistart(
    X_gmm  = prep$X_gmm,
    X_cox  = prep$X_cox,
    time   = time,
    status = status,
    K      = K_best,
    lambda = lambda_final,
    alpha  = alpha,
    max_iter = max_iter * 2,
    surv_weight = gamma_best,
    verbose = verbose,
    n_starts = n_starts_final
  )
  
  list(
    fit        = final_fit,
    cv_result  = cv_res,
    K_selected = K_best,
    gamma_selected = gamma_best,
    preprocessing = prep
  )
}

cat("GeMCox_improved_v3.R loaded.\n")
cat("  New functions: preprocess_log_transform, preprocess_train, preprocess_test\n")
cat("  New functions: shared_breslow, cox_poisson_weighted\n")
cat("  New functions: cv_select_K_gamma_v3, gemcox_pipeline\n")
