###############################################################
## GeM-Cox v4 patch
##
## Purpose:
##   - Fix leakage from preprocessing-before-CV in v3
##   - Implement the shared-baseline EM fit recommended for small n
##   - Correct the 1-SE rule (use standard error, not SD)
##   - Add CVIA078 helpers that avoid PCA and default to
##     interpretable delta-only clustering features
##
## Notes:
##   - This file SOURCES GeMCox_improved_v2.R and overrides / adds
##     higher-level functions; it does not modify v2 on disk.
##   - Because this is a patch layer, utility functions from v2 are used
##     directly where possible.
###############################################################

source("GeMCox_improved_v2.R")

###############################################################
## 1. Fold-safe preprocessing
###############################################################

skewness_v4 <- function(x) {
  x <- x[is.finite(x)]
  n <- length(x)
  if (n < 3) return(0)
  m <- mean(x)
  s <- stats::sd(x)
  if (!is.finite(s) || s < 1e-12) return(0)
  mean((x - m)^3) / (s^3)
}

preprocess_train_v4 <- function(X_gmm, X_cox,
                                log_transform = TRUE,
                                skew_threshold = 1.0,
                                impute = c("median", "mean")) {
  impute <- match.arg(impute)
  X_gmm <- as.matrix(X_gmm)
  X_cox <- as.matrix(X_cox)
  
  colnames(X_gmm) <- make.names(colnames(X_gmm), unique = TRUE)
  colnames(X_cox) <- make.names(colnames(X_cox), unique = TRUE)
  
  choose_log_cols <- function(X) {
    if (!log_transform) return(integer(0))
    sk <- apply(X, 2, skewness_v4)
    ok_nonneg <- apply(X, 2, function(z) all(z[is.finite(z)] >= 0))
    which(is.finite(sk) & sk > skew_threshold & ok_nonneg)
  }
  
  gmm_log_cols <- choose_log_cols(X_gmm)
  cox_log_cols <- choose_log_cols(X_cox)
  
  if (length(gmm_log_cols) > 0) {
    X_gmm[, gmm_log_cols] <- log1p(pmax(X_gmm[, gmm_log_cols, drop = FALSE], 0))
  }
  if (length(cox_log_cols) > 0) {
    X_cox[, cox_log_cols] <- log1p(pmax(X_cox[, cox_log_cols, drop = FALSE], 0))
  }
  
  get_fill <- function(X) {
    if (impute == "median") {
      out <- apply(X, 2, function(z) stats::median(z, na.rm = TRUE))
    } else {
      out <- colMeans(X, na.rm = TRUE)
    }
    out[!is.finite(out)] <- 0
    out
  }
  
  gmm_fill <- get_fill(X_gmm)
  cox_fill <- get_fill(X_cox)
  
  fill_na <- function(X, fill) {
    for (j in seq_len(ncol(X))) {
      idx <- which(!is.finite(X[, j]))
      if (length(idx) > 0) X[idx, j] <- fill[j]
    }
    X
  }
  
  X_gmm_imp <- fill_na(X_gmm, gmm_fill)
  X_cox_imp <- fill_na(X_cox, cox_fill)
  
  list(
    X_gmm = X_gmm_imp,
    X_cox = X_cox_imp,
    gmm_log_cols = gmm_log_cols,
    cox_log_cols = cox_log_cols,
    gmm_fill = gmm_fill,
    cox_fill = cox_fill,
    log_transform = log_transform,
    skew_threshold = skew_threshold,
    impute = impute,
    gmm_names = colnames(X_gmm_imp),
    cox_names = colnames(X_cox_imp)
  )
}

preprocess_test_v4 <- function(X_gmm_new, X_cox_new, prep) {
  X_gmm_new <- as.matrix(X_gmm_new)
  X_cox_new <- as.matrix(X_cox_new)
  
  colnames(X_gmm_new) <- make.names(colnames(X_gmm_new), unique = TRUE)
  colnames(X_cox_new) <- make.names(colnames(X_cox_new), unique = TRUE)
  
  if (!identical(colnames(X_gmm_new), prep$gmm_names)) {
    stop("X_gmm_new columns do not match training columns.")
  }
  if (!identical(colnames(X_cox_new), prep$cox_names)) {
    stop("X_cox_new columns do not match training columns.")
  }
  
  if (length(prep$gmm_log_cols) > 0) {
    X_gmm_new[, prep$gmm_log_cols] <- log1p(pmax(X_gmm_new[, prep$gmm_log_cols, drop = FALSE], 0))
  }
  if (length(prep$cox_log_cols) > 0) {
    X_cox_new[, prep$cox_log_cols] <- log1p(pmax(X_cox_new[, prep$cox_log_cols, drop = FALSE], 0))
  }
  
  for (j in seq_len(ncol(X_gmm_new))) {
    idx <- which(!is.finite(X_gmm_new[, j]))
    if (length(idx) > 0) X_gmm_new[idx, j] <- prep$gmm_fill[j]
  }
  for (j in seq_len(ncol(X_cox_new))) {
    idx <- which(!is.finite(X_cox_new[, j]))
    if (length(idx) > 0) X_cox_new[idx, j] <- prep$cox_fill[j]
  }
  
  list(X_gmm = X_gmm_new, X_cox = X_cox_new)
}

###############################################################
## 2. Shared-baseline pieces
###############################################################

shared_breslow_v4 <- function(time, status, eta_list, tau,
                              denom_floor = 1e-8,
                              max_jump = 10,
                              max_cumhaz = 500) {
  status <- as.integer(status)
  uniqT <- sort(unique(time[status == 1]))
  if (length(uniqT) == 0) {
    return(list(time = numeric(0), jump = numeric(0), cumhaz = numeric(0)))
  }
  
  C <- ncol(tau)
  exp_eta_list <- lapply(eta_list, function(eta) exp(clamp(as.numeric(eta), -30, 30)))
  
  jump <- numeric(length(uniqT))
  for (m in seq_along(uniqT)) {
    t_m <- uniqT[m]
    num <- sum(status[time == t_m])
    at_risk <- (time >= t_m)
    denom <- 0
    for (c in seq_len(C)) {
      denom <- denom + sum(tau[at_risk, c] * exp_eta_list[[c]][at_risk])
    }
    denom <- max(denom, denom_floor)
    jump[m] <- min(num / denom, max_jump)
  }
  
  list(time = uniqT, jump = jump, cumhaz = pmin(cumsum(jump), max_cumhaz))
}

cox_poisson_weighted_v4 <- function(time, status, X, w, baseline,
                                    lambda = 0.1,
                                    alpha = 0.5,
                                    min_eff_events = 0.5,
                                    fallback_fit = NULL) {
  X <- as.matrix(X)
  status <- as.integer(status)
  w <- as.numeric(w)
  w[!is.finite(w) | w < 0] <- 0
  
  scaler <- make_x_scaler(X)
  Xs <- apply_x_scaler(X, scaler)
  eff_events <- sum(w * status)
  
  if (!is.finite(eff_events) || eff_events < min_eff_events) {
    if (!is.null(fallback_fit)) {
      out <- fallback_fit
      out$baseline <- baseline
      out$ok <- FALSE
      out$reason <- "fallback_low_events"
      out$eff_events <- eff_events
      return(out)
    }
    return(list(
      beta_s = rep(0, ncol(X)),
      beta_orig = rep(0, ncol(X)),
      scaler = scaler,
      baseline = baseline,
      lambda = NA_real_,
      ok = FALSE,
      reason = "low_events",
      eff_events = eff_events
    ))
  }
  
  Lambda_t <- get_cumhaz_at(baseline, time)
  Lambda_t <- pmax(Lambda_t, 1e-10)
  offset <- log(Lambda_t)
  
  fit <- tryCatch({
    gfit <- glmnet::glmnet(
      x = Xs,
      y = status,
      family = "poisson",
      weights = w,
      offset = offset,
      alpha = alpha,
      lambda = lambda,
      standardize = FALSE,
      maxit = 1e5
    )
    beta_s <- as.numeric(stats::coef(gfit, s = lambda))[-1]
    list(beta_s = beta_s, ok = TRUE, reason = "ok", lambda = lambda)
  }, error = function(e) NULL)
  
  if (is.null(fit)) {
    if (!is.null(fallback_fit)) {
      out <- fallback_fit
      out$baseline <- baseline
      out$ok <- FALSE
      out$reason <- "fallback_poisson_failed"
      out$eff_events <- eff_events
      return(out)
    }
    fit <- list(beta_s = rep(0, ncol(X)), ok = FALSE, reason = "poisson_failed", lambda = lambda)
  }
  
  beta_s <- fit$beta_s
  beta_s[!is.finite(beta_s)] <- 0
  
  list(
    beta_s = beta_s,
    beta_orig = beta_s / scaler$scale,
    scaler = scaler,
    baseline = baseline,
    lambda = fit$lambda,
    ok = isTRUE(fit$ok),
    reason = fit$reason,
    eff_events = eff_events
  )
}

###############################################################
## 3. Shared-baseline EM fitter
###############################################################

gemcox_full_shared_v4 <- function(
    X_gmm, X_cox, time, status,
    K = 2,
    lambda = NULL,
    alpha = 0.5,
    max_iter = 60,
    tol = 1e-4,
    verbose = TRUE,
    gmm_ridge = 1e-6,
    gmm_diag = TRUE,
    normalize_gmm_by_dim = FALSE,
    surv_weight = 1.0,
    temp = 1.0,
    pi_floor = 0.02,
    min_eff_events = 1.0,
    init_method = c("supervised", "kmeans"),
    init_eta_weight = 0.75,
    init_seed = 1) {
  
  init_method <- match.arg(init_method)
  X_gmm <- as.matrix(X_gmm)
  X_cox <- as.matrix(X_cox)
  time <- as.numeric(time)
  status <- as.integer(status)
  
  n <- nrow(X_gmm)
  stopifnot(nrow(X_cox) == n, length(time) == n, length(status) == n)
  
  q <- ncol(X_gmm)
  p <- ncol(X_cox)
  gmm_log_density_scale <- if (isTRUE(normalize_gmm_by_dim)) 1 / max(q, 1) else 1
  
  colnames(X_gmm) <- make.names(colnames(X_gmm), unique = TRUE)
  colnames(X_cox) <- make.names(colnames(X_cox), unique = TRUE)
  
  gmm_scaler <- make_x_scaler(X_gmm)
  X_gmm_s <- apply_x_scaler(X_gmm, gmm_scaler)
  
  if (is.null(lambda) || is.character(lambda)) {
    lambda_use <- select_lambda(time, status, X_cox, alpha = alpha)
  } else {
    lambda_use <- as.numeric(lambda)
  }
  
  global_coxfit <- cox_glmnet_weighted(
    time = time,
    status = status,
    X = X_cox,
    w = rep(1, n),
    lambda = lambda_use,
    alpha = alpha,
    min_eff_events = 0
  )
  
  tau <- initialize_tau(
    X_gmm_s = X_gmm_s,
    X_cox = X_cox,
    global_coxfit = global_coxfit,
    K = K,
    init_method = init_method,
    init_eta_weight = init_eta_weight,
    init_seed = init_seed
  )
  
  pi_c <- pmax(colMeans(tau), 1e-8)
  pi_c <- pi_c / sum(pi_c)
  
  mu_gmm <- matrix(0, q, K)
  Sigma_list <- vector("list", K)
  coxfit <- vector("list", K)
  beta <- matrix(0, p, K, dimnames = list(colnames(X_cox), paste0("C", seq_len(K))))
  
  ## Initial GMM moments
  for (k in seq_len(K)) {
    w <- tau[, k]
    wsum <- sum(w)
    if (!is.finite(wsum) || wsum < 1e-8) {
      mu_gmm[, k] <- rep(0, q)
      Sigma_list[[k]] <- diag(1, q)
    } else {
      mu_gmm[, k] <- colSums(X_gmm_s * w) / wsum
      xc <- sweep(X_gmm_s, 2, mu_gmm[, k], "-")
      S <- (t(xc * sqrt(w)) %*% (xc * sqrt(w))) / wsum
      S <- (S + t(S)) / 2
      if (gmm_diag) S <- diag(diag(S), q)
      Sigma_list[[k]] <- as.matrix(S) + diag(gmm_ridge, q)
    }
  }
  
  ## Initial Cox betas from weighted Cox, then convert to shared baseline
  for (k in seq_len(K)) {
    coxfit[[k]] <- cox_glmnet_weighted(
      time = time,
      status = status,
      X = X_cox,
      w = tau[, k],
      lambda = lambda_use,
      alpha = alpha,
      min_eff_events = min_eff_events,
      fallback_fit = global_coxfit
    )
  }
  eta_list <- lapply(seq_len(K), function(k) {
    eta_from_scaled(X_cox, coxfit[[k]]$beta_s, coxfit[[k]]$scaler, eta_cap = 12)
  })
  shared_base <- shared_breslow_v4(time, status, eta_list, tau)
  for (k in seq_len(K)) {
    global_fb <- global_coxfit
    global_fb$baseline <- shared_base
    coxfit[[k]] <- cox_poisson_weighted_v4(
      time = time,
      status = status,
      X = X_cox,
      w = tau[, k],
      baseline = shared_base,
      lambda = lambda_use,
      alpha = alpha,
      min_eff_events = min_eff_events,
      fallback_fit = global_fb
    )
    beta[, k] <- coxfit[[k]]$beta_orig
  }
  
  loglik_trace <- numeric(max_iter)
  temp_use <- if (!is.finite(temp) || temp <= 0) 1.0 else temp
  
  for (iter in seq_len(max_iter)) {
    pi_c <- pmax(pi_c, 1e-8)
    pi_c <- pi_c / sum(pi_c)
    
    ## E-step
    log_resp <- matrix(NA_real_, n, K)
    for (k in seq_len(K)) {
      lgmm <- dmvnorm_log(X_gmm_s, mu_gmm[, k], Sigma_list[[k]], ridge = gmm_ridge)
      lsurv <- log_surv_density(time, status, X_cox, coxfit[[k]])
      log_resp[, k] <- log(pi_c[k]) + gmm_log_density_scale * lgmm + surv_weight * lsurv
    }
    log_resp[!is.finite(log_resp)] <- -700
    
    if (abs(temp_use - 1) < 1e-9) {
      loglik_i <- rowLogSumExp(log_resp)
    } else {
      loglik_i <- rowLogSumExp(log_resp / temp_use) * temp_use
    }
    loglik <- sum(loglik_i)
    loglik_trace[iter] <- loglik
    
    log_resp_t <- log_resp / temp_use
    lse_t <- rowLogSumExp(log_resp_t)
    tau <- exp(sweep(log_resp_t, 1, lse_t, "-"))
    tau <- pmax(tau, 1e-12)
    tau <- tau / rowSums(tau)
    
    if (verbose) cat(sprintf("Iter %d: loglik = %.3f\n", iter, loglik))
    
    if (iter > 1) {
      rel <- abs(loglik - loglik_trace[iter - 1]) / (abs(loglik_trace[iter - 1]) + 1e-8)
      if (rel < tol) {
        if (verbose) cat("Converged.\n")
        loglik_trace <- loglik_trace[seq_len(iter)]
        break
      }
    }
    
    ## M-step: pi
    pi_c <- colMeans(tau)
    pi_c <- pmax(pi_c, pi_floor)
    pi_c <- pi_c / sum(pi_c)
    
    ## M-step: GMM
    for (k in seq_len(K)) {
      w <- tau[, k]
      wsum <- sum(w)
      if (!is.finite(wsum) || wsum < 1e-8) next
      mu_gmm[, k] <- colSums(X_gmm_s * w) / wsum
      xc <- sweep(X_gmm_s, 2, mu_gmm[, k], "-")
      S <- (t(xc * sqrt(w)) %*% (xc * sqrt(w))) / wsum
      S <- (S + t(S)) / 2
      if (gmm_diag) S <- diag(diag(S), q)
      Sigma_list[[k]] <- as.matrix(S) + diag(gmm_ridge, q)
    }
    
    ## M-step: shared baseline then cluster-wise Poisson updates
    eta_list <- lapply(seq_len(K), function(k) {
      eta_from_scaled(X_cox, coxfit[[k]]$beta_s, coxfit[[k]]$scaler, eta_cap = 12)
    })
    shared_base <- shared_breslow_v4(time, status, eta_list, tau)
    
    for (k in seq_len(K)) {
      global_fb <- global_coxfit
      global_fb$baseline <- shared_base
      coxfit[[k]] <- cox_poisson_weighted_v4(
        time = time,
        status = status,
        X = X_cox,
        w = tau[, k],
        baseline = shared_base,
        lambda = lambda_use,
        alpha = alpha,
        min_eff_events = min_eff_events,
        fallback_fit = global_fb
      )
      beta[, k] <- coxfit[[k]]$beta_orig
    }
  }
  
  clusterid <- max.col(tau, ties.method = "first")
  eff_events <- colSums(tau * status)
  
  structure(list(
    tau = tau,
    clusterid = clusterid,
    pi = pi_c,
    gmm_scaler = gmm_scaler,
    mu_gmm = mu_gmm,
    Sigma_list = Sigma_list,
    beta = beta,
    beta_list = lapply(seq_len(K), function(k) beta[, k]),
    coxfit = coxfit,
    global_coxfit = global_coxfit,
    loglik = loglik_trace,
    eff_events = eff_events,
    lambda_used = lambda_use,
    K = K,
    gmm_diag = gmm_diag,
    gmm_ridge = gmm_ridge,
    gmm_log_density_scale = gmm_log_density_scale,
    normalize_gmm_by_dim = normalize_gmm_by_dim,
    surv_weight = surv_weight,
    temp = temp_use,
    min_eff_events = min_eff_events,
    init_method = init_method,
    init_eta_weight = init_eta_weight,
    baseline_mode = "shared",
    call = match.call()
  ), class = "gemcox_fit")
}

gemcox_full_multistart_shared_v4 <- function(...,
                                             n_starts = 10,
                                             init_seeds = NULL,
                                             verbose = FALSE) {
  if (is.null(init_seeds)) init_seeds <- seq_len(n_starts)
  if (length(init_seeds) < n_starts) init_seeds <- rep(init_seeds, length.out = n_starts)
  
  fits <- vector("list", n_starts)
  final_ll <- rep(NA_real_, n_starts)
  
  for (s in seq_len(n_starts)) {
    fit_s <- tryCatch(
      gemcox_full_shared_v4(..., init_seed = init_seeds[s], verbose = FALSE),
      error = function(e) NULL
    )
    fits[[s]] <- fit_s
    final_ll[s] <- if (is.null(fit_s)) NA_real_ else last_finite(fit_s$loglik)
  }
  
  ok <- is.finite(final_ll)
  if (!any(ok)) stop("All multistart fits failed.")
  best_idx <- which.max(final_ll)
  best_fit <- fits[[best_idx]]
  best_fit$multistart <- list(
    n_starts = n_starts,
    init_seeds = init_seeds,
    final_loglik = final_ll,
    best_start = best_idx
  )
  if (verbose) {
    cat(sprintf("Selected shared-baseline multistart solution with final loglik = %.3f (start %d)\n",
                final_ll[best_idx], best_idx))
  }
  best_fit
}

###############################################################
## 4. CV with fold-safe preprocessing + correct 1-SE rule
###############################################################

cv_select_K_gamma_v4 <- function(X_gmm, X_cox, time, status,
                                 K_grid = 1:3,
                                 gamma_grid = c(0, 0.25, 0.5, 1.0),
                                 nfolds = 5,
                                 alpha = 0.5,
                                 max_iter = 60,
                                 tol = 1e-4,
                                 metric = c("cindex", "heldout_loglik"),
                                 risk_type = c("event_prob", "soft_lp", "hard_lp"),
                                 t_eval = NULL,
                                 n_starts = 5,
                                 gmm_diag = TRUE,
                                 normalize_gmm_by_dim = FALSE,
                                 pi_floor = 0.02,
                                 min_eff_events = 1.0,
                                 init_method = c("supervised", "kmeans"),
                                 init_eta_weight = 0.75,
                                 log_transform = TRUE,
                                 skew_threshold = 1.0,
                                 impute = c("median", "mean"),
                                 seed = 1,
                                 verbose = FALSE,
                                 use_1se_rule = TRUE,
                                 baseline_mode = c("shared", "cluster")) {
  
  metric <- match.arg(metric)
  risk_type <- match.arg(risk_type)
  init_method <- match.arg(init_method)
  impute <- match.arg(impute)
  baseline_mode <- match.arg(baseline_mode)
  
  X_gmm <- as.matrix(X_gmm)
  X_cox <- as.matrix(X_cox)
  fold <- make_stratified_folds(status, nfolds = nfolds, seed = seed)
  
  grid <- expand.grid(K = K_grid, gamma = gamma_grid, stringsAsFactors = FALSE)
  fold_records <- vector("list", nrow(grid))
  summary_df <- data.frame(
    K = grid$K,
    gamma = grid$gamma,
    mean_score = NA_real_,
    sd_score = NA_real_,
    se_score = NA_real_,
    n_ok_folds = NA_integer_,
    n_fail_folds = NA_integer_,
    stringsAsFactors = FALSE
  )
  
  for (g in seq_len(nrow(grid))) {
    K_cur <- grid$K[g]
    gamma_cur <- grid$gamma[g]
    if (verbose) message("CV for K=", K_cur, ", gamma=", gamma_cur)
    
    scores <- rep(NA_real_, nfolds)
    reasons <- rep(NA_character_, nfolds)
    
    for (v in seq_len(nfolds)) {
      tr <- fold != v
      te <- fold == v
      
      prep_tr <- preprocess_train_v4(
        X_gmm[tr, , drop = FALSE],
        X_cox[tr, , drop = FALSE],
        log_transform = log_transform,
        skew_threshold = skew_threshold,
        impute = impute
      )
      prep_te <- preprocess_test_v4(
        X_gmm[te, , drop = FALSE],
        X_cox[te, , drop = FALSE],
        prep_tr
      )
      
      fitter <- if (baseline_mode == "shared") gemcox_full_multistart_shared_v4 else gemcox_full_multistart
      
      fit <- tryCatch(
        fitter(
          X_gmm = prep_tr$X_gmm,
          X_cox = prep_tr$X_cox,
          time = time[tr],
          status = status[tr],
          K = K_cur,
          lambda = NULL,
          alpha = alpha,
          max_iter = max_iter,
          tol = tol,
          verbose = FALSE,
          surv_weight = gamma_cur,
          gmm_diag = gmm_diag,
          normalize_gmm_by_dim = normalize_gmm_by_dim,
          pi_floor = pi_floor,
          min_eff_events = min_eff_events,
          init_method = init_method,
          init_eta_weight = init_eta_weight,
          n_starts = n_starts,
          init_seeds = seed + seq_len(n_starts) + 100 * v + 1000 * g
        ),
        error = function(e) e
      )
      
      if (inherits(fit, "error") || is.null(fit)) {
        reasons[v] <- if (inherits(fit, "error")) paste("fit_error:", fit$message) else "fit_null"
        next
      }
      
      tau_te <- tryCatch(
        predict_tau(fit, prep_te$X_gmm),
        error = function(e) e
      )
      if (inherits(tau_te, "error") || is.null(tau_te)) {
        reasons[v] <- if (inherits(tau_te, "error")) paste("tau_error:", tau_te$message) else "tau_null"
        next
      }
      
      t_eval_fold <- t_eval
      if (is.null(t_eval_fold)) {
        ev_train <- time[tr & status == 1]
        t_eval_fold <- if (length(ev_train) > 0) stats::median(ev_train) else stats::median(time[tr])
      }
      
      sc <- if (metric == "heldout_loglik") {
        tryCatch(
          heldout_loglik(
            fit = fit,
            X_gmm_new = prep_te$X_gmm,
            X_cox_new = prep_te$X_cox,
            time_new = time[te],
            status_new = status[te],
            tau_new = tau_te
          ),
          error = function(e) NA_real_
        )
      } else {
        risk_te <- tryCatch(
          predict_risk(
            fit,
            X_cox_new = prep_te$X_cox,
            tau_new = tau_te,
            type = risk_type,
            t_eval = t_eval_fold
          ),
          error = function(e) e
        )
        if (inherits(risk_te, "error") || is.null(risk_te)) {
          reasons[v] <- if (inherits(risk_te, "error")) paste("risk_error:", risk_te$message) else "risk_null"
          NA_real_
        } else {
          risk_te <- as.numeric(risk_te) + stats::rnorm(length(risk_te), sd = 1e-8)
          tryCatch(c_index(time[te], status[te], risk_te), error = function(e) NA_real_)
        }
      }
      
      scores[v] <- sc
      reasons[v] <- ifelse(is.finite(sc), "ok", paste0(metric, "_na"))
    }
    
    ok <- is.finite(scores)
    summary_df$mean_score[g] <- if (any(ok)) mean(scores[ok]) else NA_real_
    summary_df$sd_score[g] <- if (sum(ok) >= 2) stats::sd(scores[ok]) else NA_real_
    summary_df$se_score[g] <- if (sum(ok) >= 2) stats::sd(scores[ok]) / sqrt(sum(ok)) else 0
    summary_df$n_ok_folds[g] <- sum(ok)
    summary_df$n_fail_folds[g] <- nfolds - sum(ok)
    
    fold_records[[g]] <- data.frame(
      K = K_cur,
      gamma = gamma_cur,
      fold = seq_len(nfolds),
      score = scores,
      reason = reasons,
      stringsAsFactors = FALSE
    )
  }
  
  ok_rows <- which(is.finite(summary_df$mean_score))
  best <- NULL
  best_method <- NULL
  if (length(ok_rows) > 0) {
    best_idx <- ok_rows[which.max(summary_df$mean_score[ok_rows])]
    best <- summary_df[best_idx, , drop = FALSE]
    best_method <- "max_mean"
    
    if (use_1se_rule) {
      threshold <- summary_df$mean_score[best_idx] - summary_df$se_score[best_idx]
      eligible <- summary_df[ok_rows, , drop = FALSE]
      eligible <- eligible[eligible$mean_score >= threshold, , drop = FALSE]
      eligible <- eligible[order(eligible$K, eligible$gamma), , drop = FALSE]
      if (nrow(eligible) > 0) {
        best <- eligible[1, , drop = FALSE]
        best_method <- "1se_rule"
      }
    }
  }
  
  list(
    metric = metric,
    summary = summary_df,
    folds = fold_records,
    best = best,
    best_method = best_method,
    baseline_mode = baseline_mode
  )
}

###############################################################
## 5. Final pipeline
###############################################################

gemcox_pipeline_v4 <- function(
    X_gmm, X_cox, time, status,
    K_grid = 1:3,
    gamma_grid = c(0, 0.25, 0.5, 1.0),
    alpha = 0.5,
    log_transform = TRUE,
    skew_threshold = 1.0,
    impute = c("median", "mean"),
    use_1se_rule = TRUE,
    nfolds = 5,
    n_starts = 5,
    n_starts_final = 10,
    max_iter = 60,
    seed = 1,
    verbose = TRUE,
    baseline_mode = c("shared", "cluster")) {
  
  baseline_mode <- match.arg(baseline_mode)
  impute <- match.arg(impute)
  
  if (verbose) cat("Cross-validating (fold-safe preprocessing)...\n")
  cv_res <- cv_select_K_gamma_v4(
    X_gmm = X_gmm,
    X_cox = X_cox,
    time = time,
    status = status,
    K_grid = K_grid,
    gamma_grid = gamma_grid,
    nfolds = nfolds,
    alpha = alpha,
    max_iter = max_iter,
    n_starts = n_starts,
    log_transform = log_transform,
    skew_threshold = skew_threshold,
    impute = impute,
    use_1se_rule = use_1se_rule,
    seed = seed,
    verbose = verbose,
    baseline_mode = baseline_mode
  )
  
  if (is.null(cv_res$best) || nrow(cv_res$best) != 1) {
    stop("CV failed to select a valid model.")
  }
  
  prep_full <- preprocess_train_v4(
    X_gmm, X_cox,
    log_transform = log_transform,
    skew_threshold = skew_threshold,
    impute = impute
  )
  
  fitter <- if (baseline_mode == "shared") gemcox_full_multistart_shared_v4 else gemcox_full_multistart
  final_fit <- fitter(
    X_gmm = prep_full$X_gmm,
    X_cox = prep_full$X_cox,
    time = time,
    status = status,
    K = cv_res$best$K[1],
    lambda = NULL,
    alpha = alpha,
    max_iter = max_iter * 2,
    surv_weight = cv_res$best$gamma[1],
    verbose = verbose,
    n_starts = n_starts_final,
    init_seeds = seed + seq_len(n_starts_final)
  )
  
  list(
    fit = final_fit,
    cv_result = cv_res,
    preprocessing = prep_full,
    K_selected = cv_res$best$K[1],
    gamma_selected = cv_res$best$gamma[1],
    baseline_mode = baseline_mode
  )
}

###############################################################
## 6. CVIA078 helpers (no PCA)
###############################################################

prepare_cvia078_data_v4 <- function(dat_long,
                                    gmm_feature_mode = c("delta_only", "v20_only", "v20_plus_delta"),
                                    cox_feature_mode = c("v20_plus_delta_covars", "delta_covars", "selected_plus_covars"),
                                    missing_threshold = 0.05,
                                    corr_prune = 0.85,
                                    add_covars = TRUE) {
  gmm_feature_mode <- match.arg(gmm_feature_mode)
  cox_feature_mode <- match.arg(cox_feature_mode)
  
  needed <- c("SUBJID", "AGE", "SEX", "SICKLE", "daysfollowup", "infected", "visitid", "value", "variable")
  miss_needed <- setdiff(needed, names(dat_long))
  if (length(miss_needed) > 0) stop("Missing required columns: ", paste(miss_needed, collapse = ", "))
  
  suppressPackageStartupMessages({
    require(dplyr)
    require(tidyr)
  })
  
  dat_wide <- dat_long |>
    dplyr::mutate(feature_visit = paste0(make.names(variable), "_", visitid)) |>
    dplyr::select(SUBJID, AGE, SEX, SICKLE, daysfollowup, infected, feature_visit, value) |>
    tidyr::pivot_wider(names_from = feature_visit, values_from = value)
  
  ## Delta features
  base_vars <- unique(make.names(dat_long$variable))
  for (v in base_vars) {
    c5 <- paste0(v, "_V5")
    c20 <- paste0(v, "_V20")
    if (c5 %in% names(dat_wide) && c20 %in% names(dat_wide)) {
      dat_wide[[paste0(v, "_dV20_V5")]] <- dat_wide[[c20]] - dat_wide[[c5]]
    }
  }
  
  ## Numeric covariates
  dat_wide$SEX_M <- as.integer(toupper(substr(dat_wide$SEX, 1, 1)) == "M")
  dat_wide$SICKLE_ABNORMAL <- as.integer(!grepl("NORMAL", dat_wide$SICKLE, ignore.case = TRUE))
  
  meta_cols <- c("SUBJID", "AGE", "SEX", "SICKLE", "daysfollowup", "infected", "SEX_M", "SICKLE_ABNORMAL")
  feat_cols <- setdiff(names(dat_wide), meta_cols)
  miss_rate <- sapply(dat_wide[, feat_cols, drop = FALSE], function(z) mean(!is.finite(z)))
  
  v20_cols <- names(miss_rate)[grepl("_V20$", names(miss_rate)) & miss_rate <= missing_threshold]
  delta_cols <- names(miss_rate)[grepl("_dV20_V5$", names(miss_rate)) & miss_rate <= missing_threshold]
  sel_cols <- switch(gmm_feature_mode,
                     delta_only = delta_cols,
                     v20_only = v20_cols,
                     v20_plus_delta = c(v20_cols, delta_cols))
  
  ## Optional correlation pruning to keep interpretable representatives
  if (length(sel_cols) > 1 && is.finite(corr_prune) && corr_prune < 1) {
    X0 <- dat_wide[, sel_cols, drop = FALSE]
    cor0 <- suppressWarnings(stats::cor(X0, use = "pairwise.complete.obs"))
    keep <- character(0)
    for (nm in sel_cols) {
      if (length(keep) == 0) {
        keep <- c(keep, nm)
      } else {
        if (all(abs(cor0[nm, keep]) < corr_prune | !is.finite(cor0[nm, keep]))) {
          keep <- c(keep, nm)
        }
      }
    }
    sel_cols <- keep
  }
  
  cox_cols <- switch(cox_feature_mode,
                     v20_plus_delta_covars = unique(c(v20_cols, delta_cols,
                                                      if (add_covars) c("AGE", "SEX_M", "SICKLE_ABNORMAL") else NULL)),
                     delta_covars = unique(c(delta_cols,
                                             if (add_covars) c("AGE", "SEX_M", "SICKLE_ABNORMAL") else NULL)),
                     selected_plus_covars = unique(c(sel_cols,
                                                     if (add_covars) c("AGE", "SEX_M", "SICKLE_ABNORMAL") else NULL)))
  
  X_gmm <- as.matrix(dat_wide[, sel_cols, drop = FALSE])
  X_cox <- as.matrix(dat_wide[, cox_cols, drop = FALSE])
  storage.mode(X_gmm) <- "numeric"
  storage.mode(X_cox) <- "numeric"
  
  list(
    data_wide = dat_wide,
    X_gmm = X_gmm,
    X_cox = X_cox,
    time = dat_wide$daysfollowup,
    status = dat_wide$infected,
    ids = dat_wide$SUBJID,
    gmm_features = sel_cols,
    cox_features = cox_cols,
    missing_rate = miss_rate
  )
}

fit_cvia078_gemcox_v4 <- function(dat_long,
                                  K_grid = 1:2,
                                  gamma_grid = c(0, 0.25, 0.5, 1.0),
                                  gmm_feature_mode = "delta_only",
                                  cox_feature_mode = "v20_plus_delta_covars",
                                  missing_threshold = 0.05,
                                  corr_prune = 0.85,
                                  alpha = 0.5,
                                  nfolds = 5,
                                  n_starts = 5,
                                  n_starts_final = 10,
                                  max_iter = 60,
                                  seed = 1,
                                  verbose = TRUE,
                                  baseline_mode = "shared") {
  prep <- prepare_cvia078_data_v4(
    dat_long,
    gmm_feature_mode = gmm_feature_mode,
    cox_feature_mode = cox_feature_mode,
    missing_threshold = missing_threshold,
    corr_prune = corr_prune
  )
  
  fit <- gemcox_pipeline_v4(
    X_gmm = prep$X_gmm,
    X_cox = prep$X_cox,
    time = prep$time,
    status = prep$status,
    K_grid = K_grid,
    gamma_grid = gamma_grid,
    alpha = alpha,
    log_transform = TRUE,
    use_1se_rule = TRUE,
    nfolds = nfolds,
    n_starts = n_starts,
    n_starts_final = n_starts_final,
    max_iter = max_iter,
    seed = seed,
    verbose = verbose,
    baseline_mode = baseline_mode
  )
  
  list(data = prep, model = fit)
}

cat("GeMCox_improved_v4_patch.R loaded.\n")
cat("  Main fixes: fold-safe preprocessing, shared baseline fit, correct 1-SE rule.\n")
cat("  Main helper: fit_cvia078_gemcox_v4(dat_long, ...) with default delta-only clustering.\n")




