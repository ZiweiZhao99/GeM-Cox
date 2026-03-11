###############################################################
## GeM-Cox: Generative Mixture of Cox Regression Models
## Next improved version
##
## Main changes relative to the uploaded version:
##   1) Fixes static code issues (undefined min_eff_events,
##      unsupported x_scale/clamp_eta arguments).
##   2) Makes training-time and prediction-time tau consistent:
##      the GMM log-density scale used in the E-step is stored
##      and reused in predict_tau().
##   3) Uses model-based prediction for CV by default:
##      predicted event probability at t_eval from the mixture
##      survival S(t|x) = sum_c tau_c(x) S_c(t|x).
##   4) Adds held-out mixture log-likelihood as an alternative
##      CV criterion.
##   5) Uses diagonal GMM covariance by default, which is much
##      more stable when K > 1 or effective cluster size is small.
##   6) Adds multistart EM to reduce local-optimum / empty-cluster
##      collapse.
##   7) Keeps lambda selection outside EM and allows joint tuning
##      of K and gamma (= survival weight in E-step).
##
## Notes:
##   - This is still a penalized / generalized EM implementation.
##   - surv_weight != 1 or temp != 1 makes the fitting objective a
##     tempered / pseudo-likelihood rather than the exact stated
##     generative likelihood.
###############################################################

library(survival)
library(glmnet)
library(MASS)

###############################################################
## Utilities
###############################################################
clamp <- function(x, lo = -30, hi = 30) pmin(pmax(x, lo), hi)

make_x_scaler <- function(X) {
  X <- as.matrix(X)
  center <- colMeans(X, na.rm = TRUE)
  scalev <- apply(X, 2, stats::sd, na.rm = TRUE)
  scalev[!is.finite(scalev) | scalev <= 0] <- 1
  list(center = center, scale = scalev)
}

apply_x_scaler <- function(X, scaler) {
  X <- as.matrix(X)
  sweep(sweep(X, 2, scaler$center, "-"), 2, scaler$scale, "/")
}

eta_from_scaled <- function(X, beta_s, scaler, eta_cap = 12) {
  Xs <- apply_x_scaler(X, scaler)
  eta <- as.numeric(Xs %*% beta_s)
  clamp(eta, -eta_cap, eta_cap)
}

rowLogSumExp <- function(A) {
  m <- apply(A, 1, max)
  m[!is.finite(m)] <- 0
  s <- rowSums(exp(A - m))
  s <- pmax(s, 1e-300)
  m + log(s)
}

last_finite <- function(x) {
  xf <- x[is.finite(x)]
  if (length(xf) == 0) NA_real_ else tail(xf, 1)
}

make_stratified_folds <- function(status, nfolds = 5, seed = 1) {
  set.seed(seed)
  fold <- integer(length(status))
  for (grp in c(0, 1)) {
    idx <- which(status == grp)
    fold[idx] <- sample(rep(seq_len(nfolds), length.out = length(idx)))
  }
  fold
}

###############################################################
## Multivariate normal log-density
###############################################################
dmvnorm_log <- function(X, mu, Sigma, ridge = 1e-6) {
  X <- as.matrix(X)
  mu <- as.numeric(mu)
  d <- length(mu)
  xc <- sweep(X, 2, mu, "-")

  S2 <- as.matrix(Sigma) + diag(ridge, d)
  cholS <- tryCatch(chol(S2), error = function(e) NULL)
  if (is.null(cholS)) {
    S2 <- as.matrix(Sigma) + diag(1e-3, d)
    cholS <- chol(S2)
  }

  Sinv <- chol2inv(cholS)
  logdet <- 2 * sum(log(diag(cholS)))
  quad <- rowSums((xc %*% Sinv) * xc)
  -0.5 * (d * log(2 * pi) + logdet + quad)
}

###############################################################
## Weighted Breslow estimator (cluster-specific baseline)
###############################################################
breslow_weighted <- function(time, status, eta, w,
                             denom_floor = 1e-8,
                             max_jump = 10,
                             max_cumhaz = 500) {
  w <- as.numeric(w)
  w[!is.finite(w) | w < 0] <- 0
  status <- as.integer(status)
  eta <- clamp(as.numeric(eta), -30, 30)
  exp_eta <- exp(eta)

  if (sum(w * status) < 1e-10) {
    return(list(time = numeric(0), jump = numeric(0), cumhaz = numeric(0)))
  }

  uniqT <- sort(unique(time[status == 1]))
  jump <- numeric(length(uniqT))

  for (m in seq_along(uniqT)) {
    t_m <- uniqT[m]
    num <- sum(w[time == t_m & status == 1])
    at_risk <- (time >= t_m)
    denom <- sum(w[at_risk] * exp_eta[at_risk])
    denom <- max(denom, denom_floor)
    jump[m] <- min(num / denom, max_jump)
  }

  cumhaz <- pmin(cumsum(jump), max_cumhaz)
  list(time = uniqT, jump = jump, cumhaz = cumhaz)
}

get_cumhaz_at <- function(baseline, times) {
  tt <- as.numeric(times)
  if (length(baseline$time) == 0) return(rep(0, length(tt)))
  idx <- findInterval(tt, baseline$time)
  ifelse(idx > 0, baseline$cumhaz[idx], 0)
}

###############################################################
## Lambda pre-selection (once per training fold)
###############################################################
select_lambda <- function(time, status, X,
                          alpha = 0,
                          lambda_choice = c("lambda.1se", "lambda.min"),
                          nfolds_inner = 3,
                          fallback_lambda = 0.1) {
  lambda_choice <- match.arg(lambda_choice)
  X <- as.matrix(X)
  scaler <- make_x_scaler(X)
  Xs <- apply_x_scaler(X, scaler)

  lam <- tryCatch({
    cvfit <- cv.glmnet(
      x = Xs,
      y = Surv(time, status),
      family = "cox",
      alpha = alpha,
      standardize = FALSE,
      nfolds = nfolds_inner,
      maxit = 1e5
    )
    if (lambda_choice == "lambda.1se") cvfit$lambda.1se else cvfit$lambda.min
  }, error = function(e) fallback_lambda)

  if (!is.finite(lam) || lam <= 0) lam <- fallback_lambda
  lam
}

###############################################################
## Weighted Cox fit via glmnet + Breslow baseline
###############################################################
cox_glmnet_weighted <- function(time, status, X, w,
                                lambda = 0.1,
                                alpha = 0,
                                min_eff_events = 0.5,
                                denom_floor = 1e-8,
                                max_jump = 10,
                                max_cumhaz = 500,
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
      out <- fallback_fit
      out$ok <- FALSE
      out$reason <- "fallback_global_low_events"
      out$eff_events <- n_eff_event
      return(out)
    }
    beta_s <- rep(0, ncol(X))
    baseline <- breslow_weighted(
      time, status, eta = rep(0, length(time)), w = w,
      denom_floor = denom_floor,
      max_jump = max_jump,
      max_cumhaz = max_cumhaz
    )
    return(list(
      beta_s = beta_s,
      beta_orig = beta_s / scaler$scale,
      scaler = scaler,
      baseline = baseline,
      lambda = NA_real_,
      ok = FALSE,
      reason = "zero_effective_events",
      eff_events = n_eff_event
    ))
  }

  if (is.character(lambda)) {
    warning("cox_glmnet_weighted: lambda should be a numeric scalar. Falling back to 0.1.")
    lambda <- 0.1
  }

  fit <- tryCatch({
    gfit <- glmnet(
      x = Xs,
      y = Surv(time, status),
      family = "cox",
      weights = w,
      alpha = alpha,
      lambda = lambda,
      standardize = FALSE,
      maxit = 1e5
    )
    beta_s <- as.numeric(coef(gfit, s = lambda))
    list(beta_s = beta_s, lambda = lambda, ok = TRUE, reason = "ok")
  }, error = function(e) NULL)

  if (is.null(fit)) {
    gfit <- tryCatch(glmnet(
      x = Xs,
      y = Surv(time, status),
      family = "cox",
      weights = w,
      alpha = 0,
      lambda = lambda * 2,
      standardize = FALSE,
      maxit = 1e5
    ), error = function(e) NULL)

    if (!is.null(gfit)) {
      fit <- list(
        beta_s = as.numeric(coef(gfit, s = lambda * 2)),
        lambda = lambda * 2,
        ok = TRUE,
        reason = "ridge_fallback"
      )
    } else if (!is.null(fallback_fit)) {
      out <- fallback_fit
      out$ok <- FALSE
      out$reason <- "fallback_global_glmnet_failed"
      out$eff_events <- n_eff_event
      return(out)
    } else {
      fit <- list(
        beta_s = rep(0, ncol(X)),
        lambda = NA_real_,
        ok = FALSE,
        reason = "glmnet_failed"
      )
    }
  }

  beta_s <- fit$beta_s
  beta_s[!is.finite(beta_s)] <- 0
  beta_orig <- beta_s / scaler$scale

  eta <- eta_from_scaled(X, beta_s, scaler, eta_cap = 12)
  baseline <- breslow_weighted(
    time, status, eta = eta, w = w,
    denom_floor = denom_floor,
    max_jump = max_jump,
    max_cumhaz = max_cumhaz
  )

  list(
    beta_s = beta_s,
    beta_orig = beta_orig,
    scaler = scaler,
    baseline = baseline,
    lambda = fit$lambda,
    ok = isTRUE(fit$ok),
    reason = fit$reason,
    eff_events = n_eff_event
  )
}

###############################################################
## Survival contribution for one cluster
###############################################################
log_surv_density <- function(time, status, X, coxfit) {
  X <- as.matrix(X)
  status <- as.integer(status)
  n <- nrow(X)

  eta <- eta_from_scaled(X, coxfit$beta_s, coxfit$scaler, eta_cap = 12)
  exp_eta <- exp(eta)

  ev_t <- coxfit$baseline$time
  jump <- coxfit$baseline$jump
  cumhaz <- coxfit$baseline$cumhaz

  if (length(ev_t) == 0) {
    return(rep(0, n))
  }

  idx <- findInterval(time, ev_t)
  Lambda <- ifelse(idx > 0, cumhaz[idx], 0)

  haz0 <- rep(1, n)
  ev_idx <- which(status == 1)
  if (length(ev_idx) > 0) {
    safe_idx <- pmax(idx[ev_idx], 1L)    # guard: idx=0 if test time < first training event
    haz0[ev_idx] <- pmax(jump[safe_idx], 1e-12)
  }

  status * (log(haz0) + eta) - Lambda * exp_eta
}

###############################################################
## Internal helper: initialize responsibilities
###############################################################
initialize_tau <- function(X_gmm_s, X_cox, global_coxfit,
                           K = 2,
                           init_method = c("kmeans", "supervised"),
                           init_eta_weight = 0.75,
                           init_seed = 1) {
  init_method <- match.arg(init_method)
  set.seed(init_seed)

  if (K == 1) {
    tau <- matrix(1, nrow(X_gmm_s), 1)
    colnames(tau) <- "C1"
    return(tau)
  }

  if (init_method == "supervised") {
    eta0 <- eta_from_scaled(X_cox, global_coxfit$beta_s, global_coxfit$scaler, eta_cap = 12)
    eta0 <- as.numeric(scale(eta0))
    eta0[!is.finite(eta0)] <- 0
    X_init <- cbind(X_gmm_s, init_eta_weight * eta0)
  } else {
    X_init <- X_gmm_s
  }

  km <- kmeans(X_init, centers = K, nstart = 50)
  Z0 <- km$cluster

  tau <- matrix(0, nrow(X_gmm_s), K)
  tau[cbind(seq_len(nrow(X_gmm_s)), Z0)] <- 1
  colnames(tau) <- paste0("C", seq_len(K))
  tau
}

###############################################################
## Main EM routine
###############################################################
gemcox_full <- function(
    X_gmm, X_cox, time, status,
    K = 2,
    lambda = NULL,
    alpha = 0,
    max_iter = 50,
    tol = 1e-4,
    verbose = TRUE,
    gmm_ridge = 1e-6,
    gmm_diag = TRUE,
    normalize_gmm_by_dim = TRUE,
    surv_weight = 1.0,
    temp = 1.0,
    pi_floor = 0.02,
    min_eff_events = 0.5,
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
  stopifnot(all(status %in% c(0, 1)))

  q <- ncol(X_gmm)
  p <- ncol(X_cox)
  gmm_log_density_scale <- if (isTRUE(normalize_gmm_by_dim)) 1 / max(q, 1) else 1

  colnames(X_gmm) <- make.names(colnames(X_gmm), unique = TRUE)
  colnames(X_cox) <- make.names(colnames(X_cox), unique = TRUE)

  gmm_scaler <- make_x_scaler(X_gmm)
  X_gmm_s <- apply_x_scaler(X_gmm, gmm_scaler)

  if (is.null(lambda) || is.character(lambda)) {
    if (verbose) cat("Pre-selecting lambda via cv.glmnet (K=1, unweighted)...\n")
    lambda_use <- select_lambda(time, status, X_cox, alpha = alpha)
    if (verbose) cat(sprintf("  Selected lambda = %.5f\n", lambda_use))
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

  ## Initial M-step
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

    coxfit[[k]] <- cox_glmnet_weighted(
      time = time,
      status = status,
      X = X_cox,
      w = w,
      lambda = lambda_use,
      alpha = alpha,
      min_eff_events = min_eff_events,
      fallback_fit = global_coxfit
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

    ## M-step: mixing proportions
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

    ## M-step: Cox
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
    call = match.call()
  ), class = "gemcox_fit")
}

###############################################################
## Multi-start wrapper
###############################################################
gemcox_full_multistart <- function(...,
                                   n_starts = 10,
                                   init_seeds = NULL,
                                   verbose = FALSE) {
  if (is.null(init_seeds)) init_seeds <- seq_len(n_starts)
  if (length(init_seeds) < n_starts) {
    init_seeds <- rep(init_seeds, length.out = n_starts)
  }

  fits <- vector("list", n_starts)
  final_ll <- rep(NA_real_, n_starts)

  for (s in seq_len(n_starts)) {
    fit_s <- tryCatch(
      gemcox_full(..., init_seed = init_seeds[s], verbose = FALSE),
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
    cat(sprintf("Selected multistart solution with final loglik = %.3f (start %d)\n",
                final_ll[best_idx], best_idx))
  }
  best_fit
}

###############################################################
## Prediction utilities
###############################################################
predict_tau <- function(fit, X_gmm_new) {
  X_gmm_new <- as.matrix(X_gmm_new)
  Xg <- apply_x_scaler(X_gmm_new, fit$gmm_scaler)

  K <- length(fit$pi)
  log_resp <- matrix(NA_real_, nrow(Xg), K)
  gmm_scale <- if (!is.null(fit$gmm_log_density_scale)) fit$gmm_log_density_scale else 1

  for (k in seq_len(K)) {
    log_resp[, k] <- log(pmax(fit$pi[k], 1e-12)) +
      gmm_scale * dmvnorm_log(Xg, fit$mu_gmm[, k], fit$Sigma_list[[k]], ridge = fit$gmm_ridge)
  }

  log_resp[!is.finite(log_resp)] <- -700
  lse <- rowLogSumExp(log_resp)
  tau <- exp(sweep(log_resp, 1, lse, "-"))
  tau <- pmax(tau, 1e-12)
  tau / rowSums(tau)
}

predict_survival_component <- function(fit, X_cox_new, k, times) {
  X_cox_new <- as.matrix(X_cox_new)
  tt <- as.numeric(times)

  eta_k <- eta_from_scaled(X_cox_new, fit$coxfit[[k]]$beta_s, fit$coxfit[[k]]$scaler, eta_cap = 12)
  exp_eta_k <- exp(eta_k)
  Lambda_k <- get_cumhaz_at(fit$coxfit[[k]]$baseline, tt)

  Sk <- exp(-tcrossprod(exp_eta_k, Lambda_k))
  if (length(tt) == 1) return(as.numeric(Sk[, 1]))
  colnames(Sk) <- paste0("t_", signif(tt, 4))
  Sk
}

predict_survival_mixture <- function(fit, X_cox_new, X_gmm_new = NULL,
                                     tau_new = NULL, times) {
  X_cox_new <- as.matrix(X_cox_new)
  if (is.null(tau_new)) {
    if (is.null(X_gmm_new)) stop("Provide X_gmm_new or tau_new.")
    tau_new <- predict_tau(fit, X_gmm_new)
  }

  tt <- as.numeric(times)
  K <- ncol(tau_new)
  n <- nrow(X_cox_new)
  m <- length(tt)
  Smix <- matrix(0, n, m)

  for (k in seq_len(K)) {
    Sk <- predict_survival_component(fit, X_cox_new, k = k, times = tt)
    if (is.vector(Sk)) Sk <- matrix(Sk, nrow = n, ncol = m)
    Smix <- Smix + matrix(tau_new[, k], nrow = n, ncol = m) * Sk
  }

  if (m == 1) return(as.numeric(Smix[, 1]))
  colnames(Smix) <- paste0("t_", signif(tt, 4))
  Smix
}

predict_event_prob_horizon <- function(fit, X_cox_new, X_gmm_new = NULL,
                                       tau_new = NULL, t_eval) {
  1 - predict_survival_mixture(
    fit = fit,
    X_cox_new = X_cox_new,
    X_gmm_new = X_gmm_new,
    tau_new = tau_new,
    times = t_eval
  )
}

predict_risk <- function(fit, X_cox_new, X_gmm_new = NULL, tau_new = NULL,
                         type = c("event_prob", "soft_lp", "hard_lp"),
                         t_eval = NULL) {
  type <- match.arg(type)
  X_cox_new <- as.matrix(X_cox_new)

  if (is.null(tau_new)) {
    if (is.null(X_gmm_new)) stop("Provide X_gmm_new or tau_new for prediction.")
    tau_new <- predict_tau(fit, X_gmm_new)
  }

  if (type == "event_prob") {
    if (is.null(t_eval)) stop("t_eval is required when type = 'event_prob'.")
    return(predict_event_prob_horizon(
      fit = fit,
      X_cox_new = X_cox_new,
      tau_new = tau_new,
      t_eval = t_eval
    ))
  }

  K <- ncol(tau_new)
  eta_mat <- sapply(seq_len(K), function(k) {
    eta_from_scaled(X_cox_new, fit$coxfit[[k]]$beta_s, fit$coxfit[[k]]$scaler, eta_cap = 12)
  })
  if (is.vector(eta_mat)) eta_mat <- matrix(eta_mat, ncol = K)

  if (type == "soft_lp") {
    return(as.numeric(rowSums(tau_new * eta_mat)))
  }

  k_hat <- max.col(tau_new, ties.method = "first")
  eta <- numeric(nrow(X_cox_new))
  for (k in seq_len(K)) {
    idx <- which(k_hat == k)
    if (length(idx) > 0) eta[idx] <- eta_mat[idx, k]
  }
  eta
}

###############################################################
## Held-out mixture log-likelihood
## log p(t,delta | x) = log sum_c tau_c(x) f_c(t,delta | x)
###############################################################
heldout_loglik <- function(fit, X_gmm_new, X_cox_new, time_new, status_new,
                           tau_new = NULL) {
  X_gmm_new  <- as.matrix(X_gmm_new)
  X_cox_new  <- as.matrix(X_cox_new)
  time_new   <- as.numeric(time_new)
  status_new <- as.integer(status_new)

  K <- fit$K

  # ── WHY tau-based, not full joint density ─────────────────────────────────
  # Correct CV metric: marginal survival log-likelihood weighted by GMM-informed tau.
  #
  #   heldout_loglik = sum_i log[ sum_k tau_k(x_i) * f_Cox_k(t_i, d_i | x_i) ]
  #
  # where tau_k(x_i) = predict_tau(fit, X_gmm_new) uses ONLY the GMM posterior
  # (no Cox information at test time -- correct for prediction setting).
  #
  # Why not full joint density log(pi_k) + (1/q)*log f_GMM_k + log f_Cox_k?
  # With q=33 and only 2 features carrying separation signal, the 31 noise
  # features produce a GMM term of magnitude ~31*(1/33)*log f_noise ≈ -30,
  # which swamps the Cox term (~-2 per subject). The K=2 model cannot overcome
  # this noise and K=1 always wins -- K selection fails entirely.
  #
  # The tau-based formula avoids this: the GMM information enters only through
  # tau (the soft cluster assignment), which concentrates the signal from the
  # truly informative features. Cross-validation then penalises overfitted K
  # through worse held-out Cox fit on the test fold.
  #
  # Historical note: the original code had this correct formula but was broken
  # because lambda=2.9 shrunk all Cox betas to ~0, making f_Cox_1 ≈ f_Cox_2
  # for all K. Fixed by setting lambda=0.05. With correct lambda AND GMM-only
  # predict_tau, this formula correctly discriminates K.
  # ──────────────────────────────────────────────────────────────────────────

  if (is.null(tau_new)) {
    tau_new <- predict_tau(fit, X_gmm_new)   # GMM-only posterior (correct)
  }

  K <- ncol(tau_new)
  lmat <- matrix(NA_real_, nrow(X_cox_new), K)

  for (k in seq_len(K)) {
    lsurv_k <- log_surv_density(time_new, status_new, X_cox_new, fit$coxfit[[k]])
    lmat[, k] <- log(pmax(tau_new[, k], 1e-300)) + lsurv_k
  }

  sum(rowLogSumExp(lmat))
}

c_index <- function(time, status, risk_score) {
  # survival::concordance(Surv ~ LP)$concordance returns:
  #   P(LP_i > LP_j | T_i > T_j)  — "larger predictor → longer survival"
  # For Cox LP (higher LP = higher hazard = shorter time), a GOOD model
  # gives this value close to 0.  Harrell's C is 1 − concordance.
  cc <- survival::concordance(Surv(time, status) ~ risk_score)
  1 - as.numeric(cc$concordance)
}

###############################################################
## Joint CV for K and gamma (= surv_weight)
## metric = "cindex" uses predicted event probability by t_eval
## metric = "heldout_loglik" uses held-out mixture log-likelihood
###############################################################
cv_select_K_gamma <- function(X_gmm, X_cox, time, status,
                              K_grid = 1:3,
                              gamma_grid = c(0.5, 1),
                              nfolds = 3,
                              lambda = NULL,
                              alpha = 0,
                              max_iter = 60,
                              tol = 1e-4,
                              metric = c("cindex", "heldout_loglik"),
                              risk_type = c("event_prob", "soft_lp", "hard_lp"),
                              t_eval = NULL,
                              n_starts = 10,
                              gmm_diag = TRUE,
                              normalize_gmm_by_dim = TRUE,
                              pi_floor = 0.02,
                              min_eff_events = 0.5,
                              init_method = c("supervised", "kmeans"),
                              init_eta_weight = 0.75,
                              seed = 1,
                              verbose = FALSE) {

  metric <- match.arg(metric)
  risk_type <- match.arg(risk_type)
  init_method <- match.arg(init_method)

  X_gmm <- as.matrix(X_gmm)
  X_cox <- as.matrix(X_cox)
  n <- nrow(X_gmm)
  stopifnot(nrow(X_cox) == n, length(time) == n, length(status) == n)

  fold <- make_stratified_folds(status, nfolds = nfolds, seed = seed)

  grid <- expand.grid(
    K = K_grid,
    gamma = gamma_grid,
    stringsAsFactors = FALSE
  )

  fold_records <- vector("list", nrow(grid))
  summary_df <- data.frame(
    K = grid$K,
    gamma = grid$gamma,
    mean_score = NA_real_,
    sd_score = NA_real_,
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

      lambda_fold <- if (is.null(lambda) || is.character(lambda)) {
        tryCatch(
          select_lambda(time[tr], status[tr], X_cox[tr, , drop = FALSE], alpha = alpha),
          error = function(e) 0.1
        )
      } else {
        as.numeric(lambda)
      }

      t_eval_fold <- t_eval
      if (is.null(t_eval_fold)) {
        ev_train <- time[tr & status == 1]
        t_eval_fold <- if (length(ev_train) > 0) stats::median(ev_train) else stats::median(time[tr])
      }

      fit <- tryCatch(
        gemcox_full_multistart(
          X_gmm = X_gmm[tr, , drop = FALSE],
          X_cox = X_cox[tr, , drop = FALSE],
          time = time[tr],
          status = status[tr],
          K = K_cur,
          lambda = lambda_fold,
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
        predict_tau(fit, X_gmm[te, , drop = FALSE]),
        error = function(e) e
      )
      if (inherits(tau_te, "error") || is.null(tau_te)) {
        reasons[v] <- if (inherits(tau_te, "error")) paste("tau_error:", tau_te$message) else "tau_null"
        next
      }

      if (metric == "heldout_loglik") {
        sc <- tryCatch(
          heldout_loglik(
            fit = fit,
            X_gmm_new = X_gmm[te, , drop = FALSE],
            X_cox_new = X_cox[te, , drop = FALSE],
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
            X_cox_new = X_cox[te, , drop = FALSE],
            tau_new = tau_te,
            type = risk_type,
            t_eval = t_eval_fold
          ),
          error = function(e) e
        )
        if (inherits(risk_te, "error") || is.null(risk_te)) {
          reasons[v] <- if (inherits(risk_te, "error")) paste("risk_error:", risk_te$message) else "risk_null"
          next
        }
        risk_te <- as.numeric(risk_te) + stats::rnorm(length(risk_te), sd = 1e-8)
        sc <- tryCatch(c_index(time[te], status[te], risk_te), error = function(e) NA_real_)
      }

      scores[v] <- sc
      reasons[v] <- ifelse(is.finite(sc), "ok", paste0(metric, "_na"))
    }

    ok <- is.finite(scores)
    summary_df$mean_score[g] <- if (any(ok)) mean(scores[ok]) else NA_real_
    summary_df$sd_score[g] <- if (sum(ok) >= 2) stats::sd(scores[ok]) else NA_real_
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

  best_idx <- which.max(ifelse(is.finite(summary_df$mean_score), summary_df$mean_score, -Inf))
  best <- if (length(best_idx) == 1 && is.finite(summary_df$mean_score[best_idx])) {
    summary_df[best_idx, , drop = FALSE]
  } else {
    NULL
  }

  list(
    metric = metric,
    summary = summary_df,
    folds = fold_records,
    best = best
  )
}

###############################################################
## Backward-compatible wrapper: select K only by C-index
###############################################################
cv_select_K_cindex <- function(X_gmm, X_cox, time, status,
                               K_grid = 1:3,
                               nfolds = 3,
                               lambda = NULL,
                               alpha = 0,
                               max_iter = 60,
                               tol = 1e-4,
                               surv_weight = 1.0,
                               risk_type = c("event_prob", "soft_lp", "hard_lp"),
                               t_eval = NULL,
                               n_starts = 10,
                               gmm_diag = TRUE,
                               normalize_gmm_by_dim = TRUE,
                               pi_floor = 0.02,
                               min_eff_events = 0.5,
                               init_method = c("supervised", "kmeans"),
                               init_eta_weight = 0.75,
                               seed = 1,
                               verbose = FALSE) {
  cv_select_K_gamma(
    X_gmm = X_gmm,
    X_cox = X_cox,
    time = time,
    status = status,
    K_grid = K_grid,
    gamma_grid = surv_weight,
    nfolds = nfolds,
    lambda = lambda,
    alpha = alpha,
    max_iter = max_iter,
    tol = tol,
    metric = "cindex",
    risk_type = match.arg(risk_type),
    t_eval = t_eval,
    n_starts = n_starts,
    gmm_diag = gmm_diag,
    normalize_gmm_by_dim = normalize_gmm_by_dim,
    pi_floor = pi_floor,
    min_eff_events = min_eff_events,
    init_method = match.arg(init_method),
    init_eta_weight = init_eta_weight,
    seed = seed,
    verbose = verbose
  )
}

###############################################################
## Final fit helper after CV
###############################################################
fit_best_gemcox <- function(X_gmm, X_cox, time, status,
                            cv_result,
                            lambda = NULL,
                            alpha = 0,
                            max_iter = 60,
                            tol = 1e-4,
                            n_starts = 10,
                            gmm_diag = TRUE,
                            normalize_gmm_by_dim = TRUE,
                            pi_floor = 0.02,
                            min_eff_events = 0.5,
                            init_method = c("supervised", "kmeans"),
                            init_eta_weight = 0.75,
                            seed = 1,
                            verbose = TRUE) {
  init_method <- match.arg(init_method)

  if (is.null(cv_result$best) || nrow(cv_result$best) != 1) {
    stop("cv_result$best is missing; cannot fit final model.")
  }

  K_best <- cv_result$best$K[1]
  gamma_best <- cv_result$best$gamma[1]

  lambda_use <- if (is.null(lambda) || is.character(lambda)) {
    select_lambda(time, status, X_cox, alpha = alpha)
  } else {
    as.numeric(lambda)
  }

  gemcox_full_multistart(
    X_gmm = X_gmm,
    X_cox = X_cox,
    time = time,
    status = status,
    K = K_best,
    lambda = lambda_use,
    alpha = alpha,
    max_iter = max_iter,
    tol = tol,
    verbose = verbose,
    surv_weight = gamma_best,
    gmm_diag = gmm_diag,
    normalize_gmm_by_dim = normalize_gmm_by_dim,
    pi_floor = pi_floor,
    min_eff_events = min_eff_events,
    init_method = init_method,
    init_eta_weight = init_eta_weight,
    n_starts = n_starts,
    init_seeds = seed + seq_len(n_starts)
  )
}
