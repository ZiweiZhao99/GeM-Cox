###############################################################
## GeM-Cox v4 revisions patch
##
## Append this AFTER sourcing GeMCox_improved_v4.R, or replace
## the functions with the same names in that file.
##
## Main changes:
##   1) normalize_gmm_by_dim defaults to FALSE
##   2) kmeans becomes the default initializer
##   3) CV screening rejects fits with tiny clusters / tiny effective events
##   4) selection policy can use screened 1-SE or screened max-mean
##   5) pipeline forwards the new controls
###############################################################

pairwise_beta_sep_v4 <- function(beta_mat) {
  beta_mat <- as.matrix(beta_mat)
  K <- ncol(beta_mat)
  if (K <= 1) return(Inf)
  dmin <- Inf
  for (i in 1:(K - 1)) {
    for (j in (i + 1):K) {
      dij <- sqrt(mean((beta_mat[, i] - beta_mat[, j])^2, na.rm = TRUE))
      if (is.finite(dij)) dmin <- min(dmin, dij)
    }
  }
  dmin
}

fit_screen_v4 <- function(fit, status,
                          min_cluster_prop = 0.08,
                          min_eff_events_per_cluster = 6,
                          min_beta_sep = 0.02,
                          require_nontrivial_K = FALSE) {
  tau <- fit$tau
  K <- ncol(tau)
  n <- nrow(tau)
  prop <- colMeans(tau)
  eff_events <- colSums(tau * as.integer(status))
  beta_sep <- if (!is.null(fit$beta)) pairwise_beta_sep_v4(fit$beta) else Inf

  ok <- TRUE
  reason <- "ok"
  if (any(prop < min_cluster_prop)) {
    ok <- FALSE
    reason <- "small_cluster_prop"
  }
  if (ok && any(eff_events < min_eff_events_per_cluster)) {
    ok <- FALSE
    reason <- "small_eff_events"
  }
  if (ok && require_nontrivial_K && K > 1 && is.finite(beta_sep) && beta_sep < min_beta_sep) {
    ok <- FALSE
    reason <- "beta_not_separated"
  }

  list(
    ok = ok,
    reason = reason,
    min_cluster_prop = min(prop),
    min_eff_events = min(eff_events),
    beta_sep = beta_sep,
    cluster_prop = prop,
    eff_events = eff_events
  )
}

cv_select_K_gamma_v4 <- function(X_gmm, X_cox, time, status,
                                 K_grid = 1:3,
                                 gamma_grid = c(0, 0.25, 0.5, 1.0, 2.0),
                                 nfolds = 5,
                                 alpha = 0.5,
                                 max_iter = 60,
                                 tol = 1e-4,
                                 metric = c("cindex", "heldout_loglik"),
                                 risk_type = c("event_prob", "soft_lp", "hard_lp"),
                                 t_eval = NULL,
                                 n_starts = 8,
                                 gmm_diag = TRUE,
                                 normalize_gmm_by_dim = FALSE,
                                 pi_floor = 0.02,
                                 min_eff_events = 1.0,
                                 init_method = c("kmeans", "supervised"),
                                 init_eta_weight = 0.75,
                                 log_transform = TRUE,
                                 skew_threshold = 1.0,
                                 impute = c("median", "mean"),
                                 seed = 1,
                                 verbose = FALSE,
                                 use_1se_rule = TRUE,
                                 baseline_mode = c("shared", "cluster"),
                                 selection_policy = c("screened_1se", "screened_max", "1se", "max_mean"),
                                 min_cluster_prop = 0.08,
                                 min_eff_events_per_cluster = 6,
                                 min_beta_sep = 0.02,
                                 require_nontrivial_K = FALSE) {

  metric <- match.arg(metric)
  risk_type <- match.arg(risk_type)
  init_method <- match.arg(init_method)
  impute <- match.arg(impute)
  baseline_mode <- match.arg(baseline_mode)
  selection_policy <- match.arg(selection_policy)

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
    min_cluster_prop_mean = NA_real_,
    min_eff_events_mean = NA_real_,
    beta_sep_mean = NA_real_,
    stringsAsFactors = FALSE
  )

  for (g in seq_len(nrow(grid))) {
    K_cur <- grid$K[g]
    gamma_cur <- grid$gamma[g]
    if (verbose) message("CV for K=", K_cur, ", gamma=", gamma_cur)

    scores <- rep(NA_real_, nfolds)
    reasons <- rep(NA_character_, nfolds)
    min_prop_vec <- rep(NA_real_, nfolds)
    min_eff_vec <- rep(NA_real_, nfolds)
    beta_sep_vec <- rep(NA_real_, nfolds)

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

      scr <- fit_screen_v4(
        fit = fit,
        status = status[tr],
        min_cluster_prop = min_cluster_prop,
        min_eff_events_per_cluster = min_eff_events_per_cluster,
        min_beta_sep = min_beta_sep,
        require_nontrivial_K = require_nontrivial_K
      )
      min_prop_vec[v] <- scr$min_cluster_prop
      min_eff_vec[v] <- scr$min_eff_events
      beta_sep_vec[v] <- scr$beta_sep
      if (!scr$ok) {
        reasons[v] <- paste0("screen_fail:", scr$reason)
        next
      }

      tau_te <- tryCatch(predict_tau(fit, prep_te$X_gmm), error = function(e) e)
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
          predict_risk(fit, X_cox_new = prep_te$X_cox, tau_new = tau_te,
                       type = risk_type, t_eval = t_eval_fold),
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
    summary_df$min_cluster_prop_mean[g] <- if (any(is.finite(min_prop_vec))) mean(min_prop_vec, na.rm = TRUE) else NA_real_
    summary_df$min_eff_events_mean[g] <- if (any(is.finite(min_eff_vec))) mean(min_eff_vec, na.rm = TRUE) else NA_real_
    summary_df$beta_sep_mean[g] <- if (any(is.finite(beta_sep_vec))) mean(beta_sep_vec, na.rm = TRUE) else NA_real_

    fold_records[[g]] <- data.frame(
      K = K_cur,
      gamma = gamma_cur,
      fold = seq_len(nfolds),
      score = scores,
      reason = reasons,
      min_cluster_prop = min_prop_vec,
      min_eff_events = min_eff_vec,
      beta_sep = beta_sep_vec,
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
        if (selection_policy %in% c("screened_1se", "1se")) {
          best <- eligible[1, , drop = FALSE]
          best_method <- if (selection_policy == "screened_1se") "screened_1se" else "1se_rule"
        }
      }
    }

    if (selection_policy %in% c("screened_max", "max_mean")) {
      best <- summary_df[best_idx, , drop = FALSE]
      best_method <- if (selection_policy == "screened_max") "screened_max" else "max_mean"
    }
  }

  list(
    metric = metric,
    summary = summary_df,
    folds = fold_records,
    best = best,
    best_method = best_method,
    baseline_mode = baseline_mode,
    selection_policy = selection_policy
  )
}

gemcox_pipeline_v4 <- function(
    X_gmm, X_cox, time, status,
    K_grid = 1:3,
    gamma_grid = c(0, 0.25, 0.5, 1.0, 2.0),
    alpha = 0.5,
    log_transform = TRUE,
    skew_threshold = 1.0,
    impute = c("median", "mean"),
    use_1se_rule = TRUE,
    nfolds = 5,
    n_starts = 8,
    n_starts_final = 12,
    max_iter = 60,
    seed = 1,
    verbose = TRUE,
    baseline_mode = c("shared", "cluster"),
    gmm_diag = TRUE,
    normalize_gmm_by_dim = FALSE,
    init_method = c("kmeans", "supervised"),
    selection_policy = c("screened_1se", "screened_max", "1se", "max_mean"),
    min_cluster_prop = 0.08,
    min_eff_events_per_cluster = 6,
    min_beta_sep = 0.02,
    require_nontrivial_K = FALSE) {

  baseline_mode <- match.arg(baseline_mode)
  impute <- match.arg(impute)
  init_method <- match.arg(init_method)
  selection_policy <- match.arg(selection_policy)

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
    gmm_diag = gmm_diag,
    normalize_gmm_by_dim = normalize_gmm_by_dim,
    init_method = init_method,
    log_transform = log_transform,
    skew_threshold = skew_threshold,
    impute = impute,
    use_1se_rule = use_1se_rule,
    seed = seed,
    verbose = verbose,
    baseline_mode = baseline_mode,
    selection_policy = selection_policy,
    min_cluster_prop = min_cluster_prop,
    min_eff_events_per_cluster = min_eff_events_per_cluster,
    min_beta_sep = min_beta_sep,
    require_nontrivial_K = require_nontrivial_K
  )

  if (is.null(cv_res$best) || nrow(cv_res$best) != 1) stop("CV failed to select a valid model.")

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
    gmm_diag = gmm_diag,
    normalize_gmm_by_dim = normalize_gmm_by_dim,
    init_method = init_method,
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
    baseline_mode = baseline_mode,
    selection_policy = selection_policy
  )
}
