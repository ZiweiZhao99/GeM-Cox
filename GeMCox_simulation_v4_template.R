###############################################################
## GeM-Cox Simulation v4 template
##
## Purpose:
##   - Move the simulation closer to CVIA078
##   - Avoid K=4 by default for n = 117
##   - Separate clustering features from Cox features
##   - Optionally impose CVIA078-like missingness and tied follow-up times
##   - Use the fold-safe / shared-baseline pipeline in
##     GeMCox_improved_v4_patch.R
##
## Typical workflow:
##   1) Read the real CVIA078 long-format data
##   2) Build an empirical template from interpretable features
##   3) Simulate transformed-scale datasets from that template
##   4) Fit GeM-Cox with log_transform = FALSE inside simulation,
##      because the features are already on the transformed scale
###############################################################

required_pkgs <- c("survival", "glmnet", "mclust", "MASS", "dplyr", "tidyr")
for (pkg in required_pkgs) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg, repos = "https://cloud.r-project.org")
  }
}

library(survival)
library(glmnet)
library(mclust)
library(MASS)
library(dplyr)
library(tidyr)

source("GeMCox_improved_v4_patch.R")

###############################################################
## 1. Build an empirical CVIA078 template
###############################################################

build_cvia078_template_v4 <- function(dat_long,
                                      gmm_feature_mode = "delta_only",
                                      cox_feature_mode = "v20_plus_delta_covars",
                                      missing_threshold = 0.05,
                                      corr_prune = 0.85) {
  real_dat <- prepare_cvia078_data_v4(
    dat_long,
    gmm_feature_mode = gmm_feature_mode,
    cox_feature_mode = cox_feature_mode,
    missing_threshold = missing_threshold,
    corr_prune = corr_prune
  )
  
  ## Put data on the modeling scale used inside CV folds, then keep that scale
  ## for simulation to avoid an unnecessary transform/inverse-transform cycle.
  prep_scale <- preprocess_train_v4(
    real_dat$X_gmm,
    real_dat$X_cox,
    log_transform = TRUE,
    skew_threshold = 1.0,
    impute = "median"
  )
  
  Xg <- prep_scale$X_gmm
  Xc <- prep_scale$X_cox
  
  cov_ridge <- function(S, ridge = 1e-4) {
    S <- as.matrix(S)
    S[!is.finite(S)] <- 0
    S <- (S + t(S)) / 2
    S + diag(ridge, nrow(S))
  }
  
  list(
    X_gmm_center = colMeans(Xg),
    X_gmm_cov = cov_ridge(stats::cov(Xg)),
    X_cox_center = colMeans(Xc),
    X_cox_cov = cov_ridge(stats::cov(Xc)),
    gmm_features = colnames(Xg),
    cox_features = colnames(Xc),
    gmm_missing_rate = colMeans(!is.finite(as.matrix(real_dat$X_gmm))),
    cox_missing_rate = colMeans(!is.finite(as.matrix(real_dat$X_cox))),
    time_support = sort(unique(real_dat$time)),
    time_event_support = sort(unique(real_dat$time[real_dat$status == 1])),
    n = nrow(Xg),
    event_rate = mean(real_dat$status),
    real_data = real_dat,
    prep_scale = prep_scale
  )
}

###############################################################
## 2. CVIA078-like simulator
###############################################################

nearest_support_v4 <- function(x, support) {
  if (length(support) == 0) return(x)
  vapply(x, function(z) support[which.min(abs(support - z))], numeric(1))
}

simulate_cvia078_like_v4 <- function(template,
                                     n = template$n,
                                     K_true = 2,
                                     pi_c = NULL,
                                     separation = c("low", "med", "high"),
                                     signal_gmm = NULL,
                                     signal_cox = NULL,
                                     delta_mu = NULL,
                                     delta_beta = NULL,
                                     event_rate_target = template$event_rate,
                                     baseline_shape = 1.5,
                                     baseline_mode = c("shared", "mild_cluster_shift"),
                                     add_missing = TRUE,
                                     seed = 1) {
  separation <- match.arg(separation)
  baseline_mode <- match.arg(baseline_mode)
  set.seed(seed)
  
  if (is.null(pi_c)) pi_c <- rep(1 / K_true, K_true)
  
  if (is.null(delta_mu)) {
    delta_mu <- switch(separation, low = 0.35, med = 0.65, high = 1.0)
  }
  if (is.null(delta_beta)) {
    delta_beta <- switch(separation, low = 0.25, med = 0.45, high = 0.70)
  }
  
  q <- length(template$gmm_features)
  p <- length(template$cox_features)
  if (is.null(signal_gmm)) signal_gmm <- template$gmm_features[seq_len(min(6, q))]
  if (is.null(signal_cox)) signal_cox <- intersect(template$cox_features, signal_gmm)
  if (length(signal_cox) == 0) signal_cox <- template$cox_features[seq_len(min(6, p))]
  
  idx_gmm <- match(signal_gmm, template$gmm_features)
  idx_gmm <- idx_gmm[is.finite(idx_gmm)]
  idx_cox <- match(signal_cox, template$cox_features)
  idx_cox <- idx_cox[is.finite(idx_cox)]
  
  mu0_g <- template$X_gmm_center
  mu0_c <- template$X_cox_center
  Sg <- template$X_gmm_cov
  Sc <- template$X_cox_cov
  
  positions <- seq(-(K_true - 1) / 2, (K_true - 1) / 2, length.out = K_true)
  mu_g_list <- lapply(seq_len(K_true), function(k) {
    mu <- mu0_g
    if (length(idx_gmm) > 0) {
      mu[idx_gmm] <- mu[idx_gmm] + positions[k] * delta_mu / sqrt(seq_along(idx_gmm))
    }
    mu
  })
  
  beta_base <- rep(0, p)
  if (length(idx_cox) > 0) {
    signs <- rep(c(1, -1, 1, -1, 1, -1), length.out = length(idx_cox))
    beta_base[idx_cox] <- signs * seq(0.6, 0.25, length.out = length(idx_cox))
  }
  beta_list <- lapply(seq_len(K_true), function(k) {
    b <- beta_base
    if (length(idx_cox) > 0) {
      b[idx_cox] <- b[idx_cox] + positions[k] * delta_beta / sqrt(seq_along(idx_cox))
    }
    b
  })
  
  Z <- sample(seq_len(K_true), size = n, replace = TRUE, prob = pi_c)
  X_gmm <- matrix(NA_real_, n, q, dimnames = list(NULL, template$gmm_features))
  X_cox <- matrix(NA_real_, n, p, dimnames = list(NULL, template$cox_features))
  
  for (k in seq_len(K_true)) {
    idx <- which(Z == k)
    if (length(idx) == 0) next
    X_gmm[idx, ] <- MASS::mvrnorm(length(idx), mu = mu_g_list[[k]], Sigma = Sg)
    
    ## Generate Cox features around the empirical center, but align shared signal
    ## columns with the same cluster location whenever names overlap.
    Xc_k <- MASS::mvrnorm(length(idx), mu = mu0_c, Sigma = Sc)
    common_nm <- intersect(template$gmm_features, template$cox_features)
    if (length(common_nm) > 0) {
      jg <- match(common_nm, template$gmm_features)
      jc <- match(common_nm, template$cox_features)
      Xc_k[, jc] <- X_gmm[idx, jg, drop = FALSE]
    }
    X_cox[idx, ] <- Xc_k
  }
  
  ## Optional missingness mimicking the real data on the preprocessed scale.
  if (add_missing) {
    for (j in seq_len(q)) {
      miss_j <- stats::rbinom(n, 1, template$gmm_missing_rate[j]) == 1
      if (any(miss_j)) X_gmm[miss_j, j] <- NA_real_
    }
    for (j in seq_len(p)) {
      miss_j <- stats::rbinom(n, 1, template$cox_missing_rate[j]) == 1
      if (any(miss_j)) X_cox[miss_j, j] <- NA_real_
    }
  }
  
  ## Event times on the Cox-feature scale.
  scale_pos <- if (baseline_mode == "shared") rep(0, K_true) else positions * 0.25
  baseline_scale <- 150 * exp(scale_pos)
  
  T_true <- numeric(n)
  for (k in seq_len(K_true)) {
    idx <- which(Z == k)
    if (length(idx) == 0) next
    Xc_obs <- X_cox[idx, , drop = FALSE]
    Xc_imp <- Xc_obs
    for (j in seq_len(ncol(Xc_imp))) {
      z <- Xc_imp[, j]
      z[!is.finite(z)] <- mu0_c[j]
      Xc_imp[, j] <- z
    }
    eta <- as.numeric(Xc_imp %*% beta_list[[k]])
    eta <- pmin(pmax(eta, -4), 4)
    U <- stats::runif(length(idx))
    T_true[idx] <- baseline_scale[k] * (-log(U))^(1 / baseline_shape) * exp(-eta / baseline_shape)
  }
  
  calibrate_cens_v4 <- function(target, Tlat, tol = 0.005) {
    lo <- 1e-6
    hi <- 2
    for (it in 1:60) {
      mid <- (lo + hi) / 2
      C_try <- stats::rexp(length(Tlat), rate = mid)
      er <- mean(Tlat <= C_try)
      if (abs(er - target) < tol) break
      if (er < target) hi <- mid else lo <- mid
    }
    mid
  }
  
  cens_rate <- calibrate_cens_v4(event_rate_target, T_true)
  C_time <- stats::rexp(n, rate = cens_rate)
  time <- pmin(T_true, C_time)
  status <- as.integer(T_true <= C_time)
  
  ## Impose CVIA078-like ties by snapping to observed follow-up support.
  support <- if (length(template$time_support) > 0) template$time_support else sort(unique(round(time)))
  time <- nearest_support_v4(time, support)
  time <- pmax(time, min(support))
  
  list(
    X_gmm = X_gmm,
    X_cox = X_cox,
    time = time,
    status = status,
    Z_true = Z,
    params = list(
      K_true = K_true,
      pi_c = pi_c,
      beta_list = beta_list,
      mu_g_list = mu_g_list,
      delta_mu = delta_mu,
      delta_beta = delta_beta,
      separation = separation,
      event_rate_actual = mean(status),
      baseline_mode = baseline_mode,
      signal_gmm = signal_gmm,
      signal_cox = signal_cox
    )
  )
}

###############################################################
## 3. Single replicate and study runner
###############################################################

run_replicate_cvia078_v4 <- function(sim_data,
                                     K_grid = 1:3,
                                     gamma_grid = c(0, 0.25, 0.5, 1.0),
                                     nfolds = 5,
                                     alpha = 0.5,
                                     max_iter = 60,
                                     n_starts = 5,
                                     n_starts_final = 10,
                                     baseline_mode = "shared",
                                     verbose = FALSE) {
  out <- list(
    K_true = sim_data$params$K_true,
    K_selected = NA_integer_,
    K_correct = NA_integer_,
    gamma_selected = NA_real_,
    ARI = NA_real_,
    cindex_cv = NA_real_,
    cindex_full = NA_real_,
    n_events = sum(sim_data$status),
    event_rate_actual = mean(sim_data$status),
    error_msg = NA_character_
  )
  
  fit <- tryCatch(
    gemcox_pipeline_v4(
      X_gmm = sim_data$X_gmm,
      X_cox = sim_data$X_cox,
      time = sim_data$time,
      status = sim_data$status,
      K_grid = K_grid,
      gamma_grid = gamma_grid,
      alpha = alpha,
      log_transform = FALSE,
      use_1se_rule = TRUE,
      nfolds = nfolds,
      n_starts = n_starts,
      n_starts_final = n_starts_final,
      max_iter = max_iter,
      verbose = verbose,
      baseline_mode = baseline_mode
    ),
    error = function(e) e
  )
  
  if (inherits(fit, "error")) {
    out$error_msg <- fit$message
    return(out)
  }
  
  out$K_selected <- fit$K_selected
  out$K_correct <- as.integer(out$K_selected == out$K_true)
  out$gamma_selected <- fit$gamma_selected
  out$cindex_cv <- fit$cv_result$best$mean_score[1]
  
  tau <- fit$fit$tau
  K_hat <- fit$K_selected
  eta_mat <- sapply(seq_len(K_hat), function(k) {
    eta_from_scaled(fit$preprocessing$X_cox,
                    fit$fit$coxfit[[k]]$beta_s,
                    fit$fit$coxfit[[k]]$scaler)
  })
  if (is.vector(eta_mat)) eta_mat <- matrix(eta_mat, ncol = K_hat)
  out$cindex_full <- c_index(sim_data$time, sim_data$status, rowSums(tau * eta_mat))
  
  if (out$K_correct == 1 && out$K_true > 1) {
    out$ARI <- mclust::adjustedRandIndex(sim_data$Z_true, fit$fit$clusterid)
  }
  
  out
}

run_simulation_study_cvia078_v4 <- function(template,
                                            scenario_grid = NULL,
                                            nsim = 50,
                                            seed = 1,
                                            baseline_mode = "shared") {
  if (is.null(scenario_grid)) {
    scenario_grid <- expand.grid(
      K_true = c(1, 2, 3),
      separation = c("low", "med", "high"),
      event_rate_target = c(template$event_rate, 0.30),
      stringsAsFactors = FALSE
    )
  }
  scenario_grid$scenario_id <- seq_len(nrow(scenario_grid))
  
  all_res <- vector("list", nrow(scenario_grid))
  for (i in seq_len(nrow(scenario_grid))) {
    sc <- scenario_grid[i, ]
    message(sprintf("Scenario %d: K=%d sep=%s er=%.2f",
                    sc$scenario_id, sc$K_true, sc$separation, sc$event_rate_target))
    
    res_i <- lapply(seq_len(nsim), function(r) {
      dat <- simulate_cvia078_like_v4(
        template = template,
        n = template$n,
        K_true = sc$K_true,
        separation = sc$separation,
        event_rate_target = sc$event_rate_target,
        baseline_mode = baseline_mode,
        seed = seed + 1000 * i + r
      )
      ans <- run_replicate_cvia078_v4(dat, baseline_mode = baseline_mode, verbose = FALSE)
      data.frame(
        scenario_id = sc$scenario_id,
        rep = r,
        K_true = sc$K_true,
        separation = sc$separation,
        event_rate_target = sc$event_rate_target,
        ans,
        stringsAsFactors = FALSE
      )
    })
    
    all_res[[i]] <- dplyr::bind_rows(res_i)
  }
  
  raw_df <- dplyr::bind_rows(all_res)
  summary_df <- raw_df |>
    dplyr::group_by(K_true, separation, event_rate_target) |>
    dplyr::summarise(
      n_reps = dplyr::n(),
      K_sel_accuracy = mean(K_correct, na.rm = TRUE),
      K_sel_mean = mean(K_selected, na.rm = TRUE),
      ARI_mean_given_correctK = mean(ARI, na.rm = TRUE),
      cindex_cv_mean = mean(cindex_cv, na.rm = TRUE),
      cindex_full_mean = mean(cindex_full, na.rm = TRUE),
      gamma_mean = mean(gamma_selected, na.rm = TRUE),
      .groups = "drop"
    )
  
  list(raw = raw_df, summary = summary_df, scenario_grid = scenario_grid)
}

cat("GeMCox_simulation_v4_template.R loaded.\n")
cat("  Main idea: use build_cvia078_template_v4() then simulate_cvia078_like_v4().\n")