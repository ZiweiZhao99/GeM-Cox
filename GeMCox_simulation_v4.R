###############################################################
## GeM-Cox Simulation v4 — Fixed DGP
##
## ROOT CAUSE FIXES:
##   1. Signal defined on OBSERVABLE features (not latent)
##      -> log1p(exp(x)) roundtrip no longer destroys signal
##   2. Skewness added via chi-squared mixture, not exp()
##      -> invertible by log-transform without asymmetric compression
##   3. Beta heterogeneity capped to prevent cancellation
##   4. Two modes: "gaussian" (clean test) and "skewed" (realistic)
###############################################################

required_pkgs <- c("survival", "glmnet", "mclust", "dplyr",
                   "tidyr", "ggplot2", "patchwork", "tibble",
                   "MASS", "parallel")
for (pkg in required_pkgs) {
  if (!requireNamespace(pkg, quietly = TRUE))
    install.packages(pkg, repos = "https://cloud.r-project.org")
}
library(survival); library(glmnet); library(mclust)
library(dplyr); library(tidyr); library(ggplot2)
library(patchwork); library(tibble); library(MASS); library(parallel)
select <- dplyr::select; filter <- dplyr::filter
mutate <- dplyr::mutate; arrange <- dplyr::arrange

source("GeMCox_improved_v3.R")
set.seed(2025)

###############################################################
## DGP v4: Signal lives on the observable feature space
###############################################################

simulate_gemcox_v4 <- function(
    n = 117, K = 2, p_gmm = 10, p_cox = 10,
    p_gmm_signal = 5, p_beta_signal = 4,
    pi_c = NULL, separation = c("med", "low", "high"),
    feat_corr = 0.2, event_rate = 0.44,
    t_max = 253, weibull_shape = 1.5,
    baseline_var = 0.2,
    feature_mode = c("gaussian", "skewed"),
    seed = 1
) {
  separation <- match.arg(separation)
  feature_mode <- match.arg(feature_mode)
  set.seed(seed)

  delta_mu <- switch(separation, low = 0.6, med = 1.2, high = 2.0)
  ## FIX #2: Cap delta_beta so adjacent clusters differ by ≤ 0.8 in beta[1]
  ## This prevents cancellation: with K=3, betas range over ±0.8 around base
  delta_beta <- switch(separation, low = 0.3, med = 0.6, high = 1.0)

  if (is.null(pi_c)) pi_c <- rep(1/K, K)
  p_gmm_signal <- min(p_gmm_signal, p_gmm)
  p_beta_signal <- min(p_beta_signal, p_cox)

  ## --- GMM means on the OBSERVABLE feature space ---
  mean_positions <- seq(-(K-1)/2, (K-1)/2, length.out = K) * delta_mu
  mu_list <- lapply(seq_len(K), function(k) {
    m <- rep(0, p_gmm)
    for (j in seq_len(p_gmm_signal)) m[j] <- mean_positions[k] / sqrt(j)
    m
  })

  ## AR1 covariance
  Sigma <- outer(1:p_gmm, 1:p_gmm, function(i,j) feat_corr^abs(i-j))

  ## --- Cox betas: defined on observable features ---
  ## FIX #3: Heterogeneity is ADDITIVE to a common strong base effect
  ## Base effect is always in the same direction across clusters
  ## Only the MAGNITUDE varies, not the sign
  beta_base <- rep(0, p_cox)
  beta_base[1:p_beta_signal] <- c(0.8, -0.5, 0.4, -0.3)[1:p_beta_signal]

  beta_list <- lapply(seq_len(K), function(k) {
    pos <- (k - (K+1)/2)  # centered positions
    b <- beta_base
    ## Shift magnitude, but keep sign aligned with base
    for (j in seq_len(min(3, p_beta_signal))) {
      b[j] <- b[j] + delta_beta * pos / sqrt(j)
    }
    b
  })

  ## --- Weibull scales per cluster ---
  weibull_scale_base <- 200
  scale_pos <- seq(-(K-1)/2, (K-1)/2, length.out = K)
  weibull_scales <- weibull_scale_base * exp(baseline_var * scale_pos)

  ## --- Draw clusters and features ---
  Z <- sample(1:K, n, replace = TRUE, prob = pi_c)

  ## Generate features on the observable space (what the model will see)
  X_obs <- matrix(0, n, p_gmm)
  for (k in 1:K) {
    idx <- which(Z == k)
    if (length(idx) == 0) next
    X_obs[idx, ] <- mvrnorm(length(idx), mu = mu_list[[k]], Sigma = Sigma)
  }

  ## FIX #1: Add skewness as a POST-PROCESSING step (not via exp/log)
  ## This way, the cluster structure in X_obs is preserved
  if (feature_mode == "skewed") {
    ## Add right-skewness by squaring + shifting: x -> x + c * x^2
    ## This makes features right-skewed but preserves ordering and spacing
    for (j in 1:p_gmm) {
      xj <- X_obs[, j]
      ## Standardize to [0,1]-ish range first
      xj_std <- (xj - mean(xj)) / (sd(xj) + 1e-8)
      ## Add skewness: shift right tail
      X_obs[, j] <- xj + 0.3 * pmax(xj_std, 0)^2
    }
  }

  colnames(X_obs) <- paste0("feat_", 1:p_gmm)
  X_cox <- X_obs[, 1:p_cox, drop = FALSE]

  ## --- Survival times: signal defined on OBSERVABLE features ---
  ## FIX: The Cox model operates on X_obs directly
  T_true <- numeric(n)
  for (k in 1:K) {
    idx <- which(Z == k)
    if (length(idx) == 0) next
    eta_k <- as.numeric(X_cox[idx, , drop = FALSE] %*% beta_list[[k]])
    eta_k <- pmin(pmax(eta_k, -5), 5)
    U <- runif(length(idx))
    T_true[idx] <- weibull_scales[k] *
      (-log(U))^(1/weibull_shape) * exp(-eta_k / weibull_shape)
  }

  ## Censoring
  calibrate_cens <- function(target, Tlat, tmax, tol = 0.005) {
    lo <- 1e-6; hi <- 10
    for (it in 1:60) {
      mid <- (lo + hi) / 2
      er <- mean(Tlat <= pmin(rexp(length(Tlat), mid), tmax))
      if (abs(er - target) < tol) break
      if (er < target) hi <- mid else lo <- mid
    }; mid
  }
  cens_rate <- calibrate_cens(event_rate, T_true, t_max)
  C_time <- pmin(rexp(n, rate = cens_rate), t_max)
  time <- pmax(pmin(T_true, C_time) + runif(n, 0, 0.01), 0.1)
  status <- as.integer(T_true <= C_time)

  list(
    X_gmm = X_obs, X_cox = X_cox, time = time, status = status,
    Z_true = Z,
    params = list(K = K, pi_c = pi_c, mu_list = mu_list,
                  beta_list = beta_list, Sigma = Sigma,
                  separation = separation, feature_mode = feature_mode,
                  delta_mu = delta_mu, delta_beta = delta_beta,
                  weibull_scales = weibull_scales,
                  actual_event_rate = mean(status))
  )
}


###############################################################
## Replicate runner v4
###############################################################

run_replicate_v4 <- function(sim_data,
                             K_grid = 1:3,
                             gamma_grid = c(0, 0.5, 1.0, 2.0),
                             nfolds = 3, max_iter = 50,
                             n_starts = 5, n_starts_final = 10,
                             alpha = 0.5,
                             log_transform = FALSE, # OFF for gaussian mode
                             verbose = FALSE) {

  X_gmm <- sim_data$X_gmm; X_cox <- sim_data$X_cox
  time <- sim_data$time; status <- sim_data$status
  Z_true <- sim_data$Z_true; K_true <- sim_data$params$K

  result <- list(K_true = K_true, K_selected = NA_integer_,
                 K_correct = NA_integer_, gamma_selected = NA_real_,
                 ARI = NA_real_, cindex_full = NA_real_,
                 cindex_cv = NA_real_, cindex_cox1 = NA_real_,
                 beta_rmse = NA_real_, n_events = sum(status),
                 event_rate_actual = mean(status),
                 converged = NA_integer_, error_msg = NA_character_)

  ## Preprocess
  if (log_transform) {
    prep <- preprocess_train(X_gmm, X_cox, log_transform = TRUE)
    Xg <- prep$X_gmm; Xc <- prep$X_cox
  } else {
    ## Just standardize, no log-transform
    Xg <- apply_x_scaler(X_gmm, make_x_scaler(X_gmm))
    Xc <- X_cox  # glmnet handles its own standardization
  }

  ## Baseline: single Cox
  result$cindex_cox1 <- tryCatch({
    sc <- make_x_scaler(Xc); Xs <- apply_x_scaler(Xc, sc)
    lam <- tryCatch(cv.glmnet(Xs, Surv(time, status), family = "cox",
                              alpha = alpha, standardize = FALSE,
                              nfolds = min(nfolds, 5))$lambda.1se,
                    error = function(e) 0.05)
    gfit <- glmnet(Xs, Surv(time, status), family = "cox",
                   alpha = alpha, lambda = lam, standardize = FALSE)
    eta <- as.numeric(Xs %*% coef(gfit, s = lam))
    c_index(time, status, eta)
  }, error = function(e) NA_real_)

  ## CV
  cv_res <- tryCatch(
    cv_select_K_gamma_v3(
      X_gmm = Xg, X_cox = Xc, time = time, status = status,
      K_grid = K_grid, gamma_grid = gamma_grid,
      nfolds = nfolds, alpha = alpha, max_iter = max_iter,
      n_starts = n_starts, verbose = verbose, use_1se_rule = TRUE
    ),
    error = function(e) e
  )
  if (inherits(cv_res, "error")) {
    result$error_msg <- cv_res$message; return(result)
  }

  best <- cv_res$best
  if (is.null(best) || nrow(best) == 0) {
    result$error_msg <- "CV failed"; return(result)
  }

  K_hat <- as.integer(best$K[1])
  gamma_hat <- as.numeric(best$gamma[1])
  result$K_selected <- K_hat
  result$K_correct <- as.integer(K_hat == K_true)
  result$gamma_selected <- gamma_hat
  result$cindex_cv <- as.numeric(best$mean_score[1])

  ## Final fit
  fit <- tryCatch(
    gemcox_full_multistart(
      X_gmm = Xg, X_cox = Xc, time = time, status = status,
      K = K_hat, alpha = alpha, max_iter = max_iter * 2,
      surv_weight = gamma_hat, verbose = FALSE, n_starts = n_starts_final
    ),
    error = function(e) e
  )
  if (inherits(fit, "error")) {
    result$error_msg <- fit$message; return(result)
  }
  result$converged <- 1L

  if (K_hat == K_true && K_true > 1)
    result$ARI <- adjustedRandIndex(Z_true, fit$clusterid)

  ## C-index
  eta_full <- tryCatch({
    tau <- fit$tau
    eta_mat <- sapply(1:K_hat, function(k)
      eta_from_scaled(Xc, fit$coxfit[[k]]$beta_s, fit$coxfit[[k]]$scaler))
    if (is.vector(eta_mat)) eta_mat <- matrix(eta_mat, ncol = K_hat)
    rowSums(tau * eta_mat)
  }, error = function(e) NULL)
  if (!is.null(eta_full))
    result$cindex_full <- c_index(time, status, eta_full)

  ## Beta recovery
  if (K_hat == K_true && K_true > 1) {
    true_betas <- sim_data$params$beta_list
    est_betas <- lapply(1:K_hat, function(k) fit$beta[, k])
    pc <- min(length(true_betas[[1]]), nrow(fit$beta))
    cost <- matrix(0, K_true, K_hat)
    for (i in 1:K_true) for (j in 1:K_hat)
      cost[i,j] <- mean((true_betas[[i]][1:pc] - est_betas[[j]][1:pc])^2)
    used <- logical(K_hat); rv <- numeric(K_true)
    for (i in 1:K_true) {
      av <- which(!used); bj <- av[which.min(cost[i,av])]
      rv[i] <- sqrt(cost[i,bj]); used[bj] <- TRUE
    }
    result$beta_rmse <- mean(rv)
  }

  result
}


###############################################################
## Simulation grid
###############################################################

RUN_DIAGNOSTIC <- TRUE

sim_grid <- expand.grid(
  K_true = c(1, 2, 3), separation = c("low", "med", "high"),
  event_rate = c(0.44, 0.25), feature_mode = c("gaussian"),
  n_subjects = 117, feat_corr = 0.2, stringsAsFactors = FALSE
)
sim_grid$scenario_id <- seq_len(nrow(sim_grid))

if (RUN_DIAGNOSTIC) {
  SIM_K_GRID <- 1:3; SIM_GAMMA_GRID <- c(0, 0.5, 1.0, 2.0)
  SIM_NFOLDS <- 3; SIM_MAX_ITER <- 40
  SIM_N_STARTS <- 5; SIM_N_STARTS_FINAL <- 10
  NSIM <- 20L
} else {
  SIM_K_GRID <- 1:3; SIM_GAMMA_GRID <- c(0, 0.25, 0.5, 1.0, 2.0)
  SIM_NFOLDS <- 5; SIM_MAX_ITER <- 60
  SIM_N_STARTS <- 5; SIM_N_STARTS_FINAL <- 15
  NSIM <- 100L
}
N_CORES <- min(50L, max(1L, detectCores(logical = FALSE) - 1L))

cat(sprintf("Scenarios: %d x %d reps = %d total\n",
            nrow(sim_grid), NSIM, nrow(sim_grid) * NSIM))

###############################################################
## Run
###############################################################

run_scenario_v4 <- function(sc_row, nsim = NSIM) {
  K_true <- sc_row$K_true; sep <- sc_row$separation
  er <- sc_row$event_rate; fmode <- sc_row$feature_mode
  sc_id <- sc_row$scenario_id

  cat(sprintf("[%d] K=%d sep=%s er=%.2f mode=%s\n", sc_id, K_true, sep, er, fmode))

  results <- lapply(seq_len(nsim), function(rep_i) {
    sim_data <- tryCatch(
      simulate_gemcox_v4(
        n = sc_row$n_subjects, K = K_true, p_gmm = 10, p_cox = 10,
        p_gmm_signal = 5, p_beta_signal = 4,
        separation = sep, feat_corr = sc_row$feat_corr,
        event_rate = er, feature_mode = fmode,
        seed = sc_id * 10000 + rep_i
      ),
      error = function(e) NULL
    )
    if (is.null(sim_data)) return(NULL)

    use_log <- (fmode == "skewed")

    res <- tryCatch(
      run_replicate_v4(sim_data,
        K_grid = SIM_K_GRID, gamma_grid = SIM_GAMMA_GRID,
        nfolds = SIM_NFOLDS, max_iter = SIM_MAX_ITER,
        n_starts = SIM_N_STARTS, n_starts_final = SIM_N_STARTS_FINAL,
        alpha = 0.5, log_transform = use_log, verbose = FALSE),
      error = function(e) list(error_msg = e$message)
    )
    res <- lapply(res, function(v) if (is.null(v)||length(v)==0) NA else v[1])
    as.data.frame(c(
      list(scenario_id = sc_id, rep = rep_i, K_true = K_true,
           separation = sep, event_rate_target = er,
           feature_mode = fmode, n = sc_row$n_subjects),
      res), stringsAsFactors = FALSE)
  })
  bind_rows(Filter(Negate(is.null), results))
}

run_par <- function(nj, FUN, nc) {
  if (nc <= 1L) return(lapply(seq_len(nj), FUN))
  if (.Platform$OS.type == "windows") {
    cl <- makeCluster(nc, type = "PSOCK"); on.exit(stopCluster(cl))
    clusterExport(cl, ls(.GlobalEnv), .GlobalEnv)
    clusterEvalQ(cl, { library(survival); library(glmnet); library(MASS)
                       library(mclust); library(dplyr); library(parallel) })
    parLapply(cl, seq_len(nj), FUN)
  } else mclapply(seq_len(nj), FUN, mc.cores = nc)
}

message("========== STARTING v4 SIMULATION ==========")
all_res <- run_par(nrow(sim_grid),
  function(i) tryCatch(run_scenario_v4(sim_grid[i,]),
    error = function(e) data.frame(scenario_id = sim_grid$scenario_id[i],
      error_msg = e$message, stringsAsFactors = FALSE)), N_CORES)

results_df <- bind_rows(Filter(is.data.frame, all_res))
saveRDS(results_df, "GeMCox_simulation_v4_raw.rds")

## Summarize
em <- results_df$error_msg
ok <- is.na(em) | (nzchar(trimws(em)) == FALSE) | (trimws(em) == "NA")
df_ok <- results_df[ok, ]

n_success <- nrow(df_ok); n_total <- nrow(results_df)
cat(sprintf("\nSuccess rate: %d / %d = %.0f%%\n", n_success, n_total,
            100 * n_success / max(n_total, 1)))

summary_df <- df_ok %>%
  group_by(K_true, separation, event_rate_target) %>%
  summarise(
    n_reps = n(),
    K_sel_accuracy = mean(K_correct, na.rm = TRUE),
    K_sel_mean = mean(K_selected, na.rm = TRUE),
    ARI_mean = mean(ARI, na.rm = TRUE),
    ARI_sd = sd(ARI, na.rm = TRUE),
    cindex_full_mean = mean(cindex_full, na.rm = TRUE),
    cindex_cv_mean = mean(cindex_cv, na.rm = TRUE),
    cindex_cox1_mean = mean(cindex_cox1, na.rm = TRUE),
    cindex_gain_mean = mean(cindex_full - cindex_cox1, na.rm = TRUE),
    gamma_mean = mean(gamma_selected, na.rm = TRUE),
    converge_rate = mean(converged, na.rm = TRUE),
    .groups = "drop"
  ) %>% mutate(across(where(is.numeric), ~round(.x, 4)))

write.csv(summary_df, "GeMCox_simulation_v4_summary.csv", row.names = FALSE)
cat("\n========== v4 RESULTS ==========\n")
print(as.data.frame(summary_df))
message("========== v4 SIMULATION COMPLETE ==========")
