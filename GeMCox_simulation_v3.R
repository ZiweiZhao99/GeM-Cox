###############################################################
## GeM-Cox Simulation v3
## Uses GeMCox_improved_v3.R (sources v2 internally)
##
## Key changes from v2:
##   - Sources v3 (log-transform, shared Breslow, 1-SE rule built in)
##   - DGP generates log-normal features (realistic)
##   - alpha = 0.5 (elastic net) by default
##   - gamma grid includes 0 (unsupervised baseline)
##   - 1-SE rule via cv_select_K_gamma_v3
##   - Preprocessing applied before fitting
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

source("GeMCox_improved_v3.R")   # loads v2 + v3 additions
set.seed(2025)

###############################################################
## DGP v3: log-normal features, block correlation, covariates
###############################################################

simulate_gemcox_v3 <- function(
    n = 117, K = 2, p_gmm = 10, p_cox = 10,
    p_gmm_signal = 5, p_beta_signal = 4,
    pi_c = NULL, separation = c("med", "low", "high"),
    feat_corr = 0.2, event_rate = 0.44,
    t_max = 253, weibull_shape = 1.5,
    baseline_var = 0.3,
    lognormal = TRUE,   # NEW: generate skewed features
    seed = 1
) {
  separation <- match.arg(separation)
  set.seed(seed)

  delta_mu <- switch(separation, low = 0.8, med = 1.5, high = 2.5)
  delta_beta <- switch(separation, low = 0.5, med = 1.0, high = 1.8)

  if (is.null(pi_c)) pi_c <- rep(1/K, K)
  p_gmm_signal <- min(p_gmm_signal, p_gmm)

  ## GMM means: spread across p_gmm_signal features
  mean_positions <- seq(-(K-1)/2, (K-1)/2, length.out = K) * delta_mu
  mu_list <- lapply(seq_len(K), function(k) {
    m <- rep(0, p_gmm)
    for (j in seq_len(p_gmm_signal)) {
      m[j] <- mean_positions[k] / sqrt(j)
    }
    m
  })

  ## AR1 covariance
  ar1_cov <- function(p, rho) {
    S <- outer(1:p, 1:p, function(i,j) rho^abs(i-j)); S
  }
  Sigma <- ar1_cov(p_gmm, feat_corr)

  ## Cox coefficients: NO division by (K-1)
  p_beta_signal <- min(p_beta_signal, p_cox)
  beta_base <- rep(0, p_cox)
  beta_base[1:p_beta_signal] <- c(1.2, -0.8, 0.6, -0.5)[1:p_beta_signal]

  beta_list <- lapply(seq_len(K), function(k) {
    pos <- (k - (K+1)/2)
    b <- beta_base
    for (j in seq_len(min(3, p_beta_signal))) {
      direction <- ifelse(j %% 2 == 1, 1, -1) * sign(beta_base[j])
      b[j] <- b[j] + delta_beta * pos * direction / sqrt(j)
    }
    b
  })

  ## Cluster-specific Weibull scales
  weibull_scale_base <- 200
  scale_positions <- seq(-(K-1)/2, (K-1)/2, length.out = K)
  weibull_scales <- weibull_scale_base * exp(baseline_var * scale_positions)

  ## Draw data
  Z <- sample(1:K, n, replace = TRUE, prob = pi_c)
  X_latent <- matrix(0, n, p_gmm)
  for (k in 1:K) {
    idx <- which(Z == k)
    if (length(idx) == 0) next
    X_latent[idx, ] <- mvrnorm(length(idx), mu = mu_list[[k]], Sigma = Sigma)
  }

  ## Log-normal transform (observable features are skewed)
  if (lognormal) {
    X_gmm <- exp(X_latent)     # right-skewed
  } else {
    X_gmm <- X_latent
  }
  colnames(X_gmm) <- paste0("feat_", 1:p_gmm)
  X_cox <- X_gmm[, 1:p_cox, drop = FALSE]

  ## Survival times (use LATENT Gaussian for true signal)
  T_true <- numeric(n)
  for (k in 1:K) {
    idx <- which(Z == k)
    if (length(idx) == 0) next
    eta_k <- as.numeric(X_latent[idx, 1:p_cox, drop = FALSE] %*% beta_list[[k]])
    eta_k <- pmin(pmax(eta_k, -5), 5)
    U <- runif(length(idx))
    T_true[idx] <- weibull_scales[k] * (-log(U))^(1/weibull_shape) * exp(-eta_k / weibull_shape)
  }

  ## Censoring calibration
  calibrate_cens <- function(target, Tlat, tmax, tol = 0.005) {
    lo <- 1e-6; hi <- 10
    for (it in 1:60) {
      mid <- (lo + hi) / 2
      er <- mean(Tlat <= pmin(rexp(length(Tlat), mid), tmax))
      if (abs(er - target) < tol) break
      if (er < target) hi <- mid else lo <- mid
    }
    mid
  }
  cens_rate <- calibrate_cens(event_rate, T_true, t_max)
  C_time <- pmin(rexp(n, rate = cens_rate), t_max)
  time <- pmin(T_true, C_time) + runif(n, 0, 0.01)
  time <- pmax(time, 0.1)
  status <- as.integer(T_true <= C_time)

  list(
    X_gmm = X_gmm, X_cox = X_cox, time = time, status = status,
    Z_true = Z,
    params = list(K = K, pi_c = pi_c, mu_list = mu_list,
                  beta_list = beta_list, Sigma = Sigma,
                  separation = separation, lognormal = lognormal,
                  delta_mu = delta_mu, delta_beta = delta_beta,
                  weibull_scales = weibull_scales,
                  actual_event_rate = mean(status))
  )
}


###############################################################
## Single-replicate evaluation using v3 pipeline
###############################################################

run_replicate_v3 <- function(sim_data,
                             K_grid = 1:3,
                             gamma_grid = c(0, 0.25, 0.5, 1.0, 2.0),
                             nfolds = 5, max_iter = 50,
                             n_starts = 5, n_starts_final = 10,
                             alpha = 0.5,     # elastic net
                             lambda = NULL,    # NULL = cv.glmnet selects
                             log_transform = TRUE,
                             verbose = FALSE) {

  X_gmm <- sim_data$X_gmm; X_cox <- sim_data$X_cox
  time <- sim_data$time; status <- sim_data$status
  Z_true <- sim_data$Z_true; K_true <- sim_data$params$K
  n <- nrow(X_gmm)

  result <- list(K_true = K_true, K_selected = NA_integer_,
                 K_correct = NA_integer_, gamma_selected = NA_real_,
                 ARI = NA_real_, cindex_full = NA_real_,
                 cindex_cv = NA_real_, cindex_cox1 = NA_real_,
                 beta_rmse = NA_real_, n_events = sum(status),
                 event_rate_actual = mean(status),
                 converged = NA_integer_, error_msg = NA_character_)

  ## Preprocess (log-transform + impute + standardize)
  prep <- preprocess_train(X_gmm, X_cox, log_transform = log_transform)

  ## Baseline: single Cox
  result$cindex_cox1 <- tryCatch({
    sc <- make_x_scaler(prep$X_cox)
    Xs <- apply_x_scaler(prep$X_cox, sc)
    gfit <- glmnet(Xs, Surv(time, status), family = "cox",
                   alpha = alpha, standardize = FALSE)
    lam <- tryCatch(cv.glmnet(Xs, Surv(time, status), family = "cox",
                              alpha = alpha, standardize = FALSE,
                              nfolds = 3)$lambda.1se,
                    error = function(e) 0.05)
    eta <- as.numeric(Xs %*% coef(gfit, s = lam))
    c_index(time, status, eta)
  }, error = function(e) NA_real_)

  ## CV with 1-SE rule (built into v3)
  cv_res <- tryCatch(
    cv_select_K_gamma_v3(
      X_gmm = prep$X_gmm, X_cox = prep$X_cox,
      time = time, status = status,
      K_grid = K_grid, gamma_grid = gamma_grid,
      nfolds = nfolds, alpha = alpha,
      max_iter = max_iter, n_starts = n_starts,
      verbose = verbose,
      use_1se_rule = TRUE
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
      X_gmm = prep$X_gmm, X_cox = prep$X_cox,
      time = time, status = status, K = K_hat,
      alpha = alpha, max_iter = max_iter * 2,
      surv_weight = gamma_hat, verbose = FALSE,
      n_starts = n_starts_final
    ),
    error = function(e) e
  )
  if (inherits(fit, "error")) {
    result$error_msg <- fit$message; return(result)
  }
  result$converged <- 1L

  ## ARI
  if (K_hat == K_true && K_true > 1) {
    result$ARI <- adjustedRandIndex(Z_true, fit$clusterid)
  }

  ## C-index (training)
  eta_full <- tryCatch({
    tau <- fit$tau
    eta_mat <- sapply(1:K_hat, function(k)
      eta_from_scaled(prep$X_cox, fit$coxfit[[k]]$beta_s, fit$coxfit[[k]]$scaler))
    if (is.vector(eta_mat)) eta_mat <- matrix(eta_mat, ncol = K_hat)
    rowSums(tau * eta_mat)
  }, error = function(e) NULL)
  if (!is.null(eta_full))
    result$cindex_full <- c_index(time, status, eta_full)

  ## Beta recovery
  if (K_hat == K_true && K_true > 1) {
    true_betas <- sim_data$params$beta_list
    est_betas <- lapply(1:K_hat, function(k) fit$beta[, k])
    p_common <- min(length(true_betas[[1]]), nrow(fit$beta))
    cost <- matrix(0, K_true, K_hat)
    for (i in 1:K_true) for (j in 1:K_hat)
      cost[i,j] <- mean((true_betas[[i]][1:p_common] - est_betas[[j]][1:p_common])^2)
    used <- logical(K_hat); rmse_vals <- numeric(K_true)
    for (i in 1:K_true) {
      avail <- which(!used)
      bj <- avail[which.min(cost[i, avail])]
      rmse_vals[i] <- sqrt(cost[i, bj]); used[bj] <- TRUE
    }
    result$beta_rmse <- mean(rmse_vals)
  }

  result
}


###############################################################
## Simulation grid and execution
###############################################################

RUN_DIAGNOSTIC <- TRUE

sim_grid <- expand.grid(
  K_true = c(1, 2, 3, 4), separation = c("low", "med", "high"),
  event_rate = c(0.44, 0.25), n_subjects = 117,
  feat_corr = 0.2, stringsAsFactors = FALSE
)
sim_grid$scenario_id <- seq_len(nrow(sim_grid))

if (RUN_DIAGNOSTIC) {
  SIM_K_GRID <- 1:3;  SIM_GAMMA_GRID <- c(0, 0.5, 1.0, 2.0)
  SIM_NFOLDS <- 3;    SIM_MAX_ITER <- 30
  SIM_N_STARTS <- 3;  SIM_N_STARTS_FINAL <- 5
  NSIM <- 15L
} else {
  SIM_K_GRID <- 1:3;  SIM_GAMMA_GRID <- c(0, 0.25, 0.5, 1.0, 2.0)
  SIM_NFOLDS <- 5;    SIM_MAX_ITER <- 50
  SIM_N_STARTS <- 5;  SIM_N_STARTS_FINAL <- 10
  NSIM <- 100L
}

SIM_ALPHA <- 0.5
N_CORES <- min(50L, max(1L, detectCores(logical = FALSE) - 1L))

run_scenario_v3 <- function(sc_row, nsim = NSIM) {
  K_true <- sc_row$K_true; sep <- sc_row$separation
  er <- sc_row$event_rate; nn <- sc_row$n_subjects
  rho <- sc_row$feat_corr; sc_id <- sc_row$scenario_id

  cat(sprintf("[Scenario %d] K=%d sep=%s er=%.2f\n", sc_id, K_true, sep, er))

  results <- lapply(seq_len(nsim), function(rep_i) {
    sim_data <- tryCatch(
      simulate_gemcox_v3(
        n = nn, K = K_true, p_gmm = 10, p_cox = 10,
        p_gmm_signal = 5, p_beta_signal = 4,
        separation = sep, feat_corr = rho, event_rate = er,
        t_max = 253, lognormal = TRUE,
        seed = sc_id * 10000 + rep_i
      ),
      error = function(e) NULL
    )
    if (is.null(sim_data)) return(NULL)

    res <- tryCatch(
      run_replicate_v3(sim_data,
                       K_grid = SIM_K_GRID, gamma_grid = SIM_GAMMA_GRID,
                       nfolds = SIM_NFOLDS, max_iter = SIM_MAX_ITER,
                       n_starts = SIM_N_STARTS,
                       n_starts_final = SIM_N_STARTS_FINAL,
                       alpha = SIM_ALPHA, log_transform = TRUE,
                       verbose = FALSE),
      error = function(e) list(error_msg = e$message)
    )
    res <- lapply(res, function(v) if (is.null(v) || length(v) == 0) NA else v[1])
    as.data.frame(c(
      list(scenario_id = sc_id, rep = rep_i, K_true = K_true,
           separation = sep, event_rate_target = er, n = nn, feat_corr = rho),
      res), stringsAsFactors = FALSE)
  })
  bind_rows(Filter(Negate(is.null), results))
}

run_parallel_v3 <- function(n_jobs, FUN, nc) {
  if (nc <= 1L) return(lapply(seq_len(n_jobs), FUN))
  if (.Platform$OS.type == "windows") {
    cl <- makeCluster(nc, type = "PSOCK")
    on.exit(stopCluster(cl)); clusterExport(cl, ls(.GlobalEnv), .GlobalEnv)
    clusterEvalQ(cl, { library(survival); library(glmnet); library(MASS)
                       library(mclust); library(dplyr); library(parallel) })
    parLapply(cl, seq_len(n_jobs), FUN)
  } else {
    mclapply(seq_len(n_jobs), FUN, mc.cores = nc)
  }
}

message("========== STARTING v3 SIMULATION ==========")
all_res <- run_parallel_v3(nrow(sim_grid),
  function(i) tryCatch(run_scenario_v3(sim_grid[i,]),
    error = function(e) data.frame(scenario_id = sim_grid$scenario_id[i],
      error_msg = e$message, stringsAsFactors = FALSE)), N_CORES)

results_df <- bind_rows(Filter(is.data.frame, all_res))
saveRDS(results_df, "GeMCox_simulation_v3_raw.rds")

## Summarize
em <- results_df$error_msg
ok <- is.na(em) | (nzchar(trimws(em)) == FALSE) | (trimws(em) == "NA")
df_ok <- results_df[ok, ]

summary_df <- df_ok %>%
  group_by(K_true, separation, event_rate_target) %>%
  summarise(
    n_reps = n(),
    K_sel_accuracy = mean(K_correct, na.rm = TRUE),
    K_sel_mean = mean(K_selected, na.rm = TRUE),
    ARI_mean = mean(ARI, na.rm = TRUE),
    cindex_full_mean = mean(cindex_full, na.rm = TRUE),
    cindex_cv_mean = mean(cindex_cv, na.rm = TRUE),
    cindex_cox1_mean = mean(cindex_cox1, na.rm = TRUE),
    cindex_gain_mean = mean(cindex_full - cindex_cox1, na.rm = TRUE),
    gamma_mean = mean(gamma_selected, na.rm = TRUE),
    converge_rate = mean(converged, na.rm = TRUE),
    .groups = "drop"
  ) %>% mutate(across(where(is.numeric), ~round(.x, 4)))

write.csv(summary_df, "GeMCox_simulation_v3_summary.csv", row.names = FALSE)

cat("\n========== v3 RESULTS ==========\n")
print(as.data.frame(summary_df))
message("========== v3 SIMULATION COMPLETE ==========")
