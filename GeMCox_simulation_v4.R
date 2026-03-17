###############################################################
## GeM-Cox Simulation v4 — Full Sweep over n and K
##
## Main goals:
##   1) Sweep sample size: n = 117, 250, 500
##   2) Sweep true K: K_true = 1, 2, 3
##   3) Keep model selection grid K_GRID = 1:3
##   4) Retain v4 fix: normalize_gmm_by_dim = FALSE
##   5) Summarize performance by n, K_true, separation, event_rate
##
## Notes:
##   - This is much larger than the previous fast run.
##   - Recommended:
##       * NSIM = 20 for a quick check
##       * NSIM = 50 for a good draft result
##       * NSIM = 100 for final manuscript-level run
###############################################################

required_pkgs <- c(
  "survival", "glmnet", "mclust", "dplyr",
  "tidyr", "tibble", "MASS", "parallel"
)
for (pkg in required_pkgs) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg, repos = "https://cloud.r-project.org")
  }
}

library(survival)
library(glmnet)
library(mclust)
library(dplyr)
library(tidyr)
library(tibble)
library(MASS)
library(parallel)

select  <- dplyr::select
filter  <- dplyr::filter
mutate  <- dplyr::mutate
arrange <- dplyr::arrange

source("GeMCox_improved_v4.R")
set.seed(2025)

###############################################################
## 0. User settings
###############################################################

## Full scenario sweep
N_GRID       <- c(117L, 250L, 500L)
K_TRUE_GRID  <- 1:3
SEP_GRID     <- c("low", "med", "high")
ER_GRID      <- c(0.44, 0.25)

## Model selection grid
K_GRID       <- 1:3
GAMMA_GRID   <- c(0, 0.25, 0.5, 1.0, 2.0)

## Simulation replicates
NSIM         <- 50L     # set to 20 for quick check; 100 for final run

## Fitting controls
CV_FOLDS      <- 3L
CV_MAX_ITER   <- 30L
CV_N_STARTS   <- 3L
FINAL_STARTS  <- 5L
ALPHA         <- 0.5

## Parallel
N_CORES <- min(24L, max(1L, parallel::detectCores(logical = FALSE) - 1L))

###############################################################
## 1. Data-generating mechanism
###############################################################
simulate_v4 <- function(
    n = 117,
    K = 2,
    p = 10,
    p_sig = 5,
    p_beta = 4,
    separation = c("med", "low", "high"),
    rho = 0.2,
    event_rate = 0.44,
    seed = 1
) {
  separation <- match.arg(separation)
  set.seed(seed)
  
  dmu <- switch(separation,
                low  = 0.6,
                med  = 1.2,
                high = 2.0)
  
  dbeta <- switch(separation,
                  low  = 0.3,
                  med  = 0.6,
                  high = 1.0)
  
  pi_c <- rep(1 / K, K)
  pos  <- seq(-(K - 1) / 2, (K - 1) / 2, length.out = K)
  
  ## Cluster means
  mu_list <- lapply(seq_len(K), function(k) {
    m <- rep(0, p)
    for (j in seq_len(min(p_sig, p))) {
      m[j] <- pos[k] * dmu / sqrt(j)
    }
    m
  })
  
  ## Correlated Gaussian covariance
  Sigma <- outer(seq_len(p), seq_len(p), function(i, j) rho^abs(i - j))
  
  ## Cluster-specific Cox coefficients
  beta_base <- rep(0, p)
  beta_base[seq_len(min(p_beta, p))] <- c(0.8, -0.5, 0.4, -0.3)[seq_len(min(p_beta, p))]
  
  beta_list <- lapply(seq_len(K), function(k) {
    b <- beta_base
    for (j in seq_len(min(3, p_beta))) {
      b[j] <- b[j] + dbeta * pos[k] / sqrt(j)
    }
    b
  })
  
  ## Optional mild cluster difference in baseline timescale
  wscales <- 200 * exp(0.2 * pos)
  
  ## Draw latent class + features
  Z <- sample(seq_len(K), n, replace = TRUE, prob = pi_c)
  X <- matrix(0, nrow = n, ncol = p)
  
  for (k in seq_len(K)) {
    idx <- which(Z == k)
    if (!length(idx)) next
    X[idx, ] <- MASS::mvrnorm(n = length(idx), mu = mu_list[[k]], Sigma = Sigma)
  }
  colnames(X) <- paste0("f", seq_len(p))
  
  ## Generate event times from observable X and class-specific betas
  T_true <- numeric(n)
  for (k in seq_len(K)) {
    idx <- which(Z == k)
    if (!length(idx)) next
    eta <- as.numeric(X[idx, , drop = FALSE] %*% beta_list[[k]])
    eta <- pmin(pmax(eta, -5), 5)
    T_true[idx] <- wscales[k] *
      (-log(stats::runif(length(idx))))^(1 / 1.5) *
      exp(-eta / 1.5)
  }
  
  ## Calibrate censoring to target event rate
  lo <- 1e-6
  hi <- 10
  mid <- NA_real_
  for (it in 1:50) {
    mid <- (lo + hi) / 2
    er <- mean(T_true <= pmin(stats::rexp(n, rate = mid), 253))
    if (abs(er - event_rate) < 0.005) break
    if (er < event_rate) {
      hi <- mid
    } else {
      lo <- mid
    }
  }
  
  C <- pmin(stats::rexp(n, rate = mid), 253)
  time <- pmax(pmin(T_true, C) + stats::runif(n, 0, 0.01), 0.1)
  status <- as.integer(T_true <= C)
  
  list(
    X_gmm = X,
    X_cox = X,
    time = time,
    status = status,
    Z_true = Z,
    params = list(
      n = n,
      K = K,
      separation = separation,
      dmu = dmu,
      dbeta = dbeta,
      event_rate_target = event_rate,
      event_rate_actual = mean(status),
      beta_list = beta_list
    )
  )
}

###############################################################
## 2. Single replicate
###############################################################
run_rep <- function(
    sd,
    K_grid = 1:3,
    gamma_grid = c(0, 0.25, 0.5, 1.0, 2.0),
    nfolds = 3,
    max_iter = 30,
    n_starts = 3,
    n_starts_final = 5,
    alpha = 0.5
) {
  out <- list(
    n = sd$params$n,
    K_true = sd$params$K,
    K_selected = NA_integer_,
    K_correct = NA_integer_,
    K_over = NA_integer_,
    K_under = NA_integer_,
    gamma_selected = NA_real_,
    ARI = NA_real_,
    ARI_defined = NA_integer_,
    cindex_full = NA_real_,
    cindex_cv = NA_real_,
    cindex_cox1 = NA_real_,
    cindex_gain = NA_real_,
    beta_rmse = NA_real_,
    event_rate_actual = mean(sd$status),
    n_events = sum(sd$status),
    converged = NA_integer_,
    error_msg = NA_character_
  )
  
  ## Reference single Cox model
  out$cindex_cox1 <- tryCatch({
    sc <- make_x_scaler(sd$X_cox)
    Xs <- apply_x_scaler(sd$X_cox, sc)
    
    lam <- tryCatch(
      cv.glmnet(
        x = Xs,
        y = survival::Surv(sd$time, sd$status),
        family = "cox",
        alpha = alpha,
        standardize = FALSE,
        nfolds = 3
      )$lambda.1se,
      error = function(e) 0.05
    )
    
    gf <- glmnet(
      x = Xs,
      y = survival::Surv(sd$time, sd$status),
      family = "cox",
      alpha = alpha,
      lambda = lam,
      standardize = FALSE
    )
    
    c_index(sd$time, sd$status, as.numeric(Xs %*% coef(gf, s = lam)))
  }, error = function(e) NA_real_)
  
  ## CV model selection
  cv <- tryCatch(
    cv_select_K_gamma_v4(
      X_gmm = sd$X_gmm,
      X_cox = sd$X_cox,
      time = sd$time,
      status = sd$status,
      K_grid = K_grid,
      gamma_grid = gamma_grid,
      nfolds = nfolds,
      alpha = alpha,
      max_iter = max_iter,
      n_starts = n_starts,
      log_transform = FALSE,
      normalize_gmm_by_dim = FALSE,   ## critical fix
      verbose = FALSE,
      use_1se_rule = TRUE,
      baseline_mode = "shared"
    ),
    error = function(e) e
  )
  
  if (inherits(cv, "error")) {
    out$error_msg <- cv$message
    return(out)
  }
  
  best <- cv$best
  if (is.null(best) || nrow(best) == 0) {
    out$error_msg <- "no model selected"
    return(out)
  }
  
  K_hat <- as.integer(best$K[1])
  gamma_hat <- as.numeric(best$gamma[1])
  
  out$K_selected <- K_hat
  out$K_correct  <- as.integer(K_hat == sd$params$K)
  out$K_over     <- as.integer(K_hat > sd$params$K)
  out$K_under    <- as.integer(K_hat < sd$params$K)
  out$gamma_selected <- gamma_hat
  out$cindex_cv  <- as.numeric(best$mean_score[1])
  
  ## Final fit
  prep <- preprocess_train_v4(sd$X_gmm, sd$X_cox, log_transform = FALSE)
  
  fit <- tryCatch(
    gemcox_full_multistart_shared_v4(
      X_gmm = prep$X_gmm,
      X_cox = prep$X_cox,
      time = sd$time,
      status = sd$status,
      K = K_hat,
      alpha = alpha,
      max_iter = max_iter * 2,
      surv_weight = gamma_hat,
      normalize_gmm_by_dim = FALSE,   ## critical fix
      verbose = FALSE,
      n_starts = n_starts_final
    ),
    error = function(e) e
  )
  
  if (inherits(fit, "error")) {
    out$error_msg <- fit$message
    return(out)
  }
  
  out$converged <- 1L
  
  ## ARI only defined when selected K matches true K and K>1
  out$ARI_defined <- as.integer(K_hat == sd$params$K && sd$params$K > 1)
  if (isTRUE(out$ARI_defined == 1L)) {
    out$ARI <- mclust::adjustedRandIndex(sd$Z_true, fit$clusterid)
  }
  
  ## Full-data discrimination from fitted mixture
  eta_mix <- tryCatch({
    tau <- fit$tau
    eta_mat <- sapply(seq_len(K_hat), function(k) {
      eta_from_scaled(prep$X_cox, fit$coxfit[[k]]$beta_s, fit$coxfit[[k]]$scaler)
    })
    if (is.vector(eta_mat)) eta_mat <- matrix(eta_mat, ncol = K_hat)
    rowSums(tau * eta_mat)
  }, error = function(e) NULL)
  
  if (!is.null(eta_mix)) {
    out$cindex_full <- c_index(sd$time, sd$status, eta_mix)
    out$cindex_gain <- out$cindex_full - out$cindex_cox1
  }
  
  ## Beta RMSE when K is correctly selected and K>1
  if (K_hat == sd$params$K && sd$params$K > 1) {
    tb <- sd$params$beta_list
    eb <- lapply(seq_len(K_hat), function(k) fit$beta[, k])
    
    pcomp <- min(length(tb[[1]]), nrow(fit$beta))
    cost <- matrix(0, nrow = sd$params$K, ncol = K_hat)
    
    for (i in seq_len(sd$params$K)) {
      for (j in seq_len(K_hat)) {
        cost[i, j] <- mean((tb[[i]][1:pcomp] - eb[[j]][1:pcomp])^2)
      }
    }
    
    used <- logical(K_hat)
    rmse_vec <- numeric(sd$params$K)
    for (i in seq_len(sd$params$K)) {
      avail <- which(!used)
      best_j <- avail[which.min(cost[i, avail])]
      rmse_vec[i] <- sqrt(cost[i, best_j])
      used[best_j] <- TRUE
    }
    out$beta_rmse <- mean(rmse_vec)
  }
  
  out
}

###############################################################
## 3. Scenario grid
###############################################################
sim_grid <- expand.grid(
  n = N_GRID,
  K_true = K_TRUE_GRID,
  separation = SEP_GRID,
  event_rate = ER_GRID,
  stringsAsFactors = FALSE
)
sim_grid <- sim_grid %>%
  mutate(sid = dplyr::row_number()) %>%
  arrange(n, K_true, separation, event_rate)

###############################################################
## 4. Run one scenario
###############################################################
run_scenario <- function(sc) {
  cat(sprintf(
    "[sid=%02d] n=%d  K=%d  sep=%s  er=%.2f\n",
    sc$sid, sc$n, sc$K_true, sc$separation, sc$event_rate
  ))
  
  out_list <- lapply(seq_len(NSIM), function(ri) {
    dat <- tryCatch(
      simulate_v4(
        n = sc$n,
        K = sc$K_true,
        separation = sc$separation,
        event_rate = sc$event_rate,
        seed = sc$sid * 10000 + ri
      ),
      error = function(e) NULL
    )
    
    if (is.null(dat)) return(NULL)
    
    res <- tryCatch(
      run_rep(
        sd = dat,
        K_grid = K_GRID,
        gamma_grid = GAMMA_GRID,
        nfolds = CV_FOLDS,
        max_iter = CV_MAX_ITER,
        n_starts = CV_N_STARTS,
        n_starts_final = FINAL_STARTS,
        alpha = ALPHA
      ),
      error = function(e) list(error_msg = e$message)
    )
    
    res <- lapply(res, function(v) {
      if (is.null(v) || length(v) == 0) NA else v[1]
    })
    
    as.data.frame(
      c(
        list(
          sid = sc$sid,
          rep = ri,
          n = sc$n,
          K_true = sc$K_true,
          separation = sc$separation,
          event_rate = sc$event_rate
        ),
        res
      ),
      stringsAsFactors = FALSE
    )
  })
  
  dplyr::bind_rows(Filter(Negate(is.null), out_list))
}

###############################################################
## 5. Execute simulation
###############################################################
message("========== GeM-Cox v4 FULL SWEEP START ==========")

cat(sprintf("Scenario grid: %d rows\n", nrow(sim_grid)))
cat(sprintf("Replicates per scenario: %d\n", NSIM))
cat(sprintf("Total replicate jobs: %d\n", nrow(sim_grid) * NSIM))
cat(sprintf("K selection grid: %s\n", paste(K_GRID, collapse = ",")))
cat(sprintf("Gamma grid: %s\n", paste(GAMMA_GRID, collapse = ",")))
cat(sprintf("Using up to %d cores\n", N_CORES))

t0 <- proc.time()

if (N_CORES > 1 && .Platform$OS.type == "windows") {
  cl <- makeCluster(N_CORES, type = "PSOCK")
  on.exit(stopCluster(cl), add = TRUE)
  
  clusterExport(cl, ls(.GlobalEnv), .GlobalEnv)
  clusterEvalQ(cl, {
    library(survival)
    library(glmnet)
    library(MASS)
    library(mclust)
    library(dplyr)
    library(parallel)
    NULL
  })
  
  all_res <- parLapply(
    cl,
    X = seq_len(nrow(sim_grid)),
    fun = function(i) {
      tryCatch(
        run_scenario(sim_grid[i, ]),
        error = function(e) data.frame(
          sid = sim_grid$sid[i],
          n = sim_grid$n[i],
          K_true = sim_grid$K_true[i],
          separation = sim_grid$separation[i],
          event_rate = sim_grid$event_rate[i],
          error_msg = e$message,
          stringsAsFactors = FALSE
        )
      )
    }
  )
} else if (N_CORES > 1) {
  all_res <- mclapply(
    X = seq_len(nrow(sim_grid)),
    FUN = function(i) {
      tryCatch(
        run_scenario(sim_grid[i, ]),
        error = function(e) data.frame(
          sid = sim_grid$sid[i],
          n = sim_grid$n[i],
          K_true = sim_grid$K_true[i],
          separation = sim_grid$separation[i],
          event_rate = sim_grid$event_rate[i],
          error_msg = e$message,
          stringsAsFactors = FALSE
        )
      )
    },
    mc.cores = N_CORES
  )
} else {
  all_res <- lapply(seq_len(nrow(sim_grid)), function(i) run_scenario(sim_grid[i, ]))
}

elapsed <- (proc.time() - t0)[3]

results_df <- dplyr::bind_rows(Filter(is.data.frame, all_res))
saveRDS(results_df, "GeMCox_simulation_v4_raw.rds")

###############################################################
## 6. Summaries
###############################################################
err <- results_df$error_msg
ok <- is.na(err) | !nzchar(trimws(err)) | trimws(err) == "NA"
df_ok <- results_df[ok, , drop = FALSE]

cat(sprintf(
  "\nSuccessful reps: %d / %d (%.1f%%)\n",
  nrow(df_ok), nrow(results_df),
  100 * nrow(df_ok) / max(1, nrow(results_df))
))
cat(sprintf("Wall time: %.1f sec (%.1f min)\n", elapsed, elapsed / 60))

summary_df <- df_ok %>%
  group_by(n, K_true, separation, event_rate) %>%
  summarise(
    n_reps                 = n(),
    n_events_mean          = round(mean(n_events, na.rm = TRUE), 1),
    event_rate_mean        = round(mean(event_rate_actual, na.rm = TRUE), 3),
    
    K_sel_accuracy         = round(mean(K_correct, na.rm = TRUE), 3),
    K_sel_mean             = round(mean(K_selected, na.rm = TRUE), 2),
    K_over_rate            = round(mean(K_over, na.rm = TRUE), 3),
    K_under_rate           = round(mean(K_under, na.rm = TRUE), 3),
    
    ARI_n_defined          = sum(ARI_defined %in% 1, na.rm = TRUE),
    ARI_mean               = round(mean(ARI, na.rm = TRUE), 3),
    ARI_sd                 = round(sd(ARI, na.rm = TRUE), 3),
    
    beta_rmse_mean         = round(mean(beta_rmse, na.rm = TRUE), 3),
    
    cindex_full_mean       = round(mean(cindex_full, na.rm = TRUE), 3),
    cindex_cv_mean         = round(mean(cindex_cv, na.rm = TRUE), 3),
    cindex_cox1_mean       = round(mean(cindex_cox1, na.rm = TRUE), 3),
    cindex_gain_mean       = round(mean(cindex_gain, na.rm = TRUE), 3),
    
    gamma_mean             = round(mean(gamma_selected, na.rm = TRUE), 2),
    converge_rate          = round(mean(converged, na.rm = TRUE), 3),
    .groups = "drop"
  ) %>%
  arrange(n, K_true, separation, event_rate)

write.csv(summary_df, "GeMCox_simulation_v4_summary.csv", row.names = FALSE)

## Additional compact summary by n and K_true only
summary_nk_df <- df_ok %>%
  group_by(n, K_true) %>%
  summarise(
    n_reps                 = n(),
    K_sel_accuracy         = round(mean(K_correct, na.rm = TRUE), 3),
    K_sel_mean             = round(mean(K_selected, na.rm = TRUE), 2),
    K_over_rate            = round(mean(K_over, na.rm = TRUE), 3),
    K_under_rate           = round(mean(K_under, na.rm = TRUE), 3),
    ARI_n_defined          = sum(ARI_defined %in% 1, na.rm = TRUE),
    ARI_mean               = round(mean(ARI, na.rm = TRUE), 3),
    beta_rmse_mean         = round(mean(beta_rmse, na.rm = TRUE), 3),
    cindex_full_mean       = round(mean(cindex_full, na.rm = TRUE), 3),
    cindex_cox1_mean       = round(mean(cindex_cox1, na.rm = TRUE), 3),
    cindex_gain_mean       = round(mean(cindex_gain, na.rm = TRUE), 3),
    gamma_mean             = round(mean(gamma_selected, na.rm = TRUE), 2),
    converge_rate          = round(mean(converged, na.rm = TRUE), 3),
    .groups = "drop"
  ) %>%
  arrange(n, K_true)

write.csv(summary_nk_df, "GeMCox_simulation_v4_summary_by_nK.csv", row.names = FALSE)

###############################################################
## 7. Print results
###############################################################
cat("\n========== FULL SUMMARY ==========\n")
print(as.data.frame(summary_df), right = FALSE, row.names = FALSE)

cat("\n========== COMPACT SUMMARY BY n AND K_true ==========\n")
print(as.data.frame(summary_nk_df), right = FALSE, row.names = FALSE)

cat(sprintf("\nTotal wall time: %.1f min\n", elapsed / 60))
cat("Saved:\n")
cat("  - GeMCox_simulation_v4_raw.rds\n")
cat("  - GeMCox_simulation_v4_summary.csv\n")
cat("  - GeMCox_simulation_v4_summary_by_nK.csv\n")