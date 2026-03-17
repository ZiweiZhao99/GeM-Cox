###############################################################
## GeM-Cox Simulation v5
##
## Main changes vs v4:
##   1) n sweep: 117, 250, 500
##   2) K_true sweep: 1, 2, 3
##   3) DGP uses X_gmm != X_cox
##   4) Cluster-specific Cox effects use sign-flip/disjoint patterns,
##      so a single Cox model is a less competitive baseline
##   5) CV uses cluster-validity screening
##   6) Save both unconditional ARI and conditional ARI
###############################################################

required_pkgs <- c("survival", "glmnet", "mclust", "dplyr",
                   "tidyr", "tibble", "MASS", "parallel")
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

source("GeMCox_improved_v4.R")
source("GeMCox_improved_v4_revisions_patch.R")

set.seed(2025)

###############################################################
## 0. User settings
###############################################################
N_GRID       <- c(117L, 250L, 500L)
K_TRUE_GRID  <- 1:3
SEP_GRID     <- c("low", "med", "high")
ER_GRID      <- c(0.44, 0.25)
PATTERN_GRID <- c("signflip", "disjoint")

K_GRID       <- 1:3
GAMMA_GRID   <- c(0, 0.25, 0.5, 1.0, 2.0)

NSIM            <- 20L
CV_FOLDS        <- 3L
CV_MAX_ITER     <- 35L
CV_N_STARTS     <- 6L
FINAL_STARTS    <- 12L
ALPHA           <- 0.5
SELECTION_POLICY <- "screened_1se"
INIT_METHOD      <- "kmeans"
N_CORES <- min(24L, max(1L, parallel::detectCores(logical = FALSE) - 1L))

###############################################################
## 1. Data-generating mechanism
###############################################################
build_beta_patterns_v5 <- function(K, p, sep_scale = 1, pattern = c("signflip", "disjoint")) {
  pattern <- match.arg(pattern)
  B <- matrix(0, nrow = p, ncol = K)

  if (K == 1) {
    B[1:4, 1] <- sep_scale * c(0.8, -0.6, 0.5, -0.3)
    return(B)
  }

  if (pattern == "signflip") {
    if (K >= 2) {
      B[1:4, 1] <- sep_scale * c(1.2, -0.9, 0.7, 0.0)
      B[1:4, 2] <- sep_scale * c(-1.2, 0.9, 0.0, 0.7)
    }
    if (K >= 3) {
      B[1:4, 3] <- sep_scale * c(0.0, 0.0, -1.0, 1.0)
      B[5:6, 3] <- sep_scale * c(0.8, -0.8)
    }
  } else {
    if (K >= 1) B[c(1, 2, 7), 1] <- sep_scale * c(1.1, -0.8, 0.6)
    if (K >= 2) B[c(3, 4, 8), 2] <- sep_scale * c(-1.0, 0.8, 0.6)
    if (K >= 3) B[c(5, 6, 9), 3] <- sep_scale * c(1.0, -0.9, 0.5)
  }

  B
}

simulate_v5 <- function(n = 117,
                        K = 2,
                        q_gmm = 6,
                        p_cox = 10,
                        separation = c("low", "med", "high"),
                        pattern = c("signflip", "disjoint"),
                        rho_gmm = 0.25,
                        rho_cox = 0.20,
                        event_rate = 0.44,
                        seed = 1) {
  separation <- match.arg(separation)
  pattern <- match.arg(pattern)
  set.seed(seed)

  sep_mu <- switch(separation, low = 0.5, med = 1.0, high = 1.6)
  sep_beta <- switch(separation, low = 0.7, med = 1.0, high = 1.4)

  pos <- seq(-(K - 1) / 2, (K - 1) / 2, length.out = K)
  pi_c <- rep(1 / K, K)

  ## Clustering features
  Sigma_gmm <- outer(seq_len(q_gmm), seq_len(q_gmm), function(i, j) rho_gmm^abs(i - j))
  mu_list <- lapply(seq_len(K), function(k) {
    mu <- rep(0, q_gmm)
    mu[1:min(4, q_gmm)] <- pos[k] * sep_mu * c(1.0, 0.8, 0.6, 0.5)[1:min(4, q_gmm)]
    mu
  })

  Z <- sample(seq_len(K), size = n, replace = TRUE, prob = pi_c)
  X_gmm <- matrix(0, nrow = n, ncol = q_gmm)
  for (k in seq_len(K)) {
    idx <- which(Z == k)
    if (length(idx) == 0) next
    X_gmm[idx, ] <- MASS::mvrnorm(length(idx), mu = mu_list[[k]], Sigma = Sigma_gmm)
  }
  colnames(X_gmm) <- paste0("g", seq_len(q_gmm))

  ## Cox features: first few are noisy copies of GMM signal, rest are extra covariates
  Sigma_cox_extra <- outer(seq_len(p_cox - min(4, q_gmm)), seq_len(p_cox - min(4, q_gmm)), function(i, j) rho_cox^abs(i - j))
  X_cox <- matrix(0, nrow = n, ncol = p_cox)
  n_link <- min(4, q_gmm, p_cox)
  X_cox[, 1:n_link] <- X_gmm[, 1:n_link, drop = FALSE] + matrix(rnorm(n * n_link, sd = 0.35), n, n_link)
  if (p_cox > n_link) {
    X_cox[, (n_link + 1):p_cox] <- MASS::mvrnorm(n, mu = rep(0, p_cox - n_link), Sigma = Sigma_cox_extra)
  }
  colnames(X_cox) <- paste0("x", seq_len(p_cox))

  ## Make the latent class identifiable from X_gmm, but the Cox relation class-specific
  B <- build_beta_patterns_v5(K = K, p = p_cox, sep_scale = sep_beta, pattern = pattern)

  T_true <- numeric(n)
  base_scale <- 180
  for (k in seq_len(K)) {
    idx <- which(Z == k)
    if (length(idx) == 0) next
    eta <- as.numeric(X_cox[idx, , drop = FALSE] %*% B[, k])
    eta <- pmin(pmax(eta, -4.5), 4.5)
    T_true[idx] <- base_scale * (-log(runif(length(idx))))^(1 / 1.4) * exp(-eta / 1.4)
  }

  ## Calibrate censoring to target event rate
  lo <- 1e-6
  hi <- 10
  mid <- NA_real_
  for (it in 1:60) {
    mid <- (lo + hi) / 2
    er <- mean(T_true <= pmin(rexp(n, rate = mid), 253))
    if (abs(er - event_rate) < 0.005) break
    if (er < event_rate) hi <- mid else lo <- mid
  }

  C <- pmin(rexp(n, rate = mid), 253)
  time <- pmax(pmin(T_true, C) + runif(n, 0, 0.01), 0.1)
  status <- as.integer(T_true <= C)

  list(
    X_gmm = X_gmm,
    X_cox = X_cox,
    time = time,
    status = status,
    Z_true = Z,
    params = list(
      n = n,
      K = K,
      separation = separation,
      pattern = pattern,
      event_rate_target = event_rate,
      event_rate_actual = mean(status),
      beta_mat = B
    )
  )
}

###############################################################
## 2. Single replicate
###############################################################
fit_single_cox_v5 <- function(X, time, status, alpha = 0.5) {
  sc <- make_x_scaler(X)
  Xs <- apply_x_scaler(X, sc)
  lam <- tryCatch(
    cv.glmnet(x = Xs, y = survival::Surv(time, status), family = "cox",
              alpha = alpha, standardize = FALSE, nfolds = 3)$lambda.1se,
    error = function(e) 0.05
  )
  gf <- glmnet(x = Xs, y = survival::Surv(time, status), family = "cox",
               alpha = alpha, lambda = lam, standardize = FALSE)
  eta <- as.numeric(Xs %*% coef(gf, s = lam))
  list(cindex = c_index(time, status, eta), lambda = lam)
}

run_rep_v5 <- function(sd,
                       K_grid = 1:3,
                       gamma_grid = c(0, 0.25, 0.5, 1.0, 2.0),
                       nfolds = 3,
                       max_iter = 35,
                       n_starts = 6,
                       n_starts_final = 12,
                       alpha = 0.5,
                       selection_policy = "screened_1se",
                       init_method = "kmeans") {
  out <- list(
    n = sd$params$n,
    K_true = sd$params$K,
    separation = sd$params$separation,
    pattern = sd$params$pattern,
    event_rate_target = sd$params$event_rate_target,
    event_rate_actual = mean(sd$status),
    n_events = sum(sd$status),
    K_selected = NA_integer_,
    K_correct = NA_integer_,
    K_over = NA_integer_,
    K_under = NA_integer_,
    gamma_selected = NA_real_,
    selection_method = NA_character_,
    ARI_all = NA_real_,
    ARI_conditional = NA_real_,
    ARI_conditional_defined = NA_integer_,
    cindex_full = NA_real_,
    cindex_cv = NA_real_,
    cindex_cox1 = NA_real_,
    cindex_gain = NA_real_,
    beta_rmse_conditional = NA_real_,
    min_cluster_prop = NA_real_,
    min_eff_events = NA_real_,
    beta_sep = NA_real_,
    converged = NA_integer_,
    error_msg = NA_character_
  )

  out$cindex_cox1 <- tryCatch(fit_single_cox_v5(sd$X_cox, sd$time, sd$status, alpha = alpha)$cindex,
                              error = function(e) NA_real_)

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
      gmm_diag = TRUE,
      normalize_gmm_by_dim = FALSE,
      init_method = init_method,
      log_transform = FALSE,
      verbose = FALSE,
      use_1se_rule = TRUE,
      baseline_mode = "shared",
      selection_policy = selection_policy,
      min_cluster_prop = 0.08,
      min_eff_events_per_cluster = ifelse(sd$params$K <= 2, 8, 6),
      min_beta_sep = 0.03,
      require_nontrivial_K = FALSE
    ),
    error = function(e) e
  )
  if (inherits(cv, "error")) {
    out$error_msg <- cv$message
    return(out)
  }
  if (is.null(cv$best) || nrow(cv$best) == 0) {
    out$error_msg <- "no model selected"
    return(out)
  }

  K_hat <- as.integer(cv$best$K[1])
  gamma_hat <- as.numeric(cv$best$gamma[1])
  out$K_selected <- K_hat
  out$K_correct <- as.integer(K_hat == sd$params$K)
  out$K_over <- as.integer(K_hat > sd$params$K)
  out$K_under <- as.integer(K_hat < sd$params$K)
  out$gamma_selected <- gamma_hat
  out$cindex_cv <- as.numeric(cv$best$mean_score[1])
  out$selection_method <- cv$best_method

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
      gmm_diag = TRUE,
      normalize_gmm_by_dim = FALSE,
      init_method = init_method,
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

  scr <- fit_screen_v4(fit, sd$status,
                       min_cluster_prop = 0,
                       min_eff_events_per_cluster = 0,
                       min_beta_sep = 0,
                       require_nontrivial_K = FALSE)
  out$min_cluster_prop <- scr$min_cluster_prop
  out$min_eff_events <- scr$min_eff_events
  out$beta_sep <- scr$beta_sep

  if (sd$params$K > 1) {
    out$ARI_all <- adjustedRandIndex(sd$Z_true, fit$clusterid)
  }
  out$ARI_conditional_defined <- as.integer(K_hat == sd$params$K && sd$params$K > 1)
  if (isTRUE(out$ARI_conditional_defined == 1L)) {
    out$ARI_conditional <- adjustedRandIndex(sd$Z_true, fit$clusterid)
  }

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

  if (K_hat == sd$params$K && sd$params$K > 1) {
    tb <- sd$params$beta_mat
    eb <- as.matrix(fit$beta)
    pcomp <- min(nrow(tb), nrow(eb))
    cost <- matrix(0, nrow = ncol(tb), ncol = ncol(eb))
    for (i in seq_len(ncol(tb))) {
      for (j in seq_len(ncol(eb))) {
        cost[i, j] <- mean((tb[1:pcomp, i] - eb[1:pcomp, j])^2)
      }
    }
    used <- logical(ncol(eb))
    rmse_vec <- numeric(ncol(tb))
    for (i in seq_len(ncol(tb))) {
      avail <- which(!used)
      bj <- avail[which.min(cost[i, avail])]
      rmse_vec[i] <- sqrt(cost[i, bj])
      used[bj] <- TRUE
    }
    out$beta_rmse_conditional <- mean(rmse_vec)
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
  pattern = PATTERN_GRID,
  stringsAsFactors = FALSE
) |>
  dplyr::mutate(sid = dplyr::row_number()) |>
  dplyr::arrange(n, K_true, separation, event_rate, pattern)

run_scenario_v5 <- function(sc) {
  cat(sprintf("[sid=%03d] n=%d K=%d sep=%s er=%.2f patt=%s\n",
              sc$sid, sc$n, sc$K_true, sc$separation, sc$event_rate, sc$pattern))
  out_list <- lapply(seq_len(NSIM), function(ri) {
    dat <- tryCatch(
      simulate_v5(
        n = sc$n,
        K = sc$K_true,
        separation = sc$separation,
        pattern = sc$pattern,
        event_rate = sc$event_rate,
        seed = sc$sid * 10000 + ri
      ),
      error = function(e) NULL
    )
    if (is.null(dat)) return(NULL)

    res <- tryCatch(
      run_rep_v5(
        sd = dat,
        K_grid = K_GRID,
        gamma_grid = GAMMA_GRID,
        nfolds = CV_FOLDS,
        max_iter = CV_MAX_ITER,
        n_starts = CV_N_STARTS,
        n_starts_final = FINAL_STARTS,
        alpha = ALPHA,
        selection_policy = SELECTION_POLICY,
        init_method = INIT_METHOD
      ),
      error = function(e) list(error_msg = e$message)
    )

    res <- lapply(res, function(v) if (is.null(v) || length(v) == 0) NA else v[1])
    as.data.frame(c(list(sid = sc$sid, rep = ri), res), stringsAsFactors = FALSE)
  })
  dplyr::bind_rows(Filter(Negate(is.null), out_list))
}

###############################################################
## 4. Run simulation
###############################################################
message("========== GeM-Cox v5 simulation start ==========")
cat(sprintf("Scenarios: %d\n", nrow(sim_grid)))
cat(sprintf("Replicates/scenario: %d\n", NSIM))
cat(sprintf("Total replicate jobs: %d\n", nrow(sim_grid) * NSIM))
cat(sprintf("K grid: %s\n", paste(K_GRID, collapse = ",")))
cat(sprintf("Gamma grid: %s\n", paste(GAMMA_GRID, collapse = ",")))
cat(sprintf("Selection policy: %s\n", SELECTION_POLICY))
cat(sprintf("Init method: %s\n", INIT_METHOD))
cat(sprintf("Cores: %d\n", N_CORES))

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
    NULL
  })
  all_res <- parLapply(cl, seq_len(nrow(sim_grid)), function(i) {
    tryCatch(run_scenario_v5(sim_grid[i, ]),
             error = function(e) data.frame(sid = sim_grid$sid[i], error_msg = e$message,
                                            stringsAsFactors = FALSE))
  })
} else if (N_CORES > 1) {
  all_res <- mclapply(seq_len(nrow(sim_grid)), function(i) {
    tryCatch(run_scenario_v5(sim_grid[i, ]),
             error = function(e) data.frame(sid = sim_grid$sid[i], error_msg = e$message,
                                            stringsAsFactors = FALSE))
  }, mc.cores = N_CORES)
} else {
  all_res <- lapply(seq_len(nrow(sim_grid)), function(i) run_scenario_v5(sim_grid[i, ]))
}

elapsed <- (proc.time() - t0)[3]
results_df <- dplyr::bind_rows(Filter(is.data.frame, all_res))
saveRDS(results_df, "GeMCox_simulation_v5_raw.rds")

###############################################################
## 5. Summaries
###############################################################
em <- results_df$error_msg
ok <- is.na(em) | !nzchar(trimws(em)) | trimws(em) == "NA"
df_ok <- results_df[ok, , drop = FALSE]

summary_df <- df_ok |>
  dplyr::group_by(n, K_true, separation, event_rate_target, pattern) |>
  dplyr::summarise(
    n_reps                    = dplyr::n(),
    n_events_mean             = round(mean(n_events, na.rm = TRUE), 1),
    event_rate_mean           = round(mean(event_rate_actual, na.rm = TRUE), 3),
    K_sel_accuracy            = round(mean(K_correct, na.rm = TRUE), 3),
    K_sel_mean                = round(mean(K_selected, na.rm = TRUE), 2),
    K_over_rate               = round(mean(K_over, na.rm = TRUE), 3),
    K_under_rate              = round(mean(K_under, na.rm = TRUE), 3),
    ARI_all_mean              = round(mean(ARI_all, na.rm = TRUE), 3),
    ARI_all_sd                = round(sd(ARI_all, na.rm = TRUE), 3),
    ARI_conditional_n         = sum(ARI_conditional_defined %in% 1, na.rm = TRUE),
    ARI_conditional_mean      = round(mean(ARI_conditional, na.rm = TRUE), 3),
    ARI_conditional_sd        = round(sd(ARI_conditional, na.rm = TRUE), 3),
    beta_rmse_conditional_mean= round(mean(beta_rmse_conditional, na.rm = TRUE), 3),
    cindex_full_mean          = round(mean(cindex_full, na.rm = TRUE), 3),
    cindex_cv_mean            = round(mean(cindex_cv, na.rm = TRUE), 3),
    cindex_cox1_mean          = round(mean(cindex_cox1, na.rm = TRUE), 3),
    cindex_gain_mean          = round(mean(cindex_gain, na.rm = TRUE), 3),
    gamma_mean                = round(mean(gamma_selected, na.rm = TRUE), 2),
    min_cluster_prop_mean     = round(mean(min_cluster_prop, na.rm = TRUE), 3),
    min_eff_events_mean       = round(mean(min_eff_events, na.rm = TRUE), 2),
    beta_sep_mean             = round(mean(beta_sep, na.rm = TRUE), 3),
    converge_rate             = round(mean(converged, na.rm = TRUE), 3),
    .groups = "drop"
  ) |>
  dplyr::arrange(n, K_true, separation, event_rate_target, pattern)

summary_nk_df <- df_ok |>
  dplyr::group_by(n, K_true, pattern) |>
  dplyr::summarise(
    n_reps                    = dplyr::n(),
    K_sel_accuracy            = round(mean(K_correct, na.rm = TRUE), 3),
    K_sel_mean                = round(mean(K_selected, na.rm = TRUE), 2),
    ARI_all_mean              = round(mean(ARI_all, na.rm = TRUE), 3),
    ARI_conditional_n         = sum(ARI_conditional_defined %in% 1, na.rm = TRUE),
    ARI_conditional_mean      = round(mean(ARI_conditional, na.rm = TRUE), 3),
    cindex_full_mean          = round(mean(cindex_full, na.rm = TRUE), 3),
    cindex_cox1_mean          = round(mean(cindex_cox1, na.rm = TRUE), 3),
    cindex_gain_mean          = round(mean(cindex_gain, na.rm = TRUE), 3),
    min_cluster_prop_mean     = round(mean(min_cluster_prop, na.rm = TRUE), 3),
    min_eff_events_mean       = round(mean(min_eff_events, na.rm = TRUE), 2),
    beta_sep_mean             = round(mean(beta_sep, na.rm = TRUE), 3),
    converge_rate             = round(mean(converged, na.rm = TRUE), 3),
    .groups = "drop"
  ) |>
  dplyr::arrange(n, K_true, pattern)

write.csv(summary_df, "GeMCox_simulation_v5_summary.csv", row.names = FALSE)
write.csv(summary_nk_df, "GeMCox_simulation_v5_summary_by_nK.csv", row.names = FALSE)

cat(sprintf("\nSuccess reps: %d / %d (%.1f%%)\n",
            nrow(df_ok), nrow(results_df), 100 * nrow(df_ok) / max(1, nrow(results_df))))
cat(sprintf("Wall time: %.1f sec (%.1f min)\n", elapsed, elapsed / 60))
cat("Saved:\n")
cat("  - GeMCox_simulation_v5_raw.rds\n")
cat("  - GeMCox_simulation_v5_summary.csv\n")
cat("  - GeMCox_simulation_v5_summary_by_nK.csv\n")
