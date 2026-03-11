###############################################################
## Quick diagnostic: Does v4 shared-baseline EM work on clean data?
## Run this BEFORE the full simulation to verify the core algorithm
###############################################################

library(MASS); library(survival); library(glmnet); library(mclust)
source("GeMCox_improved_v4.R")

set.seed(2025)

## ============================================================
## TEST 1: Can GMM alone recover clusters? (no survival)
## ============================================================
cat("=== TEST 1: Pure GMM cluster recovery (no survival) ===\n")

n <- 200; K <- 2
mu1 <- c(1.5, 1.0, rep(0, 8))
mu2 <- c(-1.5, -1.0, rep(0, 8))
Sigma <- diag(10)
Z <- rep(1:2, each = 100)
X <- rbind(mvrnorm(100, mu1, Sigma), mvrnorm(100, mu2, Sigma))

km <- kmeans(X, 2, nstart = 50)
cat(sprintf("  K-means ARI: %.4f\n", adjustedRandIndex(Z, km$cluster)))

## ============================================================
## TEST 2: v4 shared-baseline EM on cleanest possible data
## ============================================================
cat("\n=== TEST 2: v4 shared-baseline EM, K=2, Gaussian, n=300 ===\n")

simulate_clean <- function(n = 300, K = 2, delta_mu = 2.0, delta_beta = 0.8, seed = 1) {
  set.seed(seed)
  p <- 10; p_sig <- 5
  pi_c <- rep(1/K, K)
  positions <- seq(-(K-1)/2, (K-1)/2, length.out = K) * delta_mu

  mu_list <- lapply(1:K, function(k) {
    m <- rep(0, p); for (j in 1:p_sig) m[j] <- positions[k]/sqrt(j); m
  })
  Sigma <- outer(1:p, 1:p, function(i,j) 0.2^abs(i-j))

  ## Betas: same sign, different magnitude
  beta_base <- c(0.8, -0.5, 0.4, -0.3, rep(0, p-4))
  beta_list <- lapply(1:K, function(k) {
    pos <- (k - (K+1)/2)
    b <- beta_base
    for (j in 1:3) b[j] <- b[j] + delta_beta * pos / sqrt(j)
    b
  })

  Z <- sample(1:K, n, replace = TRUE, prob = pi_c)
  X <- matrix(0, n, p)
  for (k in 1:K) {
    idx <- which(Z == k)
    X[idx, ] <- mvrnorm(length(idx), mu_list[[k]], Sigma)
  }
  colnames(X) <- paste0("f", 1:p)

  ## Survival on OBSERVABLE features
  T_true <- numeric(n)
  for (k in 1:K) {
    idx <- which(Z == k)
    eta <- as.numeric(X[idx, ] %*% beta_list[[k]])
    eta <- pmin(pmax(eta, -5), 5)
    T_true[idx] <- 200 * (-log(runif(length(idx))))^(1/1.5) * exp(-eta/1.5)
  }
  crate <- 0.005  # gives ~40-50% events
  C <- pmin(rexp(n, crate), 253)
  time <- pmax(pmin(T_true, C) + runif(n, 0, 0.01), 0.1)
  status <- as.integer(T_true <= C)

  list(X_gmm = X, X_cox = X, time = time, status = status,
       Z_true = Z, beta_list = beta_list, mu_list = mu_list)
}

dat <- simulate_clean(n = 300, K = 2, delta_mu = 2.0, delta_beta = 0.8, seed = 42)
cat(sprintf("  Events: %d/%d = %.0f%%\n", sum(dat$status), length(dat$status),
            100*mean(dat$status)))
cat(sprintf("  True beta[1]: cluster1=%.3f cluster2=%.3f\n",
            dat$beta_list[[1]][1], dat$beta_list[[2]][1]))

## Fit v4 shared-baseline EM
fit <- tryCatch(
  gemcox_full_multistart_shared_v4(
    X_gmm = dat$X_gmm, X_cox = dat$X_cox,
    time = dat$time, status = dat$status,
    K = 2, alpha = 0.5, max_iter = 100,
    surv_weight = 1.0, verbose = FALSE, n_starts = 20
  ),
  error = function(e) { cat("  ERROR:", e$message, "\n"); NULL }
)

if (!is.null(fit)) {
  ari <- adjustedRandIndex(dat$Z_true, fit$clusterid)
  cat(sprintf("  v4 shared EM ARI: %.4f\n", ari))
  cat(sprintf("  Cluster sizes: %s (true: %s)\n",
              paste(table(fit$clusterid), collapse="/"),
              paste(table(dat$Z_true), collapse="/")))
  cat(sprintf("  Pi: [%.3f, %.3f]\n", fit$pi[1], fit$pi[2]))
  cat(sprintf("  Est beta[1]: c1=%.3f c2=%.3f\n", fit$beta[1,1], fit$beta[1,2]))
  cat(sprintf("  Est beta[2]: c1=%.3f c2=%.3f\n", fit$beta[2,1], fit$beta[2,2]))

  ## Compare v2 separate-baseline EM
  fit_v2 <- tryCatch(
    gemcox_full_multistart(
      X_gmm = dat$X_gmm, X_cox = dat$X_cox,
      time = dat$time, status = dat$status,
      K = 2, lambda = 0.05, alpha = 0.5, max_iter = 100,
      surv_weight = 1.0, verbose = FALSE, n_starts = 20
    ),
    error = function(e) { cat("  v2 ERROR:", e$message, "\n"); NULL }
  )
  if (!is.null(fit_v2)) {
    ari_v2 <- adjustedRandIndex(dat$Z_true, fit_v2$clusterid)
    cat(sprintf("\n  v2 separate EM ARI: %.4f\n", ari_v2))
    cat(sprintf("  v2 Est beta[1]: c1=%.3f c2=%.3f\n", fit_v2$beta[1,1], fit_v2$beta[1,2]))
  }
}

## ============================================================
## TEST 3: Sweep gamma values to find optimal
## ============================================================
cat("\n=== TEST 3: Gamma sweep on K=2 clean data ===\n")

for (gam in c(0, 0.25, 0.5, 1.0, 2.0, 5.0)) {
  fit_g <- tryCatch(
    gemcox_full_multistart_shared_v4(
      X_gmm = dat$X_gmm, X_cox = dat$X_cox,
      time = dat$time, status = dat$status,
      K = 2, alpha = 0.5, max_iter = 80,
      surv_weight = gam, verbose = FALSE, n_starts = 10
    ),
    error = function(e) NULL
  )
  if (!is.null(fit_g)) {
    ari_g <- adjustedRandIndex(dat$Z_true, fit_g$clusterid)
    cat(sprintf("  gamma=%.2f -> ARI=%.4f  pi=[%.2f,%.2f]  beta1=[%.3f,%.3f]\n",
                gam, ari_g, fit_g$pi[1], fit_g$pi[2],
                fit_g$beta[1,1], fit_g$beta[1,2]))
  } else {
    cat(sprintf("  gamma=%.2f -> FAILED\n", gam))
  }
}

## ============================================================
## TEST 4: Replicate the full CV pipeline on one clean dataset
## ============================================================
cat("\n=== TEST 4: Full v4 CV pipeline on K=2 clean data, n=117 ===\n")

dat117 <- simulate_clean(n = 117, K = 2, delta_mu = 1.5, delta_beta = 0.6, seed = 42)
cat(sprintf("  Events: %d/%d = %.0f%%\n", sum(dat117$status), length(dat117$status),
            100*mean(dat117$status)))

cv_res <- tryCatch(
  cv_select_K_gamma_v4(
    X_gmm = dat117$X_gmm, X_cox = dat117$X_cox,
    time = dat117$time, status = dat117$status,
    K_grid = 1:3, gamma_grid = c(0, 0.5, 1.0, 2.0),
    nfolds = 3, alpha = 0.5, max_iter = 40,
    n_starts = 5, log_transform = FALSE, verbose = FALSE,
    use_1se_rule = TRUE, baseline_mode = "shared"
  ),
  error = function(e) e
)

if (!inherits(cv_res, "error")) {
  s <- cv_res$summary[order(-cv_res$summary$mean_score), ]
  s <- s[is.finite(s$mean_score), ]
  cat("  CV scores (top 10):\n")
  print(head(s[, c("K", "gamma", "mean_score", "se_score", "n_ok_folds")], 10))
  cat(sprintf("\n  Selected: K=%d gamma=%.2f (method: %s)\n",
              cv_res$best$K, cv_res$best$gamma, cv_res$best_method))
  cat(sprintf("  Correct K? %s\n", ifelse(cv_res$best$K == 2, "YES", "NO")))
} else {
  cat("  CV FAILED:", cv_res$message, "\n")
}

cat("\n=== DIAGNOSTIC COMPLETE ===\n")
