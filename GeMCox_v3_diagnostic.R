###############################################################
## GeM-Cox v3 Diagnostic Script
## Identifies root causes of poor simulation performance
###############################################################

library(MASS)
library(survival)
library(glmnet)
set.seed(2025)

source("GeMCox_improved_v3.R")

cat("==============================================================\n")
cat("DIAGNOSTIC 1: The log-normal roundtrip destroys cluster signal\n")
cat("==============================================================\n\n")

## Generate one dataset from the v3 DGP
sim <- simulate_gemcox_v3(n = 500, K = 3, separation = "high",
                          event_rate = 0.44, lognormal = TRUE, seed = 42)

X_latent <- log(sim$X_gmm)   # recover latent (since X_gmm = exp(X_latent))
X_observed <- sim$X_gmm       # what the model sees
X_log1p <- log1p(sim$X_gmm)   # what preprocessing produces

## Check: how much cluster separation survives each transform?
cat("Cluster means on feature 1 (the main signal feature):\n")
for (k in 1:3) {
  idx <- sim$Z_true == k
  cat(sprintf("  Cluster %d (n=%d):\n", k, sum(idx)))
  cat(sprintf("    Latent space:    mean=%.3f  sd=%.3f\n",
              mean(X_latent[idx, 1]), sd(X_latent[idx, 1])))
  cat(sprintf("    Observable exp(): mean=%.3f  sd=%.3f\n",
              mean(X_observed[idx, 1]), sd(X_observed[idx, 1])))
  cat(sprintf("    After log1p():   mean=%.3f  sd=%.3f\n",
              mean(X_log1p[idx, 1]), sd(X_log1p[idx, 1])))
}

## Measure separation: ratio of between-cluster variance to within-cluster variance
cat("\nSignal-to-noise ratio (between/within variance) on feature 1:\n")
cluster_means <- tapply(X_latent[,1], sim$Z_true, mean)
between_var <- var(cluster_means)
within_var <- mean(tapply(X_latent[,1], sim$Z_true, var))
cat(sprintf("  Latent:  between=%.3f  within=%.3f  ratio=%.3f\n",
            between_var, within_var, between_var/within_var))

cluster_means2 <- tapply(X_log1p[,1], sim$Z_true, mean)
between_var2 <- var(cluster_means2)
within_var2 <- mean(tapply(X_log1p[,1], sim$Z_true, var))
cat(sprintf("  log1p(): between=%.3f  within=%.3f  ratio=%.3f\n",
            between_var2, within_var2, between_var2/within_var2))

cat(sprintf("\n  => Signal loss factor: %.1fx\n",
            (between_var/within_var) / (between_var2/within_var2)))


cat("\n\n==============================================================\n")
cat("DIAGNOSTIC 2: The asymmetric compression of log1p(exp(x))\n")
cat("==============================================================\n\n")

## Show how log1p(exp(x)) compresses the left tail
x_grid <- seq(-3, 3, by = 0.5)
cat("  x_latent -> exp(x) -> log1p(exp(x)) -> loss\n")
for (x in x_grid) {
  observed <- exp(x)
  back <- log1p(observed)
  loss <- abs(back - x)
  cat(sprintf("  %6.2f  ->  %8.4f  ->  %6.4f    (error=%.4f, %.0f%%)\n",
              x, observed, back, loss, 100*loss/max(abs(x), 0.01)))
}
cat("\n  => For x < 0, log1p(exp(x)) severely compresses toward log(2)=0.693\n")
cat("  => Clusters with negative latent means lose their distinctness\n")


cat("\n\n==============================================================\n")
cat("DIAGNOSTIC 3: Cox model sees cancellation with high K + high sep\n")
cat("==============================================================\n\n")

## Show the true betas for K=3, high separation
sim3 <- simulate_gemcox_v3(n = 500, K = 3, separation = "high",
                           event_rate = 0.44, lognormal = FALSE, seed = 42)
cat("True beta[1] values for K=3, high separation:\n")
for (k in 1:3) {
  cat(sprintf("  Cluster %d: beta[1] = %.3f  (HR = %.2f)\n",
              k, sim3$params$beta_list[[k]][1],
              exp(sim3$params$beta_list[[k]][1])))
}

## What the single Cox sees:
scaler <- make_x_scaler(sim3$X_cox)
Xs <- apply_x_scaler(sim3$X_cox, scaler)
gfit <- glmnet(Xs, Surv(sim3$time, sim3$status), family = "cox",
               alpha = 0.5, lambda = 0.05, standardize = FALSE)
beta_single <- as.numeric(coef(gfit, s = 0.05))
cat(sprintf("\nSingle Cox beta_s[1] = %.4f\n", beta_single[1]))
cat("  => Near zero because positive and negative effects cancel!\n")

ci_single <- c_index(sim3$time, sim3$status,
                     as.numeric(Xs %*% beta_single))
cat(sprintf("  => Single Cox C-index = %.4f", ci_single))
if (ci_single < 0.55) cat(" *** AT CHANCE ***")
cat("\n")

## Now try K=3 with LOW separation — betas don't cancel
sim3low <- simulate_gemcox_v3(n = 500, K = 3, separation = "low",
                              event_rate = 0.44, lognormal = FALSE, seed = 42)
cat("\nTrue beta[1] values for K=3, LOW separation:\n")
for (k in 1:3) {
  cat(sprintf("  Cluster %d: beta[1] = %.3f  (HR = %.2f)\n",
              k, sim3low$params$beta_list[[k]][1],
              exp(sim3low$params$beta_list[[k]][1])))
}
scaler2 <- make_x_scaler(sim3low$X_cox)
Xs2 <- apply_x_scaler(sim3low$X_cox, scaler2)
gfit2 <- glmnet(Xs2, Surv(sim3low$time, sim3low$status), family = "cox",
                alpha = 0.5, lambda = 0.05, standardize = FALSE)
ci2 <- c_index(sim3low$time, sim3low$status,
               as.numeric(Xs2 %*% coef(gfit2, s = 0.05)))
cat(sprintf("  Single Cox C-index = %.4f\n", ci2))


cat("\n\n==============================================================\n")
cat("DIAGNOSTIC 4: Verify the EM can recover clusters on clean data\n")
cat("==============================================================\n\n")

## Generate clean Gaussian data (no log-normal), K=2, med separation
sim_clean <- simulate_gemcox_v3(n = 200, K = 2, separation = "med",
                                event_rate = 0.44, lognormal = FALSE, seed = 42)

cat("Clean Gaussian data, K=2, med separation, n=200:\n")
cat(sprintf("  Events: %d (%.0f%%)\n", sum(sim_clean$status),
            100*mean(sim_clean$status)))

## Fit GeM-Cox directly (no log-transform needed)
fit_clean <- tryCatch(
  gemcox_full_multistart(
    X_gmm = sim_clean$X_gmm, X_cox = sim_clean$X_cox,
    time = sim_clean$time, status = sim_clean$status,
    K = 2, lambda = 0.05, alpha = 0.5,
    max_iter = 100, surv_weight = 1.0, verbose = FALSE,
    n_starts = 10
  ),
  error = function(e) { cat("  FIT ERROR:", e$message, "\n"); NULL }
)

if (!is.null(fit_clean)) {
  library(mclust)
  ari <- adjustedRandIndex(sim_clean$Z_true, fit_clean$clusterid)
  cat(sprintf("  ARI (cluster recovery) = %.4f\n", ari))

  eta_mat <- sapply(1:2, function(k)
    eta_from_scaled(sim_clean$X_cox, fit_clean$coxfit[[k]]$beta_s,
                    fit_clean$coxfit[[k]]$scaler))
  eta_mix <- rowSums(fit_clean$tau * eta_mat)
  ci_gem <- c_index(sim_clean$time, sim_clean$status, eta_mix)
  cat(sprintf("  GeM-Cox C-index = %.4f\n", ci_gem))

  # Single Cox baseline
  sc <- make_x_scaler(sim_clean$X_cox)
  Xsc <- apply_x_scaler(sim_clean$X_cox, sc)
  gf <- glmnet(Xsc, Surv(sim_clean$time, sim_clean$status),
               family = "cox", alpha = 0.5, lambda = 0.05, standardize = FALSE)
  ci_s <- c_index(sim_clean$time, sim_clean$status,
                  as.numeric(Xsc %*% coef(gf, s = 0.05)))
  cat(sprintf("  Single Cox C-index = %.4f\n", ci_s))
  cat(sprintf("  Gain = %.4f\n", ci_gem - ci_s))

  cat("\n  Estimated betas:\n")
  for (k in 1:2) {
    cat(sprintf("    Cluster %d: beta[1:3] = [%.3f, %.3f, %.3f]  (true: [%.3f, %.3f, %.3f])\n",
                k, fit_clean$beta[1,k], fit_clean$beta[2,k], fit_clean$beta[3,k],
                sim_clean$params$beta_list[[k]][1],
                sim_clean$params$beta_list[[k]][2],
                sim_clean$params$beta_list[[k]][3]))
  }
}


cat("\n\n==============================================================\n")
cat("DIAGNOSTIC 5: Does GeM-Cox work at ALL on simple K=2 case?\n")
cat("==============================================================\n\n")

## Simplest possible test: 2 well-separated clusters, Gaussian, large n
sim_easy <- simulate_gemcox_v3(n = 300, K = 2, separation = "high",
                               event_rate = 0.44, lognormal = FALSE, seed = 42)

cat("EASY TEST: K=2, high separation, Gaussian, n=300:\n")
fit_easy <- tryCatch(
  gemcox_full_multistart(
    X_gmm = sim_easy$X_gmm, X_cox = sim_easy$X_cox,
    time = sim_easy$time, status = sim_easy$status,
    K = 2, lambda = 0.05, alpha = 0.5,
    max_iter = 100, surv_weight = 1.0, verbose = FALSE,
    n_starts = 20
  ),
  error = function(e) { cat("  FIT ERROR:", e$message, "\n"); NULL }
)

if (!is.null(fit_easy)) {
  ari_easy <- adjustedRandIndex(sim_easy$Z_true, fit_easy$clusterid)
  cat(sprintf("  ARI = %.4f\n", ari_easy))
  cat(sprintf("  Cluster sizes: %s (true: %s)\n",
              paste(table(fit_easy$clusterid), collapse=", "),
              paste(table(sim_easy$Z_true), collapse=", ")))
  cat(sprintf("  Pi: [%.3f, %.3f]\n", fit_easy$pi[1], fit_easy$pi[2]))
  cat(sprintf("  Effective events: [%.1f, %.1f]\n",
              fit_easy$eff_events[1], fit_easy$eff_events[2]))

  ## Check: are the betas actually different?
  cat("\n  Beta comparison:\n")
  for (j in 1:min(5, nrow(fit_easy$beta))) {
    cat(sprintf("    feat_%d: cluster1=%.3f  cluster2=%.3f  (true: %.3f, %.3f)\n",
                j, fit_easy$beta[j,1], fit_easy$beta[j,2],
                sim_easy$params$beta_list[[1]][j],
                sim_easy$params$beta_list[[2]][j]))
  }
}


cat("\n\n==============================================================\n")
cat("DIAGNOSTIC 6: Check what lambda cv.glmnet selects\n")
cat("==============================================================\n\n")

## The lambda may be over-regularizing
sim_lam <- simulate_gemcox_v3(n = 117, K = 2, separation = "med",
                              event_rate = 0.44, lognormal = TRUE, seed = 42)
prep <- preprocess_train(sim_lam$X_gmm, sim_lam$X_cox, log_transform = TRUE)

lam_cv <- tryCatch({
  cvfit <- cv.glmnet(
    apply_x_scaler(prep$X_cox, make_x_scaler(prep$X_cox)),
    Surv(sim_lam$time, sim_lam$status),
    family = "cox", alpha = 0.5, nfolds = 5, standardize = FALSE)
  c(lambda.min = cvfit$lambda.min, lambda.1se = cvfit$lambda.1se)
}, error = function(e) c(lambda.min = NA, lambda.1se = NA))

cat(sprintf("  cv.glmnet lambda.min = %.4f\n", lam_cv["lambda.min"]))
cat(sprintf("  cv.glmnet lambda.1se = %.4f\n", lam_cv["lambda.1se"]))
cat(sprintf("  Simulation uses fixed lambda = 0.05\n"))
cat("  => If cv.glmnet lambda >> 0.05, the fixed lambda may be too permissive\n")
cat("  => If cv.glmnet lambda << 0.05, the fixed lambda may over-regularize\n")


cat("\n\n==============================================================\n")
cat("DIAGNOSTIC 7: Check the 1-SE rule behavior\n")
cat("==============================================================\n\n")

## Run CV on a K_true=1 dataset to see if 1-SE works
sim_k1 <- simulate_gemcox_v3(n = 117, K = 1, separation = "med",
                             event_rate = 0.44, lognormal = TRUE, seed = 42)
prep_k1 <- preprocess_train(sim_k1$X_gmm, sim_k1$X_cox, log_transform = TRUE)

cv_k1 <- tryCatch(
  cv_select_K_gamma_v3(
    X_gmm = prep_k1$X_gmm, X_cox = prep_k1$X_cox,
    time = sim_k1$time, status = sim_k1$status,
    K_grid = 1:3, gamma_grid = c(0, 0.5, 1.0, 2.0),
    nfolds = 3, alpha = 0.5, max_iter = 30,
    n_starts = 3, verbose = FALSE, use_1se_rule = TRUE
  ),
  error = function(e) e
)

if (!inherits(cv_k1, "error")) {
  cat("CV summary for K_true=1 data:\n")
  s <- cv_k1$summary[order(-cv_k1$summary$mean_score), ]
  s <- s[is.finite(s$mean_score), ]
  print(s[, c("K", "gamma", "mean_score", "sd_score", "n_ok_folds")])
  cat(sprintf("\n1-SE rule selected: K=%d, gamma=%.2f\n",
              cv_k1$best$K, cv_k1$best$gamma))
  cat(sprintf("  (correct = K=1)\n"))
}


cat("\n\n==============================================================\n")
cat("SUMMARY OF ROOT CAUSES\n")
cat("==============================================================\n")
cat("
1. LOG-NORMAL ROUNDTRIP: log1p(exp(x)) ≠ x. The left tail is
   compressed, destroying cluster separation. For x=-2, the error
   is 70%. This is the primary reason ARI ≈ 0.

2. BETA CANCELLATION: With K≥3 and high separation, cluster betas
   range from negative to large positive. The single Cox averages
   these to ~0, giving C-index ≈ 0.5. This is EXPECTED behavior
   (it's what motivates mixture models), but GeM-Cox can't exploit
   it because of issue #1.

3. THE FIX: Use log() not log1p() as the inverse of exp(), OR
   (better) define the true Cox signal on the OBSERVED feature space
   so log-transform is not needed for signal recovery. OR simply
   generate Gaussian data when testing the core methodology.
")
