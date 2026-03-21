###############################################################
## CVIA078 Analysis v2: Proper feature reduction + CV safeguards
##
## Root cause of v1 failure:
##   44 Cox features vs 52 events → lambda.1se shrinks all to 0
##   → single Cox C-index = 0.50 (random baseline)
##   → K=3 γ=5 selected on 1 lucky fold (SE=0 artefact)
##
## Fixes:
##   1. Reduce Cox features to ≤15 via univariate screening
##   2. Require n_ok_folds ≥ 3 for CV selection
##   3. Test multiple GMM feature sets (6, 10, 16)
##   4. Paired gamma=1 vs tuned on each configuration
###############################################################

library(survival); library(glmnet); library(MASS)
library(dplyr); library(tidyr); library(mclust)

source("GeMCox_improved_v4.R")
source("GeMCox_improved_v4_revisions_patch.R")

## ---- Load CVIA078 ----
dat_long <- readxl::read_excel("CVIA078.xlsx")
prep_full <- prepare_cvia078_data_v4(
  dat_long,
  gmm_feature_mode  = "delta_only",
  cox_feature_mode  = "v20_plus_delta_covars",
  missing_threshold = 0.10,
  corr_prune        = 0.85,
  add_covars        = TRUE
)

cat("=== FULL DATA ===\n")
cat(sprintf("  n=%d, events=%d, GMM features=%d, Cox features=%d\n",
            nrow(prep_full$X_gmm), sum(prep_full$status),
            ncol(prep_full$X_gmm), ncol(prep_full$X_cox)))

time   <- prep_full$time
status <- prep_full$status
n      <- length(time)
n_ev   <- sum(status)


###############################################################
## 1. Univariate feature screening for Cox features
###############################################################
cat("\n=== UNIVARIATE SCREENING ===\n")

## Log-transform + impute all features first
prep0 <- preprocess_train_v4(prep_full$X_gmm, prep_full$X_cox,
                             log_transform = TRUE)
X_all <- prep0$X_cox

## Univariate Cox p-values
uv_pvals <- numeric(ncol(X_all))
uv_cindex <- numeric(ncol(X_all))
names(uv_pvals) <- colnames(X_all)
names(uv_cindex) <- colnames(X_all)

for (j in seq_len(ncol(X_all))) {
  xj <- X_all[, j]
  if (sd(xj, na.rm=TRUE) < 1e-10) {
    uv_pvals[j] <- 1; uv_cindex[j] <- 0.5; next
  }
  fit_j <- tryCatch(coxph(Surv(time, status) ~ xj), error=function(e) NULL)
  if (is.null(fit_j)) {
    uv_pvals[j] <- 1; uv_cindex[j] <- 0.5
  } else {
    uv_pvals[j] <- summary(fit_j)$coefficients[, "Pr(>|z|)"]
    uv_cindex[j] <- tryCatch(c_index(time, status, xj), error=function(e) 0.5)
  }
}

## Sort by p-value
uv_df <- data.frame(
  feature = names(uv_pvals),
  pval    = uv_pvals,
  cindex  = uv_cindex,
  stringsAsFactors = FALSE
) %>% arrange(pval)

cat("  Top 20 features by univariate Cox p-value:\n")
print(head(uv_df, 20), row.names=FALSE)

## Select top features
## Rule: events-per-variable ≥ 5 → max features = 52/5 ≈ 10
## Use top 10 + 3 covariates (AGE, SEX_M, SICKLE_ABNORMAL) = 13
N_COX_TOP <- 10
covar_cols <- intersect(c("AGE", "SEX_M", "SICKLE_ABNORMAL"), colnames(X_all))
top_immune <- head(uv_df$feature[!uv_df$feature %in% covar_cols], N_COX_TOP)
cox_selected <- unique(c(top_immune, covar_cols))

cat(sprintf("\n  Selected %d Cox features (top %d immune + %d covariates):\n",
            length(cox_selected), length(top_immune), length(covar_cols)))
cat(paste("   ", cox_selected, collapse="\n"), "\n")
cat(sprintf("  Events-per-variable: %.1f\n", n_ev / length(cox_selected)))


###############################################################
## 2. Single Cox with reduced features
###############################################################
cat("\n=== SINGLE COX WITH REDUCED FEATURES ===\n")

X_cox_sel <- prep_full$X_cox[, intersect(cox_selected, colnames(prep_full$X_cox)),
                             drop=FALSE]
## Add covariates that may be in X_cox but not in top_immune
for (cv in covar_cols) {
  if (!cv %in% colnames(X_cox_sel) && cv %in% colnames(prep_full$X_cox))
    X_cox_sel <- cbind(X_cox_sel, prep_full$X_cox[, cv, drop=FALSE])
}

prep_sel <- preprocess_train_v4(prep_full$X_gmm, X_cox_sel, log_transform=TRUE)
sc <- make_x_scaler(prep_sel$X_cox)
Xs <- apply_x_scaler(prep_sel$X_cox, sc)
cvfit1 <- cv.glmnet(Xs, Surv(time, status), family="cox",
                    alpha=0.5, standardize=FALSE, nfolds=5)
eta1 <- as.numeric(predict(cvfit1, Xs, s="lambda.1se"))
ci_cox1 <- c_index(time, status, eta1)
cat(sprintf("  C-index: %.4f (was 0.5000 with 44 features)\n", ci_cox1))
cat(sprintf("  lambda.1se: %.4f (was 0.2368)\n", cvfit1$lambda.1se))
cat(sprintf("  Nonzero coefficients: %d / %d\n",
            sum(coef(cvfit1, s="lambda.1se") != 0), ncol(Xs)))


###############################################################
## 3. GeM-Cox with 3 GMM feature set sizes
###############################################################

## Helper: safe CV + refit with minimum fold requirement
run_gemcox <- function(X_gmm, X_cox, gamma_grid, label,
                       min_ok_folds = 3, K_grid = 1:2) {
  cat(sprintf("\n--- %s (q=%d, p=%d, gamma=[%s]) ---\n",
              label, ncol(X_gmm), ncol(X_cox),
              paste(gamma_grid, collapse=",")))
  
  cv <- tryCatch(cv_select_K_gamma_v4(
    X_gmm=X_gmm, X_cox=X_cox, time=time, status=status,
    K_grid=K_grid, gamma_grid=gamma_grid,
    nfolds=5, alpha=0.5, max_iter=60, n_starts=8,
    log_transform=TRUE, normalize_gmm_by_dim=FALSE,
    verbose=FALSE, use_1se_rule=TRUE,
    baseline_mode="shared", selection_policy="screened_1se"
  ), error=function(e) { cat("  CV ERROR:", e$message, "\n"); NULL })
  
  if (is.null(cv)) return(NULL)
  
  ## SAFEGUARD: require ≥ min_ok_folds
  s <- cv$summary
  s_safe <- s[!is.na(s$n_ok_folds) & s$n_ok_folds >= min_ok_folds, , drop=FALSE]
  
  if (nrow(s_safe) == 0) {
    cat("  No configuration passed screening on ≥", min_ok_folds, "folds\n")
    ## Fall back to K=1
    s_safe <- s[s$K == 1 & is.finite(s$mean_score), , drop=FALSE]
    if (nrow(s_safe) == 0) return(NULL)
  }
  
  ## Apply 1-SE rule on safe configurations only
  best_idx <- which.max(s_safe$mean_score)
  best_score <- s_safe$mean_score[best_idx]
  best_se <- s_safe$se_score[best_idx]
  if (!is.finite(best_se)) best_se <- 0
  threshold <- best_score - best_se
  
  eligible <- s_safe[s_safe$mean_score >= threshold, , drop=FALSE]
  eligible <- eligible[order(eligible$K, eligible$gamma), , drop=FALSE]
  best <- eligible[1, , drop=FALSE]
  
  cat(sprintf("  CV selected: K=%d, gamma=%.1f (CV C-index=%.4f, SE=%.4f, folds=%d)\n",
              best$K, best$gamma, best$mean_score, best$se_score, best$n_ok_folds))
  
  ## Print full CV table
  cat("  Full CV table (safe rows, sorted by score):\n")
  s_print <- s_safe[order(-s_safe$mean_score), ]
  print(s_print[, c("K","gamma","mean_score","se_score","n_ok_folds")], row.names=FALSE)
  
  ## Refit
  prep <- preprocess_train_v4(X_gmm, X_cox, log_transform=TRUE)
  fit <- tryCatch(gemcox_full_multistart_shared_v4(
    X_gmm=prep$X_gmm, X_cox=prep$X_cox,
    time=time, status=status,
    K=best$K, alpha=0.5, max_iter=120,
    surv_weight=best$gamma, normalize_gmm_by_dim=FALSE,
    verbose=FALSE, n_starts=15
  ), error=function(e) { cat("  FIT ERROR:", e$message, "\n"); NULL })
  
  if (is.null(fit)) return(NULL)
  
  eta <- tryCatch({
    tau<-fit$tau; em<-sapply(1:fit$K,function(k)
      eta_from_scaled(prep$X_cox,fit$coxfit[[k]]$beta_s,fit$coxfit[[k]]$scaler))
    if(is.vector(em))em<-matrix(em,ncol=fit$K)
    rowSums(tau*em)},error=function(e)NULL)
  
  ci <- if(!is.null(eta)) c_index(time, status, eta) else NA
  cat(sprintf("  Full-data C-index: %.4f (gain: %+.4f)\n", ci, ci - ci_cox1))
  cat(sprintf("  Cluster sizes: %s, Pi: %s\n",
              paste(table(fit$clusterid), collapse="/"),
              paste(round(fit$pi, 3), collapse=", ")))
  cat(sprintf("  Eff events: %s\n",
              paste(round(fit$eff_events, 1), collapse=", ")))
  
  ## Top betas per cluster
  if (fit$K > 1) {
    cat("  Top betas (|beta| > 0.05):\n")
    for (k in 1:fit$K) {
      bk <- fit$beta[, k]
      big <- which(abs(bk) > 0.05)
      if (length(big) > 0) {
        big <- big[order(-abs(bk[big]))]
        cat(sprintf("    Cluster %d: %s\n", k,
                    paste(sprintf("%s=%.3f", names(bk[big]), bk[big]), collapse=", ")))
      }
    }
  }
  
  list(cv=cv, best=best, fit=fit, ci=ci, label=label)
}


###############################################################
## Run configurations
###############################################################

## GMM feature subsets
gmm_all   <- prep_full$gmm_features  # 16 features
## Top 6 by univariate signal among GMM features
gmm_uv    <- uv_df$feature[uv_df$feature %in% gmm_all]
gmm_top6  <- head(gmm_uv, 6)
gmm_top10 <- head(gmm_uv, 10)

cat(sprintf("\n=== GMM FEATURE SETS ===\n"))
cat(sprintf("  Full:  %d features\n", length(gmm_all)))
cat(sprintf("  Top10: %d features [%s]\n", length(gmm_top10),
            paste(head(gmm_top10, 5), collapse=", ")))
cat(sprintf("  Top6:  %d features [%s]\n", length(gmm_top6),
            paste(gmm_top6, collapse=", ")))

results <- list()

## --- Config 1: Full GMM (q=16), reduced Cox (p~13) ---
cat("\n\n########## CONFIG 1: q=16 GMM, reduced Cox ##########\n")
res1A <- run_gemcox(prep_full$X_gmm, X_cox_sel, gamma_grid=1.0,
                    label="q16_gamma1")
res1B <- run_gemcox(prep_full$X_gmm, X_cox_sel,
                    gamma_grid=c(0.5, 1.0, 2.0, 5.0),
                    label="q16_tuned")

## --- Config 2: Top 10 GMM (q=10), reduced Cox ---
X_gmm_10 <- prep_full$X_gmm[, gmm_top10[gmm_top10 %in% colnames(prep_full$X_gmm)],
                            drop=FALSE]
cat("\n\n########## CONFIG 2: q=10 GMM, reduced Cox ##########\n")
res2A <- run_gemcox(X_gmm_10, X_cox_sel, gamma_grid=1.0,
                    label="q10_gamma1")
res2B <- run_gemcox(X_gmm_10, X_cox_sel,
                    gamma_grid=c(0.5, 1.0, 2.0, 5.0),
                    label="q10_tuned")

## --- Config 3: Top 6 GMM (q=6), reduced Cox ---
X_gmm_6 <- prep_full$X_gmm[, gmm_top6[gmm_top6 %in% colnames(prep_full$X_gmm)],
                           drop=FALSE]
cat("\n\n########## CONFIG 3: q=6 GMM, reduced Cox ##########\n")
res3A <- run_gemcox(X_gmm_6, X_cox_sel, gamma_grid=1.0,
                    label="q6_gamma1")
res3B <- run_gemcox(X_gmm_6, X_cox_sel,
                    gamma_grid=c(0.5, 1.0, 2.0, 5.0),
                    label="q6_tuned")


###############################################################
## Summary table
###############################################################
cat("\n\n========================================\n")
cat("CVIA078 SUMMARY\n")
cat("========================================\n\n")
cat(sprintf("  Single Cox (p=%d): C-index = %.4f\n", ncol(X_cox_sel), ci_cox1))

make_row <- function(res) {
  if (is.null(res)) return(data.frame(label=NA, K=NA, gamma=NA,
                                      cv_ci=NA, full_ci=NA, gain=NA))
  data.frame(
    label   = res$label,
    K       = res$best$K,
    gamma   = res$best$gamma,
    cv_ci   = round(res$best$mean_score, 4),
    full_ci = round(res$ci, 4),
    gain    = round(res$ci - ci_cox1, 4),
    stringsAsFactors = FALSE
  )
}

summary_tab <- rbind(
  make_row(res1A), make_row(res1B),
  make_row(res2A), make_row(res2B),
  make_row(res3A), make_row(res3B)
)
summary_tab <- summary_tab[!is.na(summary_tab$label), ]
print(summary_tab, row.names=FALSE)

cat("\n========================================\n")
cat("DONE\n")
cat("========================================\n")
