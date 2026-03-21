###############################################################
## Gamma Ablation v3: Signal Alignment
##
## Three conditions (same p=10, same n, same seeds):
##   ALIGNED:  features 1-5 drive BOTH GMM means and Cox betas
##   PARTIAL:  feature 1 shared; 2-3 GMM-only; 4-6 Cox-only
##   MISALIGNED: features 1-3 GMM-only; 4-6 Cox-only; no overlap
##
## Hypothesis: gamma tuning helps when GMM and Cox signals
## are on DIFFERENT features, because gamma>1 steers the EM
## toward survival-relevant clusters that the GMM alone misses.
##
## Target: ~8 min on 24 cores
###############################################################

required_pkgs <- c("survival","glmnet","mclust","dplyr",
                   "tidyr","tibble","MASS","parallel")
for (pkg in required_pkgs)
  if (!requireNamespace(pkg, quietly=TRUE))
    install.packages(pkg, repos="https://cloud.r-project.org")

library(survival); library(glmnet); library(mclust)
library(dplyr); library(tidyr); library(tibble)
library(MASS); library(parallel)
select <- dplyr::select; filter <- dplyr::filter
mutate <- dplyr::mutate; arrange <- dplyr::arrange

source("GeMCox_improved_v4.R")
set.seed(2025)

###############################################################
## DGP with configurable signal alignment
###############################################################
simulate_alignment <- function(
    n = 117, K = 2, p = 10,
    alignment = c("aligned", "partial", "misaligned"),
    separation = c("med", "high"),
    rho = 0.2, event_rate = 0.44, seed = 1
) {
  alignment  <- match.arg(alignment)
  separation <- match.arg(separation)
  set.seed(seed)
  
  dmu   <- switch(separation, low=0.6, med=1.2, high=2.0)
  dbeta <- switch(separation, low=0.3, med=0.6, high=1.0)
  pi_c  <- rep(1/K, K)
  pos   <- seq(-(K-1)/2, (K-1)/2, length.out=K)
  
  ## ---- Define which features carry which signal ----
  ##
  ## "aligned":    features 1-5 carry both GMM means and Cox heterogeneity
  ## "partial":    feature 1 shared; 2-3 GMM-only means; 4-6 Cox-only heterogeneity
  ## "misaligned": features 1-3 GMM-only means; 4-6 Cox-only heterogeneity; NO overlap
  
  ## GMM cluster means: which features get nonzero mu shifts
  gmm_signal_idx <- switch(alignment,
                           aligned    = 1:5,
                           partial    = 1:3,    # features 1,2,3 have cluster mean differences
                           misaligned = 1:3     # features 1,2,3 have cluster mean differences
  )
  
  ## Cox heterogeneity: which features get cluster-specific beta shifts
  cox_hetero_idx <- switch(alignment,
                           aligned    = 1:3,    # same as GMM signal features
                           partial    = c(1, 4, 5, 6),  # feature 1 shared; 4,5,6 Cox-only
                           misaligned = 4:6     # completely different from GMM features
  )
  
  ## Base Cox effect: features with nonzero beta across ALL clusters
  beta_base <- rep(0, p)
  ## In all conditions, features 1-6 have some baseline effect
  beta_base[1:6] <- c(0.5, -0.3, 0.2, 0.4, -0.3, 0.2)
  
  ## Cluster means
  mu_list <- lapply(1:K, function(k) {
    m <- rep(0, p)
    for (j in gmm_signal_idx) {
      rank_in_signal <- which(gmm_signal_idx == j)
      m[j] <- pos[k] * dmu / sqrt(rank_in_signal)
    }
    m
  })
  
  ## Cluster-specific betas
  beta_list <- lapply(1:K, function(k) {
    b <- beta_base
    for (j in cox_hetero_idx) {
      rank_in_hetero <- which(cox_hetero_idx == j)
      b[j] <- b[j] + dbeta * pos[k] / sqrt(rank_in_hetero)
    }
    b
  })
  
  Sigma <- outer(1:p, 1:p, function(i,j) rho^abs(i-j))
  wscales <- 200 * exp(0.2 * pos)
  
  ## Generate data
  Z <- sample(1:K, n, replace=TRUE, prob=pi_c)
  X <- matrix(0, n, p)
  for (k in 1:K) {
    i <- which(Z == k); if (!length(i)) next
    X[i,] <- mvrnorm(length(i), mu_list[[k]], Sigma)
  }
  colnames(X) <- paste0("f", 1:p)
  
  ## Survival times from observable X
  Tt <- numeric(n)
  for (k in 1:K) {
    i <- which(Z == k); if (!length(i)) next
    eta <- pmin(pmax(as.numeric(X[i,,drop=F] %*% beta_list[[k]]), -5), 5)
    Tt[i] <- wscales[k] * (-log(runif(length(i))))^(1/1.5) * exp(-eta/1.5)
  }
  
  ## Calibrate censoring
  lo<-1e-6; hi<-10
  for (it in 1:50) {
    mid <- (lo+hi)/2
    er <- mean(Tt <= pmin(rexp(n, mid), 253))
    if (abs(er - event_rate) < 0.005) break
    if (er < event_rate) hi <- mid else lo <- mid
  }
  C <- pmin(rexp(n, mid), 253)
  time <- pmax(pmin(Tt, C) + runif(n, 0, 0.01), 0.1)
  status <- as.integer(Tt <= C)
  
  list(
    X_gmm = X, X_cox = X,
    time = time, status = status, Z_true = Z,
    params = list(K=K, alignment=alignment, separation=separation,
                  gmm_signal_idx=gmm_signal_idx,
                  cox_hetero_idx=cox_hetero_idx,
                  beta_list=beta_list, mu_list=mu_list,
                  er_actual=mean(status))
  )
}


###############################################################
## Paired run: fixed gamma=1 vs tuned gamma
###############################################################
run_pair <- function(sd, K_grid=1:3, nfolds=3, max_iter=30,
                     n_starts=3, n_starts_final=5, alpha=0.5) {
  K_true <- sd$params$K
  
  cox1 <- tryCatch({
    sc <- make_x_scaler(sd$X_cox); Xs <- apply_x_scaler(sd$X_cox, sc)
    lam <- tryCatch(cv.glmnet(Xs, Surv(sd$time, sd$status), family="cox",
                              alpha=alpha, standardize=FALSE, nfolds=3)$lambda.1se, error=function(e)0.05)
    gf <- glmnet(Xs, Surv(sd$time, sd$status), family="cox",
                 alpha=alpha, lambda=lam, standardize=FALSE)
    c_index(sd$time, sd$status, as.numeric(Xs %*% coef(gf, s=lam)))
  }, error=function(e) NA_real_)
  
  run_one <- function(gamma_grid, label) {
    r <- list(method=label, K_true=K_true, K_selected=NA_integer_,
              K_correct=NA_integer_, gamma_selected=NA_real_,
              ARI=NA_real_, cindex_full=NA_real_, cindex_cv=NA_real_,
              cindex_cox1=cox1, converged=NA_integer_, error_msg=NA_character_)
    
    cv <- tryCatch(cv_select_K_gamma_v4(
      X_gmm=sd$X_gmm, X_cox=sd$X_cox,
      time=sd$time, status=sd$status,
      K_grid=K_grid, gamma_grid=gamma_grid,
      nfolds=nfolds, alpha=alpha, max_iter=max_iter, n_starts=n_starts,
      log_transform=FALSE, normalize_gmm_by_dim=FALSE,
      verbose=FALSE, use_1se_rule=TRUE, baseline_mode="shared"
    ), error=function(e) e)
    if (inherits(cv,"error")) { r$error_msg <- cv$message; return(r) }
    b <- cv$best; if (is.null(b)||nrow(b)==0) { r$error_msg <- "no model"; return(r) }
    Kh <- as.integer(b$K[1]); gh <- as.numeric(b$gamma[1])
    r$K_selected <- Kh; r$K_correct <- as.integer(Kh == K_true)
    r$gamma_selected <- gh; r$cindex_cv <- as.numeric(b$mean_score[1])
    
    prep <- preprocess_train_v4(sd$X_gmm, sd$X_cox, log_transform=FALSE)
    fit <- tryCatch(gemcox_full_multistart_shared_v4(
      X_gmm=prep$X_gmm, X_cox=prep$X_cox,
      time=sd$time, status=sd$status, K=Kh, alpha=alpha,
      max_iter=max_iter*2, surv_weight=gh, normalize_gmm_by_dim=FALSE,
      verbose=FALSE, n_starts=n_starts_final), error=function(e) e)
    if (inherits(fit,"error")) { r$error_msg <- fit$message; return(r) }
    r$converged <- 1L
    if (Kh==K_true && K_true>1) r$ARI <- adjustedRandIndex(sd$Z_true, fit$clusterid)
    
    eta <- tryCatch({
      tau <- fit$tau
      em <- sapply(1:Kh, function(k) eta_from_scaled(prep$X_cox,
                                                     fit$coxfit[[k]]$beta_s, fit$coxfit[[k]]$scaler))
      if (is.vector(em)) em <- matrix(em, ncol=Kh)
      rowSums(tau * em)
    }, error=function(e) NULL)
    if (!is.null(eta)) r$cindex_full <- c_index(sd$time, sd$status, eta)
    r
  }
  
  list(
    fixed = run_one(1.0, "gamma_fixed_1"),
    tuned = run_one(c(0.5, 1.0, 2.0, 5.0), "gamma_tuned")
  )
}


###############################################################
## Grid: alignment × separation × K_true
###############################################################
sim_grid <- expand.grid(
  K_true    = c(1, 2),
  alignment = c("aligned", "partial", "misaligned"),
  separation = c("med", "high"),
  event_rate = 0.44,
  stringsAsFactors = FALSE
)
sim_grid$sid <- seq_len(nrow(sim_grid))

NSIM    <- 30L
N_CORES <- min(24L, max(1L, detectCores(logical=FALSE)-1L))
cat(sprintf("Scenarios: %d | Reps: %d | Total pairs: %d | Cores: %d\n",
            nrow(sim_grid), NSIM, nrow(sim_grid)*NSIM, N_CORES))

###############################################################
## Execute
###############################################################
run_scenario <- function(sc) {
  cat(sprintf("[%d] K=%d align=%s sep=%s\n",
              sc$sid, sc$K_true, sc$alignment, sc$separation))
  out <- lapply(seq_len(NSIM), function(ri) {
    d <- tryCatch(simulate_alignment(
      n=117, K=sc$K_true, alignment=sc$alignment,
      separation=sc$separation, event_rate=sc$event_rate,
      seed=sc$sid*10000+ri), error=function(e) NULL)
    if (is.null(d)) return(NULL)
    
    pair <- tryCatch(run_pair(d), error=function(e)
      list(fixed=list(method="gamma_fixed_1",error_msg=e$message),
           tuned=list(method="gamma_tuned",error_msg=e$message)))
    
    mk <- function(r) {
      r <- lapply(r, function(v) if(is.null(v)||length(v)==0) NA else v[1])
      as.data.frame(c(list(sid=sc$sid, rep=ri, K_true=sc$K_true,
                           alignment=sc$alignment, separation=sc$separation), r),
                    stringsAsFactors=FALSE)
    }
    rbind(mk(pair$fixed), mk(pair$tuned))
  })
  dplyr::bind_rows(Filter(Negate(is.null), out))
}

message("========== ALIGNMENT ABLATION START ==========")
t0 <- proc.time()

if (N_CORES > 1 && .Platform$OS.type == "windows") {
  cl <- makeCluster(N_CORES, type="PSOCK"); on.exit(stopCluster(cl), add=TRUE)
  clusterExport(cl, ls(.GlobalEnv), .GlobalEnv)
  clusterEvalQ(cl, {library(survival);library(glmnet);library(MASS)
    library(mclust);library(dplyr);library(parallel)})
  all_res <- parLapply(cl, seq_len(nrow(sim_grid)),
                       function(i) tryCatch(run_scenario(sim_grid[i,]),
                                            error=function(e) data.frame(sid=sim_grid$sid[i],
                                                                         error_msg=e$message, stringsAsFactors=FALSE)))
} else if (N_CORES > 1) {
  all_res <- mclapply(seq_len(nrow(sim_grid)),
                      function(i) tryCatch(run_scenario(sim_grid[i,]),
                                           error=function(e) data.frame(sid=sim_grid$sid[i],
                                                                        error_msg=e$message, stringsAsFactors=FALSE)),
                      mc.cores=N_CORES)
} else {
  all_res <- lapply(seq_len(nrow(sim_grid)), function(i) run_scenario(sim_grid[i,]))
}

elapsed <- (proc.time()-t0)[3]
results_df <- dplyr::bind_rows(Filter(is.data.frame, all_res))
saveRDS(results_df, "GeMCox_alignment_ablation_raw.rds")


###############################################################
## Analysis
###############################################################
em <- results_df$error_msg
ok <- is.na(em) | !nzchar(trimws(em)) | trimws(em)=="NA"
df_ok <- results_df[ok,]
cat(sprintf("\nSuccess: %d/%d (%.0f%%)  Time: %.1f min\n",
            nrow(df_ok), nrow(results_df),
            100*nrow(df_ok)/max(nrow(results_df),1), elapsed/60))

## Pivot
comp <- df_ok %>%
  select(sid, rep, K_true, alignment, separation, method,
         K_correct, ARI, cindex_full, cindex_cox1) %>%
  mutate(gain = cindex_full - cindex_cox1) %>%
  pivot_wider(names_from=method,
              values_from=c(K_correct, ARI, cindex_full, gain),
              names_sep=".")

## ========== KEY TABLE ==========
cat("\n==========================================\n")
cat("GAMMA TUNING BENEFIT BY ALIGNMENT (K_true=2)\n")
cat("==========================================\n\n")

by_align <- comp %>%
  filter(K_true == 2) %>%
  group_by(alignment, separation) %>%
  summarise(
    n = n(),
    ## K-selection
    Ksel_fixed = mean(K_correct.gamma_fixed_1, na.rm=TRUE),
    Ksel_tuned = mean(K_correct.gamma_tuned, na.rm=TRUE),
    Ksel_diff  = mean(K_correct.gamma_tuned - K_correct.gamma_fixed_1, na.rm=TRUE),
    ## ARI
    ARI_fixed  = mean(ARI.gamma_fixed_1, na.rm=TRUE),
    ARI_tuned  = mean(ARI.gamma_tuned, na.rm=TRUE),
    ARI_diff   = mean(ARI.gamma_tuned - ARI.gamma_fixed_1, na.rm=TRUE),
    ## C-index gain
    gain_fixed = mean(gain.gamma_fixed_1, na.rm=TRUE),
    gain_tuned = mean(gain.gamma_tuned, na.rm=TRUE),
    gain_diff  = mean(gain.gamma_tuned - gain.gamma_fixed_1, na.rm=TRUE),
    .groups = "drop"
  ) %>%
  mutate(across(where(is.numeric) & !matches("^n$"), ~round(.x, 3)))

print(as.data.frame(by_align), row.names=FALSE)

## Aggregate by alignment
cat("\n--- Aggregated by alignment (both separations) ---\n\n")
by_align_agg <- comp %>%
  filter(K_true == 2) %>%
  group_by(alignment) %>%
  summarise(
    Ksel_fixed = round(mean(K_correct.gamma_fixed_1, na.rm=TRUE), 3),
    Ksel_tuned = round(mean(K_correct.gamma_tuned, na.rm=TRUE), 3),
    Ksel_diff  = round(mean(K_correct.gamma_tuned - K_correct.gamma_fixed_1, na.rm=TRUE), 3),
    gain_fixed = round(mean(gain.gamma_fixed_1, na.rm=TRUE), 3),
    gain_tuned = round(mean(gain.gamma_tuned, na.rm=TRUE), 3),
    gain_diff  = round(mean(gain.gamma_tuned - gain.gamma_fixed_1, na.rm=TRUE), 3),
    ## Paired Wilcoxon on C-index gain
    p_gain = round(tryCatch(
      wilcox.test(gain.gamma_tuned, gain.gamma_fixed_1, paired=TRUE)$p.value,
      error=function(e) NA), 4),
    .groups="drop"
  )
print(as.data.frame(by_align_agg), row.names=FALSE)

## Gamma values selected
cat("\n--- Selected gamma (tuned, K_true=2) ---\n\n")
gam_sel <- df_ok %>%
  filter(K_true==2, method=="gamma_tuned", !is.na(gamma_selected)) %>%
  group_by(alignment, separation) %>%
  summarise(
    gamma_mean = round(mean(gamma_selected), 2),
    gamma_med  = median(gamma_selected),
    pct_gt1    = round(mean(gamma_selected > 1), 2),
    .groups="drop"
  )
print(as.data.frame(gam_sel), row.names=FALSE)

## K=1 false positive
cat("\n--- K_true=1 false positive control ---\n\n")
fp <- comp %>% filter(K_true==1) %>%
  group_by(alignment) %>%
  summarise(
    Ksel_fixed = round(mean(K_correct.gamma_fixed_1, na.rm=TRUE), 3),
    Ksel_tuned = round(mean(K_correct.gamma_tuned, na.rm=TRUE), 3),
    .groups="drop"
  )
print(as.data.frame(fp), row.names=FALSE)

cat("\n--- Expected pattern ---\n")
cat("  ALIGNED:    gain_diff ≈ 0   (gamma doesn't matter)\n")
cat("  PARTIAL:    gain_diff > 0   (gamma helps moderately)\n")
cat("  MISALIGNED: gain_diff >> 0  (gamma essential)\n")

write.csv(by_align, "GeMCox_alignment_ablation_summary.csv", row.names=FALSE)
cat(sprintf("\nTotal time: %.1f min\n", elapsed/60))
message("========== ALIGNMENT ABLATION COMPLETE ==========")