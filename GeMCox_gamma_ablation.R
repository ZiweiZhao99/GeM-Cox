###############################################################
## GeM-Cox: Gamma Ablation Study
##
## Paired comparison on IDENTICAL data:
##   A) gamma = 1 fixed  (standard EM, no tuning)
##   B) gamma CV-selected from {0, 0.25, 0.5, 1.0, 2.0}
##
## If B outperforms A, the survival weight is justified.
## If A ≈ B, gamma tuning is unnecessary overhead.
##
## Target: < 10 min on 24-core Windows
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
## DGP (identical to simulation v4)
###############################################################
simulate_v4 <- function(n=117, K=2, p=10, p_sig=5, p_beta=4,
                         separation=c("med","low","high"),
                         rho=0.2, event_rate=0.44, seed=1) {
  separation <- match.arg(separation)
  set.seed(seed)
  dmu <- switch(separation, low=0.6, med=1.2, high=2.0)
  dbeta <- switch(separation, low=0.3, med=0.6, high=1.0)
  pi_c <- rep(1/K, K)
  pos <- seq(-(K-1)/2, (K-1)/2, length.out=K)

  mu_list <- lapply(1:K, function(k){m<-rep(0,p)
    for(j in 1:min(p_sig,p)) m[j]<-pos[k]*dmu/sqrt(j); m})
  Sigma <- outer(1:p,1:p,function(i,j) rho^abs(i-j))
  bb <- rep(0,p); bb[1:min(p_beta,p)] <- c(0.8,-0.5,0.4,-0.3)[1:min(p_beta,p)]
  beta_list <- lapply(1:K, function(k){b<-bb
    for(j in 1:min(3,p_beta)) b[j]<-b[j]+dbeta*pos[k]/sqrt(j); b})
  wscales <- 200*exp(0.2*pos)

  Z <- sample(1:K, n, replace=TRUE, prob=pi_c)
  X <- matrix(0, n, p)
  for(k in 1:K){i<-which(Z==k); if(!length(i))next
    X[i,] <- mvrnorm(length(i), mu_list[[k]], Sigma)}
  colnames(X) <- paste0("f",1:p)

  Tt <- numeric(n)
  for(k in 1:K){i<-which(Z==k); if(!length(i))next
    eta<-pmin(pmax(as.numeric(X[i,,drop=F]%*%beta_list[[k]]),-5),5)
    Tt[i]<-wscales[k]*(-log(runif(length(i))))^(1/1.5)*exp(-eta/1.5)}

  lo<-1e-6;hi<-10
  for(it in 1:50){mid<-(lo+hi)/2
    er<-mean(Tt<=pmin(rexp(n,mid),253))
    if(abs(er-event_rate)<0.005)break
    if(er<event_rate)hi<-mid else lo<-mid}
  C<-pmin(rexp(n,mid),253)
  time<-pmax(pmin(Tt,C)+runif(n,0,0.01),0.1)
  status<-as.integer(Tt<=C)

  list(X_gmm=X, X_cox=X, time=time, status=status, Z_true=Z,
       params=list(K=K,beta_list=beta_list,separation=separation,
                   er_actual=mean(status)))
}

###############################################################
## Single replicate: run BOTH conditions on same data
###############################################################
run_rep_paired <- function(sd, K_grid=1:3,
    nfolds=3, max_iter=30, n_starts=3,
    n_starts_final=5, alpha=0.5) {

  K_true <- sd$params$K

  ## Shared baseline: single Cox
  cox1 <- tryCatch({
    sc<-make_x_scaler(sd$X_cox); Xs<-apply_x_scaler(sd$X_cox,sc)
    lam<-tryCatch(cv.glmnet(Xs,Surv(sd$time,sd$status),family="cox",
      alpha=alpha,standardize=FALSE,nfolds=3)$lambda.1se,error=function(e)0.05)
    gf<-glmnet(Xs,Surv(sd$time,sd$status),family="cox",
      alpha=alpha,lambda=lam,standardize=FALSE)
    c_index(sd$time,sd$status,as.numeric(Xs%*%coef(gf,s=lam)))
  },error=function(e) NA_real_)

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
    ),error=function(e) e)
    if(inherits(cv,"error")){r$error_msg<-cv$message;return(r)}

    b<-cv$best
    if(is.null(b)||nrow(b)==0){r$error_msg<-"no model";return(r)}
    Kh<-as.integer(b$K[1]); gh<-as.numeric(b$gamma[1])
    r$K_selected<-Kh; r$K_correct<-as.integer(Kh==K_true)
    r$gamma_selected<-gh; r$cindex_cv<-as.numeric(b$mean_score[1])

    prep<-preprocess_train_v4(sd$X_gmm,sd$X_cox,log_transform=FALSE)
    fit<-tryCatch(gemcox_full_multistart_shared_v4(
      X_gmm=prep$X_gmm,X_cox=prep$X_cox,
      time=sd$time,status=sd$status,K=Kh,alpha=alpha,
      max_iter=max_iter*2,surv_weight=gh,
      normalize_gmm_by_dim=FALSE,
      verbose=FALSE,n_starts=n_starts_final),error=function(e)e)
    if(inherits(fit,"error")){r$error_msg<-fit$message;return(r)}
    r$converged<-1L

    if(Kh==K_true && K_true>1) r$ARI<-adjustedRandIndex(sd$Z_true,fit$clusterid)

    eta<-tryCatch({tau<-fit$tau
      em<-sapply(1:Kh,function(k)eta_from_scaled(prep$X_cox,
        fit$coxfit[[k]]$beta_s,fit$coxfit[[k]]$scaler))
      if(is.vector(em))em<-matrix(em,ncol=Kh)
      rowSums(tau*em)},error=function(e)NULL)
    if(!is.null(eta)) r$cindex_full<-c_index(sd$time,sd$status,eta)
    r
  }

  ## --- Condition A: gamma = 1 fixed ---
  resA <- run_one(gamma_grid = 1.0, label = "gamma_fixed_1")

  ## --- Condition B: gamma CV-tuned ---
  resB <- run_one(gamma_grid = c(0, 0.25, 0.5, 1.0, 2.0), label = "gamma_cv_tuned")

  list(A=resA, B=resB)
}

###############################################################
## Simulation grid
###############################################################
sim_grid <- expand.grid(
  K_true=c(1,2,3), separation=c("low","med","high"),
  event_rate=c(0.44, 0.25), n=117L,  ## add 250L here for full run
  stringsAsFactors=FALSE)
sim_grid$sid <- seq_len(nrow(sim_grid))

NSIM <- 25L          # 30 reps: enough for paired comparison
N_CORES <- min(24L, max(1L, detectCores(logical=FALSE)-1L))

cat(sprintf("Scenarios: %d | Reps: %d | Total pairs: %d | Cores: %d\n",
            nrow(sim_grid), NSIM, nrow(sim_grid)*NSIM, N_CORES))

###############################################################
## Run
###############################################################
run_scenario <- function(sc) {
  cat(sprintf("[%d] K=%d sep=%s er=%.2f n=%d\n",
              sc$sid,sc$K_true,sc$separation,sc$event_rate,sc$n))
  out <- lapply(seq_len(NSIM), function(ri) {
    d <- tryCatch(simulate_v4(n=sc$n, K=sc$K_true, separation=sc$separation,
      event_rate=sc$event_rate, seed=sc$sid*10000+ri),error=function(e)NULL)
    if(is.null(d)) return(NULL)

    pair <- tryCatch(run_rep_paired(d),
      error=function(e) list(
        A=list(method="gamma_fixed_1",error_msg=e$message),
        B=list(method="gamma_cv_tuned",error_msg=e$message)))

    make_row <- function(r) {
      r <- lapply(r,function(v)if(is.null(v)||length(v)==0)NA else v[1])
      as.data.frame(c(list(sid=sc$sid,rep=ri,K_true=sc$K_true,
        separation=sc$separation,event_rate=sc$event_rate,n=sc$n),r),
        stringsAsFactors=FALSE)
    }
    rbind(make_row(pair$A), make_row(pair$B))
  })
  dplyr::bind_rows(Filter(Negate(is.null), out))
}

message("========== GAMMA ABLATION START ==========")
t0 <- proc.time()

if (N_CORES > 1 && .Platform$OS.type == "windows") {
  cl <- makeCluster(N_CORES, type="PSOCK")
  on.exit(stopCluster(cl), add=TRUE)
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
saveRDS(results_df, "GeMCox_gamma_ablation_raw.rds")

###############################################################
## Summary
###############################################################
em <- results_df$error_msg
ok <- is.na(em) | !nzchar(trimws(em)) | trimws(em)=="NA"
df_ok <- results_df[ok,]

cat(sprintf("\nSuccess: %d/%d (%.0f%%)  Time: %.1f min\n",
    nrow(df_ok), nrow(results_df),
    100*nrow(df_ok)/max(nrow(results_df),1), elapsed/60))

summary_df <- df_ok %>%
  group_by(method, K_true, separation, event_rate, n) %>%
  summarise(
    n_reps           = n(),
    K_sel_accuracy   = round(mean(K_correct, na.rm=TRUE), 3),
    K_sel_mean       = round(mean(K_selected, na.rm=TRUE), 2),
    ARI_mean         = round(mean(ARI, na.rm=TRUE), 3),
    cindex_full_mean = round(mean(cindex_full, na.rm=TRUE), 3),
    cindex_cv_mean   = round(mean(cindex_cv, na.rm=TRUE), 3),
    cindex_cox1_mean = round(mean(cindex_cox1, na.rm=TRUE), 3),
    cindex_gain      = round(mean(cindex_full - cindex_cox1, na.rm=TRUE), 3),
    gamma_mean       = round(mean(gamma_selected, na.rm=TRUE), 2),
    .groups="drop") %>%
  arrange(K_true, separation, event_rate, n, method)

write.csv(summary_df, "GeMCox_gamma_ablation_summary.csv", row.names=FALSE)

## ---- Head-to-head comparison ----
cat("\n========== HEAD-TO-HEAD: gamma=1 fixed vs gamma CV-tuned ==========\n\n")

comp <- df_ok %>%
  select(sid, rep, K_true, separation, event_rate, n, method,
         K_correct, ARI, cindex_full, cindex_cox1) %>%
  mutate(cindex_gain = cindex_full - cindex_cox1) %>%
  pivot_wider(names_from=method,
              values_from=c(K_correct, ARI, cindex_full, cindex_gain, cindex_cox1),
              names_sep=".")

## Paired comparison by scenario
paired <- comp %>%
  group_by(K_true, separation, event_rate, n) %>%
  summarise(
    n_pairs = n(),
    ## K-selection
    Ksel_fixed = mean(K_correct.gamma_fixed_1, na.rm=TRUE),
    Ksel_tuned = mean(K_correct.gamma_cv_tuned, na.rm=TRUE),
    Ksel_diff  = mean(K_correct.gamma_cv_tuned - K_correct.gamma_fixed_1, na.rm=TRUE),
    ## ARI (only where both selected correct K)
    ARI_fixed  = mean(ARI.gamma_fixed_1, na.rm=TRUE),
    ARI_tuned  = mean(ARI.gamma_cv_tuned, na.rm=TRUE),
    ## C-index gain
    gain_fixed = mean(cindex_gain.gamma_fixed_1, na.rm=TRUE),
    gain_tuned = mean(cindex_gain.gamma_cv_tuned, na.rm=TRUE),
    gain_diff  = mean(cindex_gain.gamma_cv_tuned - cindex_gain.gamma_fixed_1, na.rm=TRUE),
    .groups="drop"
  ) %>%
  mutate(across(where(is.numeric) & !matches("n_pairs"), ~round(.x, 3)))

cat("By scenario:\n")
print(as.data.frame(paired %>% filter(K_true > 1)), row.names=FALSE)

## Overall summary
cat("\n\nOverall (K_true >= 2):\n")
comp_k2 <- comp %>% filter(K_true >= 2)
cat(sprintf("  K-sel accuracy:  fixed=%.3f  tuned=%.3f  diff=%+.3f  (p=%.4f)\n",
    mean(comp_k2$K_correct.gamma_fixed_1, na.rm=TRUE),
    mean(comp_k2$K_correct.gamma_cv_tuned, na.rm=TRUE),
    mean(comp_k2$K_correct.gamma_cv_tuned - comp_k2$K_correct.gamma_fixed_1, na.rm=TRUE),
    tryCatch(wilcox.test(comp_k2$K_correct.gamma_cv_tuned,
                         comp_k2$K_correct.gamma_fixed_1,
                         paired=TRUE)$p.value, error=function(e) NA)))
cat(sprintf("  C-index gain:    fixed=%.3f  tuned=%.3f  diff=%+.3f  (p=%.4f)\n",
    mean(comp_k2$cindex_gain.gamma_fixed_1, na.rm=TRUE),
    mean(comp_k2$cindex_gain.gamma_cv_tuned, na.rm=TRUE),
    mean(comp_k2$cindex_gain.gamma_cv_tuned - comp_k2$cindex_gain.gamma_fixed_1, na.rm=TRUE),
    tryCatch(wilcox.test(comp_k2$cindex_gain.gamma_cv_tuned,
                         comp_k2$cindex_gain.gamma_fixed_1,
                         paired=TRUE)$p.value, error=function(e) NA)))

## Overall for K_true = 1
cat("\nOverall (K_true = 1, false positive control):\n")
comp_k1 <- comp %>% filter(K_true == 1)
cat(sprintf("  K-sel accuracy:  fixed=%.3f  tuned=%.3f\n",
    mean(comp_k1$K_correct.gamma_fixed_1, na.rm=TRUE),
    mean(comp_k1$K_correct.gamma_cv_tuned, na.rm=TRUE)))

cat(sprintf("\nTotal time: %.1f min\n", elapsed/60))
message("========== GAMMA ABLATION COMPLETE ==========")
