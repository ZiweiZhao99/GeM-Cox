###############################################################
## GeM-Cox Simulation v4 — Production Run
##
## ROOT CAUSE FIX: normalize_gmm_by_dim = FALSE
## Target: < 10 min on 24-core Windows PC
##
## Timing budget (24 cores, PSOCK):
##   18 scenarios → 1 batch of 18 workers (fits in 24 cores)
##   20 reps/scenario, each ~15s → 20*15/1 = 300s = 5 min
##   Overhead (cluster setup, summary) ~1 min
##   Total: ~6 min
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
## DGP: Gaussian, signal on observable space
###############################################################
simulate_v4 <- function(n=117, K=2, p=10, p_sig=5, p_beta=4,
                         separation=c("med","low","high"),
                         rho=0.2, event_rate=0.44, seed=1) {
  separation <- match.arg(separation)
  set.seed(seed)
  dmu  <- switch(separation, low=0.6, med=1.2, high=2.0)
  dbeta <- switch(separation, low=0.3, med=0.6, high=1.0)
  pi_c <- rep(1/K, K)
  pos <- seq(-(K-1)/2,(K-1)/2,length.out=K)

  mu_list <- lapply(1:K, function(k){m<-rep(0,p)
    for(j in 1:min(p_sig,p)) m[j]<-pos[k]*dmu/sqrt(j); m})
  Sigma <- outer(1:p,1:p,function(i,j) rho^abs(i-j))

  bb <- rep(0,p); bb[1:min(p_beta,p)] <- c(0.8,-0.5,0.4,-0.3)[1:min(p_beta,p)]
  beta_list <- lapply(1:K, function(k){b<-bb
    for(j in 1:min(3,p_beta)) b[j]<-b[j]+dbeta*pos[k]/sqrt(j); b})

  wscales <- 200*exp(0.2*pos)
  Z <- sample(1:K,n,replace=TRUE,prob=pi_c)
  X <- matrix(0,n,p)
  for(k in 1:K){i<-which(Z==k);if(!length(i))next
    X[i,]<-mvrnorm(length(i),mu_list[[k]],Sigma)}
  colnames(X)<-paste0("f",1:p)

  Tt<-numeric(n)
  for(k in 1:K){i<-which(Z==k);if(!length(i))next
    eta<-pmin(pmax(as.numeric(X[i,,drop=F]%*%beta_list[[k]]),-5),5)
    Tt[i]<-wscales[k]*(-log(runif(length(i))))^(1/1.5)*exp(-eta/1.5)}

  # Calibrate censoring
  lo<-1e-6;hi<-10
  for(it in 1:50){mid<-(lo+hi)/2
    er<-mean(Tt<=pmin(rexp(n,mid),253))
    if(abs(er-event_rate)<0.005)break
    if(er<event_rate)hi<-mid else lo<-mid}
  C<-pmin(rexp(n,mid),253)
  time<-pmax(pmin(Tt,C)+runif(n,0,0.01),0.1)
  status<-as.integer(Tt<=C)

  list(X_gmm=X,X_cox=X,time=time,status=status,Z_true=Z,
       params=list(K=K,beta_list=beta_list,separation=separation,
                   dmu=dmu,dbeta=dbeta,er_actual=mean(status)))
}

###############################################################
## Single replicate — KEY: normalize_gmm_by_dim=FALSE
###############################################################
run_rep <- function(sd, K_grid=1:3, gamma_grid=c(0,1.0,2.0),
                    nfolds=3, max_iter=30, n_starts=3,
                    n_starts_final=5, alpha=0.5) {
  r <- list(K_true=sd$params$K, K_selected=NA_integer_,
    K_correct=NA_integer_, gamma_selected=NA_real_,
    ARI=NA_real_, cindex_full=NA_real_, cindex_cv=NA_real_,
    cindex_cox1=NA_real_, beta_rmse=NA_real_,
    event_rate_actual=mean(sd$status),
    converged=NA_integer_, error_msg=NA_character_)

  # Single Cox baseline
  r$cindex_cox1 <- tryCatch({
    sc<-make_x_scaler(sd$X_cox); Xs<-apply_x_scaler(sd$X_cox,sc)
    lam<-tryCatch(cv.glmnet(Xs,Surv(sd$time,sd$status),family="cox",
      alpha=alpha,standardize=FALSE,nfolds=3)$lambda.1se,error=function(e)0.05)
    gf<-glmnet(Xs,Surv(sd$time,sd$status),family="cox",
      alpha=alpha,lambda=lam,standardize=FALSE)
    c_index(sd$time,sd$status,as.numeric(Xs%*%coef(gf,s=lam)))
  },error=function(e)NA_real_)

  # CV with THE FIX
  cv <- tryCatch(cv_select_K_gamma_v4(
    X_gmm=sd$X_gmm, X_cox=sd$X_cox,
    time=sd$time, status=sd$status,
    K_grid=K_grid, gamma_grid=gamma_grid,
    nfolds=nfolds, alpha=alpha,
    max_iter=max_iter, n_starts=n_starts,
    log_transform=FALSE,
    normalize_gmm_by_dim=FALSE,      ## THE FIX
    verbose=FALSE, use_1se_rule=TRUE,
    baseline_mode="shared"
  ),error=function(e)e)
  if(inherits(cv,"error")){r$error_msg<-cv$message;return(r)}
  b<-cv$best; if(is.null(b)||nrow(b)==0){r$error_msg<-"no model";return(r)}

  Kh<-as.integer(b$K[1]); gh<-as.numeric(b$gamma[1])
  r$K_selected<-Kh; r$K_correct<-as.integer(Kh==sd$params$K)
  r$gamma_selected<-gh; r$cindex_cv<-as.numeric(b$mean_score[1])

  # Final fit
  prep<-preprocess_train_v4(sd$X_gmm,sd$X_cox,log_transform=FALSE)
  fit<-tryCatch(gemcox_full_multistart_shared_v4(
    X_gmm=prep$X_gmm, X_cox=prep$X_cox,
    time=sd$time, status=sd$status,
    K=Kh, alpha=alpha, max_iter=max_iter*2,
    surv_weight=gh,
    normalize_gmm_by_dim=FALSE,      ## THE FIX
    verbose=FALSE, n_starts=n_starts_final
  ),error=function(e)e)
  if(inherits(fit,"error")){r$error_msg<-fit$message;return(r)}
  r$converged<-1L

  if(Kh==sd$params$K && sd$params$K>1)
    r$ARI<-adjustedRandIndex(sd$Z_true,fit$clusterid)

  eta<-tryCatch({tau<-fit$tau
    em<-sapply(1:Kh,function(k)eta_from_scaled(prep$X_cox,
      fit$coxfit[[k]]$beta_s,fit$coxfit[[k]]$scaler))
    if(is.vector(em))em<-matrix(em,ncol=Kh)
    rowSums(tau*em)},error=function(e)NULL)
  if(!is.null(eta)) r$cindex_full<-c_index(sd$time,sd$status,eta)

  if(Kh==sd$params$K && sd$params$K>1){
    tb<-sd$params$beta_list; eb<-lapply(1:Kh,function(k)fit$beta[,k])
    pc<-min(length(tb[[1]]),nrow(fit$beta))
    cost<-matrix(0,sd$params$K,Kh)
    for(i in 1:sd$params$K)for(j in 1:Kh)
      cost[i,j]<-mean((tb[[i]][1:pc]-eb[[j]][1:pc])^2)
    used<-logical(Kh);rv<-numeric(sd$params$K)
    for(i in 1:sd$params$K){av<-which(!used);bj<-av[which.min(cost[i,av])]
      rv[i]<-sqrt(cost[i,bj]);used[bj]<-TRUE}
    r$beta_rmse<-mean(rv)
  }
  r
}

###############################################################
## Grid + execution
###############################################################
sim_grid <- expand.grid(
  K_true=c(1,2,3), separation=c("low","med","high"),
  event_rate=c(0.44,0.25), stringsAsFactors=FALSE)
sim_grid$sid <- seq_len(nrow(sim_grid))

## Tuned for < 10 min on 24-core Windows
NSIM          <- 20L
K_GRID        <- 1:3
GAMMA_GRID    <- c(0, 1.0, 2.0)   # 3 values (was 4)
CV_FOLDS      <- 3
CV_MAX_ITER   <- 30                # 30 iter (was 40)
CV_N_STARTS   <- 3                 # 3 starts (was 5)
FINAL_STARTS  <- 5                 # 5 starts (was 10)
ALPHA         <- 0.5

N_CORES <- min(24L, max(1L, detectCores(logical=FALSE)-1L))

cat(sprintf("Grid: %d scenarios x %d reps = %d fits\n",
            nrow(sim_grid), NSIM, nrow(sim_grid)*NSIM))
cat(sprintf("CV grid: %d K x %d gamma x %d folds x %d starts = %d EM fits/rep\n",
            length(K_GRID), length(GAMMA_GRID), CV_FOLDS, CV_N_STARTS,
            length(K_GRID)*length(GAMMA_GRID)*CV_FOLDS*CV_N_STARTS))
cat(sprintf("Cores: %d\n", N_CORES))

run_scenario <- function(sc) {
  cat(sprintf("[%d] K=%d sep=%s er=%.2f\n",sc$sid,sc$K_true,sc$separation,sc$event_rate))
  out <- lapply(seq_len(NSIM), function(ri) {
    d <- tryCatch(simulate_v4(n=117, K=sc$K_true, separation=sc$separation,
      event_rate=sc$event_rate, seed=sc$sid*10000+ri), error=function(e)NULL)
    if(is.null(d)) return(NULL)
    res <- tryCatch(run_rep(d, K_grid=K_GRID, gamma_grid=GAMMA_GRID,
      nfolds=CV_FOLDS, max_iter=CV_MAX_ITER, n_starts=CV_N_STARTS,
      n_starts_final=FINAL_STARTS, alpha=ALPHA),
      error=function(e)list(error_msg=e$message))
    res <- lapply(res,function(v)if(is.null(v)||length(v)==0)NA else v[1])
    as.data.frame(c(list(sid=sc$sid,rep=ri,K_true=sc$K_true,
      separation=sc$separation,event_rate=sc$event_rate),res),
      stringsAsFactors=FALSE)
  })
  dplyr::bind_rows(Filter(Negate(is.null),out))
}

## Windows PSOCK parallel
message("========== v4 SIMULATION START ==========")
t0 <- proc.time()

if (N_CORES > 1 && .Platform$OS.type == "windows") {
  cl <- makeCluster(N_CORES, type="PSOCK")
  on.exit(stopCluster(cl), add=TRUE)
  clusterExport(cl, ls(.GlobalEnv), .GlobalEnv)
  clusterEvalQ(cl, {
    library(survival); library(glmnet); library(MASS)
    library(mclust); library(dplyr); library(parallel)
  })
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
  all_res <- lapply(seq_len(nrow(sim_grid)),
    function(i) run_scenario(sim_grid[i,]))
}

elapsed <- (proc.time()-t0)[3]
results_df <- dplyr::bind_rows(Filter(is.data.frame, all_res))
saveRDS(results_df, "GeMCox_simulation_v4_raw.rds")

## Summarize
em <- results_df$error_msg
ok <- is.na(em) | !nzchar(trimws(em)) | trimws(em)=="NA"
df_ok <- results_df[ok,]
cat(sprintf("\nSuccess: %d/%d (%.0f%%)  Time: %.0f sec (%.1f min)\n",
    nrow(df_ok), nrow(results_df),
    100*nrow(df_ok)/max(nrow(results_df),1), elapsed, elapsed/60))

summary_df <- df_ok %>%
  group_by(K_true, separation, event_rate) %>%
  summarise(
    n_reps           = n(),
    K_sel_accuracy   = round(mean(K_correct, na.rm=TRUE), 3),
    K_sel_mean       = round(mean(K_selected, na.rm=TRUE), 2),
    ARI_mean         = round(mean(ARI, na.rm=TRUE), 3),
    ARI_sd           = round(sd(ARI, na.rm=TRUE), 3),
    cindex_full_mean = round(mean(cindex_full, na.rm=TRUE), 3),
    cindex_cv_mean   = round(mean(cindex_cv, na.rm=TRUE), 3),
    cindex_cox1_mean = round(mean(cindex_cox1, na.rm=TRUE), 3),
    cindex_gain      = round(mean(cindex_full - cindex_cox1, na.rm=TRUE), 3),
    gamma_mean       = round(mean(gamma_selected, na.rm=TRUE), 2),
    converge_rate    = round(mean(converged, na.rm=TRUE), 2),
    .groups="drop")

write.csv(summary_df, "GeMCox_simulation_v4_summary.csv", row.names=FALSE)

cat("\n========== v4 RESULTS ==========\n")
print(as.data.frame(summary_df), right=FALSE)
cat(sprintf("\nTotal wall time: %.1f min\n", elapsed/60))
message("========== v4 COMPLETE ==========")
