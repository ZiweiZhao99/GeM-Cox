###############################################################
## GeM-Cox Simulation v4
##
## KEY FIX: normalize_gmm_by_dim = FALSE
##   The 1/q normalization with q=10 flattens GMM posteriors
##   to near-uniform, preventing cluster recovery even on clean
##   well-separated data. Turning it off restores full GMM
##   discriminative power. Gamma still balances GMM vs survival.
##
## Uses GeMCox_improved_v4.R (shared-baseline EM, fold-safe
## preprocessing, correct 1-SE rule)
##
## DGP: signal on OBSERVABLE features, Gaussian
###############################################################

required_pkgs <- c("survival","glmnet","mclust","dplyr","tidyr",
                   "ggplot2","patchwork","tibble","MASS","parallel")
for (pkg in required_pkgs)
  if (!requireNamespace(pkg, quietly=TRUE))
    install.packages(pkg, repos="https://cloud.r-project.org")

library(survival); library(glmnet); library(mclust)
library(dplyr); library(tidyr); library(ggplot2)
library(patchwork); library(tibble); library(MASS); library(parallel)
select <- dplyr::select; filter <- dplyr::filter
mutate <- dplyr::mutate; arrange <- dplyr::arrange

source("GeMCox_improved_v4.R")
set.seed(2025)

## ---- DGP: Gaussian features, signal on observable space ----
simulate_gemcox_v4 <- function(
    n=117, K=2, p_gmm=10, p_cox=10,
    p_gmm_signal=5, p_beta_signal=4,
    pi_c=NULL, separation=c("med","low","high"),
    feat_corr=0.2, event_rate=0.44,
    t_max=253, weibull_shape=1.5,
    baseline_var=0.2, seed=1) {

  separation <- match.arg(separation)
  set.seed(seed)
  delta_mu <- switch(separation, low=0.6, med=1.2, high=2.0)
  delta_beta <- switch(separation, low=0.3, med=0.6, high=1.0)
  if (is.null(pi_c)) pi_c <- rep(1/K, K)
  p_gmm_signal <- min(p_gmm_signal, p_gmm)
  p_beta_signal <- min(p_beta_signal, p_cox)

  positions <- seq(-(K-1)/2, (K-1)/2, length.out=K) * delta_mu
  mu_list <- lapply(1:K, function(k) {
    m <- rep(0, p_gmm)
    for (j in seq_len(p_gmm_signal)) m[j] <- positions[k]/sqrt(j)
    m
  })
  Sigma <- outer(1:p_gmm, 1:p_gmm, function(i,j) feat_corr^abs(i-j))

  beta_base <- rep(0, p_cox)
  beta_base[1:p_beta_signal] <- c(0.8,-0.5,0.4,-0.3)[1:p_beta_signal]
  beta_list <- lapply(1:K, function(k) {
    pos <- (k-(K+1)/2); b <- beta_base
    for (j in seq_len(min(3,p_beta_signal)))
      b[j] <- b[j] + delta_beta*pos/sqrt(j)
    b
  })

  scale_pos <- seq(-(K-1)/2,(K-1)/2,length.out=K)
  weibull_scales <- 200*exp(baseline_var*scale_pos)

  Z <- sample(1:K, n, replace=TRUE, prob=pi_c)
  X <- matrix(0, n, p_gmm)
  for (k in 1:K) {
    idx <- which(Z==k); if(!length(idx)) next
    X[idx,] <- mvrnorm(length(idx), mu_list[[k]], Sigma)
  }
  colnames(X) <- paste0("feat_",1:p_gmm)
  X_cox <- X[,1:p_cox,drop=FALSE]

  T_true <- numeric(n)
  for (k in 1:K) {
    idx <- which(Z==k); if(!length(idx)) next
    eta <- pmin(pmax(as.numeric(X_cox[idx,,drop=FALSE] %*% beta_list[[k]]),-5),5)
    T_true[idx] <- weibull_scales[k]*(-log(runif(length(idx))))^(1/weibull_shape)*exp(-eta/weibull_shape)
  }

  cal <- function(tgt,Tl,tm){lo<-1e-6;hi<-10
    for(it in 1:60){mid<-(lo+hi)/2;er<-mean(Tl<=pmin(rexp(length(Tl),mid),tm))
      if(abs(er-tgt)<0.005)break;if(er<tgt)hi<-mid else lo<-mid};mid}
  cr <- cal(event_rate,T_true,t_max)
  C <- pmin(rexp(n,cr),t_max)
  time <- pmax(pmin(T_true,C)+runif(n,0,0.01),0.1)
  status <- as.integer(T_true<=C)

  list(X_gmm=X, X_cox=X_cox, time=time, status=status, Z_true=Z,
       params=list(K=K,pi_c=pi_c,mu_list=mu_list,beta_list=beta_list,
                   Sigma=Sigma,separation=separation,
                   delta_mu=delta_mu,delta_beta=delta_beta,
                   weibull_scales=weibull_scales,
                   actual_event_rate=mean(status)))
}

## ---- Single replicate using v4 shared-baseline EM ----
## KEY: normalize_gmm_by_dim = FALSE everywhere
run_replicate_v4 <- function(sim_data, K_grid=1:3,
    gamma_grid=c(0, 0.5, 1.0, 2.0),
    nfolds=3, max_iter=50,
    n_starts=5, n_starts_final=10,
    alpha=0.5, baseline_mode="shared",
    verbose=FALSE) {

  X_gmm<-sim_data$X_gmm; X_cox<-sim_data$X_cox
  time<-sim_data$time; status<-sim_data$status
  Z_true<-sim_data$Z_true; K_true<-sim_data$params$K

  result <- list(K_true=K_true, K_selected=NA_integer_,
    K_correct=NA_integer_, gamma_selected=NA_real_,
    ARI=NA_real_, cindex_full=NA_real_, cindex_cv=NA_real_,
    cindex_cox1=NA_real_, beta_rmse=NA_real_,
    n_events=sum(status), event_rate_actual=mean(status),
    converged=NA_integer_, error_msg=NA_character_)

  ## Single Cox baseline
  result$cindex_cox1 <- tryCatch({
    sc <- make_x_scaler(X_cox); Xs <- apply_x_scaler(X_cox,sc)
    lam <- tryCatch(cv.glmnet(Xs,Surv(time,status),family="cox",
      alpha=alpha,standardize=FALSE,nfolds=3)$lambda.1se,error=function(e)0.05)
    gfit <- glmnet(Xs,Surv(time,status),family="cox",
      alpha=alpha,lambda=lam,standardize=FALSE)
    c_index(time,status,as.numeric(Xs%*%coef(gfit,s=lam)))
  }, error=function(e) NA_real_)

  ## CV — KEY FIX: normalize_gmm_by_dim=FALSE passed through
  cv_res <- tryCatch(cv_select_K_gamma_v4(
    X_gmm=X_gmm, X_cox=X_cox, time=time, status=status,
    K_grid=K_grid, gamma_grid=gamma_grid,
    nfolds=nfolds, alpha=alpha, max_iter=max_iter, n_starts=n_starts,
    log_transform=FALSE,
    normalize_gmm_by_dim=FALSE,   ## <-- THE FIX
    verbose=verbose, use_1se_rule=TRUE,
    baseline_mode=baseline_mode
  ), error=function(e) e)

  if (inherits(cv_res,"error")){result$error_msg<-cv_res$message;return(result)}
  best <- cv_res$best
  if (is.null(best)||nrow(best)==0){result$error_msg<-"CV no model";return(result)}

  K_hat<-as.integer(best$K[1]); gamma_hat<-as.numeric(best$gamma[1])
  result$K_selected<-K_hat; result$K_correct<-as.integer(K_hat==K_true)
  result$gamma_selected<-gamma_hat
  result$cindex_cv<-as.numeric(best$mean_score[1])

  ## Final fit
  prep <- preprocess_train_v4(X_gmm,X_cox,log_transform=FALSE)
  fitter <- if(baseline_mode=="shared") gemcox_full_multistart_shared_v4 else gemcox_full_multistart
  fit <- tryCatch(fitter(X_gmm=prep$X_gmm,X_cox=prep$X_cox,
    time=time,status=status,K=K_hat,alpha=alpha,max_iter=max_iter*2,
    surv_weight=gamma_hat,
    normalize_gmm_by_dim=FALSE,   ## <-- THE FIX
    verbose=FALSE,n_starts=n_starts_final),
    error=function(e) e)
  if (inherits(fit,"error")){result$error_msg<-fit$message;return(result)}
  result$converged <- 1L

  if (K_hat==K_true && K_true>1)
    result$ARI <- adjustedRandIndex(Z_true,fit$clusterid)

  eta_full <- tryCatch({
    tau<-fit$tau
    em<-sapply(1:K_hat,function(k) eta_from_scaled(prep$X_cox,
      fit$coxfit[[k]]$beta_s,fit$coxfit[[k]]$scaler))
    if(is.vector(em)) em<-matrix(em,ncol=K_hat)
    rowSums(tau*em)},error=function(e)NULL)
  if (!is.null(eta_full)) result$cindex_full<-c_index(time,status,eta_full)

  if (K_hat==K_true && K_true>1) {
    tb<-sim_data$params$beta_list; eb<-lapply(1:K_hat,function(k)fit$beta[,k])
    pc<-min(length(tb[[1]]),nrow(fit$beta))
    cost<-matrix(0,K_true,K_hat)
    for(i in 1:K_true)for(j in 1:K_hat)
      cost[i,j]<-mean((tb[[i]][1:pc]-eb[[j]][1:pc])^2)
    used<-logical(K_hat);rv<-numeric(K_true)
    for(i in 1:K_true){av<-which(!used);bj<-av[which.min(cost[i,av])]
      rv[i]<-sqrt(cost[i,bj]);used[bj]<-TRUE}
    result$beta_rmse<-mean(rv)
  }
  result
}

## ---- Grid ----
RUN_DIAGNOSTIC <- TRUE
sim_grid <- expand.grid(K_true=c(1,2,3), separation=c("low","med","high"),
  event_rate=c(0.44,0.25), n_subjects=117, feat_corr=0.2,
  stringsAsFactors=FALSE)
sim_grid$scenario_id <- seq_len(nrow(sim_grid))

if (RUN_DIAGNOSTIC) {
  SIM_K_GRID<-1:3; SIM_GAMMA_GRID<-c(0,0.5,1.0,2.0)
  SIM_NFOLDS<-3; SIM_MAX_ITER<-40
  SIM_N_STARTS<-5; SIM_N_STARTS_FINAL<-10
  NSIM<-25L
} else {
  SIM_K_GRID<-1:3; SIM_GAMMA_GRID<-c(0,0.25,0.5,1.0,2.0)
  SIM_NFOLDS<-5; SIM_MAX_ITER<-60
  SIM_N_STARTS<-5; SIM_N_STARTS_FINAL<-15
  NSIM<-100L
}
N_CORES <- min(50L, max(1L, detectCores(logical=FALSE)-1L))
cat(sprintf("Scenarios:%d Reps:%d Total:%d Cores:%d\n",
            nrow(sim_grid),NSIM,nrow(sim_grid)*NSIM,N_CORES))

run_scenario <- function(sc, nsim=NSIM) {
  cat(sprintf("[%d] K=%d sep=%s er=%.2f\n",
              sc$scenario_id,sc$K_true,sc$separation,sc$event_rate))
  results <- lapply(seq_len(nsim), function(ri) {
    sd <- tryCatch(simulate_gemcox_v4(n=sc$n_subjects,K=sc$K_true,
      p_gmm=10,p_cox=10,p_gmm_signal=5,p_beta_signal=4,
      separation=sc$separation,feat_corr=sc$feat_corr,
      event_rate=sc$event_rate,seed=sc$scenario_id*10000+ri),
      error=function(e)NULL)
    if(is.null(sd))return(NULL)
    res <- tryCatch(run_replicate_v4(sd,K_grid=SIM_K_GRID,
      gamma_grid=SIM_GAMMA_GRID,nfolds=SIM_NFOLDS,max_iter=SIM_MAX_ITER,
      n_starts=SIM_N_STARTS,n_starts_final=SIM_N_STARTS_FINAL,
      alpha=0.5,baseline_mode="shared",verbose=FALSE),
      error=function(e)list(error_msg=e$message))
    res<-lapply(res,function(v)if(is.null(v)||length(v)==0)NA else v[1])
    as.data.frame(c(list(scenario_id=sc$scenario_id,rep=ri,
      K_true=sc$K_true,separation=sc$separation,
      event_rate_target=sc$event_rate,n=sc$n_subjects),
      res),stringsAsFactors=FALSE)
  })
  bind_rows(Filter(Negate(is.null),results))
}

run_par <- function(nj,FUN,nc){
  if(nc<=1L)return(lapply(seq_len(nj),FUN))
  if(.Platform$OS.type=="windows"){
    cl<-makeCluster(nc,type="PSOCK");on.exit(stopCluster(cl))
    clusterExport(cl,ls(.GlobalEnv),.GlobalEnv)
    clusterEvalQ(cl,{library(survival);library(glmnet);library(MASS)
                      library(mclust);library(dplyr);library(parallel)})
    parLapply(cl,seq_len(nj),FUN)
  } else mclapply(seq_len(nj),FUN,mc.cores=nc)
}

message("========== STARTING v4 SIMULATION (normalize_gmm_by_dim=FALSE) ==========")
t0<-proc.time()
all_res <- run_par(nrow(sim_grid),function(i)tryCatch(run_scenario(sim_grid[i,]),
  error=function(e)data.frame(scenario_id=sim_grid$scenario_id[i],
    error_msg=e$message,stringsAsFactors=FALSE)),N_CORES)
elapsed<-(proc.time()-t0)[3]
results_df <- bind_rows(Filter(is.data.frame,all_res))
saveRDS(results_df,"GeMCox_simulation_v4_raw.rds")

em<-results_df$error_msg
ok<-is.na(em)|(nzchar(trimws(em))==FALSE)|(trimws(em)=="NA")
df_ok<-results_df[ok,]
cat(sprintf("\nSuccess:%d/%d(%.0f%%) in %.0fs\n",nrow(df_ok),nrow(results_df),
    100*nrow(df_ok)/max(nrow(results_df),1),elapsed))

summary_df <- df_ok %>%
  group_by(K_true,separation,event_rate_target) %>%
  summarise(n_reps=n(),
    K_sel_accuracy=mean(K_correct,na.rm=TRUE),
    K_sel_mean=mean(K_selected,na.rm=TRUE),
    ARI_mean=mean(ARI,na.rm=TRUE), ARI_sd=sd(ARI,na.rm=TRUE),
    cindex_full_mean=mean(cindex_full,na.rm=TRUE),
    cindex_cv_mean=mean(cindex_cv,na.rm=TRUE),
    cindex_cox1_mean=mean(cindex_cox1,na.rm=TRUE),
    cindex_gain_mean=mean(cindex_full-cindex_cox1,na.rm=TRUE),
    gamma_mean=mean(gamma_selected,na.rm=TRUE),
    converge_rate=mean(converged,na.rm=TRUE),.groups="drop") %>%
  mutate(across(where(is.numeric),~round(.x,4)))

write.csv(summary_df,"GeMCox_simulation_v4_summary.csv",row.names=FALSE)
cat("\n========== v4 RESULTS ==========\n")
print(as.data.frame(summary_df))
message(sprintf("========== v4 COMPLETE (%.0fs) ==========",elapsed))
