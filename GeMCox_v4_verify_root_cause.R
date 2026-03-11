## Verify: turning off 1/q normalization fixes the EM
library(MASS); library(survival); library(glmnet); library(mclust)
source("GeMCox_improved_v4.R")
set.seed(2025)

simulate_clean <- function(n=300, K=2, delta_mu=2.0, delta_beta=0.8, seed=42) {
  set.seed(seed); p <- 10; p_sig <- 5
  positions <- seq(-(K-1)/2,(K-1)/2,length.out=K)*delta_mu
  mu_list <- lapply(1:K, function(k) {
    m <- rep(0,p); for(j in 1:p_sig) m[j] <- positions[k]/sqrt(j); m })
  Sigma <- outer(1:p,1:p,function(i,j) 0.2^abs(i-j))
  beta_base <- c(0.8,-0.5,0.4,-0.3,rep(0,p-4))
  beta_list <- lapply(1:K, function(k) {
    pos <- (k-(K+1)/2); b <- beta_base
    for(j in 1:3) b[j] <- b[j]+delta_beta*pos/sqrt(j); b })
  Z <- sample(1:K,n,replace=TRUE)
  X <- matrix(0,n,p)
  for(k in 1:K){idx<-which(Z==k); X[idx,]<-mvrnorm(length(idx),mu_list[[k]],Sigma)}
  colnames(X) <- paste0("f",1:p)
  T_true <- numeric(n)
  for(k in 1:K){idx<-which(Z==k); eta<-pmin(pmax(as.numeric(X[idx,]%*%beta_list[[k]]),-5),5)
    T_true[idx]<-200*(-log(runif(length(idx))))^(1/1.5)*exp(-eta/1.5)}
  C <- pmin(rexp(n,0.005),253); time <- pmax(pmin(T_true,C)+runif(n,0,0.01),0.1)
  status <- as.integer(T_true<=C)
  list(X_gmm=X,X_cox=X,time=time,status=status,Z_true=Z,beta_list=beta_list)
}

dat <- simulate_clean()
cat(sprintf("Events: %d/%d = %.0f%%\n\n", sum(dat$status), length(dat$status),
            100*mean(dat$status)))

cat("=== WITH normalize_gmm_by_dim = TRUE (default, q_eff=10) ===\n")
cat("  This divides GMM log-density by 10, flattening posteriors\n\n")
for (gam in c(0, 0.5, 1.0, 2.0, 5.0)) {
  fit <- tryCatch(gemcox_full_multistart_shared_v4(
    X_gmm=dat$X_gmm, X_cox=dat$X_cox, time=dat$time, status=dat$status,
    K=2, alpha=0.5, max_iter=80, surv_weight=gam,
    normalize_gmm_by_dim=TRUE,
    verbose=FALSE, n_starts=10), error=function(e)NULL)
  if(!is.null(fit)) cat(sprintf("  gamma=%4.2f  ARI=%.4f  pi=[%.2f,%.2f]  sizes=%d/%d\n",
    gam, adjustedRandIndex(dat$Z_true,fit$clusterid), fit$pi[1], fit$pi[2],
    sum(fit$clusterid==1), sum(fit$clusterid==2)))
}

cat("\n=== WITH normalize_gmm_by_dim = FALSE (fix: q_eff=1) ===\n")
cat("  Full GMM density used — posteriors are sharp\n\n")
for (gam in c(0, 0.5, 1.0, 2.0, 5.0)) {
  fit <- tryCatch(gemcox_full_multistart_shared_v4(
    X_gmm=dat$X_gmm, X_cox=dat$X_cox, time=dat$time, status=dat$status,
    K=2, alpha=0.5, max_iter=80, surv_weight=gam,
    normalize_gmm_by_dim=FALSE,
    verbose=FALSE, n_starts=10), error=function(e)NULL)
  if(!is.null(fit)) {
    ari <- adjustedRandIndex(dat$Z_true,fit$clusterid)
    cat(sprintf("  gamma=%4.2f  ARI=%.4f  pi=[%.2f,%.2f]  beta1=[%.3f,%.3f]",
      gam, ari, fit$pi[1], fit$pi[2], fit$beta[1,1], fit$beta[1,2]))
    if (ari > 0.3) cat("  ***RECOVERED***")
    cat("\n")
  }
}

cat("\n=== ALSO: fewer GMM features (q=5 signal only, normalize=TRUE) ===\n")
cat("  1/5 = 0.20 is less damaging than 1/10 = 0.10\n\n")
X_gmm_5 <- dat$X_gmm[, 1:5, drop=FALSE]
for (gam in c(0, 0.5, 1.0, 2.0)) {
  fit <- tryCatch(gemcox_full_multistart_shared_v4(
    X_gmm=X_gmm_5, X_cox=dat$X_cox, time=dat$time, status=dat$status,
    K=2, alpha=0.5, max_iter=80, surv_weight=gam,
    normalize_gmm_by_dim=TRUE,
    verbose=FALSE, n_starts=10), error=function(e)NULL)
  if(!is.null(fit)) {
    ari <- adjustedRandIndex(dat$Z_true,fit$clusterid)
    cat(sprintf("  gamma=%4.2f  ARI=%.4f  pi=[%.2f,%.2f]",
      gam, ari, fit$pi[1], fit$pi[2]))
    if (ari > 0.3) cat("  ***RECOVERED***")
    cat("\n")
  }
}

cat("\n=== ALSO: q=3 signal features, normalize=TRUE ===\n")
X_gmm_3 <- dat$X_gmm[, 1:3, drop=FALSE]
for (gam in c(0, 0.5, 1.0, 2.0)) {
  fit <- tryCatch(gemcox_full_multistart_shared_v4(
    X_gmm=X_gmm_3, X_cox=dat$X_cox, time=dat$time, status=dat$status,
    K=2, alpha=0.5, max_iter=80, surv_weight=gam,
    normalize_gmm_by_dim=TRUE,
    verbose=FALSE, n_starts=10), error=function(e)NULL)
  if(!is.null(fit)) {
    ari <- adjustedRandIndex(dat$Z_true,fit$clusterid)
    cat(sprintf("  gamma=%4.2f  ARI=%.4f  pi=[%.2f,%.2f]",
      gam, ari, fit$pi[1], fit$pi[2]))
    if (ari > 0.3) cat("  ***RECOVERED***")
    cat("\n")
  }
}

cat("\n=== REFERENCE: k-means on same data ===\n")
km10 <- kmeans(dat$X_gmm, 2, nstart=50)
km5  <- kmeans(dat$X_gmm[,1:5], 2, nstart=50)
km3  <- kmeans(dat$X_gmm[,1:3], 2, nstart=50)
cat(sprintf("  k-means q=10: ARI=%.4f\n", adjustedRandIndex(dat$Z_true, km10$cluster)))
cat(sprintf("  k-means q=5:  ARI=%.4f\n", adjustedRandIndex(dat$Z_true, km5$cluster)))
cat(sprintf("  k-means q=3:  ARI=%.4f\n", adjustedRandIndex(dat$Z_true, km3$cluster)))
