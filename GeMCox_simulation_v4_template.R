###############################################################
## GeM-Cox Simulation v4 template
##
## Purpose:
##   - Move the simulation closer to CVIA078
##   - Avoid K=4 by default for n = 117
##   - Separate clustering features from Cox features
##   - Optionally impose CVIA078-like missingness and tied follow-up times
##   - Use the fold-safe / shared-baseline pipeline in
##     GeMCox_improved_v4_patch.R
##
## Typical workflow:
##   1) Read the real CVIA078 long-format data
##   2) Build an empirical template from interpretable features
##   3) Simulate transformed-scale datasets from that template
##   4) Fit GeM-Cox with log_transform = FALSE inside simulation,
##      because the features are already on the transformed scale
###############################################################

required_pkgs <- c("survival", "glmnet", "mclust", "MASS", "dplyr", "tidyr",
                   "ggplot2", "patchwork", "tibble", "scales", "parallel")
user_r_lib <- Sys.getenv("R_LIBS_USER")
if (!nzchar(user_r_lib)) {
  user_r_lib <- file.path(Sys.getenv("LOCALAPPDATA"), "R", "win-library",
                          paste(R.version$major, sub("\\..*$", "", R.version$minor), sep = "."))
}
if (!dir.exists(user_r_lib)) dir.create(user_r_lib, recursive = TRUE, showWarnings = FALSE)
.libPaths(c(user_r_lib, .libPaths()))
for (pkg in required_pkgs) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg, repos = "https://cloud.r-project.org", lib = user_r_lib)
  }
}

library(survival)
library(glmnet)
library(mclust)
library(MASS)
library(dplyr)
library(tidyr)
library(ggplot2)
library(patchwork)
library(tibble)
library(scales)
library(parallel)

engine_candidates <- c("GeMCox_improved_v4_patch.R", "GeMCox_improved_v4.R")
engine_path <- engine_candidates[file.exists(engine_candidates)][1]
if (!nzchar(engine_path)) stop("Could not find GeMCox_improved_v4_patch.R or GeMCox_improved_v4.R")
source(engine_path)

###############################################################
## 1. Build an empirical CVIA078 template
###############################################################

build_cvia078_template_v4 <- function(dat_long,
                                      gmm_feature_mode = "delta_only",
                                      cox_feature_mode = "v20_plus_delta_covars",
                                      missing_threshold = 0.05,
                                      corr_prune = 0.85) {
  real_dat <- prepare_cvia078_data_v4(
    dat_long,
    gmm_feature_mode = gmm_feature_mode,
    cox_feature_mode = cox_feature_mode,
    missing_threshold = missing_threshold,
    corr_prune = corr_prune
  )
  
  ## Put data on the modeling scale used inside CV folds, then keep that scale
  ## for simulation to avoid an unnecessary transform/inverse-transform cycle.
  prep_scale <- preprocess_train_v4(
    real_dat$X_gmm,
    real_dat$X_cox,
    log_transform = TRUE,
    skew_threshold = 1.0,
    impute = "median"
  )
  
  Xg <- prep_scale$X_gmm
  Xc <- prep_scale$X_cox
  
  cov_ridge <- function(S, ridge = 1e-4) {
    S <- as.matrix(S)
    S[!is.finite(S)] <- 0
    S <- (S + t(S)) / 2
    S + diag(ridge, nrow(S))
  }
  
  list(
    X_gmm_center = colMeans(Xg),
    X_gmm_cov = cov_ridge(stats::cov(Xg)),
    X_cox_center = colMeans(Xc),
    X_cox_cov = cov_ridge(stats::cov(Xc)),
    gmm_features = colnames(Xg),
    cox_features = colnames(Xc),
    gmm_missing_rate = colMeans(!is.finite(as.matrix(real_dat$X_gmm))),
    cox_missing_rate = colMeans(!is.finite(as.matrix(real_dat$X_cox))),
    time_support = sort(unique(real_dat$time)),
    time_event_support = sort(unique(real_dat$time[real_dat$status == 1])),
    n = nrow(Xg),
    event_rate = mean(real_dat$status),
    real_data = real_dat,
    prep_scale = prep_scale
  )
}

###############################################################
## 2. CVIA078-like simulator
###############################################################

nearest_support_v4 <- function(x, support) {
  if (length(support) == 0) return(x)
  vapply(x, function(z) support[which.min(abs(support - z))], numeric(1))
}

simulate_cvia078_like_v4 <- function(template,
                                     n = template$n,
                                     K_true = 2,
                                     pi_c = NULL,
                                     separation = c("low", "med", "high"),
                                     signal_gmm = NULL,
                                     signal_cox = NULL,
                                     delta_mu = NULL,
                                     delta_beta = NULL,
                                     event_rate_target = template$event_rate,
                                     baseline_shape = 1.5,
                                     baseline_mode = c("shared", "mild_cluster_shift"),
                                     add_missing = TRUE,
                                     seed = 1) {
  separation <- match.arg(separation)
  baseline_mode <- match.arg(baseline_mode)
  set.seed(seed)
  
  if (is.null(pi_c)) pi_c <- rep(1 / K_true, K_true)
  
  if (is.null(delta_mu)) {
    delta_mu <- switch(separation, low = 0.35, med = 0.65, high = 1.0)
  }
  if (is.null(delta_beta)) {
    delta_beta <- switch(separation, low = 0.25, med = 0.45, high = 0.70)
  }
  
  q <- length(template$gmm_features)
  p <- length(template$cox_features)
  if (is.null(signal_gmm)) signal_gmm <- template$gmm_features[seq_len(min(6, q))]
  if (is.null(signal_cox)) signal_cox <- intersect(template$cox_features, signal_gmm)
  if (length(signal_cox) == 0) signal_cox <- template$cox_features[seq_len(min(6, p))]
  
  idx_gmm <- match(signal_gmm, template$gmm_features)
  idx_gmm <- idx_gmm[is.finite(idx_gmm)]
  idx_cox <- match(signal_cox, template$cox_features)
  idx_cox <- idx_cox[is.finite(idx_cox)]
  
  mu0_g <- template$X_gmm_center
  mu0_c <- template$X_cox_center
  Sg <- template$X_gmm_cov
  Sc <- template$X_cox_cov
  
  positions <- seq(-(K_true - 1) / 2, (K_true - 1) / 2, length.out = K_true)
  mu_g_list <- lapply(seq_len(K_true), function(k) {
    mu <- mu0_g
    if (length(idx_gmm) > 0) {
      mu[idx_gmm] <- mu[idx_gmm] + positions[k] * delta_mu / sqrt(seq_along(idx_gmm))
    }
    mu
  })
  
  beta_base <- rep(0, p)
  if (length(idx_cox) > 0) {
    signs <- rep(c(1, -1, 1, -1, 1, -1), length.out = length(idx_cox))
    beta_base[idx_cox] <- signs * seq(0.6, 0.25, length.out = length(idx_cox))
  }
  beta_list <- lapply(seq_len(K_true), function(k) {
    b <- beta_base
    if (length(idx_cox) > 0) {
      b[idx_cox] <- b[idx_cox] + positions[k] * delta_beta / sqrt(seq_along(idx_cox))
    }
    b
  })
  
  Z <- sample(seq_len(K_true), size = n, replace = TRUE, prob = pi_c)
  X_gmm <- matrix(NA_real_, n, q, dimnames = list(NULL, template$gmm_features))
  X_cox <- matrix(NA_real_, n, p, dimnames = list(NULL, template$cox_features))
  
  for (k in seq_len(K_true)) {
    idx <- which(Z == k)
    if (length(idx) == 0) next
    X_gmm[idx, ] <- MASS::mvrnorm(length(idx), mu = mu_g_list[[k]], Sigma = Sg)
    
    ## Generate Cox features around the empirical center, but align shared signal
    ## columns with the same cluster location whenever names overlap.
    Xc_k <- MASS::mvrnorm(length(idx), mu = mu0_c, Sigma = Sc)
    common_nm <- intersect(template$gmm_features, template$cox_features)
    if (length(common_nm) > 0) {
      jg <- match(common_nm, template$gmm_features)
      jc <- match(common_nm, template$cox_features)
      Xc_k[, jc] <- X_gmm[idx, jg, drop = FALSE]
    }
    X_cox[idx, ] <- Xc_k
  }
  
  ## Optional missingness mimicking the real data on the preprocessed scale.
  if (add_missing) {
    for (j in seq_len(q)) {
      miss_j <- stats::rbinom(n, 1, template$gmm_missing_rate[j]) == 1
      if (any(miss_j)) X_gmm[miss_j, j] <- NA_real_
    }
    for (j in seq_len(p)) {
      miss_j <- stats::rbinom(n, 1, template$cox_missing_rate[j]) == 1
      if (any(miss_j)) X_cox[miss_j, j] <- NA_real_
    }
  }
  
  ## Event times on the Cox-feature scale.
  scale_pos <- if (baseline_mode == "shared") rep(0, K_true) else positions * 0.25
  baseline_scale <- 150 * exp(scale_pos)
  
  T_true <- numeric(n)
  for (k in seq_len(K_true)) {
    idx <- which(Z == k)
    if (length(idx) == 0) next
    Xc_obs <- X_cox[idx, , drop = FALSE]
    Xc_imp <- Xc_obs
    for (j in seq_len(ncol(Xc_imp))) {
      z <- Xc_imp[, j]
      z[!is.finite(z)] <- mu0_c[j]
      Xc_imp[, j] <- z
    }
    eta <- as.numeric(Xc_imp %*% beta_list[[k]])
    eta <- pmin(pmax(eta, -4), 4)
    U <- stats::runif(length(idx))
    T_true[idx] <- baseline_scale[k] * (-log(U))^(1 / baseline_shape) * exp(-eta / baseline_shape)
  }
  
  calibrate_cens_v4 <- function(target, Tlat, tol = 0.005) {
    lo <- 1e-6
    hi <- 2
    for (it in 1:60) {
      mid <- (lo + hi) / 2
      C_try <- stats::rexp(length(Tlat), rate = mid)
      er <- mean(Tlat <= C_try)
      if (abs(er - target) < tol) break
      if (er < target) hi <- mid else lo <- mid
    }
    mid
  }
  
  cens_rate <- calibrate_cens_v4(event_rate_target, T_true)
  C_time <- stats::rexp(n, rate = cens_rate)
  time <- pmin(T_true, C_time)
  status <- as.integer(T_true <= C_time)
  
  ## Impose CVIA078-like ties by snapping to observed follow-up support.
  support <- if (length(template$time_support) > 0) template$time_support else sort(unique(round(time)))
  time <- nearest_support_v4(time, support)
  time <- pmax(time, min(support))
  
  list(
    X_gmm = X_gmm,
    X_cox = X_cox,
    time = time,
    status = status,
    Z_true = Z,
    params = list(
      K_true = K_true,
      pi_c = pi_c,
      beta_list = beta_list,
      mu_g_list = mu_g_list,
      delta_mu = delta_mu,
      delta_beta = delta_beta,
      separation = separation,
      event_rate_actual = mean(status),
      baseline_mode = baseline_mode,
      signal_gmm = signal_gmm,
      signal_cox = signal_cox
    )
  )
}

###############################################################
## 3. Single replicate and study runner
###############################################################

run_replicate_cvia078_v4 <- function(sim_data,
                                     K_grid = 1:3,
                                     gamma_grid = c(0, 0.25, 0.5, 1.0),
                                     nfolds = 5,
                                     alpha = 0.5,
                                     max_iter = 60,
                                     n_starts = 5,
                                     n_starts_final = 10,
                                     baseline_mode = "shared",
                                     verbose = FALSE) {
  out <- list(
    K_true = sim_data$params$K_true,
    K_selected = NA_integer_,
    K_correct = NA_integer_,
    gamma_selected = NA_real_,
    ARI = NA_real_,
    cindex_cv = NA_real_,
    cindex_full = NA_real_,
    n_events = sum(sim_data$status),
    event_rate_actual = mean(sim_data$status),
    error_msg = NA_character_
  )
  
  fit <- tryCatch(
    gemcox_pipeline_v4(
      X_gmm = sim_data$X_gmm,
      X_cox = sim_data$X_cox,
      time = sim_data$time,
      status = sim_data$status,
      K_grid = K_grid,
      gamma_grid = gamma_grid,
      alpha = alpha,
      log_transform = FALSE,
      use_1se_rule = TRUE,
      nfolds = nfolds,
      n_starts = n_starts,
      n_starts_final = n_starts_final,
      max_iter = max_iter,
      verbose = verbose,
      baseline_mode = baseline_mode
    ),
    error = function(e) e
  )
  
  if (inherits(fit, "error")) {
    out$error_msg <- fit$message
    return(out)
  }
  
  out$K_selected <- fit$K_selected
  out$K_correct <- as.integer(out$K_selected == out$K_true)
  out$gamma_selected <- fit$gamma_selected
  out$cindex_cv <- fit$cv_result$best$mean_score[1]
  
  tau <- fit$fit$tau
  K_hat <- fit$K_selected
  eta_mat <- sapply(seq_len(K_hat), function(k) {
    eta_from_scaled(fit$preprocessing$X_cox,
                    fit$fit$coxfit[[k]]$beta_s,
                    fit$fit$coxfit[[k]]$scaler)
  })
  if (is.vector(eta_mat)) eta_mat <- matrix(eta_mat, ncol = K_hat)
  out$cindex_full <- c_index(sim_data$time, sim_data$status, rowSums(tau * eta_mat))
  
  if (out$K_correct == 1 && out$K_true > 1) {
    out$ARI <- mclust::adjustedRandIndex(sim_data$Z_true, fit$fit$clusterid)
  }
  
  out
}

resolve_study_profile_v4 <- function(profile = c("fast_10min", "diagnostic", "full"),
                                     nsim = NULL) {
  profile <- match.arg(profile)
  cfg <- switch(
    profile,
    fast_10min = list(
      profile = profile,
      nsim = 12L,
      K_grid = 1:3,
      gamma_grid = c(0, 0.5, 1.0),
      nfolds = 2L,
      alpha = 0.5,
      max_iter = 25L,
      n_starts = 2L,
      n_starts_final = 3L
    ),
    diagnostic = list(
      profile = profile,
      nsim = 24L,
      K_grid = 1:3,
      gamma_grid = c(0, 0.25, 0.5, 1.0),
      nfolds = 3L,
      alpha = 0.5,
      max_iter = 35L,
      n_starts = 3L,
      n_starts_final = 4L
    ),
    full = list(
      profile = profile,
      nsim = 50L,
      K_grid = 1:3,
      gamma_grid = c(0, 0.25, 0.5, 1.0),
      nfolds = 5L,
      alpha = 0.5,
      max_iter = 60L,
      n_starts = 5L,
      n_starts_final = 10L
    )
  )
  if (!is.null(nsim)) cfg$nsim <- as.integer(nsim)
  cfg
}

run_parallel_jobs_v4 <- function(job_grid, worker_fun, n_cores) {
  n_jobs <- nrow(job_grid)
  if (n_jobs == 0) return(list())
  if (n_cores <= 1L) return(lapply(seq_len(n_jobs), worker_fun))
  if (.Platform$OS.type == "windows") {
    cl <- makeCluster(n_cores, type = "PSOCK")
    on.exit(stopCluster(cl))
    clusterExport(cl, ls(.GlobalEnv), .GlobalEnv)
    clusterEvalQ(cl, {
      library(survival)
      library(glmnet)
      library(mclust)
      library(MASS)
      library(dplyr)
      library(tidyr)
      library(ggplot2)
      library(patchwork)
      library(tibble)
      library(scales)
      library(parallel)
    })
    parLapplyLB(cl, seq_len(n_jobs), worker_fun)
  } else {
    mclapply(seq_len(n_jobs), worker_fun, mc.cores = n_cores)
  }
}

estimate_wallclock_v4 <- function(n_jobs, cfg, n_cores) {
  fits_per_rep <- length(cfg$K_grid) * length(cfg$gamma_grid) * cfg$nfolds * cfg$n_starts +
    cfg$n_starts_final
  sec_per_fit <- switch(cfg$profile,
                        fast_10min = 0.45,
                        diagnostic = 0.60,
                        full = 0.80,
                        0.60)
  ceiling(n_jobs * fits_per_rep * sec_per_fit / max(1, n_cores) / 60)
}

run_simulation_study_cvia078_v4 <- function(template,
                                            scenario_grid = NULL,
                                            nsim = NULL,
                                            seed = 1,
                                            baseline_mode = "shared",
                                            profile = c("fast_10min", "diagnostic", "full"),
                                            n_cores = min(24L, max(1L, parallel::detectCores(logical = FALSE))),
                                            save_outputs = TRUE,
                                            output_prefix = "GeMCox_simulation_v4_template") {
  cfg <- resolve_study_profile_v4(profile = profile, nsim = nsim)
  if (is.null(scenario_grid)) {
    scenario_grid <- expand.grid(
      K_true = c(1, 2, 3),
      separation = c("low", "med", "high"),
      event_rate_target = c(template$event_rate, 0.30),
      n = template$n,
      feat_corr = 0.2,
      stringsAsFactors = FALSE
    )
  }
  if (!"n" %in% names(scenario_grid)) scenario_grid$n <- template$n
  if (!"feat_corr" %in% names(scenario_grid)) scenario_grid$feat_corr <- 0.2
  scenario_grid$scenario_id <- seq_len(nrow(scenario_grid))

  job_grid <- tidyr::crossing(scenario_grid, rep = seq_len(cfg$nsim))
  n_cores <- as.integer(max(1L, min(n_cores, nrow(job_grid))))

  message(sprintf("Running CVIA078-like study: profile=%s | scenarios=%d | reps=%d | jobs=%d | cores=%d",
                  cfg$profile, nrow(scenario_grid), cfg$nsim, nrow(job_grid), n_cores))
  message(sprintf("Settings: folds=%d | starts=%d/%d | gamma=%s | max_iter=%d",
                  cfg$nfolds, cfg$n_starts, cfg$n_starts_final,
                  paste(cfg$gamma_grid, collapse = ","), cfg$max_iter))
  message(sprintf("Estimated wall-clock (heuristic): ~%d min",
                  estimate_wallclock_v4(nrow(job_grid), cfg, n_cores)))

  all_res <- run_parallel_jobs_v4(job_grid, function(i) {
    sc <- job_grid[i, ]
    dat <- tryCatch(
      simulate_cvia078_like_v4(
        template = template,
        n = sc$n,
        K_true = sc$K_true,
        separation = sc$separation,
        event_rate_target = sc$event_rate_target,
        baseline_mode = baseline_mode,
        seed = seed + 1000 * sc$scenario_id + sc$rep
      ),
      error = function(e) e
    )

    if (inherits(dat, "error")) {
      return(data.frame(
        scenario_id = sc$scenario_id,
        rep = sc$rep,
        K_true = sc$K_true,
        separation = sc$separation,
        event_rate_target = sc$event_rate_target,
        n = sc$n,
        feat_corr = sc$feat_corr,
        baseline_mode = baseline_mode,
        error_msg = dat$message,
        stringsAsFactors = FALSE
      ))
    }

    ans <- tryCatch(
      run_replicate_cvia078_v4(
        dat,
        K_grid = cfg$K_grid,
        gamma_grid = cfg$gamma_grid,
        nfolds = cfg$nfolds,
        alpha = cfg$alpha,
        max_iter = cfg$max_iter,
        n_starts = cfg$n_starts,
        n_starts_final = cfg$n_starts_final,
        baseline_mode = baseline_mode,
        verbose = FALSE
      ),
      error = function(e) list(error_msg = e$message)
    )
    ans <- lapply(ans, function(v) if (is.null(v) || length(v) == 0) NA else v[1])
    data.frame(
      scenario_id = sc$scenario_id,
      rep = sc$rep,
      K_true = sc$K_true,
      separation = sc$separation,
      event_rate_target = sc$event_rate_target,
      n = sc$n,
      feat_corr = sc$feat_corr,
      baseline_mode = baseline_mode,
      ans,
      stringsAsFactors = FALSE
    )
  }, n_cores = n_cores)

  raw_df <- dplyr::bind_rows(all_res)
  em <- raw_df$error_msg
  ok <- is.na(em) | (nzchar(trimws(em)) == FALSE) | (trimws(em) == "NA")
  df_ok <- raw_df[ok, , drop = FALSE]

  summary_df <- df_ok |>
    dplyr::group_by(K_true, separation, event_rate_target, n, feat_corr, baseline_mode) |>
    dplyr::summarise(
      n_reps = dplyr::n(),
      K_sel_accuracy = mean(K_correct, na.rm = TRUE),
      K_sel_mean = mean(K_selected, na.rm = TRUE),
      ARI_mean_given_correctK = mean(ARI, na.rm = TRUE),
      cindex_cv_mean = mean(cindex_cv, na.rm = TRUE),
      cindex_full_mean = mean(cindex_full, na.rm = TRUE),
      gamma_mean = mean(gamma_selected, na.rm = TRUE),
      .groups = "drop"
    )

  figs <- make_figures_cvia078_v4(summary_df, raw_df,
                                  fig_dir = file.path("figures", output_prefix),
                                  prefix = output_prefix)

  if (isTRUE(save_outputs)) {
    saveRDS(raw_df, paste0(output_prefix, "_raw.rds"))
    write.csv(summary_df, paste0(output_prefix, "_summary.csv"), row.names = FALSE)
  }

  list(
    raw = raw_df,
    summary = summary_df,
    figures = figs,
    ok_rows = df_ok,
    scenario_grid = scenario_grid,
    run_config = cfg
  )
}

###############################################################
## Section 6: Figures
###############################################################

theme_gemcox <- function() {
  theme_bw(base_size = 12) +
    theme(
      panel.grid.minor = element_blank(),
      strip.background = element_rect(fill = "#1a3a5c", color = NA),
      strip.text = element_text(color = "white", face = "bold"),
      legend.position = "bottom",
      plot.title = element_text(face = "bold", size = 13),
      plot.subtitle = element_text(color = "grey40", size = 10)
    )
}

sep_levels <- c("low", "med", "high")
sep_colors <- c(low = "#e08060", med = "#5b9bd5", high = "#3a7a3a")

prep_summary <- function(df, rho = 0.2, n_val = NULL) {
  out <- df %>%
    filter(feat_corr == rho)
  if (!is.null(n_val)) out <- out %>% filter(n == n_val)
  out %>%
    mutate(
      separation = factor(separation, levels = sep_levels),
      event_rate_label = paste0("Event rate: ", event_rate_target),
      K_label = paste0("K = ", K_true),
      n_label = paste0("n = ", n)
    ) %>%
    droplevels()
}

prep_results <- function(df, rho = 0.2, n_val = 117) {
  df %>%
    filter(feat_corr == rho, n == n_val) %>%
    mutate(
      separation = factor(separation, levels = sep_levels),
      event_rate_label = paste0("Event rate: ", event_rate_target),
      K_label = paste0("K = ", K_true)
    ) %>%
    droplevels()
}

make_figures_cvia078_v4 <- function(summary_df,
                                    results_df,
                                    rho = 0.2,
                                    n_val = NULL,
                                    fig_dir = "figures",
                                    prefix = "GeMCox_simulation_v4_template") {
  dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

  figs <- list()
  summary_plot_df <- prep_summary(summary_df, rho = rho, n_val = n_val)
  n_plot <- if (is.null(n_val)) sort(unique(results_df$n))[1] else n_val
  figs$results_df_plot <- prep_results(results_df, rho = rho, n_val = n_plot)

  if (nrow(summary_plot_df) == 0) {
    warning("No rows matched the plotting filters; figures were not created.")
    return(figs)
  }

  fig_ksel <- summary_plot_df %>%
    ggplot(aes(x = separation, y = K_sel_accuracy,
               color = separation, group = separation)) +
    geom_point(size = 3) +
    geom_hline(yintercept = 1 / 4, linetype = "dashed", color = "grey60",
               linewidth = 0.5) +
    facet_grid(n_label ~ K_label + event_rate_label) +
    scale_color_manual(values = sep_colors, name = "Separation") +
    scale_y_continuous(limits = c(0, 1), labels = scales::percent) +
    labs(
      title = "K-Selection Accuracy",
      subtitle = "Proportion of replicates where CV correctly identified K_true  |  dashed = chance",
      x = "Cluster Separation",
      y = "P(K-hat = K_true)"
    ) +
    theme_gemcox()

  figs$fig_ksel <- fig_ksel

  ggsave(file.path(fig_dir, paste0(prefix, "_fig1_k_selection_accuracy.png")),
         plot = fig_ksel, width = 12, height = 6.5, dpi = 300)
  ggsave(file.path(fig_dir, paste0(prefix, "_fig1_k_selection_accuracy.pdf")),
         plot = fig_ksel, width = 12, height = 6.5)

  figs
}

cat("GeMCox_simulation_v4_template.R loaded.\n")
cat("  Main idea: use build_cvia078_template_v4() then simulate_cvia078_like_v4().\n")
cat("  Use run_simulation_study_cvia078_v4(profile = 'fast_10min') for the 24-core study run.\n")
