###############################################################
## GeM-Cox v4 Simulation Figures
##
## Inputs:
##   - GeMCox_simulation_v4_raw.rds
##   - GeMCox_simulation_v4_summary.csv
##   - GeMCox_simulation_v4_summary_by_nK.csv (optional)
##
## Outputs:
##   - fig1_K_selection_v4.pdf
##   - fig2_ARI_v4.pdf
##   - fig3_cindex_gain_v4.pdf
##   - fig4_gamma_selected_v4.pdf
##   - fig5_dashboard_v4.pdf
##   - fig6_cindex_scatter_v4.pdf
###############################################################

required_pkgs <- c("ggplot2", "dplyr", "tidyr", "patchwork", "scales")
for (pkg in required_pkgs) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg, repos = "https://cloud.r-project.org")
  }
}

library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)
library(scales)

## ---- Load data ----
raw_df  <- readRDS("GeMCox_simulation_v4_raw.rds")
summ_df <- read.csv("GeMCox_simulation_v4_summary.csv", stringsAsFactors = FALSE)

summ_nk_df <- NULL
if (file.exists("GeMCox_simulation_v4_summary_by_nK.csv")) {
  summ_nk_df <- read.csv("GeMCox_simulation_v4_summary_by_nK.csv", stringsAsFactors = FALSE)
}

## ---- Clean raw data (drop error rows) ----
em <- raw_df$error_msg
ok <- is.na(em) | !nzchar(trimws(em)) | trimws(em) == "NA"
raw_df <- raw_df[ok, , drop = FALSE]

## ---- Shared theme and palettes ----
sep_levels <- c("low", "med", "high")
sep_colors <- c(low = "#e08060", med = "#5b9bd5", high = "#3a7a3a")

theme_gem <- function(base_size = 12) {
  theme_bw(base_size = base_size) +
    theme(
      panel.grid.minor = element_blank(),
      strip.background = element_rect(fill = "#1a3a5c", color = NA),
      strip.text = element_text(color = "white", face = "bold", size = rel(0.95)),
      legend.position = "bottom",
      plot.title = element_text(face = "bold", size = 14),
      plot.subtitle = element_text(color = "grey40", size = 10),
      axis.title = element_text(size = 11),
      legend.title = element_text(size = 10),
      legend.text = element_text(size = 9)
    )
}

prep_v4 <- function(df) {
  df %>%
    mutate(
      separation = factor(separation, levels = sep_levels),
      n_label = factor(paste0("n = ", n), levels = paste0("n = ", c(117, 250, 500))),
      K_label = factor(paste0("K = ", K_true), levels = paste0("K = ", 1:3)),
      event_rate_label = factor(
        paste0("Event rate: ", sprintf("%.2f", event_rate)),
        levels = paste0("Event rate: ", sprintf("%.2f", c(0.44, 0.25)))
      )
    )
}

summ <- prep_v4(summ_df)
raw  <- prep_v4(raw_df)

## For v4 raw output, ARI is conditional on K_selected == K_true and K_true > 1
raw <- raw %>%
  mutate(
    ARI_defined = ifelse(K_true > 1 & K_selected == K_true & !is.na(ARI), 1, 0),
    gain = cindex_full - cindex_cox1
  )

###############################################################
## Figure 1: K-selection accuracy
###############################################################

fig1 <- summ %>%
  ggplot(aes(x = separation, y = K_sel_accuracy, fill = separation)) +
  geom_col(width = 0.65, alpha = 0.88) +
  geom_hline(yintercept = 1/3, linetype = "dashed", color = "grey55", linewidth = 0.4) +
  facet_grid(event_rate_label ~ n_label + K_label) +
  scale_fill_manual(values = sep_colors, name = "Separation") +
  scale_y_continuous(limits = c(0, 1.02), breaks = seq(0, 1, 0.25), labels = percent) +
  labs(
    title = "K-selection accuracy",
    subtitle = "CV + 1-SE rule; dashed line = chance baseline for K grid {1,2,3}",
    x = "Cluster separation",
    y = expression(P(hat(K) == K[true]))
  ) +
  theme_gem(10)

ggsave("fig1_K_selection_v4.pdf", fig1, width = 14, height = 6.8)
cat("Saved fig1_K_selection_v4.pdf\n")

###############################################################
## Figure 2: ARI (conditional on correct K selection)
###############################################################

ari_data <- raw %>%
  filter(K_true > 1, K_selected == K_true, !is.na(ARI))

if (nrow(ari_data) > 5) {
  fig2 <- ari_data %>%
    ggplot(aes(x = separation, y = ARI, fill = separation)) +
    geom_violin(alpha = 0.5, trim = FALSE, scale = "width", linewidth = 0.3) +
    geom_boxplot(width = 0.15, outlier.size = 0.8, fill = "white", alpha = 0.85) +
    facet_grid(event_rate_label ~ n_label + K_label) +
    scale_fill_manual(values = sep_colors, name = "Separation") +
    scale_y_continuous(limits = c(-0.1, 1), breaks = seq(0, 1, 0.25)) +
    labs(
      title = "Cluster recovery (ARI)",
      subtitle = "Conditional on K_selected = K_true and K_true > 1",
      x = "Cluster separation",
      y = "Adjusted Rand Index"
    ) +
    theme_gem(10)
} else {
  fig2 <- summ %>%
    filter(K_true > 1, !is.na(ARI_mean)) %>%
    ggplot(aes(x = separation, y = ARI_mean,
               ymin = pmax(ARI_mean - ARI_sd, -0.05),
               ymax = pmin(ARI_mean + ARI_sd, 1),
               color = separation)) +
    geom_pointrange(size = 0.5, linewidth = 0.65) +
    facet_grid(event_rate_label ~ n_label + K_label) +
    scale_color_manual(values = sep_colors, name = "Separation") +
    scale_y_continuous(limits = c(-0.1, 1)) +
    labs(
      title = "Cluster recovery (ARI)",
      subtitle = "Mean ± SD; conditional on K_selected = K_true and K_true > 1",
      x = "Cluster separation",
      y = "Adjusted Rand Index"
    ) +
    theme_gem(10)
}

ggsave("fig2_ARI_v4.pdf", fig2, width = 14, height = 6.8)
cat("Saved fig2_ARI_v4.pdf\n")

###############################################################
## Figure 3: C-index gain over single Cox
###############################################################

gain_data <- raw %>%
  filter(!is.na(cindex_full), !is.na(cindex_cox1), !is.na(gain))

fig3 <- gain_data %>%
  ggplot(aes(x = separation, y = gain, fill = separation)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey55", linewidth = 0.4) +
  geom_violin(alpha = 0.45, trim = FALSE, scale = "width", linewidth = 0.3) +
  geom_boxplot(width = 0.15, outlier.size = 0.8, fill = "white", alpha = 0.85) +
  facet_grid(event_rate_label ~ n_label + K_label) +
  scale_fill_manual(values = sep_colors, name = "Separation") +
  labs(
    title = "C-index gain over single Cox",
    subtitle = "Positive values indicate GeM-Cox outperformed the one-component Cox model",
    x = "Cluster separation",
    y = expression(Delta ~ "C-index")
  ) +
  theme_gem(10)

ggsave("fig3_cindex_gain_v4.pdf", fig3, width = 14, height = 6.8)
cat("Saved fig3_cindex_gain_v4.pdf\n")

###############################################################
## Figure 4: Selected gamma distribution
###############################################################

gamma_data <- raw %>%
  filter(!is.na(gamma_selected)) %>%
  mutate(
    gamma_label = factor(
      sprintf("%.2f", gamma_selected),
      levels = sprintf("%.2f", sort(unique(gamma_selected)))
    )
  )

fig4 <- gamma_data %>%
  ggplot(aes(x = gamma_label, fill = separation)) +
  geom_bar(position = "dodge", alpha = 0.88) +
  facet_grid(event_rate_label ~ n_label + K_label) +
  scale_fill_manual(values = sep_colors, name = "Separation") +
  labs(
    title = "Selected gamma distribution",
    subtitle = "How often each gamma value was selected by CV",
    x = expression(gamma ~ "(survival weight)"),
    y = "Count"
  ) +
  theme_gem(10)

ggsave("fig4_gamma_selected_v4.pdf", fig4, width = 14, height = 6.8)
cat("Saved fig4_gamma_selected_v4.pdf\n")

###############################################################
## Figure 5: Dashboard
###############################################################

pa <- summ %>%
  ggplot(aes(x = n_label, y = K_sel_accuracy,
             color = separation, group = separation)) +
  geom_line(linewidth = 0.7) +
  geom_point(size = 2.3) +
  facet_grid(event_rate_label ~ K_label) +
  geom_hline(yintercept = 1/3, linetype = "dashed", color = "grey60") +
  scale_color_manual(values = sep_colors, name = "Separation") +
  scale_y_continuous(limits = c(0, 1), labels = percent) +
  labs(
    title = "A. K-selection accuracy",
    x = "Sample size",
    y = expression(P(hat(K) == K[true]))
  ) +
  theme_gem(10)

pb <- summ %>%
  filter(K_true > 1) %>%
  ggplot(aes(x = n_label, y = ARI_mean,
             color = separation, group = separation)) +
  geom_line(linewidth = 0.7) +
  geom_point(size = 2.3) +
  facet_grid(event_rate_label ~ K_label) +
  scale_color_manual(values = sep_colors, name = "Separation") +
  scale_y_continuous(limits = c(-0.05, 1)) +
  labs(
    title = "B. ARI (conditional on correct K)",
    x = "Sample size",
    y = "Mean ARI"
  ) +
  theme_gem(10)

pc <- summ %>%
  ggplot(aes(x = n_label, y = cindex_gain_mean,
             color = separation, group = separation)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey60") +
  geom_line(linewidth = 0.7) +
  geom_point(size = 2.3) +
  facet_grid(event_rate_label ~ K_label) +
  scale_color_manual(values = sep_colors, name = "Separation") +
  labs(
    title = "C. C-index gain over single Cox",
    x = "Sample size",
    y = expression(Delta ~ "C-index")
  ) +
  theme_gem(10)

fig5 <- pa / pb / pc +
  plot_annotation(
    title = "GeM-Cox v4 simulation dashboard",
    subtitle = "n sweep + K sweep; shared baseline; elastic net alpha = 0.5; normalize_gmm_by_dim = FALSE in simulation fitting",
    theme = theme(
      plot.title = element_text(face = "bold", size = 15),
      plot.subtitle = element_text(color = "grey35", size = 11)
    )
  )

ggsave("fig5_dashboard_v4.pdf", fig5, width = 12, height = 12)
cat("Saved fig5_dashboard_v4.pdf\n")

###############################################################
## Figure 6: Per-replicate C-index comparison
###############################################################

scatter_data <- raw %>%
  filter(!is.na(cindex_full), !is.na(cindex_cox1), K_true > 1)

if (nrow(scatter_data) > 10) {
  fig6 <- scatter_data %>%
    ggplot(aes(x = cindex_cox1, y = cindex_full, color = separation)) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey60") +
    geom_point(alpha = 0.55, size = 1.3) +
    facet_grid(event_rate_label ~ n_label + K_label) +
    scale_color_manual(values = sep_colors, name = "Separation") +
    coord_equal(xlim = c(0.45, 0.95), ylim = c(0.45, 0.95)) +
    labs(
      title = "GeM-Cox vs single Cox: per-replicate C-index",
      subtitle = "Points above diagonal indicate better discrimination for GeM-Cox",
      x = "Single Cox C-index",
      y = "GeM-Cox C-index"
    ) +
    theme_gem(10)
  
  ggsave("fig6_cindex_scatter_v4.pdf", fig6, width = 14, height = 6.8)
  cat("Saved fig6_cindex_scatter_v4.pdf\n")
}

###############################################################
## Optional compact summary plot if summary_by_nK exists
###############################################################
if (!is.null(summ_nk_df)) {
  summ_nk <- summ_nk_df %>%
    mutate(
      n_label = factor(paste0("n = ", n), levels = paste0("n = ", c(117, 250, 500))),
      K_label = factor(paste0("K = ", K_true), levels = paste0("K = ", 1:3))
    )
  
  fig_extra <- summ_nk %>%
    ggplot(aes(x = n_label, y = K_sel_accuracy, group = K_label, color = K_label)) +
    geom_line(linewidth = 0.8) +
    geom_point(size = 2.5) +
    scale_y_continuous(limits = c(0, 1), labels = percent) +
    labs(
      title = "Compact summary: K-selection accuracy by sample size",
      x = "Sample size",
      y = "K-selection accuracy",
      color = "True K"
    ) +
    theme_gem(11)
  
  ggsave("fig_extra_summary_by_nK_v4.pdf", fig_extra, width = 8.5, height = 5)
  cat("Saved fig_extra_summary_by_nK_v4.pdf\n")
}

cat("\nAll v4 figures generated.\n")
