###############################################################
## GeM-Cox v5 Simulation Figures
##
## Inputs:
##   - GeMCox_simulation_v5_raw.rds
##   - GeMCox_simulation_v5_summary.csv
##
## Outputs:
##   - fig1_K_selection_v5.pdf
##   - fig2_ARI_all_v5.pdf
##   - fig3_ARI_conditional_v5.pdf
##   - fig4_cindex_gain_v5.pdf
##   - fig5_gamma_selected_v5.pdf
##   - fig6_dashboard_v5.pdf
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

raw_df  <- readRDS("GeMCox_simulation_v5_raw.rds")
summ_df <- read.csv("GeMCox_simulation_v5_summary.csv", stringsAsFactors = FALSE)

em <- raw_df$error_msg
ok <- is.na(em) | !nzchar(trimws(em)) | trimws(em) == "NA"
raw_df <- raw_df[ok, , drop = FALSE]

sep_levels <- c("low", "med", "high")
sep_colors <- c(low = "#e08060", med = "#5b9bd5", high = "#3a7a3a")
pat_colors <- c(signflip = "#7a3db8", disjoint = "#c0504d")

prep_fig <- function(df) {
  df |>
    mutate(
      separation = factor(separation, levels = sep_levels),
      n_label = factor(paste0("n = ", n), levels = paste0("n = ", c(117, 250, 500))),
      K_label = factor(paste0("K = ", K_true), levels = paste0("K = ", 1:3)),
      event_rate_label = factor(paste0("Event rate: ", sprintf("%.2f", event_rate_target %||% event_rate))),
      pattern_label = factor(paste0("Pattern: ", pattern))
    )
}
`%||%` <- function(x, y) if (!is.null(x)) x else y

summ <- prep_fig(summ_df)
raw  <- prep_fig(raw_df)

theme_gem <- function(base_size = 11) {
  theme_bw(base_size = base_size) +
    theme(
      panel.grid.minor = element_blank(),
      strip.background = element_rect(fill = "#1a3a5c", color = NA),
      strip.text = element_text(color = "white", face = "bold"),
      legend.position = "bottom",
      plot.title = element_text(face = "bold", size = 14),
      plot.subtitle = element_text(color = "grey35", size = 10)
    )
}

###############################################################
## Figure 1: K-selection accuracy by n, K_true, pattern
###############################################################
fig1 <- summ |>
  ggplot(aes(x = separation, y = K_sel_accuracy, fill = separation)) +
  geom_col(width = 0.68, alpha = 0.9) +
  geom_hline(yintercept = 1 / 3, linetype = "dashed", color = "grey55", linewidth = 0.4) +
  facet_grid(event_rate_label + pattern_label ~ n_label + K_label) +
  scale_fill_manual(values = sep_colors) +
  scale_y_continuous(limits = c(0, 1.02), labels = scales::percent) +
  labs(
    title = "K-selection accuracy",
    subtitle = "CV + screened 1-SE rule; dashed line = chance baseline for K grid {1,2,3}",
    x = "Separation",
    y = expression(P(hat(K) == K[true])),
    fill = "Separation"
  ) +
  theme_gem(10)

ggsave("fig1_K_selection_v5.pdf", fig1, width = 15, height = 8)

###############################################################
## Figure 2: Unconditional ARI
###############################################################
ari_all_data <- raw |>
  filter(K_true > 1, !is.na(ARI_all))

fig2 <- ari_all_data |>
  ggplot(aes(x = separation, y = ARI_all, fill = separation)) +
  geom_violin(trim = FALSE, scale = "width", alpha = 0.45, linewidth = 0.3) +
  geom_boxplot(width = 0.16, outlier.size = 0.7, fill = "white", alpha = 0.8) +
  facet_grid(event_rate_label + pattern_label ~ n_label + K_label) +
  scale_fill_manual(values = sep_colors) +
  scale_y_continuous(limits = c(-0.1, 1), breaks = seq(0, 1, 0.25)) +
  labs(
    title = "Cluster recovery: unconditional ARI",
    subtitle = "Computed for all successful fits when K_true > 1",
    x = "Separation", y = "Adjusted Rand Index", fill = "Separation"
  ) +
  theme_gem(10)

ggsave("fig2_ARI_all_v5.pdf", fig2, width = 15, height = 8)

###############################################################
## Figure 3: Conditional ARI given correct K selection
###############################################################
ari_cond_data <- raw |>
  filter(K_true > 1, ARI_conditional_defined == 1, !is.na(ARI_conditional))

fig3 <- ari_cond_data |>
  ggplot(aes(x = separation, y = ARI_conditional, fill = separation)) +
  geom_violin(trim = FALSE, scale = "width", alpha = 0.45, linewidth = 0.3) +
  geom_boxplot(width = 0.16, outlier.size = 0.7, fill = "white", alpha = 0.8) +
  facet_grid(event_rate_label + pattern_label ~ n_label + K_label) +
  scale_fill_manual(values = sep_colors) +
  scale_y_continuous(limits = c(-0.1, 1), breaks = seq(0, 1, 0.25)) +
  labs(
    title = "Cluster recovery: conditional ARI",
    subtitle = "Computed only when K_hat = K_true and K_true > 1",
    x = "Separation", y = "Adjusted Rand Index", fill = "Separation"
  ) +
  theme_gem(10)

ggsave("fig3_ARI_conditional_v5.pdf", fig3, width = 15, height = 8)

###############################################################
## Figure 4: C-index gain over single Cox
###############################################################
gain_data <- raw |>
  filter(!is.na(cindex_full), !is.na(cindex_cox1)) |>
  mutate(gain = cindex_full - cindex_cox1)

fig4 <- gain_data |>
  ggplot(aes(x = separation, y = gain, fill = separation)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey55") +
  geom_violin(trim = FALSE, scale = "width", alpha = 0.45, linewidth = 0.3) +
  geom_boxplot(width = 0.16, outlier.size = 0.7, fill = "white", alpha = 0.8) +
  facet_grid(event_rate_label + pattern_label ~ n_label + K_label) +
  scale_fill_manual(values = sep_colors) +
  labs(
    title = "C-index gain over single Cox",
    subtitle = "Positive values mean GeM-Cox outperformed the one-component Cox model",
    x = "Separation", y = expression(Delta ~ "C-index"), fill = "Separation"
  ) +
  theme_gem(10)

ggsave("fig4_cindex_gain_v5.pdf", fig4, width = 15, height = 8)

###############################################################
## Figure 5: Gamma selection by n and K_true
###############################################################
gamma_data <- raw |>
  filter(!is.na(gamma_selected)) |>
  mutate(gamma_label = factor(sprintf("%.2f", gamma_selected),
                              levels = sprintf("%.2f", sort(unique(gamma_selected)))))

fig5 <- gamma_data |>
  ggplot(aes(x = gamma_label, fill = separation)) +
  geom_bar(position = "dodge", alpha = 0.9) +
  facet_grid(event_rate_label + pattern_label ~ n_label + K_label) +
  scale_fill_manual(values = sep_colors) +
  labs(
    title = "Selected gamma distribution",
    subtitle = "Frequency of gamma values selected by screened CV",
    x = expression(gamma ~ "(survival weight)"), y = "Count", fill = "Separation"
  ) +
  theme_gem(10)

ggsave("fig5_gamma_selected_v5.pdf", fig5, width = 15, height = 8)

###############################################################
## Figure 6: Dashboard focused on sample-size effect
###############################################################
pa <- summ |>
  ggplot(aes(x = n_label, y = K_sel_accuracy, color = separation, group = separation)) +
  geom_line(linewidth = 0.7) +
  geom_point(size = 2.2) +
  facet_grid(pattern_label ~ K_label + event_rate_label) +
  scale_color_manual(values = sep_colors) +
  scale_y_continuous(limits = c(0, 1), labels = percent) +
  labs(title = "A. K-selection accuracy", x = NULL, y = expression(P(hat(K) == K[true])), color = "Separation") +
  theme_gem(9)

pb <- summ |>
  filter(K_true > 1) |>
  ggplot(aes(x = n_label, y = ARI_all_mean, color = separation, group = separation)) +
  geom_line(linewidth = 0.7) +
  geom_point(size = 2.2) +
  facet_grid(pattern_label ~ K_label + event_rate_label) +
  scale_color_manual(values = sep_colors) +
  scale_y_continuous(limits = c(-0.05, 1)) +
  labs(title = "B. Unconditional ARI", x = NULL, y = "Mean ARI", color = "Separation") +
  theme_gem(9)

pc <- summ |>
  ggplot(aes(x = n_label, y = cindex_gain_mean, color = separation, group = separation)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey60") +
  geom_line(linewidth = 0.7) +
  geom_point(size = 2.2) +
  facet_grid(pattern_label ~ K_label + event_rate_label) +
  scale_color_manual(values = sep_colors) +
  labs(title = "C. C-index gain over single Cox", x = "Sample size", y = expression(Delta ~ "C-index"), color = "Separation") +
  theme_gem(9)

fig6 <- pa / pb / pc +
  plot_annotation(
    title = "GeM-Cox v5 simulation dashboard",
    subtitle = "n sweep + K sweep + pattern sweep; screened CV; normalize_gmm_by_dim = FALSE; alpha = 0.5",
    theme = theme(plot.title = element_text(face = "bold", size = 15),
                  plot.subtitle = element_text(size = 11, color = "grey35"))
  )

ggsave("fig6_dashboard_v5.pdf", fig6, width = 16, height = 12)

cat("All v5 figures generated.\n")
