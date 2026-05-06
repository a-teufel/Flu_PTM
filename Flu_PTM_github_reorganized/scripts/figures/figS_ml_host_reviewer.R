# Reviewer figure: ML host prediction from PTM states
# Teufel et al. (2025) Genome Biology and Evolution — reviewer revision
#
# Reads results/reviewer_analyses/ml_results_summary.csv and
# ml_feature_importance.csv produced by scripts/analysis/ml_host_prediction.py
#
# Produces:
#   fig_ml_balanced_accuracy.pdf — sorted dot plot, proteins × feature set (RF)
#   fig_ml_rf_importance.pdf     — bar chart of RF site importances

library(dplyr)
library(tidyr)
library(ggplot2)
library(readr)
library(forcats)
library(here)

C_ALL  <- "#0072B2"
C_FAST <- "#D55E00"
C_SLOW <- "#009E73"

PROT_ORDER <- c(
  "PB2_polymerase", "PB1_polymerase", "PB1-F2_protein",
  "PA_polymerase",  "PA-X_protein",   "NP_protein",
  "NA_neuraminidase","M1_matrix_protein","M2_ion_channel",
  "NS1_protein",    "HA_full"
)
PROT_LABELS <- c(
  PB2_polymerase   = "PB2",    PB1_polymerase    = "PB1",
  `PB1-F2_protein` = "PB1-F2", PA_polymerase     = "PA",
  `PA-X_protein`   = "PA-X",   NP_protein        = "NP",
  NA_neuraminidase = "NA",     M1_matrix_protein = "M1",
  M2_ion_channel   = "M2",     NS1_protein       = "NS1",
  HA_full          = "HA"
)

out_dir <- here("results/reviewer_analyses")

# ---------------------------------------------------------------------------
# Figure 1: dot plot — balanced accuracy, RF only, sorted by all-sites perf
# ---------------------------------------------------------------------------
ml_df <- read_csv(here("results/reviewer_analyses/ml_results_summary.csv"),
                  show_col_types = FALSE) |>
  filter(model == "RandomForest", protein %in% PROT_ORDER) |>
  mutate(feature_set = factor(feature_set,
                              levels = c("all", "fast", "slow"),
                              labels = c("All sites", "Fast sites", "Slow sites")))

# Order proteins by all-sites balanced accuracy (high to low)
prot_order_ml <- ml_df |>
  filter(feature_set == "All sites") |>
  arrange(balanced_accuracy) |>
  pull(protein)

ml_df <- ml_df |>
  mutate(protein = factor(protein,
                          levels = prot_order_ml,
                          labels = PROT_LABELS[prot_order_ml]))

p1 <- ggplot(ml_df,
             aes(x = balanced_accuracy, y = protein,
                 colour = feature_set, shape = feature_set)) +
  geom_vline(xintercept = 1/3, linetype = "dotted",
             colour = "grey60", linewidth = 0.55) +
  geom_point(size = 3.5, stroke = 0.8, alpha = 0.9) +
  annotate("text", x = 1/3 + 0.012, y = 0.6,
           label = "Chance", colour = "grey55",
           size = 3.2, hjust = 0, vjust = 0) +
  scale_colour_manual(
    values = c("All sites" = C_ALL, "Fast sites" = C_FAST, "Slow sites" = C_SLOW),
    name = NULL
  ) +
  scale_shape_manual(
    values = c("All sites" = 19, "Fast sites" = 17, "Slow sites" = 15),
    name = NULL
  ) +
  scale_x_continuous(limits = c(0.20, 1.00),
                     breaks = seq(0.2, 1.0, 0.2),
                     labels = scales::percent_format(accuracy = 1)) +
  labs(x = "Balanced accuracy  (5-fold CV, Random Forest)",
       y = NULL) +
  theme_classic(base_size = 12) +
  theme(
    axis.line.y        = element_blank(),
    axis.ticks.y       = element_blank(),
    axis.text          = element_text(size = 11),
    axis.title.x       = element_text(size = 11, margin = margin(t = 8)),
    legend.position    = "top",
    legend.text        = element_text(size = 10),
    legend.key.size    = unit(0.55, "cm"),
    legend.spacing.x   = unit(0.4, "cm"),
    panel.grid.major.x = element_line(colour = "grey92", linewidth = 0.4),
    panel.grid.major.y = element_line(colour = "grey95", linewidth = 0.3),
    plot.margin        = margin(12, 18, 10, 10)
  )

ggsave(file.path(out_dir, "fig_ml_balanced_accuracy.pdf"),
       p1, width = 6.0, height = 5.0, device = cairo_pdf)
message("Saved fig_ml_balanced_accuracy.pdf")

# ---------------------------------------------------------------------------
# Figure 2: RF site importance — proteins with fast sites, faceted
# ---------------------------------------------------------------------------
imp_df <- read_csv(here("results/reviewer_analyses/ml_feature_importance.csv"),
                   show_col_types = FALSE)

prots_fast <- imp_df |>
  filter(is_fast, protein %in% PROT_ORDER) |>
  distinct(protein) |>
  pull(protein)
prots_fast <- intersect(PROT_ORDER, prots_fast)

imp_plot <- imp_df |>
  filter(protein %in% prots_fast) |>
  mutate(
    protein    = factor(protein, levels = prots_fast,
                        labels = PROT_LABELS[prots_fast]),
    rate_class = case_when(
      is_fast ~ "Fast-evolving",
      is_slow ~ "Slow-evolving",
      TRUE    ~ "Neutral"
    )
  ) |>
  group_by(protein) |>
  mutate(rank = rank(orig_pos, ties.method = "first")) |>
  ungroup()

top_labels <- imp_plot |>
  filter(rate_class == "Fast-evolving") |>
  group_by(protein) |>
  slice_max(importance, n = 1) |>
  ungroup()

p2 <- ggplot(imp_plot, aes(x = rank, y = importance, fill = rate_class)) +
  geom_col(width = 0.9, alpha = 0.85) +
  geom_text(data = top_labels,
            aes(label = orig_pos), colour = C_FAST,
            vjust = -0.4, size = 3.2, fontface = "bold") +
  facet_wrap(~protein, scales = "free_x", nrow = 1) +
  scale_fill_manual(
    values = c("Fast-evolving" = C_FAST, "Slow-evolving" = C_SLOW,
               "Neutral" = "grey75"),
    name = NULL
  ) +
  scale_x_continuous(breaks = NULL) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
  labs(x = "PTM site (ranked by alignment position)",
       y = "Random Forest importance") +
  theme_classic(base_size = 12) +
  theme(
    strip.background   = element_blank(),
    strip.text         = element_text(size = 11, face = "bold"),
    legend.position    = "top",
    legend.text        = element_text(size = 10),
    legend.key.size    = unit(0.5, "cm"),
    axis.ticks.x       = element_blank(),
    axis.text.x        = element_blank(),
    axis.text.y        = element_text(size = 9),
    axis.title         = element_text(size = 11),
    panel.grid.major.y = element_line(colour = "grey92", linewidth = 0.4),
    panel.spacing      = unit(1.2, "lines"),
    plot.margin        = margin(12, 12, 10, 10)
  )

ggsave(file.path(out_dir, "fig_ml_rf_importance.pdf"),
       p2, width = 2.8 * length(prots_fast), height = 4.2,
       device = cairo_pdf)
message("Saved fig_ml_rf_importance.pdf")
