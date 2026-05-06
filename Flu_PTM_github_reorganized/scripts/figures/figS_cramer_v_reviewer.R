# Reviewer figure: Per-protein Cramér's V stratified by site rate class
# Teufel et al. (2025) Genome Biology and Evolution — reviewer revision
#
# Reads results/reviewer_analyses/cramer_v_per_protein.csv and
# cramer_v_per_site.csv produced by scripts/analysis/cramer_v_by_protein.py
#
# Produces:
#   fig_cramer_v_per_protein.pdf  — Cleveland dot plot, all / fast / slow sites
#   fig_cramer_v_distribution.pdf — boxplot of per-site V by protein,
#                                   fast/slow sites overlaid as points

library(dplyr)
library(tidyr)
library(ggplot2)
library(readr)
library(forcats)
library(here)

# ---------------------------------------------------------------------------
# Shared constants
# ---------------------------------------------------------------------------
C_ALL  <- "#0072B2"   # blue  (Okabe-Ito)
C_FAST <- "#D55E00"   # vermillion
C_SLOW <- "#009E73"   # green

PROT_ORDER <- c(
  "PB2_polymerase", "PB1_polymerase", "PB1-F2_protein",
  "PA_polymerase",  "PA-X_protein",   "NP_protein",
  "NA_neuraminidase","M1_matrix_protein","M2_ion_channel",
  "NS1_protein",    "HA_full"          # NEP excluded (only 2 sites, V inflated)
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
# Figure 1: Cleveland dot plot — protein-level mean Cramér's V
# ---------------------------------------------------------------------------
prot_df <- read_csv(here("results/reviewer_analyses/cramer_v_per_protein.csv"),
                    show_col_types = FALSE) |>
  filter(protein %in% PROT_ORDER)

# Pivot to long; keep NAs so proteins without fast/slow just don't appear
long_prot <- prot_df |>
  select(protein,
         all  = cramers_v_all,  se_all,
         fast = cramers_v_fast, se_fast,
         slow = cramers_v_slow, se_slow) |>
  pivot_longer(c(all, fast, slow),
               names_to  = "site_class",
               values_to = "cramers_v") |>
  mutate(
    se = case_when(
      site_class == "all"  ~ se_all,
      site_class == "fast" ~ se_fast,
      site_class == "slow" ~ se_slow
    ),
    site_class = factor(site_class,
                        levels = c("all", "fast", "slow"),
                        labels = c("All sites", "Fast-evolving", "Slow-evolving")),
    protein = factor(protein,
                     levels = PROT_ORDER[order(
                       prot_df$cramers_v_all[match(PROT_ORDER, prot_df$protein)])],
                     labels = PROT_LABELS[PROT_ORDER[order(
                       prot_df$cramers_v_all[match(PROT_ORDER, prot_df$protein)])]])
  ) |>
  filter(!is.na(cramers_v))

p1 <- ggplot(long_prot,
             aes(x = cramers_v, y = protein,
                 colour = site_class, shape = site_class)) +
  geom_vline(xintercept = 0.10, linetype = "dashed",
             colour = "grey65", linewidth = 0.5) +
  geom_errorbarh(aes(xmin = pmax(cramers_v - se, 0),
                     xmax = cramers_v + se),
                 height = 0.0, linewidth = 0.55, alpha = 0.7) +
  geom_point(size = 3.2, stroke = 0.8) +
  annotate("text", x = 0.103, y = 0.6, label = "V = 0.10",
           colour = "grey55", size = 3.2, hjust = 0, vjust = 0) +
  scale_colour_manual(
    values = c("All sites" = C_ALL, "Fast-evolving" = C_FAST, "Slow-evolving" = C_SLOW),
    name = NULL
  ) +
  scale_shape_manual(
    values = c("All sites" = 19, "Fast-evolving" = 17, "Slow-evolving" = 15),
    name = NULL
  ) +
  scale_x_continuous(limits = c(0, 0.38), expand = c(0.01, 0)) +
  labs(x = "Cramér's V  (host × PTM symbol)", y = NULL) +
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

ggsave(file.path(out_dir, "fig_cramer_v_per_protein.pdf"),
       p1, width = 6.5, height = 5.0, device = cairo_pdf)
message("Saved fig_cramer_v_per_protein.pdf")

# ---------------------------------------------------------------------------
# Figure 2: Boxplot of per-site Cramér's V by protein,
#           fast/slow sites overlaid as coloured points
# ---------------------------------------------------------------------------
site_df <- read_csv(here("results/reviewer_analyses/cramer_v_per_site.csv"),
                    show_col_types = FALSE) |>
  filter(protein %in% PROT_ORDER, !is.na(cramers_v)) |>
  mutate(
    protein    = factor(protein,
                        levels = PROT_ORDER[order(
                          prot_df$cramers_v_all[match(PROT_ORDER, prot_df$protein)])],
                        labels = PROT_LABELS[PROT_ORDER[order(
                          prot_df$cramers_v_all[match(PROT_ORDER, prot_df$protein)])]]),
    rate_class = case_when(
      is_fast ~ "Fast-evolving",
      is_slow ~ "Slow-evolving",
      TRUE    ~ "Neutral"
    )
  )

highlighted <- filter(site_df, rate_class != "Neutral")

p2 <- ggplot(site_df, aes(x = cramers_v, y = protein)) +
  geom_boxplot(outlier.shape = NA, fill = "grey93",
               colour = "grey45", linewidth = 0.45, width = 0.55) +
  geom_jitter(data = highlighted,
              aes(colour = rate_class),
              height = 0.18, width = 0, size = 2.5, alpha = 0.85) +
  scale_colour_manual(
    values = c("Fast-evolving" = C_FAST, "Slow-evolving" = C_SLOW),
    name = NULL
  ) +
  scale_x_continuous(limits = c(0, NA), expand = c(0.01, 0)) +
  labs(x = "Cramér's V  (host × PTM symbol, per site)", y = NULL) +
  theme_classic(base_size = 12) +
  theme(
    axis.line.y        = element_blank(),
    axis.ticks.y       = element_blank(),
    axis.text          = element_text(size = 11),
    axis.title.x       = element_text(size = 11, margin = margin(t = 8)),
    legend.position    = "top",
    legend.text        = element_text(size = 10),
    legend.key.size    = unit(0.5, "cm"),
    panel.grid.major.x = element_line(colour = "grey92", linewidth = 0.4),
    plot.margin        = margin(12, 18, 10, 10)
  )

ggsave(file.path(out_dir, "fig_cramer_v_distribution.pdf"),
       p2, width = 6.5, height = 5.0, device = cairo_pdf)
message("Saved fig_cramer_v_distribution.pdf")
