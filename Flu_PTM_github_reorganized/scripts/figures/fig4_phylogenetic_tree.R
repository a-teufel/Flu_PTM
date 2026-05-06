# Figures 4, S1, S3, S5, S7: Circular phylogenetic trees with host icons
# Teufel et al. (2025) Genome Biology and Evolution
#
# Reads the consensus NEXUS tree for a given protein, colours branches
# by posterior mean branch rate (viridis plasma scale), annotates tip
# labels with host emoji and strain identity, and adds node pie charts
# for host-state proportions from ML ancestral reconstruction.
#
# Usage: set PROTEIN at the top, then source the script. Run once per
# protein to generate each supplemental tree figure.

library(ape)
library(dplyr)
library(readr)
library(stringr)
library(jsonlite)
library(viridis)
library(here)

# ---------------------------------------------------------------------------
# Configuration — change PROTEIN to generate a different figure
# ---------------------------------------------------------------------------
PROTEIN      <- "NP_protein"          # e.g. PB2_polymerase, HA_full, NA_neuraminidase
BURNIN       <- 0.25                   # fraction of MCMC samples to discard

tree_file    <- here("data/trees", paste0(PROTEIN, "_tree.nex"))
branch_log   <- here("results/host_shift", paste0(PROTEIN, "_branch_rates.log"))
mapping_file <- here("data/ptm/mappings", paste0(PROTEIN, "_id_mapping.json"))
meta_files   <- list(
  H1N1 = here("data/raw_sequences/H1N1/H1N1_seq_info.csv"),
  H5N1 = here("data/raw_sequences/H5N1/H5N1_seq_info.csv"),
  H7N9 = here("data/raw_sequences/H7N9/H7N9_seq_info.csv")
)
out_file <- here("results/summary_figures", paste0("fig_tree_", PROTEIN, ".pdf"))

# ---------------------------------------------------------------------------
# Host → emoji mapping
# Species names used as partial-match patterns (case-insensitive)
# ---------------------------------------------------------------------------
HOST_ICONS <- c(
  "Sus scrofa"         = "\U1F416",
  "pig"                = "\U1F416",
  "Homo sapiens"       = "\U1F9D1",
  "human"              = "\U1F9D1",
  "Gallus gallus"      = "\U1F414",
  "chicken"            = "\U1F414",
  "Anatidae"           = "\U1F986",
  "Anas"               = "\U1F986",
  "duck"               = "\U1F986",
  "Spatula"            = "\U1F986",
  "Mareca"             = "\U1F986",
  "Cygnus"             = "\U1F9A2",
  "swan"               = "\U1F9A2",
  "Meleagris"          = "\U1F983",
  "turkey"             = "\U1F983",
  "Numididae"          = "\U1F983",
  "Accipitriformes"    = "\U1F985",
  "Psittacidae"        = "\U1F99C",
  "Columbidae"         = "\U1F54A",
  "Phasianidae"        = "\U1F99A",
  "Phasianinae"        = "\U1F99A",
  "Mustela"            = "\U1F98A",
  "Panthera"           = "\U1F405"
)
DEFAULT_ICON <- "?"

STRAIN_COLOURS <- c(H1N1 = "#E41A1C", H5N1 = "#377EB8", H7N9 = "#4DAF4A", Unknown = "grey50")

# ---------------------------------------------------------------------------
# Helper: get host icon from species string
# ---------------------------------------------------------------------------
get_icon <- function(host) {
  if (is.na(host)) return(DEFAULT_ICON)
  match_idx <- which(str_detect(host, regex(names(HOST_ICONS), ignore_case = TRUE)))[1]
  if (is.na(match_idx)) DEFAULT_ICON else HOST_ICONS[[match_idx]]
}

# ---------------------------------------------------------------------------
# Load tree
# ---------------------------------------------------------------------------
raw_trees <- tryCatch(read.nexus(tree_file), error = function(e) read.tree(tree_file))

if (inherits(raw_trees, "multiPhylo")) {
  keep_idx  <- seq(floor(length(raw_trees) * BURNIN) + 1, length(raw_trees))
  cons_tree <- consensus(raw_trees[keep_idx])
} else {
  cons_tree <- raw_trees
}

# Ensure valid branch lengths
if (is.null(cons_tree$edge.length) || any(is.na(cons_tree$edge.length))) {
  cons_tree$edge.length <- rep(1e-3, nrow(cons_tree$edge))
}

# ---------------------------------------------------------------------------
# Load branch rates
# ---------------------------------------------------------------------------
n_edges    <- nrow(cons_tree$edge)
rate_vec   <- rep(1, n_edges)

if (file.exists(branch_log)) {
  log_tbl    <- read.table(branch_log, header = TRUE)
  burnin_row <- floor(nrow(log_tbl) * BURNIN)
  log_post   <- log_tbl[(burnin_row + 1):nrow(log_tbl), ]
  rate_cols  <- grep("^full_branch_rates\\.", colnames(log_post), value = TRUE)
  if (length(rate_cols) > 0) {
    means   <- colMeans(log_post[, rate_cols])
    use_n   <- min(length(means), n_edges)
    rate_vec[seq_len(use_n)] <- means[seq_len(use_n)]
  }
}
rate_vec[is.na(rate_vec)] <- mean(rate_vec, na.rm = TRUE)

# Colour edges by rate; near-zero branches drawn in black
near_zero     <- rate_vec < quantile(rate_vec, 0.05)
plasma_cols   <- plasma(10)
rate_cuts     <- cut(rate_vec, 10, include.lowest = TRUE)
edge_col      <- plasma_cols[rate_cuts]
edge_col[near_zero] <- "black"
edge_width          <- ifelse(near_zero, 1.5, 3)

# ---------------------------------------------------------------------------
# Load metadata and resolve tip labels
# ---------------------------------------------------------------------------
meta <- imap_dfr(meta_files, function(path, strain) {
  read_csv(path, show_col_types = FALSE) |>
    mutate(Strain = strain,
           Accession = str_trim(Accession) |> str_remove("\\.1$"))
})

json_map <- if (file.exists(mapping_file)) fromJSON(mapping_file) else list()

resolve_accession <- function(tip) {
  if (!is.null(json_map[[tip]])) return(json_map[[tip]])
  str_remove(tip, "^H\\dN\\d__") |> str_remove("_$")
}

orig_tips <- cons_tree$tip.label
short_tips <- map_chr(orig_tips, resolve_accession)

tip_meta <- tibble(acc = short_tips) |>
  left_join(meta, by = c("acc" = "Accession")) |>
  mutate(
    Strain = replace_na(Strain, "Unknown"),
    Host   = replace_na(Host,   "Unknown"),
    icon   = map_chr(Host, get_icon)
  )

cons_tree$tip.label <- paste0(
  tip_meta$icon, " ",
  tip_meta$Strain, " | ",
  tip_meta$Host
)

tip_colors <- STRAIN_COLOURS[tip_meta$Strain]

# ---------------------------------------------------------------------------
# Node pie charts: host proportions from ML ancestral reconstruction
# ---------------------------------------------------------------------------
host_factor  <- factor(tip_meta$Host, levels = unique(tip_meta$Host))
names(host_factor) <- short_tips

ace_out <- tryCatch(
  ace(host_factor, cons_tree, type = "discrete"),
  error = function(e) NULL
)

n_tips  <- length(cons_tree$tip.label)
n_nodes <- cons_tree$Nnode
all_hosts <- levels(host_factor)

if (!is.null(ace_out)) {
  anc_states <- ace_out$lik.anc  # n_nodes × n_states matrix
} else {
  # Fallback: assign all weight to most common host
  dominant <- names(sort(table(host_factor), decreasing = TRUE))[1]
  anc_states <- matrix(0, nrow = n_nodes, ncol = length(all_hosts),
                       dimnames = list(NULL, all_hosts))
  anc_states[, dominant] <- 1
}

# Build host colour vector for pie charts
# Birds are all shown in shades of red; mammals in distinct colours
host_col <- rep("grey80", length(all_hosts))
names(host_col) <- all_hosts

BIRD_PATTERNS <- c("Gallus","Anas","Anatidae","Aves","bird","duck","swan",
                   "Spatula","Mareca","Meleagris","Numididae","Cygnus",
                   "Accipitriformes","Psittacidae","Columbidae","Phasian")
for (h in all_hosts) {
  if (any(str_detect(h, regex(BIRD_PATTERNS, ignore_case = TRUE)))) {
    host_col[h] <- "#CC4444"
  } else if (str_detect(h, regex("Homo sapiens|human", ignore_case = TRUE))) {
    host_col[h] <- "#7570B3"
  } else if (str_detect(h, regex("Sus scrofa|pig|swine", ignore_case = TRUE))) {
    host_col[h] <- "#E7298A"
  } else if (h == "Unknown") {
    host_col[h] <- "grey40"
  }
}

# Reorder columns of anc_states to match host_col order
if (!is.null(colnames(anc_states))) {
  shared <- intersect(colnames(anc_states), names(host_col))
  anc_states <- anc_states[, shared, drop = FALSE]
  pie_col    <- host_col[shared]
} else {
  pie_col <- head(host_col, ncol(anc_states))
}

# ---------------------------------------------------------------------------
# Plot
# ---------------------------------------------------------------------------
n_tips_total <- length(cons_tree$tip.label)
fig_height   <- max(8, n_tips_total * 0.12)
pie_cex      <- max(0.08, min(0.20, 4 / log(n_nodes + 1)))

pdf(out_file, width = 16, height = fig_height)
par(mar = c(1, 1, 2, 12), xpd = TRUE)

plot.phylo(
  cons_tree,
  type        = "fan",
  cex         = 0.55,
  edge.width  = edge_width,
  edge.color  = edge_col,
  tip.color   = tip_colors,
  no.margin   = FALSE,
  x.lim       = c(-max(nodeHeights(cons_tree)) * 1.05,
                   max(nodeHeights(cons_tree)) * 1.05)
)

nodelabels(pie = anc_states, piecol = pie_col, cex = pie_cex)

title(main = str_replace_all(PROTEIN, "_", " "), cex.main = 1)

# Branch rate legend
usr <- par("usr")
lx  <- usr[2] * 1.02
ly  <- seq(usr[4], usr[3], length.out = 12)
legend(lx, ly[1],
       legend = c("Low (near zero)", paste0("Level ", 1:10)),
       col    = c("black", plasma_cols),
       lwd    = c(1.5, rep(3, 10)),
       cex    = 0.65, bty = "n", title = "Branch rate")

dev.off()
message("Saved: ", out_file)
