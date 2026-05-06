# Install packages if needed
# install.packages(c("tidyverse", "ggplot2", "ggseqlogo", "gridExtra"))
library(tidyverse)
library(ggplot2)
library(ggseqlogo)
library(gridExtra)
library(grid)  # Add grid library for textGrob
library(here)

# 1. Define your custom alphabet and colors
# Original alphabet (with problematic characters)
original_alphabet <- c("A", "C", "K", "M", "N", "O", "Q", "R", "S", "T", "Y", "U", "H", "L", "X", "-", "?")

# Map to alphanumeric equivalents for ggseqlogo
alphanumeric_alphabet <- c("A", "C", "K", "M", "N", "O", "Q", "R", "S", "T", "Y", "U", "H", "L", "X", "Z", "W")

# Create mapping between original and alphanumeric
symbol_mapping <- setNames(alphanumeric_alphabet, original_alphabet)
reverse_mapping <- setNames(original_alphabet, alphanumeric_alphabet)

# Original colors (using alphanumeric keys for ggseqlogo)
SYMBOL_COLOURS_ORIGINAL <- c(
  A = "#E69F00", C = "#1E90FF", K = "#56B4E9", M = "#009E73",
  N = "#F0E442", O = "#0072B2", Q = "#D55E00", R = "#CC79A7",
  S = "#984EA3", T = "#82B366", Y = "#8C564B", U = "#4DFFC3",
  H = "#FF6B6B", L = "#6A5ACD", X = "grey20",
  `-` = "darkgray", `?` = "lightgray"
)

# Colors mapped to alphanumeric equivalents
SYMBOL_COLOURS_MAPPED <- c(
  A = "#E69F00", C = "#1E90FF", K = "#56B4E9", M = "#009E73",
  N = "#F0E442", O = "#0072B2", Q = "#D55E00", R = "#CC79A7",
  S = "#984EA3", T = "#82B366", Y = "#8C564B", U = "#4DFFC3",
  H = "#FF6B6B", L = "#6A5ACD", X = "grey20",
  Z = "darkgray", W = "lightgray"  # Z = "-", W = "?"
)

# Create custom color scheme for ggseqlogo
custom_col_scheme <- make_col_scheme(
  chars = alphanumeric_alphabet,
  cols = SYMBOL_COLOURS_MAPPED[alphanumeric_alphabet]
)

# 2. Highlight rates from your table - now categorized by selection type
highlight_rates <- list(
  PB2_polymerase       = list(
    positive = c(`4` = 0.0469, `20` = 0.0480, `38` = 0.0449),
    purifying = c()
  ),
  `PB1-F2_protein`     = list(
    positive = c(),
    purifying = c(`7` = 0.0569, `8` = 0.0601, `12` = 0.0558)
  ),
  PB1_polymerase       = list(
    positive = c(),
    purifying = c(`1` = 0.0335, `8` = 0.0329, `13` = 0.0316, `14`= 0.0272, `23`= 0.0323)
  ),
  PA_polymerase        = list(
    positive = c(`13`= 0.0556, `25`= 0.0655),
    purifying = c()
  ),
  `PA-X_protein`       = list(
    positive = c(),
    purifying = c(`2` = 0.0690, `4` = 0.0696, `11` = 0.0699)
  ),
  HA_full           = list(
    positive = c(`36`= 0.0478, `77`=0.1602),
    purifying = c(`9` = 0.0093, `24`= 0.0089, `51`= 0.00946,`73`=0.009378)
  ),
  NP_protein           = list(
    positive = c(`13`= 0.0730, `15`= 0.0593),
    purifying = c()
  ),
  NA_neuraminidase     = list(
    positive = c(`33`= 0.1051,`53`= 0.1129),
    purifying = c()
  ),
  M1_matrix_protein    = list(
    positive = c(),
    purifying = c(`6` = 0.1063)
  ),
  M2_ion_channel       = list(
    positive = c(),
    purifying = c(`4` = 0.0779, `8` = 0.0747)
  ),
  NS1_protein          = list(
    positive = c(),
    purifying = c(`5` = 0.0303, `13`= 0.0332, `15`= 0.0350, `17`= 0.0341, `18`= 0.0351)
  )
)

# 3. Get all rates for any needed calculations (though we won't use gradient anymore)
all_positive_rates <- unlist(map(highlight_rates, "positive"))
all_purifying_rates <- unlist(map(highlight_rates, "purifying"))
all_rates <- c(all_positive_rates, all_purifying_rates)
rate_min <- min(all_rates)
rate_max <- max(all_rates)

# 4. Find nexus files only for proteins in highlight_rates list
nexus_files <- list.files(
  path       = file.path(normalizePath(file.path(here::here(), ".."), winslash = "/", mustWork = FALSE),
                        "cluster_output_final_converged", "processed_segments_fixed_final"),
  pattern    = "_multistate\\.txt$",
  recursive  = TRUE,
  full.names = TRUE
)

# Filter to only include files that match our highlight_rates keys
target_files <- c()
for (nex_file in nexus_files) {
  key <- sub("_multistate\\.txt$", "", basename(nex_file))
  if (key %in% names(highlight_rates)) {
    target_files <- c(target_files, nex_file)
  }
}

cat("Found", length(nexus_files), "total files\n")
if(length(target_files) == 0) {
  stop("No matching _multistate.txt files found for proteins in highlight_rates!")
}
cat("Processing", length(target_files), "files with highlight rates:\n")
cat(paste("  -", sub("_multistate\\.txt$", "", basename(target_files)), collapse = "\n"), "\n")

# Store individual plots
plot_list <- list()
seq_lengths <- c()

for (nex_file in target_files) {
  key   <- sub("_multistate\\.txt$", "", basename(nex_file))
  cat("Processing:", key, "\n")
  
  rates <- highlight_rates[[key]]  # This will always exist now
  
  # Parse MATRIX section
  lines <- readLines(nex_file)
  st    <- which(str_detect(lines, regex("^\\s*MATRIX", ignore_case=TRUE)))[1] + 1
  ends  <- which(str_trim(lines) == ";")
  en    <- ends[ends > st][1] - 1
  mat   <- str_trim(lines[st:en])
  mat   <- mat[mat != ""]
  
  # Extract sequences for ggseqlogo
  sequences <- map_chr(mat, function(line) {
    parts <- str_split(line, "\\s+")[[1]]
    original_seq <- paste(parts[-1], collapse = "")  # Remove taxon name, join sequence
    
    # Convert original symbols to alphanumeric equivalents
    mapped_seq <- original_seq
    for(i in seq_along(original_alphabet)) {
      mapped_seq <- gsub(paste0("\\", original_alphabet[i]), symbol_mapping[original_alphabet[i]], mapped_seq, fixed = TRUE)
    }
    
    return(mapped_seq)
  })
  seq_lengths[key] <- max(nchar(sequences))
  
  # Create individual sequence logo
  p <- ggseqlogo(
    sequences,  # Back to using original sequences with X's
    method = 'probability',
    namespace = alphanumeric_alphabet,
    col_scheme = custom_col_scheme
  ) +
    labs(
      title = key,
      x = "Position",
      y = "Probability"
    ) +
    theme_minimal(base_size = 9) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 10),
      axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size = 7),
      axis.text.y = element_text(size = 7),
      axis.title = element_text(size = 8),
      legend.position = "none"  # Remove individual legends
    )
  
  # Add selection pressure highlights - separate for positive and purifying
  selection_data <- highlight_rates[[key]]
  
  if (!is.null(selection_data)) {
    # Probability logos are bounded in [0, 1], so place markers slightly above 1.
    max_y <- 1
    
    # Add positive selection sites (red dots)
    if (length(selection_data$positive) > 0) {
      pos_df <- tibble(
        pos = as.integer(names(selection_data$positive)),
        rate = unname(selection_data$positive),
        type = "Positive"
      )
      
      p <- p + 
        geom_point(
          data = pos_df,
          aes(x = pos, y = max_y * 1.1),
          color = "red",
          size = 1.5,
          shape = 19,
          inherit.aes = FALSE
        )
    }
    
    # Add purifying selection sites (blue dots)
    if (length(selection_data$purifying) > 0) {
      pur_df <- tibble(
        pos = as.integer(names(selection_data$purifying)),
        rate = unname(selection_data$purifying),
        type = "Purifying"
      )
      
      p <- p + 
        geom_point(
          data = pur_df,
          aes(x = pos, y = max_y * 1.1),
          color = "blue",
          size = 1.5,
          shape = 19,
          inherit.aes = FALSE
        )
    }
  }
  
  plot_list[[key]] <- p
}

cat("Combining all plots with custom layout...\n")

# Ensure sequence lengths are numeric integers
seq_lengths <- round(as.numeric(seq_lengths))
names(seq_lengths) <- names(plot_list)

cat("Sequence lengths:\n")
for(i in seq_along(seq_lengths)) {
  cat(paste("  ", names(seq_lengths)[i], ":", seq_lengths[i], "positions\n"))
}

# Find the 3 longest sequences
longest_3 <- names(sort(seq_lengths, decreasing = TRUE)[1:3])
shorter_ones <- setdiff(names(plot_list), longest_3)

cat("Longest 3 (will span 2 columns):", paste(longest_3, collapse = ", "), "\n")
cat("Shorter ones (2 per row):", paste(shorter_ones, collapse = ", "), "\n")

# Create custom layout using grid
# We'll use grid.layout to create a custom arrangement

# Calculate number of rows needed
# 3 long plots (each takes 1 row) + ceiling(shorter_ones/2) rows for pairs
n_short_rows <- ceiling(length(shorter_ones) / 2)
total_rows <- 3 + n_short_rows

# Create layout matrix
layout_matrix <- matrix(0, nrow = total_rows, ncol = 2)

# First 3 rows: long plots spanning both columns
for(i in 1:3) {
  layout_matrix[i, ] <- i
}

# Remaining rows: short plots in pairs
plot_num <- 4
for(i in 4:total_rows) {
  if(plot_num <= length(plot_list)) {
    layout_matrix[i, 1] <- plot_num
    plot_num <- plot_num + 1
  }
  if(plot_num <= length(plot_list)) {
    layout_matrix[i, 2] <- plot_num
    plot_num <- plot_num + 1
  }
}

# Create the plot list in the right order (longest first, then shorter)
ordered_plots <- c(plot_list[longest_3], plot_list[shorter_ones])

# Use grid.arrange with custom layout
combined_plot <- grid.arrange(
  grobs = ordered_plots,
  layout_matrix = layout_matrix,
  heights = c(rep(1, nrow(layout_matrix)), 0.3)  # Add space for legend at bottom
)
dir.create(here("results/summary_figures"), showWarnings = FALSE, recursive = TRUE)
cowplot::save_plot(here("results/summary_figures/fig3_ptm_seqlogo.pdf"), combined_plot, ncol = 2, 
                   base_asp = 1.1,base_height	=8)

# Create color legend
color_legend_data <- data.frame(
  x = seq_along(alphanumeric_alphabet),
  y = rep(1, length(alphanumeric_alphabet)),
  symbol = alphanumeric_alphabet,
  original = c("A", "C", "K", "M", "N", "O", "Q", "R", "S", "T", "Y", "U", "H", "L", "X", "gap", "unknown")
)

color_legend <- ggplot(color_legend_data, aes(x = x, y = y)) +
  geom_tile(aes(fill = symbol), color = "white", size = 0.5) +
  geom_text(aes(label = symbol), color = "white", fontface = "bold", size = 3) +
  scale_fill_manual(values = SYMBOL_COLOURS_MAPPED, guide = "none") +
  labs(title = "Symbol Colors") +
  theme_void() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 10, face = "bold"),
    plot.margin = margin(5, 5, 5, 5)
  ) +
  coord_fixed(ratio = 1)

# Create mutation rate legend - now shows selection types
rate_legend_data <- data.frame(
  x = c(1, 2),
  y = c(1, 1),
  type = c("Positive Selection", "Purifying Selection"),
  color = c("red", "blue")
)

rate_legend <- ggplot(rate_legend_data, aes(x = x, y = y, color = type)) +
  geom_point(size = 4, shape = 19) +
  scale_color_manual(
    name = "Selection Type",
    values = c("Positive Selection" = "red", "Purifying Selection" = "blue")
  ) +
  labs(title = "Selection Pressure") +
  theme_void() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 10, face = "bold"),
    legend.position = "bottom",
    legend.title = element_text(size = 9),
    legend.text = element_text(size = 8),
    plot.margin = margin(5, 5, 5, 5)
  ) +
  guides(color = guide_legend(title.position = "top", title.hjust = 0.5))

# Combine legends
legend_combined <- grid.arrange(color_legend, rate_legend, ncol = 2, widths = c(2, 1))

# Add explanation text
grid.text("Red dots = positive selection sites, Blue dots = purifying selection sites", 
          x = 0.5, y = 0.02, 
          gp = gpar(fontsize = 9, fontface = "italic"))

print("Combined plot created!")