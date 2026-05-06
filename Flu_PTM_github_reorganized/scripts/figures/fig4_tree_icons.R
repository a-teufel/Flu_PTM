# Complete R Script for Phylogenetic Tree Visualization with Aggressive Branch Collapsing
# Author: Modified to collapse clades into tiny triangles

# Load required libraries
library(ape)
library(ggtree)
library(ggplot2)
library(dplyr)
library(stringr)
library(viridis)
library(jsonlite)
library(readr)
library(tidyr)
library(purrr)
library(here)

# Install and load phangorn for tree manipulation if needed
if (!requireNamespace("phangorn", quietly = TRUE)) {
  message("Installing phangorn package for tree manipulation...")
  install.packages("phangorn")
}
library(phangorn)

# Define host-to-emoji mapping
species_symbols <- list(
  "Sus scrofa"           = "\U1F416",
  "Homo sapiens"         = "\U1F9D1",
  "Mus musculus"         = "\U1F42D",
  "Gallus gallus"        = "\U1F414",
  "Phasianinae"          = "\U1F99A",
  "Phasianidae"          = "\U1F99A",
  "Panthera tigris"      = "\U1F405",
  "Mustela lutreola"     = "\U1F43A",
  "Anatidae"             = "\U1F986",
  "Passer montanus"      = "\U1F426",
  "Anas platyrhynchos"   = "\U1F986",
  "Arenaria interpres"   = "\U1F426",
  "Accipitriformes"      = "\U1F985",
  "Aves"                 = "\U1F426",
  "Psittacidae"          = "\U1F99C",
  "Columbidae"           = "\U1F54A",
  "Suricata suricatta"   = "\U1F9A1",
  "Mareca americana"     = "\U1F986",
  "Anas carolinensis"    = "\U1F986",
  "Numididae"            = "\U1F983",
  "Meleagris gallopavo"  = "\U1F983",
  "Spatula clypeata"     = "\U1F986",
  "Numididae sp."        = "\U1F983",
  "Spatula discors"      = "\U1F986",
  "Anas cyanoptera"      = "\U1F986",
  "Cygnus olor"          = "\U1F9A2",
  "Cairina moschata"     = "\U1F986",
  "Mustela putorius furo"= "\U1F98A",
  "Felis catus"          = "\U1F408",
  "Tachybaptus ruficollis"= "\U1F986",
  "Hirundo rustica"      = "\U1F426",
  "Accipitridae"         = "\U1F985",
  "Tadorna"              = "\U1F986",
  "Parus major"          = "\U1F426",
  "Aythya ferina"        = "\U1F986"
)
default_symbol <- "\U003F"  # question mark for anything missing

# Define file paths for CSV files (UPDATE THESE PATHS TO YOUR FILES)
csv_files <- list(
  h1n1 = here("data/raw_sequences/H1N1/H1N1_seq_info.csv"),
  h5n1 = here("data/raw_sequences/H5N1/H5N1_seq_info.csv"),
  h7n9 = here("data/raw_sequences/H7N9/H7N9_seq_info.csv")
)

# Function to read PTM metadata file
read_ptm_metadata <- function(file_path) {
  lines <- readLines(file_path)
  data_lines <- lines[!grepl("^#", lines)]
  
  ptm_data <- data.frame(
    Alignment_Position = integer(),
    PTM_Type = character(),
    Sequences = character(),
    Consensus_Position = integer(),
    Consensus_AA = character(),
    stringsAsFactors = FALSE
  )
  
  for (line in data_lines) {
    if (nchar(line) > 0) {
      parts <- unlist(strsplit(line, "\t"))
      if (length(parts) >= 5) {
        ptm_data <- rbind(ptm_data, data.frame(
          Alignment_Position = as.integer(parts[1]),
          PTM_Type = parts[2],
          Sequences = parts[3],
          Consensus_Position = as.integer(parts[4]),
          Consensus_AA = parts[5],
          stringsAsFactors = FALSE
        ))
      }
    }
  }
  
  return(ptm_data)
}

# Strain matching function
match_strains_to_meta <- function(tip_labels, meta_data, mapping) {
  result <- data.frame(
    label = tip_labels,
    Final_Accession = NA_character_,
    Matched_Strain = NA_character_,
    match_method = NA_character_,
    stringsAsFactors = FALSE
  )
  
  # STEP 1: Direct mapping via JSON
  if (length(mapping) > 0) {
    for (i in seq_len(nrow(result))) {
      curr_label <- result$label[i]
      
      # Case 1: Tip label is a value in the mapping
      if (curr_label %in% mapping) {
        matching_keys <- names(mapping)[mapping == curr_label]
        if (length(matching_keys) > 0) {
          for (key in matching_keys) {
            if (grepl("__([A-Z]{2,}\\d{5,})_", key)) {
              potential_acc <- sub(".*__([A-Z]{2,}\\d{5,})_.*", "\\1", key)
              idx <- match(potential_acc, meta_data$Accession)
              if (!is.na(idx)) {
                result$Final_Accession[i] <- potential_acc
                result$Matched_Strain[i] <- meta_data$Strain[idx]
                result$match_method[i] <- "json_value_reverse"
                break
              }
            }
          }
          if (is.na(result$Final_Accession[i])) {
            result$Final_Accession[i] <- curr_label
            idx <- match(curr_label, meta_data$Accession)
            if (!is.na(idx)) {
              result$Matched_Strain[i] <- meta_data$Strain[idx]
              result$match_method[i] <- "json_value_direct"
            }
          }
        }
      }
      
      # Case 2: Tip label is a key in the mapping
      if (is.na(result$Final_Accession[i]) && curr_label %in% names(mapping)) {
        mapped_value <- mapping[[curr_label]]
        if (grepl("^H[1-9]N[1-9]__", mapped_value)) {
          clean_value <- gsub("^(H[1-9]N[1-9]__)", "", mapped_value)
          clean_value <- gsub("_$", "", clean_value)
          idx <- match(clean_value, meta_data$Accession)
          if (!is.na(idx)) {
            result$Final_Accession[i] <- clean_value
            result$Matched_Strain[i] <- meta_data$Strain[idx]
            result$match_method[i] <- "json_key_h1n1_pattern"
          } else {
            result$Final_Accession[i] <- clean_value
            result$match_method[i] <- "json_key_h1n1_pattern_no_meta"
          }
        } else {
          idx <- match(mapped_value, meta_data$Accession)
          if (!is.na(idx)) {
            result$Final_Accession[i] <- mapped_value
            result$Matched_Strain[i] <- meta_data$Strain[idx]
            result$match_method[i] <- "json_key_direct"
          } else {
            result$Final_Accession[i] <- mapped_value
            result$match_method[i] <- "json_key_direct_no_meta"
          }
        }
      }
    }
  }
  
  # Step 2: Direct accession matching
  rem <- which(is.na(result$Final_Accession))
  if (length(rem) > 0) {
    acc_pat <- "^[A-Z]{2}[_]?[0-9]{6,}$"
    hit <- which(is.na(result$Final_Accession) &
                   grepl(acc_pat, result$label) &
                   result$label %in% meta_data$Accession)
    if (length(hit) > 0) {
      result$Final_Accession[hit] <- result$label[hit]
      result$Matched_Strain[hit] <- meta_data$Strain[match(result$label[hit], meta_data$Accession)]
      result$match_method[hit] <- "direct_accession"
    }
  }
  
  # Step 3: HxNy__ACC_ pattern
  rem <- which(is.na(result$Final_Accession))
  if (length(rem) > 0) {
    for (i in rem) {
      if (grepl("^H\\d+N\\d+__([A-Z]{2,}\\d{5,})_$", result$label[i])) {
        acc <- sub("^H\\d+N\\d+__([A-Z]{2,}\\d{5,})_$", "\\1", result$label[i])
        idx <- match(acc, meta_data$Accession)
        if (!is.na(idx)) {
          result$Final_Accession[i] <- acc
          result$Matched_Strain[i] <- meta_data$Strain[idx]
          result$match_method[i] <- "HxNy_pattern"
        }
      }
    }
  }
  
  # Step 4: ARG via JSON mapping
  rem <- which(is.na(result$Final_Accession))
  if (length(rem) > 0 && length(mapping) > 0) {
    for (i in rem) {
      if (grepl("^A[0-9]{8,}$", result$label[i])) {
        key <- names(mapping)[mapping == result$label[i]]
        if (length(key) && grepl("__([A-Z]{2,}\\d{5,})_\\d+_", key)) {
          acc <- sub(".*__([A-Z]{2,}\\d{5,})_\\d+_.*", "\\1", key)
          idx <- match(acc, meta_data$Accession)
          if (!is.na(idx)) {
            result$Final_Accession[i] <- acc
            result$Matched_Strain[i] <- meta_data$Strain[idx]
            result$match_method[i] <- "ARG_pattern"
          }
        }
      }
    }
  }
  
  # Step 5: No underscore suffix
  rem <- which(is.na(result$Final_Accession))
  if (length(rem) > 0) {
    label_no_underscore <- gsub("_$", "", result$label[rem])
    meta_no_underscore <- gsub("_$", "", meta_data$Accession)
    for (i in seq_along(rem)) {
      idx <- match(label_no_underscore[i], meta_no_underscore)
      if (!is.na(idx)) {
        result$Final_Accession[rem[i]] <- meta_data$Accession[idx]
        result$Matched_Strain[rem[i]] <- meta_data$Strain[idx]
        result$match_method[rem[i]] <- "no_underscore"
      }
    }
  }
  
  return(result)
}

# Function to get comprehensive information for sequences (strain, year, country, host)
get_comprehensive_info <- function(tip_labels, json_file) {
  h1 <- read_csv(csv_files$h1n1, show_col_types = FALSE) %>%
    mutate(Strain = "H1N1")
  h5 <- read_csv(csv_files$h5n1, show_col_types = FALSE) %>%
    mutate(Strain = "H5N1")
  h7 <- read_csv(csv_files$h7n9, show_col_types = FALSE) %>%
    mutate(Strain = "H7N9")
  
  all_meta <- bind_rows(h1, h5, h7)
  
  mapping <- list()
  if (file.exists(json_file)) {
    mapping <- fromJSON(json_file)
  }
  
  strain_matching <- match_strains_to_meta(tip_labels, all_meta, mapping)
  
  comprehensive_info <- data.frame(
    tip_label = tip_labels,
    Host = NA_character_,
    Country = NA_character_,
    Year = NA_character_,
    Strain_Type = NA_character_,
    stringsAsFactors = FALSE
  )
  
  for (i in seq_len(nrow(comprehensive_info))) {
    acc <- strain_matching$Final_Accession[i]
    if (!is.na(acc)) {
      meta_idx <- match(acc, all_meta$Accession)
      if (!is.na(meta_idx)) {
        comprehensive_info$Host[i] <- all_meta$Host[meta_idx]
        comprehensive_info$Country[i] <- all_meta$Country[meta_idx]
        comprehensive_info$Strain_Type[i] <- all_meta$Strain[meta_idx]
        
        # Extract year from Collection_Date
        if (!is.na(all_meta$Collection_Date[meta_idx])) {
          year_match <- str_extract(all_meta$Collection_Date[meta_idx], "\\d{4}")
          if (!is.na(year_match)) {
            comprehensive_info$Year[i] <- year_match
          }
        }
      }
    }
  }
  
  return(comprehensive_info)
}

# Function to read RevBayes output and remove burn-in
read_revbayes_rates <- function(file_path, burnin_fraction = 0.25) {
  data <- read.table(file_path, header = TRUE, sep = "\t")
  burnin_rows <- ceiling(nrow(data) * burnin_fraction)
  data_post_burnin <- data[(burnin_rows + 1):nrow(data), ]
  rate_cols <- grep("full_branch_rates", colnames(data_post_burnin))
  branch_rates <- data_post_burnin[, rate_cols]
  mean_rates <- apply(branch_rates, 2, mean)
  return(mean_rates)
}

# Function to extract sequences for a specific position
get_sequences_for_position <- function(ptm_data, position) {
  row_data <- ptm_data[ptm_data$Alignment_Position == position, ]
  if (nrow(row_data) == 0) {
    return(character(0))
  }
  sequences <- unlist(strsplit(row_data$Sequences, ","))
  sequences <- trimws(sequences)
  return(sequences)
}

# Function to get sequences for multiple positions and create label mapping
get_sequences_for_multiple_positions <- function(ptm_data, positions) {
  all_sequences <- character(0)
  sequence_to_positions <- list()
  
  for (position in positions) {
    sequences <- get_sequences_for_position(ptm_data, position)
    all_sequences <- c(all_sequences, sequences)
    
    for (seq in sequences) {
      if (seq %in% names(sequence_to_positions)) {
        sequence_to_positions[[seq]] <- c(sequence_to_positions[[seq]], position)
      } else {
        sequence_to_positions[[seq]] <- position
      }
    }
  }
  
  return(list(
    sequences = unique(all_sequences),
    sequence_to_positions = sequence_to_positions
  ))
}

# Function to create comprehensive tip labels with strain, year, country, positions, and host emojis
create_comprehensive_tip_labels <- function(tip_labels, sequence_to_positions, comprehensive_info) {
  modified_labels <- tip_labels
  
  for (i in seq_along(tip_labels)) {
    label <- tip_labels[i]
    
    # Clean up the label - remove H1N1__, H5N1__, H7N9__ prefixes
    clean_label <- label
    if (grepl("^H[1-9]N[1-9]__", label)) {
      clean_label <- sub("^H[1-9]N[1-9]__", "", label)
    }
    
    # Start with cleaned label
    final_label <- clean_label
    
    # Add metadata if available
    strain <- comprehensive_info$Strain_Type[i]
    year <- comprehensive_info$Year[i]
    country <- comprehensive_info$Country[i]
    host <- comprehensive_info$Host[i]
    
    # Build metadata string
    metadata_parts <- c()
    if (!is.na(strain)) metadata_parts <- c(metadata_parts, strain)
    if (!is.na(year)) metadata_parts <- c(metadata_parts, year)
    if (!is.na(country)) metadata_parts <- c(metadata_parts, country)
    if (!is.na(host)) metadata_parts <-c(metadata_parts, host)
    
    if (length(metadata_parts) > 0) {
      metadata_str <- paste(metadata_parts, collapse = "/")
      final_label <- paste0(clean_label, " [", metadata_str, "]")
    }
    
    # Add position annotations if sequence is highlighted
    if (label %in% names(sequence_to_positions)) {
      positions <- sequence_to_positions[[label]]
      positions_str <- paste(positions, collapse = ",")
      final_label <- paste0(final_label, " (", positions_str, ")")
    }
    
    # Add host emoji - ALWAYS add an emoji
    if (!is.na(host) && host != "Unknown" && host %in% names(species_symbols)) {
      # Known host with emoji
      emoji <- species_symbols[[host]]
    } else {
      # Unknown host or no host data - use question mark
      emoji <- default_symbol
    }
    final_label <- paste0(final_label, " ", emoji)
    
    modified_labels[i] <- final_label
  }
  
  return(modified_labels)
}

# Function to add comprehensive metadata to tip labels (for non-highlighted sequences)
add_comprehensive_metadata <- function(tip_labels, comprehensive_info) {
  modified_labels <- tip_labels
  
  for (i in seq_along(tip_labels)) {
    label <- tip_labels[i]
    
    # Clean up the label - remove H1N1__, H5N1__, H7N9__ prefixes
    clean_label <- label
    if (grepl("^H[1-9]N[1-9]__", label)) {
      clean_label <- sub("^H[1-9]N[1-9]__", "", label)
    }
    
    # Start with cleaned label
    final_label <- clean_label
    
    # Add metadata if available
    strain <- comprehensive_info$Strain_Type[i]
    year <- comprehensive_info$Year[i]
    country <- comprehensive_info$Country[i]
    host <- comprehensive_info$Host[i]
    
    # Build metadata string
    metadata_parts <- c()
    if (!is.na(strain)) metadata_parts <- c(metadata_parts, strain)
    if (!is.na(year)) metadata_parts <- c(metadata_parts, year)
    if (!is.na(country)) metadata_parts <- c(metadata_parts, country)
    if (!is.na(host)) metadata_parts <-c(metadata_parts, host)
    
    if (length(metadata_parts) > 0) {
      metadata_str <- paste(metadata_parts, collapse = "/")
      final_label <- paste0(clean_label, " [", metadata_str, "]")
    }
    
    # Add host emoji - ALWAYS add an emoji
    if (!is.na(host) && host != "Unknown" && host %in% names(species_symbols)) {
      # Known host with emoji
      emoji <- species_symbols[[host]]
    } else {
      # Unknown host or no host data - use question mark
      emoji <- default_symbol
    }
    final_label <- paste0(final_label, " ", emoji)
    
    modified_labels[i] <- final_label
  }
  
  return(modified_labels)
}

# MODIFIED FUNCTION: Identify clades to collapse - MORE AGGRESSIVE
identify_collapse_clades <- function(tree, interesting_tips, min_branch_length = 0.01, min_clade_size = 2, max_clades = 100) {
  collapse_groups <- list()
  collapse_nodes <- c()
  
  # Get all internal nodes (excluding root)
  n_tips <- length(tree$tip.label)
  if (is.null(tree$Nnode) || tree$Nnode == 0) {
    message("No internal nodes found in tree, skipping collapse")
    return(list(groups = collapse_groups, nodes = collapse_nodes))
  }
  
  internal_nodes <- (n_tips + 1):(n_tips + tree$Nnode)
  
  # Store candidates for collapsing
  candidates <- data.frame(
    node = integer(),
    descendant_count = integer(),
    branch_length = numeric(),
    stringsAsFactors = FALSE
  )
  
  for (node in internal_nodes) {
    # Skip invalid nodes
    if (is.na(node) || node <= 0) next
    
    tryCatch({
      # Get descendants
      descendants <- Descendants(tree, node, type = "tips")
      if (length(descendants) == 0 || is.null(descendants[[1]]) || length(descendants[[1]]) < min_clade_size) next
      
      descendant_tips <- tree$tip.label[descendants[[1]]]
      
      # Validate tips
      valid_tips <- descendant_tips[descendant_tips %in% tree$tip.label]
      if (length(valid_tips) < min_clade_size) next
      
      # Skip if any descendants are interesting
      if (any(valid_tips %in% interesting_tips)) next
      
      # Get branch length
      edge_idx <- which(tree$edge[, 2] == node)
      branch_length <- 0.001  # Default for missing branch lengths
      
      if (length(edge_idx) > 0 && !is.null(tree$edge.length) && length(tree$edge.length) >= edge_idx) {
        bl <- tree$edge.length[edge_idx[1]]
        if (!is.na(bl) && is.finite(bl)) {
          branch_length <- bl
        }
      }
      
      # Consider for collapse if branch length is short
      if (branch_length <= min_branch_length) {
        candidates <- rbind(candidates, data.frame(
          node = node,
          descendant_count = length(valid_tips),
          branch_length = branch_length,
          stringsAsFactors = FALSE
        ))
      }
    }, error = function(e) {
      message("Skipping node ", node, " due to error: ", e$message)
    })
  }
  
  # Sort candidates by branch length (shortest first) and descendant count (largest first)
  if (nrow(candidates) > 0) {
    candidates <- candidates[order(candidates$branch_length, -candidates$descendant_count), ]
    
    # Select top candidates, avoiding overlaps
    selected_nodes <- c()
    used_tips <- character(0)
    
    for (i in seq_len(min(nrow(candidates), max_clades))) {
      node <- candidates$node[i]
      
      tryCatch({
        descendants <- Descendants(tree, node, type = "tips")
        if (length(descendants) > 0 && !is.null(descendants[[1]])) {
          descendant_tips <- tree$tip.label[descendants[[1]]]
          valid_tips <- descendant_tips[descendant_tips %in% tree$tip.label]
          
          # Check for overlap with already selected tips
          if (length(valid_tips) >= min_clade_size && !any(valid_tips %in% used_tips)) {
            if (node > n_tips && node <= (n_tips + tree$Nnode)) {
              selected_nodes <- c(selected_nodes, node)
              used_tips <- c(used_tips, valid_tips)
              collapse_groups[[paste0("clade_", length(collapse_groups) + 1)]] <- valid_tips
            }
          }
        }
      }, error = function(e) {
        message("Skipping candidate node ", node, " due to error: ", e$message)
      })
    }
    
    collapse_nodes <- selected_nodes
  }
  
  return(list(groups = collapse_groups, nodes = collapse_nodes))
}

# MODIFIED MAIN FUNCTION: Create tree with aggressive collapsing, no green dots, and tiny triangles
create_collapsed_tree <- function(tree_file, ptm_file, revbayes_file, json_file,
                                  red_positions = NULL, blue_positions = NULL,
                                  collapse_uninteresting = TRUE, min_branch_length = 0.01,
                                  tree_layout = "circular") {
  
  # Read and validate tree
  tree <- read.nexus(tree_file)
  if (is.null(tree$edge) || nrow(tree$edge) == 0) {
    stop("Invalid tree: no edges found")
  }
  if (is.null(tree$Nnode) || tree$Nnode == 0) {
    stop("Invalid tree: no internal nodes found")
  }
  
  # Get comprehensive information for all tip labels
  comprehensive_info <- get_comprehensive_info(tree$tip.label, json_file)
  
  # Get PTM data
  ptm_data <- read_ptm_metadata(ptm_file)
  
  # Combine all positions
  all_positions <- c(red_positions, blue_positions)
  
  # Get sequences for all positions
  position_data <- get_sequences_for_multiple_positions(ptm_data, all_positions)
  
  # Get sequences for red positions
  red_sequences <- character(0)
  if (!is.null(red_positions)) {
    for (pos in red_positions) {
      red_sequences <- c(red_sequences, get_sequences_for_position(ptm_data, pos))
    }
    red_sequences <- unique(red_sequences)
  }
  
  # Get sequences for blue positions
  blue_sequences <- character(0)
  if (!is.null(blue_positions)) {
    for (pos in blue_positions) {
      blue_sequences <- c(blue_sequences, get_sequences_for_position(ptm_data, pos))
    }
    blue_sequences <- unique(blue_sequences)
  }
  
  # Handle overlapping sequences (red takes priority over blue)
  blue_sequences <- setdiff(blue_sequences, red_sequences)
  
  # All interesting sequences
  interesting_sequences <- c(red_sequences, blue_sequences)
  
  message(paste("Red sequences:", length(red_sequences)))
  message(paste("Blue sequences:", length(blue_sequences)))
  if (length(intersect(red_sequences, blue_sequences)) > 0) {
    message("Note: Some sequences have both red and blue positions - they will be colored RED (red takes priority)")
  }
  
  # Identify clades to collapse aggressively
  collapse_data <- list(groups = list(), nodes = c())
  if (collapse_uninteresting && length(interesting_sequences) > 0) {
    message("Identifying clades with short branches for aggressive collapsing...")
    collapse_data <- identify_collapse_clades(tree, interesting_sequences, min_branch_length,
                                              min_clade_size = 2, max_clades = 100)
    
    if (length(collapse_data$groups) > 0) {
      message(paste("Found", length(collapse_data$groups), "clades to collapse"))
      tryCatch({
        tree <- groupOTU(tree, collapse_data$groups)
      }, error = function(e) {
        message("Error in groupOTU, proceeding without collapse: ", e$message)
        collapse_data <- list(groups = list(), nodes = c())
      })
    } else {
      message("No suitable clades found for collapse")
    }
  }
  
  # Get branch rates
  rates <- read_revbayes_rates(revbayes_file)
  
  # Create edge data frame for branch rates
  n_branches <- nrow(tree$edge)
  if (length(rates) < n_branches) {
    rates <- c(rates, rep(NA, n_branches - length(rates)))
  } else if (length(rates) > n_branches) {
    rates <- rates[1:n_branches]
  }
  edge_data <- data.frame(
    node = tree$edge[, 2],
    rate = rates
  )
  
  # Create title with position information
  title_parts <- character(0)
  
  # Create base plot with thicker branches
  p <- ggtree(tree, branch.length = "none", layout = tree_layout, size = 1) %<+% edge_data +
    aes(color = rate) +
    scale_color_viridis_c(name = "Branch Rate", option = "magma", na.value = "grey50") +
    theme(legend.position = "right") 
  
  # ===========================
  # Replace collapse loop to produce very small triangles
  # ===========================
  if (length(collapse_data$nodes) > 0) {
    message(paste("Collapsing", length(collapse_data$nodes), "clades into tiny triangles..."))
    for (node in collapse_data$nodes) {
      # Set expand = 0.01 to shrink the triangle width
      p <- collapse(p, node = node, mode = "none")
    }
    message(paste("Successfully collapsed", length(collapse_data$nodes), "clades"))
  } else {
    message("No clades were collapsed")
  }
  # ===========================
  
  # Update plot title with collapse summary
  n_collapsed <- length(collapse_data$nodes)
  if (n_collapsed > 0) {
    collapse_summary <- paste0("n = ", n_collapsed, " clades collapsed")
    message(collapse_summary)
  } else {
    message("No branches were collapsed")
  }
  
  # Get sequences for positions (matching current tree tips)
  red_sequences <- intersect(red_sequences, tree$tip.label)
  blue_sequences <- intersect(blue_sequences, tree$tip.label)
  
  # Create comprehensive tip labels for all sequences
  comprehensive_tip_labels <- create_comprehensive_tip_labels(tree$tip.label,
                                                              position_data$sequence_to_positions,
                                                              comprehensive_info)
  
  # Add tip labels - non-highlighted in small gray text
  non_highlighted_indices <- !tree$tip.label %in% c(red_sequences, blue_sequences)
  non_highlighted_labels <- add_comprehensive_metadata(tree$tip.label[non_highlighted_indices],
                                                       comprehensive_info[non_highlighted_indices, , drop = FALSE])
  p <- p + geom_tiplab(data = function(x) {
    non_high <- subset(x, !label %in% c(red_sequences, blue_sequences))
    if (nrow(non_high) > 0) {
      non_high$comprehensive_label <- non_highlighted_labels[match(non_high$label,
                                                                   tree$tip.label[non_highlighted_indices])]
      return(non_high)
    } else {
      return(non_high)
    }
  },
  aes(label = comprehensive_label),
  size = 1, hjust = 0, color = "gray60")
  
  # Add tip labels - red sequences (bold, small)
  if (length(red_sequences) > 0) {
    p <- p + geom_tiplab(data = function(x) {
      red_data <- subset(x, label %in% red_sequences)
      if (nrow(red_data) > 0) {
        red_data$comprehensive_label <- comprehensive_tip_labels[match(red_data$label,
                                                                       tree$tip.label)]
        return(red_data)
      } else {
        return(red_data)
      }
    },
    aes(label = comprehensive_label),
    size = 1, hjust = 0, color = "red", fontface = "bold")
  }
  
  # Add tip labels - blue sequences (bold, small)
  if (length(blue_sequences) > 0) {
    p <- p + geom_tiplab(data = function(x) {
      blue_data <- subset(x, label %in% blue_sequences)
      if (nrow(blue_data) > 0) {
        blue_data$comprehensive_label <- comprehensive_tip_labels[match(blue_data$label,
                                                                        tree$tip.label)]
        return(blue_data)
      } else {
        return(blue_data)
      }
    },
    aes(label = comprehensive_label),
    size = 1, hjust = 0, color = "blue", fontface = "bold")
  }
  
  return(p)
}

# ==============================================================================
# USAGE EXAMPLE - AGGRESSIVE COLLAPSED TREE WITH TINY TRIANGLES
# ==============================================================================

# Create tree that aggressively collapses short branches into very small triangles
plot <- create_collapsed_tree(
  tree_file = here("data/trees/HA_full_tree.nex"),
  ptm_file = here("data/ptm/mappings/HA_full_position_metadata.txt"),
  revbayes_file = here("mcmc_logs/HA_full_branch_rates.log"),
  json_file = here("data/ptm/mappings/HA_full_id_mapping.json"),
  red_positions = c(182, 570),
  blue_positions = c(56, 143, 281, 519),
  collapse_uninteresting = TRUE,
  min_branch_length = 0.0001,  # Aggressive collapsing for short branches
  tree_layout = "circular"
)

print(plot)


# (3) Open a JPEG device at 5 in wide, calculated height inches, and a decent resolution:
dir.create(here("results/summary_figures"), showWarnings = FALSE, recursive = TRUE)
png(
  filename = here("results/summary_figures/HA_collapsed_tree.png"),
  width    =  9.5,           # 5 inches wide
  height = 9,
  units    = "in",
  res      = 600           # 300 dpi
)

# (4) Draw the plot onto that device:
print(plot)

# (5) Close the device:
dev.off()


