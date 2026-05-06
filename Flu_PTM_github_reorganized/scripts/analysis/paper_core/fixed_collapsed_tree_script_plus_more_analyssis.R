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
  h1n1 = "C:/Users/Ashley/Documents/flu_project/H1N1/H1N1_seq_info.csv",
  h5n1 = "C:/Users/Ashley/Documents/flu_project/H5N1/H5N1_seq_info.csv",
  h7n9 = "C:/Users/Ashley/Documents/flu_project/H7N9/H7N9_seq_info.csv"
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
  tree_file = "HA_full_tree.nex",
  ptm_file = "HA_full_position_metadata.txt",
  revbayes_file = "HA_full_branch_rates.log",
  json_file = "HA_full_id_mapping.json",
  red_positions = c(182, 570),
  blue_positions = c(56, 143, 281, 519),
  collapse_uninteresting = TRUE,
  min_branch_length = 0.0001,  # Aggressive collapsing for short branches
  tree_layout = "circular"
)

print(plot)



# (3) Open a JPEG device at 5 in wide, calculated height inches, and a decent resolution:
#png(
#  filename = "HA_collapsed_tree.png",
#  width    =  9.5,           # 5 inches wide
#  height = 9,
#  units    = "in",
#  res      = 600           # 300 dpi
#)

# (4) Draw the plot onto that device:
#print(plot)

# (5) Close the device:
#dev.off()



# ==============================================================================
# WRAPPER FUNCTION + DIAGNOSTICS (FULL COPY-PASTE SAFE)
# ==============================================================================

create_collapsed_tree_with_diagnostics <- function(tree_file, ptm_file, revbayes_file, json_file,
                                                   red_positions = NULL, blue_positions = NULL,
                                                   collapse_uninteresting = TRUE, min_branch_length = 0.01,
                                                   tree_layout = "circular") {
  
  # --- replicate core data prep steps here ---
  
  tree <- read.nexus(tree_file)
  comprehensive_info <- get_comprehensive_info(tree$tip.label, json_file)
  ptm_data <- read_ptm_metadata(ptm_file)
  
  all_positions <- c(red_positions, blue_positions)
  position_data <- get_sequences_for_multiple_positions(ptm_data, all_positions)
  
  red_sequences <- character(0)
  if (!is.null(red_positions)) {
    for (pos in red_positions) {
      red_sequences <- c(red_sequences, get_sequences_for_position(ptm_data, pos))
    }
    red_sequences <- unique(red_sequences)
  }
  
  blue_sequences <- character(0)
  if (!is.null(blue_positions)) {
    for (pos in blue_positions) {
      blue_sequences <- c(blue_sequences, get_sequences_for_position(ptm_data, pos))
    }
    blue_sequences <- unique(blue_sequences)
  }
  
  blue_sequences <- setdiff(blue_sequences, red_sequences)
  interesting_sequences <- c(red_sequences, blue_sequences)
  
  collapse_data <- list(groups = list(), nodes = c())
  if (collapse_uninteresting && length(interesting_sequences) > 0) {
    collapse_data <- identify_collapse_clades(tree, interesting_sequences, min_branch_length,
                                              min_clade_size = 2, max_clades = 100)
    if (length(collapse_data$groups) > 0) {
      tree <- groupOTU(tree, collapse_data$groups)
    }
  }
  
  rates <- read_revbayes_rates(revbayes_file)
  
  # Now call your original plotting function
  plot <- create_collapsed_tree(tree_file = tree_file,
                                ptm_file = ptm_file,
                                revbayes_file = revbayes_file,
                                json_file = json_file,
                                red_positions = red_positions,
                                blue_positions = blue_positions,
                                collapse_uninteresting = collapse_uninteresting,
                                min_branch_length = min_branch_length,
                                tree_layout = tree_layout)
  
  # Return everything
  return(list(
    plot = plot,
    collapse_data = collapse_data,
    rates = rates,
    tree = tree,
    comprehensive_info = comprehensive_info,
    red_sequences = intersect(red_sequences, tree$tip.label),
    blue_sequences = intersect(blue_sequences, tree$tip.label)
  ))
}

# ==============================================================================
# RUN TREE WITH DIAGNOSTICS
# ==============================================================================

tree_result <- create_collapsed_tree_with_diagnostics(
  tree_file = "HA_full_tree.nex",
  ptm_file = "HA_full_position_metadata.txt",
  revbayes_file = "HA_full_branch_rates.log",
  json_file = "HA_full_id_mapping.json",
  red_positions = c(182, 570),
  blue_positions = c(56, 143, 281, 519),
  collapse_uninteresting = TRUE,
  min_branch_length = 0.0001,
  tree_layout = "circular"
)

print(tree_result$plot)

# ==============================================================================
# DIAGNOSTICS BLOCK
# ==============================================================================

collapse_data <- tree_result$collapse_data
rates <- tree_result$rates
tree <- tree_result$tree
comprehensive_info <- tree_result$comprehensive_info
red_sequences <- tree_result$red_sequences
blue_sequences <- tree_result$blue_sequences

# (1) Summary of clades collapsed
cat("Number of clades collapsed:", length(collapse_data$nodes), "\n")
if (length(collapse_data$groups) > 0) {
  clade_sizes <- sapply(collapse_data$groups, length)
  cat("Sizes of collapsed clades:\n")
  print(clade_sizes)
}

# (2) Branch rate summary
cat("Branch rate summary:\n")
cat("Mean:", mean(rates, na.rm = TRUE), "\n")
cat("SD:", sd(rates, na.rm = TRUE), "\n")
cat("Min:", min(rates, na.rm = TRUE), "\n")
cat("Max:", max(rates, na.rm = TRUE), "\n")

# (3) Red / blue sequences
cat("Number of RED sequences (fast PTM sites):", length(red_sequences), "\n")
cat("Number of BLUE sequences (slow PTM sites):", length(blue_sequences), "\n")

# (4) Host species represented overall
cat("Host species represented (count):\n")
host_counts <- table(comprehensive_info$Host)
print(host_counts)

# (5) Hosts represented in RED sequences
red_hosts <- comprehensive_info %>%
  filter(tip_label %in% red_sequences) %>%
  pull(Host)
cat("Host species in RED sequences:\n")
print(table(red_hosts))

# (6) Hosts represented in BLUE sequences
blue_hosts <- comprehensive_info %>%
  filter(tip_label %in% blue_sequences) %>%
  pull(Host)
cat("Host species in BLUE sequences:\n")
print(table(blue_hosts))

# (7) Branch rate summary for RED and BLUE sequences
edge_df <- data.frame(parent = tree$edge[, 1], child = tree$edge[, 2])
n_tips <- length(tree$tip.label)

red_indices <- which(tree$tip.label %in% red_sequences)
blue_indices <- which(tree$tip.label %in% blue_sequences)

red_edge_indices <- which(tree$edge[, 2] %in% red_indices)
blue_edge_indices <- which(tree$edge[, 2] %in% blue_indices)

cat("Branch rates for RED sequences:\n")
if (length(red_edge_indices) > 0) {
  cat("Mean:", mean(rates[red_edge_indices], na.rm = TRUE), "\n")
  cat("SD:", sd(rates[red_edge_indices], na.rm = TRUE), "\n")
  cat("Min:", min(rates[red_edge_indices], na.rm = TRUE), "\n")
  cat("Max:", max(rates[red_edge_indices], na.rm = TRUE), "\n")
} else {
  cat("No matching branches found for RED sequences.\n")
}

cat("Branch rates for BLUE sequences:\n")
if (length(blue_edge_indices) > 0) {
  cat("Mean:", mean(rates[blue_edge_indices], na.rm = TRUE), "\n")
  cat("SD:", sd(rates[blue_edge_indices], na.rm = TRUE), "\n")
  cat("Min:", min(rates[blue_edge_indices], na.rm = TRUE), "\n")
  cat("Max:", max(rates[blue_edge_indices], na.rm = TRUE), "\n")
} else {
  cat("No matching branches found for BLUE sequences.\n")
}

# (8) Branch rate summary by host species
tip_node_indices <- 1:length(tree$tip.label)
tip_rates <- rates[match(tip_node_indices, tree$edge[, 2])]

host_rate_df <- data.frame(
  Host = comprehensive_info$Host,
  Rate = tip_rates
) %>%
  filter(!is.na(Rate))

cat("Branch rate summary by host species:\n")
host_summary <- host_rate_df %>%
  group_by(Host) %>%
  summarize(
    Mean_Rate = mean(Rate, na.rm = TRUE),
    SD_Rate = sd(Rate, na.rm = TRUE),
    Min_Rate = min(Rate, na.rm = TRUE),
    Max_Rate = max(Rate, na.rm = TRUE),
    n = n()
  )
print(host_summary, n=40)




# ============================
# Group rare hosts into "Other" and re-run host association tests for BOTH Fast (Red) and Slow (Blue) PTM sites
# ============================

# Set threshold for grouping rare hosts (e.g. <5 sequences → "Other")
min_host_count <- 5

# Count total sequences per host
host_counts_total <- table(comprehensive_info$Host)

# Mark "rare" hosts
rare_hosts <- names(host_counts_total[host_counts_total < min_host_count])

# Create new column grouping rare hosts as "Other"
comprehensive_info$Host_grouped <- ifelse(comprehensive_info$Host %in% rare_hosts, "Other", comprehensive_info$Host)

# ==================================
# Test for Fast PTM (RED) sites
# ==================================

# Mark which sequences are RED (fast PTM site) or not
comprehensive_info$Fast_PTM <- ifelse(comprehensive_info$tip_label %in% red_sequences, "Yes", "No")

# Build grouped contingency table: Host_grouped × Fast_PTM
host_fast_table_grouped <- table(comprehensive_info$Host_grouped, comprehensive_info$Fast_PTM)

cat("Grouped contingency table (Host vs Fast PTM site present):\n")
print(host_fast_table_grouped)

# Run Fisher's exact test (simulate p-value for large table)
cat("\nFisher's Exact Test result (Grouped Fast PTM sites):\n")
fisher_result_fast_grouped <- fisher.test(host_fast_table_grouped, simulate.p.value = TRUE)
print(fisher_result_fast_grouped)

# ==================================
# Test for Slow PTM (BLUE) sites
# ==================================

# Mark which sequences are BLUE (slow PTM site) or not
comprehensive_info$Slow_PTM <- ifelse(comprehensive_info$tip_label %in% blue_sequences, "Yes", "No")

# Build grouped contingency table: Host_grouped × Slow_PTM
host_slow_table_grouped <- table(comprehensive_info$Host_grouped, comprehensive_info$Slow_PTM)

cat("\nGrouped contingency table (Host vs Slow PTM site present):\n")
print(host_slow_table_grouped)

# Run Fisher's exact test (simulate p-value for large table)
cat("\nFisher's Exact Test result (Grouped Slow PTM sites):\n")
fisher_result_slow_grouped <- fisher.test(host_slow_table_grouped, simulate.p.value = TRUE)
print(fisher_result_slow_grouped)









# Host Switching Analysis for PTM Sites
# This script analyzes whether sequences with mutations at specific positions
# show patterns of host switching
# Host Switching Analysis for PTM Sites
# This script analyzes whether sequences with mutations at specific positions
# show patterns of host switching
# Host Switching Analysis for PTM Sites
# This script analyzes whether sequences with mutations at specific positions
# show patterns of host switching

library(ape)
library(dplyr)
library(tidyr)
library(ggplot2)
library(readr)
library(jsonlite)
library(stringr)
library(phytools)
library(ggtree)
library(gridExtra)

# Source the original functions from your script
# (You'll need to source or include the functions from your original script)

# Additional functions for host switching analysis

# Function to analyze host patterns in phylogenetic context
analyze_host_switching <- function(tree, tip_info, sequences_of_interest) {
  # Create a mapping of tips to hosts
  host_mapping <- setNames(tip_info$Host, tip_info$tip_label)
  
  # Clean up host names for consistency
  host_mapping <- gsub("Homo sapiens", "Human", host_mapping)
  host_mapping <- gsub("Sus scrofa", "Swine", host_mapping)
  host_mapping <- gsub("Gallus gallus", "Chicken", host_mapping)
  host_mapping <- gsub("Anas platyrhynchos", "Duck", host_mapping)
  host_mapping[is.na(host_mapping)] <- "Unknown"
  
  # Identify nodes with sequences of interest
  interest_tips <- which(tree$tip.label %in% sequences_of_interest)
  
  # Find internal nodes that are ancestors of sequences of interest
  ancestors_of_interest <- unique(unlist(lapply(interest_tips, function(tip) {
    nodepath(tree, from = tip, to = length(tree$tip.label) + 1)
  })))
  
  # Analyze host transitions
  host_transitions <- data.frame(
    node = integer(),
    parent_node = integer(),
    descendant_hosts = character(),
    n_hosts = integer(),
    has_interest_seq = logical(),
    stringsAsFactors = FALSE
  )
  
  # For each internal node
  for (node in unique(tree$edge[,1])) {
    descendants <- Descendants(tree, node, type = "tips")[[1]]
    desc_hosts <- unique(host_mapping[tree$tip.label[descendants]])
    desc_hosts <- desc_hosts[desc_hosts != "Unknown"]
    
    has_interest <- any(tree$tip.label[descendants] %in% sequences_of_interest)
    
    host_transitions <- rbind(host_transitions, data.frame(
      node = node,
      parent_node = Ancestors(tree, node, type = "parent"),
      descendant_hosts = paste(desc_hosts, collapse = ";"),
      n_hosts = length(desc_hosts),
      has_interest_seq = has_interest,
      stringsAsFactors = FALSE
    ))
  }
  
  return(list(
    host_mapping = host_mapping,
    host_transitions = host_transitions
  ))
}

# Function to calculate host switching metrics
calculate_host_switching_metrics <- function(tree, tip_info, red_sequences, blue_sequences) {
  # Combine all sequences of interest
  all_interest_sequences <- c(red_sequences, blue_sequences)
  
  # Get host information
  host_analysis <- analyze_host_switching(tree, tip_info, all_interest_sequences)
  
  # Calculate metrics for nodes with and without sequences of interest
  metrics <- host_analysis$host_transitions %>%
    group_by(has_interest_seq) %>%
    summarise(
      mean_hosts = mean(n_hosts, na.rm = TRUE),
      median_hosts = median(n_hosts, na.rm = TRUE),
      max_hosts = max(n_hosts, na.rm = TRUE),
      n_multi_host = sum(n_hosts > 1),
      n_nodes = n(),
      prop_multi_host = n_multi_host / n_nodes
    )
  
  # Statistical test
  with_interest <- host_analysis$host_transitions$n_hosts[host_analysis$host_transitions$has_interest_seq]
  without_interest <- host_analysis$host_transitions$n_hosts[!host_analysis$host_transitions$has_interest_seq]
  
  wilcox_test <- wilcox.test(with_interest, without_interest)
  
  return(list(
    metrics = metrics,
    wilcox_test = wilcox_test,
    host_analysis = host_analysis
  ))
}

# Function to analyze host switching patterns by position
analyze_host_switching_by_position <- function(tree, tip_info, ptm_data, positions, position_names = NULL) {
  if (is.null(position_names)) {
    position_names <- as.character(positions)
  }
  
  results <- list()
  
  for (i in seq_along(positions)) {
    pos <- positions[i]
    pos_name <- position_names[i]
    
    # Get sequences for this position
    sequences <- get_sequences_for_position(ptm_data, pos)
    sequences <- intersect(sequences, tree$tip.label)
    
    if (length(sequences) > 0) {
      # Get host information for these sequences
      seq_hosts <- tip_info$Host[tip_info$tip_label %in% sequences]
      seq_hosts <- seq_hosts[!is.na(seq_hosts) & seq_hosts != ""]
      
      if (length(seq_hosts) > 0) {
        # Analyze host patterns for these sequences
        host_counts <- table(seq_hosts)
        
        # Calculate diversity metrics
        shannon_diversity <- if (length(host_counts) > 1) {
          -sum((host_counts/sum(host_counts)) * log(host_counts/sum(host_counts)))
        } else {
          0
        }
        
        results[[pos_name]] <- list(
          position = pos,
          n_sequences = length(sequences),
          host_counts = host_counts,
          n_unique_hosts = length(host_counts),
          shannon_diversity = shannon_diversity,
          sequences = sequences
        )
      } else {
        # No valid host information
        results[[pos_name]] <- list(
          position = pos,
          n_sequences = length(sequences),
          host_counts = table(character(0)),
          n_unique_hosts = 0,
          shannon_diversity = 0,
          sequences = sequences
        )
      }
    }
  }
  
  return(results)
}

# Function to create host switching visualization
create_host_switching_plot <- function(switching_metrics) {
  # Create bar plot comparing nodes with/without sequences of interest
  p1 <- ggplot(switching_metrics$metrics, aes(x = factor(has_interest_seq), y = mean_hosts)) +
    geom_bar(stat = "identity", fill = c("gray70", "red")) +
    geom_errorbar(aes(ymin = mean_hosts - 0.1, ymax = mean_hosts + 0.1), width = 0.2) +
    labs(x = "Has Sequence of Interest", y = "Mean Number of Host Species per Node",
         title = "Host Diversity at Phylogenetic Nodes") +
    scale_x_discrete(labels = c("No", "Yes")) +
    theme_minimal()
  
  # Create proportion plot
  p2 <- ggplot(switching_metrics$metrics, aes(x = factor(has_interest_seq), y = prop_multi_host)) +
    geom_bar(stat = "identity", fill = c("gray70", "red")) +
    labs(x = "Has Sequence of Interest", y = "Proportion of Multi-Host Nodes",
         title = "Frequency of Host Switching") +
    scale_x_discrete(labels = c("No", "Yes")) +
    theme_minimal()
  
  return(list(p1 = p1, p2 = p2))
}

# Function to create host distribution plot by position
create_host_distribution_plot <- function(position_analysis) {
  # Prepare data for plotting
  plot_data <- data.frame()
  
  for (pos_name in names(position_analysis)) {
    pos_data <- position_analysis[[pos_name]]
    
    # Check if there are any host counts
    if (length(pos_data$host_counts) > 0) {
      host_df <- data.frame(
        Position = pos_name,
        Host = names(pos_data$host_counts),
        Count = as.numeric(pos_data$host_counts),
        stringsAsFactors = FALSE
      )
      plot_data <- rbind(plot_data, host_df)
    }
  }
  
  # Check if we have any data to plot
  if (nrow(plot_data) == 0) {
    message("No host data available for plotting")
    # Create empty plot with message
    p <- ggplot() + 
      annotate("text", x = 1, y = 1, label = "No host data available", size = 6) +
      theme_void()
    return(p)
  }
  
  # Clean up host names for better visualization
  plot_data$Host <- gsub("Homo sapiens", "Human", plot_data$Host)
  plot_data$Host <- gsub("Sus scrofa", "Swine", plot_data$Host)
  plot_data$Host <- gsub("Gallus gallus", "Chicken", plot_data$Host)
  plot_data$Host[is.na(plot_data$Host)] <- "Unknown"
  
  # Create stacked bar plot
  p <- ggplot(plot_data, aes(x = Position, y = Count, fill = Host)) +
    geom_bar(stat = "identity") +
    labs(title = "Host Distribution by PTM Position",
         x = "Position", y = "Number of Sequences") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    scale_fill_brewer(palette = "Set3")
  
  return(p)
}

# Function to perform phylogenetic analysis of host switching
phylogenetic_host_analysis <- function(tree, tip_info, sequences_of_interest) {
  # Clean host names first
  tip_info$Host[is.na(tip_info$Host)] <- "Unknown"
  tip_info$Host[tip_info$Host == ""] <- "Unknown"
  
  # Create simplified host categories
  simplified_hosts <- character(length(tip_info$Host))
  for (i in seq_along(tip_info$Host)) {
    host <- tip_info$Host[i]
    if (grepl("Homo sapiens", host, ignore.case = TRUE)) {
      simplified_hosts[i] <- "Human"
    } else if (grepl("Sus scrofa", host, ignore.case = TRUE)) {
      simplified_hosts[i] <- "Swine"
    } else if (grepl("Gallus|Anas|Anatidae|duck|avian|bird", host, ignore.case = TRUE)) {
      simplified_hosts[i] <- "Avian"
    } else if (host == "Unknown" || is.na(host)) {
      simplified_hosts[i] <- "Unknown"
    } else {
      simplified_hosts[i] <- "Other"
    }
  }
  
  # Create mapping
  host_mapping <- setNames(factor(simplified_hosts), tip_info$tip_label)
  
  # Check if tree is rooted and binary
  if (!is.rooted(tree)) {
    message("Tree is not rooted. Rooting at midpoint...")
    tree <- midpoint.root(tree)
  }
  
  # Check if tree is binary (fully dichotomous)
  if (!is.binary(tree)) {
    message("Tree is not fully dichotomous. Resolving polytomies...")
    tree <- multi2di(tree, random = FALSE)
  }
  
  # Try ancestral state reconstruction
  anc_states <- NULL
  tryCatch({
    anc_states <- ace(host_mapping, tree, type = "discrete", model = "ER")
  }, error = function(e) {
    message("Ancestral state reconstruction failed: ", e$message)
  })
  
  # If ASR failed, just use the simple analysis
  if (is.null(anc_states)) {
    return(list(
      anc_states = NULL,
      transitions_to_interest = NA,
      transitions_to_other = NA,
      total_branches_to_interest = sum(tree$tip.label %in% sequences_of_interest),
      error = "ASR failed"
    ))
  }
  
  # Count transitions along branches leading to sequences of interest
  transitions_to_interest <- 0
  transitions_to_other <- 0
  
  # Check if ancestral state reconstruction was successful
  if (is.null(anc_states)) {
    message("Ancestral state reconstruction failed. Skipping transition counting.")
    return(list(
      anc_states = NULL,
      transitions_to_interest = NA,
      transitions_to_other = NA,
      total_branches_to_interest = sum(tree$tip.label %in% sequences_of_interest)
    ))
  }
  
  for (i in 1:nrow(tree$edge)) {
    parent <- tree$edge[i, 1]
    child <- tree$edge[i, 2]
    
    # Skip if we can't determine states
    if (parent <= length(tree$tip.label) || child > length(tree$tip.label) + tree$Nnode) next
    
    # Get most likely states
    parent_state <- which.max(anc_states$lik.anc[parent - length(tree$tip.label), ])
    
    if (child <= length(tree$tip.label)) {
      # Terminal branch
      child_state <- as.numeric(host_mapping[tree$tip.label[child]])
      is_interest <- tree$tip.label[child] %in% sequences_of_interest
    } else {
      # Internal branch
      child_state <- which.max(anc_states$lik.anc[child - length(tree$tip.label), ])
      # Check if any descendants are sequences of interest
      descendants <- tryCatch({
        Descendants(tree, child, type = "tips")[[1]]
      }, error = function(e) {
        numeric(0)
      })
      is_interest <- any(tree$tip.label[descendants] %in% sequences_of_interest)
    }
    
    # Count transitions
    if (!is.na(parent_state) && !is.na(child_state) && parent_state != child_state) {
      if (is_interest) {
        transitions_to_interest <- transitions_to_interest + 1
      } else {
        transitions_to_other <- transitions_to_other + 1
      }
    }
  }
  
}
  
  # Alternative: Simpler host switching analysis without ancestral state reconstruction
  simple_host_switching_analysis <- function(tree, tip_info, sequences_of_interest) {
    # This function analyzes host switching patterns without requiring ASR
    
    host_switches <- data.frame(
      sequence = character(),
      seq_host = character(),
      sister_hosts = character(),
      n_sister_hosts = integer(),
      is_switch = logical(),
      stringsAsFactors = FALSE
    )
    
    for (seq in sequences_of_interest) {
      tip_idx <- which(tree$tip.label == seq)
      if (length(tip_idx) == 0) next
      
      # Find the parent node
      parent_edge_idx <- which(tree$edge[,2] == tip_idx)
      if (length(parent_edge_idx) == 0) next
      
      parent_node <- tree$edge[parent_edge_idx, 1]
      
      # Find sister taxa
      sister_edges <- which(tree$edge[,1] == parent_node & tree$edge[,2] != tip_idx)
      sister_nodes <- tree$edge[sister_edges, 2]
      
      # Get all descendant tips of sister nodes
      sister_tips <- c()
      for (node in sister_nodes) {
        if (node <= length(tree$tip.label)) {
          # It's a tip
          sister_tips <- c(sister_tips, tree$tip.label[node])
        } else {
          # It's an internal node
          desc <- tryCatch({
            Descendants(tree, node, type = "tips")[[1]]
          }, error = function(e) numeric(0))
          if (length(desc) > 0) {
            sister_tips <- c(sister_tips, tree$tip.label[desc])
          }
        }
      }
      
      # Get host information
      seq_host <- tip_info$Host[tip_info$tip_label == seq]
      sister_hosts <- unique(tip_info$Host[tip_info$tip_label %in% sister_tips])
      sister_hosts <- sister_hosts[!is.na(sister_hosts) & sister_hosts != "Unknown"]
      
      # Determine if this represents a host switch
      is_switch <- FALSE
      if (!is.na(seq_host) && seq_host != "Unknown" && length(sister_hosts) > 0) {
        is_switch <- !seq_host %in% sister_hosts
      }
      
      host_switches <- rbind(host_switches, data.frame(
        sequence = seq,
        seq_host = ifelse(is.na(seq_host), "Unknown", seq_host),
        sister_hosts = paste(sister_hosts, collapse = ";"),
        n_sister_hosts = length(sister_hosts),
        is_switch = is_switch,
        stringsAsFactors = FALSE
      ))
    }
    
    return(host_switches)
  }
  
  # Main analysis function
  perform_host_switching_analysis <- function(tree_file, ptm_file, json_file,
                                              red_positions, blue_positions,
                                              output_prefix = "host_switching") {
    
    # Read tree
    tree <- read.nexus(tree_file)
    
    # Get comprehensive information
    comprehensive_info <- get_comprehensive_info(tree$tip.label, json_file)
    
    # Read PTM data
    ptm_data <- read_ptm_metadata(ptm_file)
    
    # Get sequences for positions
    red_sequences <- character(0)
    if (!is.null(red_positions)) {
      for (pos in red_positions) {
        red_sequences <- c(red_sequences, get_sequences_for_position(ptm_data, pos))
      }
      red_sequences <- unique(intersect(red_sequences, tree$tip.label))
    }
    
    blue_sequences <- character(0)
    if (!is.null(blue_positions)) {
      for (pos in blue_positions) {
        blue_sequences <- c(blue_sequences, get_sequences_for_position(ptm_data, pos))
      }
      blue_sequences <- unique(intersect(blue_sequences, tree$tip.label))
    }
    
    # Remove overlaps (red takes priority)
    blue_sequences <- setdiff(blue_sequences, red_sequences)
    
    # Calculate host switching metrics
    message("Calculating host switching metrics...")
    switching_metrics <- calculate_host_switching_metrics(tree, comprehensive_info, 
                                                          red_sequences, blue_sequences)
    
    # Analyze by position
    message("Analyzing host distribution by position...")
    all_positions <- c(red_positions, blue_positions)
    position_names <- c(paste0("Pos", red_positions, "_red"), 
                        paste0("Pos", blue_positions, "_blue"))
    position_analysis <- analyze_host_switching_by_position(tree, comprehensive_info, 
                                                            ptm_data, all_positions, position_names)
    
    # Phylogenetic analysis
    message("Performing phylogenetic analysis of host switching...")
    all_sequences <- c(red_sequences, blue_sequences)
    phylo_analysis <- tryCatch({
      phylogenetic_host_analysis(tree, comprehensive_info, all_sequences)
    }, error = function(e) {
      message("Standard phylogenetic analysis failed: ", e$message)
      message("Using simple host switching analysis instead...")
      
      # Use the simpler analysis
      simple_results <- simple_host_switching_analysis(tree, comprehensive_info, all_sequences)
      
      # Convert to compatible format
      list(
        anc_states = NULL,
        transitions_to_interest = sum(simple_results$is_switch),
        transitions_to_other = NA,
        total_branches_to_interest = nrow(simple_results),
        simple_analysis = simple_results
      )
    })
    
    # Create visualizations
    message("Creating visualizations...")
    switching_plots <- create_host_switching_plot(switching_metrics)
    host_dist_plot <- create_host_distribution_plot(position_analysis)
    
    # Save plots
    pdf(paste0(output_prefix, "_plots.pdf"), width = 12, height = 10)
    
    # Page 1: Host switching metrics
    grid.arrange(switching_plots$p1, switching_plots$p2, ncol = 2)
    
    # Page 2: Host distribution by position
    print(host_dist_plot)
    
    dev.off()
    
    # Generate summary report
    report <- list(
      summary = data.frame(
        Metric = c("Total sequences with red positions",
                   "Total sequences with blue positions",
                   "Mean hosts per node (with interest seq)",
                   "Mean hosts per node (without interest seq)",
                   "Wilcoxon test p-value",
                   "Host transitions/switches (interest sequences)",
                   "Host transitions/switches (other sequences)"),
        Value = c(length(red_sequences),
                  length(blue_sequences),
                  switching_metrics$metrics$mean_hosts[switching_metrics$metrics$has_interest_seq],
                  switching_metrics$metrics$mean_hosts[!switching_metrics$metrics$has_interest_seq],
                  switching_metrics$wilcox_test$p.value,
                  ifelse(is.na(phylo_analysis$transitions_to_interest), 
                         "Not calculated", phylo_analysis$transitions_to_interest),
                  ifelse(is.na(phylo_analysis$transitions_to_other), 
                         "Not calculated", phylo_analysis$transitions_to_other))
      ),
      position_summary = do.call(rbind, lapply(names(position_analysis), function(pos) {
        data.frame(
          Position = pos,
          N_Sequences = position_analysis[[pos]]$n_sequences,
          N_Hosts = position_analysis[[pos]]$n_unique_hosts,
          Shannon_Diversity = position_analysis[[pos]]$shannon_diversity,
          stringsAsFactors = FALSE
        )
      })),
      detailed_metrics = switching_metrics,
      position_analysis = position_analysis,
      phylo_analysis = phylo_analysis
    )
    
    # Add simple analysis results if available
    if (!is.null(phylo_analysis$simple_analysis)) {
      report$host_switching_sequences <- phylo_analysis$simple_analysis
      
      # Print additional information about host switches
      cat("\n=== DETAILED HOST SWITCHING ANALYSIS ===\n")
      cat(sprintf("Total sequences analyzed: %d\n", nrow(phylo_analysis$simple_analysis)))
      cat(sprintf("Sequences with apparent host switches: %d (%.1f%%)\n", 
                  sum(phylo_analysis$simple_analysis$is_switch),
                  100 * sum(phylo_analysis$simple_analysis$is_switch) / nrow(phylo_analysis$simple_analysis)))
      
      if (sum(phylo_analysis$simple_analysis$is_switch) > 0) {
        cat("\nSequences with host switches:\n")
        switches <- phylo_analysis$simple_analysis[phylo_analysis$simple_analysis$is_switch, ]
        for (i in 1:nrow(switches)) {
          cat(sprintf("  %s (%s) - sisters in: %s\n", 
                      switches$sequence[i], switches$seq_host[i], switches$sister_hosts[i]))
        }
      }
    }
    
    # Save report
    saveRDS(report, paste0(output_prefix, "_analysis_results.rds"))
    
    # Print summary
    cat("\n=== HOST SWITCHING ANALYSIS SUMMARY ===\n")
    print(report$summary)
    cat("\n=== POSITION-SPECIFIC HOST DIVERSITY ===\n")
    print(report$position_summary)
    
    return(report)
  }
  
  # SIMPLIFIED VERSION - More robust host switching analysis
  simplified_host_switching_analysis <- function(tree_file, ptm_file, json_file,
                                                 red_positions, blue_positions,
                                                 output_prefix = "host_switching_simple") {
    
    # Read tree
    tree <- read.nexus(tree_file)
    
    # Get comprehensive information
    comprehensive_info <- get_comprehensive_info(tree$tip.label, json_file)
    
    # Clean up host information
    comprehensive_info$Host[is.na(comprehensive_info$Host) | comprehensive_info$Host == ""] <- "Unknown"
    
    # Read PTM data
    ptm_data <- read_ptm_metadata(ptm_file)
    
    # Get sequences for each position
    position_summary <- data.frame()
    all_sequences_by_position <- list()
    
    # Analyze red positions
    for (pos in red_positions) {
      seqs <- get_sequences_for_position(ptm_data, pos)
      seqs <- intersect(seqs, tree$tip.label)
      
      if (length(seqs) > 0) {
        hosts <- comprehensive_info$Host[comprehensive_info$tip_label %in% seqs]
        hosts <- hosts[hosts != "Unknown"]
        
        position_summary <- rbind(position_summary, data.frame(
          Position = paste0("Pos", pos, "_red"),
          Type = "red",
          N_Sequences = length(seqs),
          N_Hosts = length(unique(hosts)),
          Host_List = paste(sort(unique(hosts)), collapse = "; "),
          stringsAsFactors = FALSE
        ))
        
        all_sequences_by_position[[paste0("Pos", pos, "_red")]] <- seqs
      }
    }
    
    # Analyze blue positions
    for (pos in blue_positions) {
      seqs <- get_sequences_for_position(ptm_data, pos)
      seqs <- intersect(seqs, tree$tip.label)
      
      if (length(seqs) > 0) {
        hosts <- comprehensive_info$Host[comprehensive_info$tip_label %in% seqs]
        hosts <- hosts[hosts != "Unknown"]
        
        position_summary <- rbind(position_summary, data.frame(
          Position = paste0("Pos", pos, "_blue"),
          Type = "blue",
          N_Sequences = length(seqs),
          N_Hosts = length(unique(hosts)),
          Host_List = paste(sort(unique(hosts)), collapse = "; "),
          stringsAsFactors = FALSE
        ))
        
        all_sequences_by_position[[paste0("Pos", pos, "_blue")]] <- seqs
      }
    }
    
    # Get all unique sequences of interest
    all_interest_sequences <- unique(unlist(all_sequences_by_position))
    
    # Simple host switching detection
    message("Analyzing potential host switches...")
    host_switch_data <- simple_host_switching_analysis(tree, comprehensive_info, all_interest_sequences)
    
    # Create visualization
    pdf(paste0(output_prefix, "_simple_analysis.pdf"), width = 12, height = 10)
    
    # Plot 1: Host distribution by position
    if (nrow(position_summary) > 0) {
      p1 <- ggplot(position_summary, aes(x = Position, y = N_Hosts, fill = Type)) +
        geom_bar(stat = "identity") +
        scale_fill_manual(values = c("red" = "red", "blue" = "blue")) +
        labs(title = "Number of Host Species by PTM Position",
             x = "Position", y = "Number of Unique Host Species") +
        theme_minimal() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1))
      print(p1)
    }
    
    # Plot 2: Sequences vs host switches
    if (nrow(host_switch_data) > 0) {
      switch_summary <- data.frame(
        Category = c("No Host Switch", "Host Switch"),
        Count = c(sum(!host_switch_data$is_switch), sum(host_switch_data$is_switch))
      )
      
      p2 <- ggplot(switch_summary, aes(x = Category, y = Count, fill = Category)) +
        geom_bar(stat = "identity") +
        scale_fill_manual(values = c("No Host Switch" = "gray70", "Host Switch" = "red")) +
        labs(title = "Host Switching in Sequences with PTM Sites",
             y = "Number of Sequences") +
        theme_minimal()
      print(p2)
    }
    
    dev.off()
    
    # Print summary
    cat("\n=== SIMPLIFIED HOST SWITCHING ANALYSIS ===\n")
    cat("\nPosition Summary:\n")
    print(position_summary)
    
    if (nrow(host_switch_data) > 0) {
      cat(sprintf("\n\nTotal sequences analyzed: %d\n", nrow(host_switch_data)))
      cat(sprintf("Sequences with apparent host switches: %d (%.1f%%)\n",
                  sum(host_switch_data$is_switch),
                  100 * sum(host_switch_data$is_switch) / nrow(host_switch_data)))
      
      # Show examples of host switches
      if (sum(host_switch_data$is_switch) > 0) {
        cat("\nExamples of host switches:\n")
        switches <- host_switch_data[host_switch_data$is_switch, ]
        for (i in 1:min(10, nrow(switches))) {
          cat(sprintf("  %s (%s) - closest relatives in: %s\n",
                      switches$sequence[i], switches$seq_host[i], switches$sister_hosts[i]))
        }
      }
    }
    
    # Save results
    results <- list(
      position_summary = position_summary,
      host_switch_data = host_switch_data,
      sequences_by_position = all_sequences_by_position
    )
    
    saveRDS(results, paste0(output_prefix, "_simple_results.rds"))
    
    return(results)
  }
  
  # ==============================================================================
  # USAGE - Run the simplified analysis
  # ==============================================================================
  
  # Run the simplified analysis (more robust, doesn't require ASR)
  simple_results <- simplified_host_switching_analysis(
    tree_file = "HA_full_tree.nex",
    ptm_file = "HA_full_position_metadata.txt",
    json_file = "HA_full_id_mapping.json",
    red_positions = c(182, 570),
    blue_positions = c(56, 143, 281, 519),
    output_prefix = "HA_host_switching_simple"
  )
 
  
  
  # Additional analyses for simplified results
  # These functions work with the output from simplified_host_switching_analysis
  
  # 1. Analyze host diversity by position
  analyze_position_host_diversity <- function(simple_results, tree, comprehensive_info) {
    cat("\n=== HOST DIVERSITY BY POSITION ===\n")
    
    # Get sequences by position from simple results
    sequences_by_position <- simple_results$sequences_by_position
    
    for (pos_name in names(sequences_by_position)) {
      sequences <- sequences_by_position[[pos_name]]
      
      # Get host information for these sequences
      hosts <- comprehensive_info$Host[comprehensive_info$tip_label %in% sequences]
      hosts <- hosts[!is.na(hosts) & hosts != "Unknown"]
      
      if (length(hosts) > 0) {
        host_counts <- table(hosts)
        
        # Calculate Shannon diversity
        shannon <- if (length(host_counts) > 1) {
          -sum((host_counts/sum(host_counts)) * log(host_counts/sum(host_counts)))
        } else {
          0
        }
        
        cat(sprintf("\n%s:\n", pos_name))
        cat(sprintf("  Total sequences: %d\n", length(sequences)))
        cat(sprintf("  Sequences with known hosts: %d\n", length(hosts)))
        cat(sprintf("  Unique host species: %d\n", length(host_counts)))
        cat(sprintf("  Shannon diversity: %.3f\n", shannon))
        cat("  Host distribution:\n")
        
        # Sort host counts by frequency
        host_counts_sorted <- sort(host_counts, decreasing = TRUE)
        for (host in names(host_counts_sorted)) {
          cat(sprintf("    %s: %d (%.1f%%)\n", 
                      host, 
                      host_counts_sorted[host],
                      100 * host_counts_sorted[host] / sum(host_counts)))
        }
      }
    }
  }
  
  # 2. Compare host switching rates between red and blue positions
  compare_position_types <- function(simple_results) {
    cat("\n=== HOST SWITCHING COMPARISON: RED vs BLUE POSITIONS ===\n")
    
    host_switch_data <- simple_results$host_switch_data
    position_data <- simple_results$sequences_by_position
    
    # Create mapping of sequences to position types
    seq_to_type <- data.frame(
      sequence = character(),
      position_type = character(),
      stringsAsFactors = FALSE
    )
    
    for (pos_name in names(position_data)) {
      type <- ifelse(grepl("_red$", pos_name), "red", "blue")
      sequences <- position_data[[pos_name]]
      
      seq_to_type <- rbind(seq_to_type, data.frame(
        sequence = sequences,
        position_type = type,
        stringsAsFactors = FALSE
      ))
    }
    
    # Remove duplicates (keep first occurrence)
    seq_to_type <- seq_to_type[!duplicated(seq_to_type$sequence), ]
    
    # Merge with host switch data
    switch_analysis <- merge(host_switch_data, seq_to_type, by = "sequence", all.x = TRUE)
    
    # Calculate statistics by position type
    for (type in c("red", "blue")) {
      type_data <- switch_analysis[switch_analysis$position_type == type, ]
      
      if (nrow(type_data) > 0) {
        n_switches <- sum(type_data$is_switch)
        n_total <- nrow(type_data)
        percent_switch <- 100 * n_switches / n_total
        
        cat(sprintf("\n%s positions:\n", toupper(type)))
        cat(sprintf("  Total sequences: %d\n", n_total))
        cat(sprintf("  Host switches: %d (%.1f%%)\n", n_switches, percent_switch))
        
        # Show most common host switches
        if (n_switches > 0) {
          cat("  Examples of host switches:\n")
          switches <- type_data[type_data$is_switch, ]
          for (i in 1:min(3, nrow(switches))) {
            cat(sprintf("    %s: %s → (sisters in %s)\n",
                        switches$sequence[i],
                        switches$seq_host[i],
                        switches$sister_hosts[i]))
          }
        }
      }
    }
    
    # Statistical test between red and blue
    red_switches <- switch_analysis$is_switch[switch_analysis$position_type == "red"]
    blue_switches <- switch_analysis$is_switch[switch_analysis$position_type == "blue"]
    
    if (length(red_switches) > 0 && length(blue_switches) > 0) {
      fisher_test <- fisher.test(matrix(c(
        sum(red_switches), length(red_switches) - sum(red_switches),
        sum(blue_switches), length(blue_switches) - sum(blue_switches)
      ), nrow = 2))
      
      cat(sprintf("\nFisher's exact test p-value: %.4f\n", fisher_test$p.value))
      if (fisher_test$p.value < 0.05) {
        cat("Significant difference in host switching rates between red and blue positions\n")
      } else {
        cat("No significant difference in host switching rates between red and blue positions\n")
      }
    }
  }
  
  # 3. Find sequences that appear in multiple positions
  find_multi_position_sequences <- function(simple_results) {
    cat("\n=== SEQUENCES WITH MULTIPLE PTM POSITIONS ===\n")
    
    sequences_by_position <- simple_results$sequences_by_position
    
    # Create a table of sequence occurrences
    all_sequences <- unlist(sequences_by_position)
    seq_counts <- table(all_sequences)
    
    # Find sequences appearing in multiple positions
    multi_pos_seqs <- names(seq_counts)[seq_counts > 1]
    
    if (length(multi_pos_seqs) > 0) {
      cat(sprintf("Found %d sequences with multiple PTM positions:\n", length(multi_pos_seqs)))
      
      for (seq in multi_pos_seqs[1:min(10, length(multi_pos_seqs))]) {
        positions <- c()
        for (pos_name in names(sequences_by_position)) {
          if (seq %in% sequences_by_position[[pos_name]]) {
            positions <- c(positions, pos_name)
          }
        }
        cat(sprintf("  %s: %s\n", seq, paste(positions, collapse = ", ")))
      }
      
      if (length(multi_pos_seqs) > 10) {
        cat(sprintf("  ... and %d more\n", length(multi_pos_seqs) - 10))
      }
    } else {
      cat("No sequences found with multiple PTM positions\n")
    }
  }
  
  # 4. Analyze host switching by host pair
  analyze_host_pairs <- function(simple_results) {
    cat("\n=== HOST SWITCHING PATTERNS BY HOST PAIR ===\n")
    
    host_switch_data <- simple_results$host_switch_data
    switches <- host_switch_data[host_switch_data$is_switch, ]
    
    if (nrow(switches) > 0) {
      # Create host pair table
      host_pairs <- data.frame()
      
      for (i in 1:nrow(switches)) {
        from_host <- switches$seq_host[i]
        to_hosts <- unlist(strsplit(switches$sister_hosts[i], ";"))
        
        for (to_host in to_hosts) {
          host_pairs <- rbind(host_pairs, data.frame(
            From = from_host,
            To = trimws(to_host),
            stringsAsFactors = FALSE
          ))
        }
      }
      
      # Count host pair frequencies
      pair_counts <- table(paste(host_pairs$From, "→", host_pairs$To))
      pair_counts <- sort(pair_counts, decreasing = TRUE)
      
      cat("Most common host switching patterns:\n")
      for (i in 1:min(10, length(pair_counts))) {
        cat(sprintf("  %s: %d occurrences\n", names(pair_counts)[i], pair_counts[i]))
      }
    } else {
      cat("No host switches detected\n")
    }
  }
  
  # 5. Create a summary report
  create_host_switch_summary <- function(simple_results, output_file = "host_switch_summary.txt") {
    sink(output_file)
    
    cat("HOST SWITCHING ANALYSIS SUMMARY REPORT\n")
    cat("=====================================\n")
    cat(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")
    
    # Position summary
    cat("POSITION SUMMARY:\n")
    print(simple_results$position_summary)
    
    # Overall statistics
    host_switch_data <- simple_results$host_switch_data
    cat(sprintf("\n\nOVERALL STATISTICS:\n"))
    cat(sprintf("Total sequences analyzed: %d\n", nrow(host_switch_data)))
    cat(sprintf("Sequences with host switches: %d (%.1f%%)\n",
                sum(host_switch_data$is_switch),
                100 * sum(host_switch_data$is_switch) / nrow(host_switch_data)))
    
    # Host switch examples
    if (sum(host_switch_data$is_switch) > 0) {
      cat("\n\nHOST SWITCH EXAMPLES:\n")
      switches <- host_switch_data[host_switch_data$is_switch, ]
      for (i in 1:min(20, nrow(switches))) {
        cat(sprintf("%s (%s) - sisters in: %s\n",
                    switches$sequence[i], switches$seq_host[i], switches$sister_hosts[i]))
      }
    }
    
    sink()
    cat(sprintf("\nSummary report saved to: %s\n", output_file))
  }
  
  # ==============================================================================
  # USAGE EXAMPLE
  # ==============================================================================
  
  # Assuming you have already run:
  # simple_results <- simplified_host_switching_analysis(...)
  
  # Load necessary data if not already loaded
  tree <- read.nexus("HA_full_tree.nex")
  comprehensive_info <- get_comprehensive_info(tree$tip.label, "HA_full_id_mapping.json")
  
  # Run all additional analyses
  analyze_position_host_diversity(simple_results, tree, comprehensive_info)
  compare_position_types(simple_results)
  find_multi_position_sequences(simple_results)
  analyze_host_pairs(simple_results)
  create_host_switch_summary(simple_results, "HA_host_switch_summary.txt")
  
  # Create additional visualization
  library(ggplot2)
  library(dplyr)
  
  # Visualize host switching by position
  plot_host_switches_by_position <- function(simple_results) {
    host_switch_data <- simple_results$host_switch_data
    position_data <- simple_results$sequences_by_position
    
    # Map sequences to positions
    seq_to_pos <- data.frame()
    for (pos_name in names(position_data)) {
      sequences <- position_data[[pos_name]]
      seq_to_pos <- rbind(seq_to_pos, data.frame(
        sequence = sequences,
        position = pos_name,
        stringsAsFactors = FALSE
      ))
    }
    
    # Merge with switch data
    plot_data <- merge(host_switch_data[, c("sequence", "is_switch")], 
                       seq_to_pos, by = "sequence", all.x = TRUE)
    
    # Summarize by position
    summary_data <- plot_data %>%
      group_by(position) %>%
      summarise(
        total = n(),
        switches = sum(is_switch),
        percent_switch = 100 * switches / total
      ) %>%
      filter(!is.na(position))
    
    # Create plot
    p <- ggplot(summary_data, aes(x = position, y = percent_switch, fill = position)) +
      geom_bar(stat = "identity") +
      geom_text(aes(label = sprintf("%d/%d", switches, total)), 
                vjust = -0.5, size = 3) +
      labs(title = "Host Switching Rate by PTM Position",
           x = "Position",
           y = "Percentage of Sequences with Host Switches") +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1),
            legend.position = "none") +
      scale_fill_manual(values = c(rep("red", sum(grepl("_red", summary_data$position))),
                                   rep("blue", sum(grepl("_blue", summary_data$position)))))
    
    return(p)
  }
  
  # Generate the plot
  p <- plot_host_switches_by_position(simple_results)
  ggsave("host_switching_by_position.pdf", p, width = 10, height = 6)  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  # Deeper analysis of host switching patterns
  # Looking for more subtle signals in the data
  
  # 1. Analyze phylogenetic distance to host switches
  analyze_phylogenetic_context <- function(tree, comprehensive_info, simple_results) {
    cat("\n=== PHYLOGENETIC CONTEXT OF PTM SEQUENCES ===\n")
    
    host_switch_data <- simple_results$host_switch_data
    sequences_by_position <- simple_results$sequences_by_position
    
    # For each position, analyze the phylogenetic distribution
    for (pos_name in names(sequences_by_position)) {
      sequences <- sequences_by_position[[pos_name]]
      
      cat(sprintf("\n%s (%d sequences):\n", pos_name, length(sequences)))
      
      # Calculate pairwise distances between sequences
      if (length(sequences) > 1) {
        # Get indices of these sequences in the tree
        seq_indices <- which(tree$tip.label %in% sequences)
        
        if (length(seq_indices) > 1) {
          # Calculate cophenetic distances
          distances <- cophenetic(tree)[seq_indices, seq_indices]
          
          # Get host information
          hosts <- comprehensive_info$Host[match(sequences, comprehensive_info$tip_label)]
          
          # Calculate within-host vs between-host distances
          within_host_dist <- c()
          between_host_dist <- c()
          
          for (i in 1:(length(sequences)-1)) {
            for (j in (i+1):length(sequences)) {
              dist <- distances[i, j]
              if (!is.na(hosts[i]) && !is.na(hosts[j]) && 
                  hosts[i] != "Unknown" && hosts[j] != "Unknown") {
                if (hosts[i] == hosts[j]) {
                  within_host_dist <- c(within_host_dist, dist)
                } else {
                  between_host_dist <- c(between_host_dist, dist)
                }
              }
            }
          }
          
          if (length(within_host_dist) > 0 && length(between_host_dist) > 0) {
            cat(sprintf("  Mean within-host distance: %.4f\n", mean(within_host_dist)))
            cat(sprintf("  Mean between-host distance: %.4f\n", mean(between_host_dist)))
            
            # Test for difference
            if (length(within_host_dist) > 1 && length(between_host_dist) > 1) {
              test <- wilcox.test(within_host_dist, between_host_dist)
              cat(sprintf("  Wilcoxon test p-value: %.4f\n", test$p.value))
              
              if (test$p.value < 0.05) {
                cat("  ** Significant difference in phylogenetic distances **\n")
              }
            }
          }
        }
      }
    }
  }
  
  # 2. Look for host-associated clades with PTM sequences
  find_host_enriched_clades <- function(tree, comprehensive_info, sequences_of_interest) {
    cat("\n=== HOST-ENRICHED CLADES CONTAINING PTM SEQUENCES ===\n")
    
    # Find all nodes
    n_tips <- length(tree$tip.label)
    internal_nodes <- (n_tips + 1):(n_tips + tree$Nnode)
    
    clade_info <- data.frame()
    
    for (node in internal_nodes) {
      descendants <- Descendants(tree, node, type = "tips")[[1]]
      
      if (length(descendants) >= 5) {  # Only consider clades with at least 5 sequences
        desc_labels <- tree$tip.label[descendants]
        
        # Count PTM sequences in this clade
        ptm_count <- sum(desc_labels %in% sequences_of_interest)
        
        if (ptm_count > 0) {
          # Get host composition
          hosts <- comprehensive_info$Host[comprehensive_info$tip_label %in% desc_labels]
          hosts <- hosts[!is.na(hosts) & hosts != "Unknown"]
          
          if (length(hosts) > 0) {
            host_table <- table(hosts)
            dominant_host <- names(which.max(host_table))
            host_purity <- max(host_table) / sum(host_table)
            
            clade_info <- rbind(clade_info, data.frame(
              node = node,
              size = length(descendants),
              ptm_sequences = ptm_count,
              ptm_percent = 100 * ptm_count / length(descendants),
              dominant_host = dominant_host,
              host_purity = host_purity,
              n_hosts = length(unique(hosts)),
              stringsAsFactors = FALSE
            ))
          }
        }
      }
    }
    
    if (nrow(clade_info) > 0) {
      # Sort by PTM percentage
      clade_info <- clade_info[order(-clade_info$ptm_percent), ]
      
      cat("\nTop clades enriched for PTM sequences:\n")
      for (i in 1:min(10, nrow(clade_info))) {
        cat(sprintf("Clade %d: %d sequences (%.1f%% with PTMs), %s-dominated (%.1f%% purity), %d hosts\n",
                    clade_info$node[i],
                    clade_info$size[i],
                    clade_info$ptm_percent[i],
                    clade_info$dominant_host[i],
                    100 * clade_info$host_purity[i],
                    clade_info$n_hosts[i]))
      }
      
      # Look for pattern: are mixed-host clades more likely to have PTMs?
      multi_host_clades <- clade_info[clade_info$n_hosts > 1, ]
      single_host_clades <- clade_info[clade_info$n_hosts == 1, ]
      
      if (nrow(multi_host_clades) > 0 && nrow(single_host_clades) > 0) {
        cat(sprintf("\n\nMulti-host clades: mean %.1f%% PTM sequences\n", 
                    mean(multi_host_clades$ptm_percent)))
        cat(sprintf("Single-host clades: mean %.1f%% PTM sequences\n", 
                    mean(single_host_clades$ptm_percent)))
      }
    }
    
    return(clade_info)
  }
  
  # 3. Temporal analysis - are host switches associated with specific time periods?
  temporal_host_analysis <- function(comprehensive_info, simple_results) {
    cat("\n=== TEMPORAL PATTERNS IN HOST SWITCHING ===\n")
    
    host_switch_data <- simple_results$host_switch_data
    
    # Add year information
    switch_years <- data.frame()
    
    for (i in 1:nrow(host_switch_data)) {
      seq <- host_switch_data$sequence[i]
      year_info <- comprehensive_info$Year[comprehensive_info$tip_label == seq]
      
      if (!is.na(year_info) && year_info != "") {
        switch_years <- rbind(switch_years, data.frame(
          sequence = seq,
          year = as.numeric(year_info),
          is_switch = host_switch_data$is_switch[i],
          host = host_switch_data$seq_host[i],
          stringsAsFactors = FALSE
        ))
      }
    }
    
    if (nrow(switch_years) > 0) {
      # Group by decade
      switch_years$decade <- floor(switch_years$year / 10) * 10
      
      decade_summary <- switch_years %>%
        group_by(decade) %>%
        summarise(
          total = n(),
          switches = sum(is_switch),
          switch_rate = 100 * switches / total
        )
      
      cat("\nHost switching rates by decade:\n")
      print(decade_summary)
      
      # Look for temporal clustering of PTM sequences
      cat("\n\nTemporal distribution of PTM sequences:\n")
      year_range <- range(switch_years$year, na.rm = TRUE)
      cat(sprintf("Year range: %d - %d\n", year_range[1], year_range[2]))
      
      # Check if switches are temporally clustered
      if (sum(switch_years$is_switch) > 2) {
        switch_years_only <- switch_years$year[switch_years$is_switch]
        cat(sprintf("Years with host switches: %s\n", 
                    paste(sort(switch_years_only), collapse = ", ")))
      }
    }
  }
  
  # 4. Network analysis of host relationships
  create_host_network_analysis <- function(tree, comprehensive_info, sequences_of_interest) {
    cat("\n=== HOST RELATIONSHIP NETWORK ===\n")
    
    # Build a network of host relationships based on phylogenetic proximity
    host_connections <- data.frame()
    
    for (seq in sequences_of_interest) {
      seq_idx <- which(tree$tip.label == seq)
      if (length(seq_idx) == 0) next
      
      # Find closely related sequences (within a certain distance)
      all_distances <- cophenetic(tree)[seq_idx, ]
      close_relatives <- which(all_distances < quantile(all_distances, 0.05) & all_distances > 0)
      
      if (length(close_relatives) > 0) {
        seq_host <- comprehensive_info$Host[comprehensive_info$tip_label == seq]
        
        for (rel_idx in close_relatives) {
          rel_host <- comprehensive_info$Host[comprehensive_info$tip_label == tree$tip.label[rel_idx]]
          
          if (!is.na(seq_host) && !is.na(rel_host) && 
              seq_host != "Unknown" && rel_host != "Unknown") {
            host_connections <- rbind(host_connections, data.frame(
              from = seq_host,
              to = rel_host,
              distance = all_distances[rel_idx],
              stringsAsFactors = FALSE
            ))
          }
        }
      }
    }
    
    if (nrow(host_connections) > 0) {
      # Count connections
      connection_counts <- table(paste(host_connections$from, "↔", host_connections$to))
      connection_counts <- sort(connection_counts, decreasing = TRUE)
      
      cat("Most common host proximity relationships:\n")
      for (i in 1:min(10, length(connection_counts))) {
        cat(sprintf("  %s: %d connections\n", names(connection_counts)[i], connection_counts[i]))
      }
      
      # Calculate host connectivity
      all_hosts <- unique(c(host_connections$from, host_connections$to))
      host_connectivity <- data.frame()
      
      for (host in all_hosts) {
        n_connections <- sum(host_connections$from == host | host_connections$to == host)
        n_other_hosts <- length(unique(c(
          host_connections$to[host_connections$from == host],
          host_connections$from[host_connections$to == host]
        )))
        
        host_connectivity <- rbind(host_connectivity, data.frame(
          host = host,
          total_connections = n_connections,
          connected_hosts = n_other_hosts,
          stringsAsFactors = FALSE
        ))
      }
      
      host_connectivity <- host_connectivity[order(-host_connectivity$connected_hosts), ]
      
      cat("\n\nHost connectivity (PTM sequences only):\n")
      print(host_connectivity)
    }
  }
  
  # 5. Statistical comparison with random expectation
  compare_to_random_expectation <- function(tree, comprehensive_info, simple_results, n_permutations = 100) {
    cat("\n=== COMPARISON TO RANDOM EXPECTATION ===\n")
    
    actual_switches <- sum(simple_results$host_switch_data$is_switch)
    n_sequences <- nrow(simple_results$host_switch_data)
    
    # Perform permutation test
    random_switches <- numeric(n_permutations)
    
    for (i in 1:n_permutations) {
      # Randomly select same number of sequences
      random_seqs <- sample(tree$tip.label, n_sequences)
      
      # Calculate host switches for random set
      random_analysis <- simple_host_switching_analysis(tree, comprehensive_info, random_seqs)
      random_switches[i] <- sum(random_analysis$is_switch)
    }
    
    # Calculate p-value
    p_value <- sum(random_switches >= actual_switches) / n_permutations
    
    cat(sprintf("Observed host switches: %d\n", actual_switches))
    cat(sprintf("Expected (random): %.1f ± %.1f\n", mean(random_switches), sd(random_switches)))
    cat(sprintf("Permutation test p-value: %.3f\n", p_value))
    
    if (p_value < 0.05) {
      cat("** PTM sequences show significantly different host switching than random **\n")
    } else {
      cat("PTM sequences show similar host switching to random expectation\n")
    }
    
    return(list(
      observed = actual_switches,
      expected = mean(random_switches),
      p_value = p_value
    ))
  }
  
  # ==============================================================================
  # RUN ALL DEEPER ANALYSES
  # ==============================================================================
  
  # Make sure you have loaded the necessary data
  tree <- read.nexus("HA_full_tree.nex")
  comprehensive_info <- get_comprehensive_info(tree$tip.label, "HA_full_id_mapping.json")
  
  # Get all sequences of interest
  all_interest_sequences <- unique(unlist(simple_results$sequences_by_position))
  
  # Run the deeper analyses
  analyze_phylogenetic_context(tree, comprehensive_info, simple_results)
  clade_info <- find_host_enriched_clades(tree, comprehensive_info, all_interest_sequences)
  temporal_host_analysis(comprehensive_info, simple_results)
  create_host_network_analysis(tree, comprehensive_info, all_interest_sequences)
  random_comparison <- compare_to_random_expectation(tree, comprehensive_info, simple_results)
  
  # Create a comprehensive visualization
  library(ggplot2)
  library(gridExtra)
  
  # Plot 1: Phylogenetic distance comparison
  plot_phylogenetic_distances <- function(tree, comprehensive_info, sequences_by_position) {
    plot_list <- list()
    
    for (pos_name in names(sequences_by_position)) {
      sequences <- sequences_by_position[[pos_name]]
      
      if (length(sequences) > 1) {
        seq_indices <- which(tree$tip.label %in% sequences)
        
        if (length(seq_indices) > 1) {
          distances <- cophenetic(tree)[seq_indices, seq_indices]
          hosts <- comprehensive_info$Host[match(sequences, comprehensive_info$tip_label)]
          
          # Create distance data
          dist_data <- data.frame()
          for (i in 1:(length(sequences)-1)) {
            for (j in (i+1):length(sequences)) {
              if (!is.na(hosts[i]) && !is.na(hosts[j])) {
                dist_data <- rbind(dist_data, data.frame(
                  distance = distances[i, j],
                  type = ifelse(hosts[i] == hosts[j], "Within-host", "Between-host"),
                  stringsAsFactors = FALSE
                ))
              }
            }
          }
          
          if (nrow(dist_data) > 0) {
            p <- ggplot(dist_data, aes(x = type, y = distance, fill = type)) +
              geom_boxplot() +
              labs(title = pos_name, y = "Phylogenetic distance") +
              theme_minimal() +
              theme(legend.position = "none")
            
            plot_list[[pos_name]] <- p
          }
        }
      }
    }
    
    if (length(plot_list) > 0) {
      do.call(grid.arrange, c(plot_list, ncol = 2))
    }
  }
  
  # Generate the plot
  pdf("deeper_host_analysis_plots.pdf", width = 12, height = 10)
  plot_phylogenetic_distances(tree, comprehensive_info, simple_results$sequences_by_position)
  dev.off()
  
  cat("\n\n=== SUMMARY ===\n")
  cat("While the overall host switching rate is low (2.7%), these deeper analyses reveal:\n")
  cat("1. Whether PTM sequences are phylogenetically clustered within or between hosts\n")
  cat("2. If PTM sequences are enriched in mixed-host clades\n")
  cat("3. Any temporal patterns in host switching\n")
  cat("4. The network structure of host relationships for PTM sequences\n")
  cat("5. Whether the observed pattern differs from random expectation\n") 
  
  
  
  # Control Analysis: Compare phylogenetic patterns between PTM and non-PTM sequences
  # This tests whether the within/between host distance pattern is specific to PTM sites
  
  library(ape)
  library(dplyr)
  library(ggplot2)
  library(gridExtra)
  
  # Function to calculate within/between host distances for any set of sequences
  calculate_host_distances <- function(tree, comprehensive_info, sequences, label = "sequences") {
    # Get indices of these sequences in the tree
    seq_indices <- which(tree$tip.label %in% sequences)
    
    if (length(seq_indices) < 2) {
      return(NULL)
    }
    
    # Calculate cophenetic distances
    distances <- cophenetic(tree)[seq_indices, seq_indices]
    
    # Get host information
    hosts <- comprehensive_info$Host[match(tree$tip.label[seq_indices], comprehensive_info$tip_label)]
    
    # Calculate within-host vs between-host distances
    within_host_dist <- c()
    between_host_dist <- c()
    
    for (i in 1:(length(seq_indices)-1)) {
      for (j in (i+1):length(seq_indices)) {
        dist <- distances[i, j]
        if (!is.na(hosts[i]) && !is.na(hosts[j]) && 
            hosts[i] != "Unknown" && hosts[j] != "Unknown") {
          if (hosts[i] == hosts[j]) {
            within_host_dist <- c(within_host_dist, dist)
          } else {
            between_host_dist <- c(between_host_dist, dist)
          }
        }
      }
    }
    
    result <- list(
      label = label,
      n_sequences = length(sequences),
      within_host_dist = within_host_dist,
      between_host_dist = between_host_dist,
      mean_within = mean(within_host_dist, na.rm = TRUE),
      mean_between = mean(between_host_dist, na.rm = TRUE),
      ratio = mean(between_host_dist, na.rm = TRUE) / mean(within_host_dist, na.rm = TRUE)
    )
    
    # Test for difference if both have sufficient data
    if (length(within_host_dist) > 1 && length(between_host_dist) > 1) {
      result$wilcox_test <- wilcox.test(within_host_dist, between_host_dist)
      result$p_value <- result$wilcox_test$p.value
    } else {
      result$p_value <- NA
    }
    
    return(result)
  }
  
  # Main control analysis function
  perform_control_analysis <- function(tree, comprehensive_info, ptm_sequences, 
                                       n_random_sets = 10, set_size = NULL) {
    
    cat("=== CONTROL ANALYSIS: PHYLOGENETIC DISTANCE PATTERNS ===\n\n")
    
    # 1. Calculate distances for ALL sequences in the tree
    cat("1. BASELINE - All sequences in tree:\n")
    all_seq_results <- calculate_host_distances(tree, comprehensive_info, 
                                                tree$tip.label, "All sequences")
    
    if (!is.null(all_seq_results)) {
      cat(sprintf("   N = %d sequences\n", all_seq_results$n_sequences))
      cat(sprintf("   Mean within-host distance: %.4f\n", all_seq_results$mean_within))
      cat(sprintf("   Mean between-host distance: %.4f\n", all_seq_results$mean_between))
      cat(sprintf("   Ratio (between/within): %.2f\n", all_seq_results$ratio))
      if (!is.na(all_seq_results$p_value)) {
        cat(sprintf("   Wilcoxon p-value: %.4f %s\n", 
                    all_seq_results$p_value,
                    ifelse(all_seq_results$p_value < 0.05, "**", "")))
      }
    }
    
    # 2. Calculate distances for PTM sequences
    cat("\n2. PTM sequences:\n")
    ptm_results <- calculate_host_distances(tree, comprehensive_info, 
                                            ptm_sequences, "PTM sequences")
    
    if (!is.null(ptm_results)) {
      cat(sprintf("   N = %d sequences\n", ptm_results$n_sequences))
      cat(sprintf("   Mean within-host distance: %.4f\n", ptm_results$mean_within))
      cat(sprintf("   Mean between-host distance: %.4f\n", ptm_results$mean_between))
      cat(sprintf("   Ratio (between/within): %.2f\n", ptm_results$ratio))
      if (!is.na(ptm_results$p_value)) {
        cat(sprintf("   Wilcoxon p-value: %.4f %s\n", 
                    ptm_results$p_value,
                    ifelse(ptm_results$p_value < 0.05, "**", "")))
      }
    }
    
    # 3. Calculate distances for non-PTM sequences
    cat("\n3. Non-PTM sequences:\n")
    non_ptm_sequences <- setdiff(tree$tip.label, ptm_sequences)
    non_ptm_results <- calculate_host_distances(tree, comprehensive_info, 
                                                non_ptm_sequences, "Non-PTM sequences")
    
    if (!is.null(non_ptm_results)) {
      cat(sprintf("   N = %d sequences\n", non_ptm_results$n_sequences))
      cat(sprintf("   Mean within-host distance: %.4f\n", non_ptm_results$mean_within))
      cat(sprintf("   Mean between-host distance: %.4f\n", non_ptm_results$mean_between))
      cat(sprintf("   Ratio (between/within): %.2f\n", non_ptm_results$ratio))
      if (!is.na(non_ptm_results$p_value)) {
        cat(sprintf("   Wilcoxon p-value: %.4f %s\n", 
                    non_ptm_results$p_value,
                    ifelse(non_ptm_results$p_value < 0.05, "**", "")))
      }
    }
    
    # 4. Random sampling control
    cat("\n4. Random sample controls:\n")
    if (is.null(set_size)) {
      set_size <- length(ptm_sequences)
    }
    
    random_results <- list()
    for (i in 1:n_random_sets) {
      random_seqs <- sample(tree$tip.label, min(set_size, length(tree$tip.label)))
      random_results[[i]] <- calculate_host_distances(tree, comprehensive_info, 
                                                      random_seqs, paste("Random set", i))
    }
    
    # Summarize random results
    random_ratios <- sapply(random_results, function(x) ifelse(is.null(x), NA, x$ratio))
    random_ratios <- random_ratios[!is.na(random_ratios)]
    
    cat(sprintf("   Generated %d random sets of %d sequences each\n", 
                n_random_sets, set_size))
    cat(sprintf("   Mean ratio (between/within): %.2f ± %.2f\n", 
                mean(random_ratios), sd(random_ratios)))
    cat(sprintf("   Range: %.2f - %.2f\n", 
                min(random_ratios), max(random_ratios)))
    
    # Store all results
    all_results <- list(
      all_sequences = all_seq_results,
      ptm_sequences = ptm_results,
      non_ptm_sequences = non_ptm_results,
      random_samples = random_results
    )
    
    return(all_results)
  }
  
  # Function to create comparison plots
  create_control_comparison_plots <- function(control_results) {
    # Extract data for plotting
    plot_data <- data.frame()
    
    # Add main categories
    categories <- c("all_sequences", "ptm_sequences", "non_ptm_sequences")
    labels <- c("All", "PTM", "Non-PTM")
    
    for (i in seq_along(categories)) {
      cat_data <- control_results[[categories[i]]]
      if (!is.null(cat_data) && length(cat_data$within_host_dist) > 0) {
        # Within-host distances
        plot_data <- rbind(plot_data, data.frame(
          Category = labels[i],
          Type = "Within-host",
          Distance = cat_data$within_host_dist,
          stringsAsFactors = FALSE
        ))
        
        # Between-host distances
        if (length(cat_data$between_host_dist) > 0) {
          plot_data <- rbind(plot_data, data.frame(
            Category = labels[i],
            Type = "Between-host",
            Distance = cat_data$between_host_dist,
            stringsAsFactors = FALSE
          ))
        }
      }
    }
    
    # Add random samples
    for (i in seq_along(control_results$random_samples)) {
      rand_data <- control_results$random_samples[[i]]
      if (!is.null(rand_data) && length(rand_data$within_host_dist) > 0) {
        plot_data <- rbind(plot_data, data.frame(
          Category = "Random",
          Type = "Within-host",
          Distance = rand_data$within_host_dist,
          stringsAsFactors = FALSE
        ))
        
        if (length(rand_data$between_host_dist) > 0) {
          plot_data <- rbind(plot_data, data.frame(
            Category = "Random",
            Type = "Between-host",
            Distance = rand_data$between_host_dist,
            stringsAsFactors = FALSE
          ))
        }
      }
    }
    
    # Create boxplot
    p1 <- ggplot(plot_data, aes(x = Category, y = Distance, fill = Type)) +
      geom_boxplot(outlier.shape = NA) +
      scale_fill_manual(values = c("Within-host" = "lightblue", "Between-host" = "coral")) +
      labs(title = "Phylogenetic Distance Patterns: PTM vs Control",
           y = "Phylogenetic Distance",
           x = "Sequence Category") +
      theme_minimal() +
      theme(legend.position = "bottom")
    
    # Create ratio comparison
    ratio_data <- data.frame()
    
    # Calculate ratios for each category
    for (i in seq_along(categories)) {
      cat_data <- control_results[[categories[i]]]
      if (!is.null(cat_data)) {
        ratio_data <- rbind(ratio_data, data.frame(
          Category = labels[i],
          Ratio = cat_data$ratio,
          Type = "Observed",
          stringsAsFactors = FALSE
        ))
      }
    }
    
    # Add random ratios
    for (rand_data in control_results$random_samples) {
      if (!is.null(rand_data)) {
        ratio_data <- rbind(ratio_data, data.frame(
          Category = "Random",
          Ratio = rand_data$ratio,
          Type = "Random",
          stringsAsFactors = FALSE
        ))
      }
    }
    
    p2 <- ggplot(ratio_data, aes(x = Category, y = Ratio, fill = Category)) +
      geom_bar(data = ratio_data[ratio_data$Type == "Observed",], 
               stat = "identity", alpha = 0.7) +
      geom_point(data = ratio_data[ratio_data$Type == "Random",], 
                 position = position_jitter(width = 0.2), alpha = 0.5) +
      geom_hline(yintercept = 1, linetype = "dashed", color = "red") +
      labs(title = "Between/Within Host Distance Ratios",
           y = "Ratio (Between/Within)",
           x = "Sequence Category") +
      theme_minimal() +
      theme(legend.position = "none")
    
    return(list(p1 = p1, p2 = p2))
  }
  
  # Function to test if PTM pattern is significantly different from background
  test_ptm_vs_background <- function(control_results) {
    cat("\n\n=== STATISTICAL COMPARISON: PTM vs BACKGROUND ===\n")
    
    # Compare PTM ratio to random ratios
    ptm_ratio <- control_results$ptm_sequences$ratio
    random_ratios <- sapply(control_results$random_samples, 
                            function(x) ifelse(is.null(x), NA, x$ratio))
    random_ratios <- random_ratios[!is.na(random_ratios)]
    
    # Calculate empirical p-value
    p_value <- sum(random_ratios >= ptm_ratio) / length(random_ratios)
    
    cat(sprintf("PTM sequences ratio: %.2f\n", ptm_ratio))
    cat(sprintf("Random samples mean ratio: %.2f ± %.2f\n", 
                mean(random_ratios), sd(random_ratios)))
    cat(sprintf("Empirical p-value: %.3f\n", p_value))
    
    if (p_value < 0.05) {
      cat("** PTM sequences show significantly different pattern than random **\n")
    } else {
      cat("PTM sequences show similar pattern to random background\n")
    }
    
    # Compare to non-PTM sequences
    if (!is.null(control_results$non_ptm_sequences)) {
      cat("\n\nDirect comparison PTM vs non-PTM:\n")
      cat(sprintf("PTM ratio: %.2f\n", ptm_ratio))
      cat(sprintf("Non-PTM ratio: %.2f\n", control_results$non_ptm_sequences$ratio))
      
      # Can also do a permutation test here
      all_within <- c(control_results$ptm_sequences$within_host_dist,
                      control_results$non_ptm_sequences$within_host_dist)
      all_between <- c(control_results$ptm_sequences$between_host_dist,
                       control_results$non_ptm_sequences$between_host_dist)
      
      n_ptm_within <- length(control_results$ptm_sequences$within_host_dist)
      n_ptm_between <- length(control_results$ptm_sequences$between_host_dist)
      
      # Permutation test
      n_perm <- 1000
      perm_diffs <- numeric(n_perm)
      
      obs_diff <- ptm_ratio - control_results$non_ptm_sequences$ratio
      
      for (i in 1:n_perm) {
        # Randomly assign distances to PTM or non-PTM
        perm_within <- sample(all_within)
        perm_between <- sample(all_between)
        
        perm_ptm_ratio <- mean(perm_between[1:n_ptm_between]) / 
          mean(perm_within[1:n_ptm_within])
        perm_nonptm_ratio <- mean(perm_between[-(1:n_ptm_between)]) / 
          mean(perm_within[-(1:n_ptm_within)])
        
        perm_diffs[i] <- perm_ptm_ratio - perm_nonptm_ratio
      }
      
      perm_p <- sum(abs(perm_diffs) >= abs(obs_diff)) / n_perm
      cat(sprintf("Permutation test p-value: %.3f\n", perm_p))
    }
  }
  
  # ==============================================================================
  # RUN THE CONTROL ANALYSIS
  # ==============================================================================
  
  # Load your data
  tree <- read.nexus("HA_full_tree.nex")
  comprehensive_info <- get_comprehensive_info(tree$tip.label, "HA_full_id_mapping.json")
  
  # Get all PTM sequences
  all_ptm_sequences <- unique(unlist(simple_results$sequences_by_position))
  
  # Run the control analysis
  control_results <- perform_control_analysis(tree, comprehensive_info, 
                                              all_ptm_sequences, 
                                              n_random_sets = 20)
  
  # Test if PTM pattern is different from background
  test_ptm_vs_background(control_results)
  
  # Create visualizations
  plots <- create_control_comparison_plots(control_results)
  
  # Save plots
  pdf("control_analysis_plots.pdf", width = 12, height = 6)
  grid.arrange(plots$p1, plots$p2, ncol = 2)
  dev.off()
  
  # Additional analysis: Check individual positions against background
  cat("\n\n=== INDIVIDUAL POSITION ANALYSIS ===\n")
  for (pos_name in names(simple_results$sequences_by_position)) {
    sequences <- simple_results$sequences_by_position[[pos_name]]
    
    if (length(sequences) > 5) {  # Only analyze positions with sufficient sequences
      cat(sprintf("\n%s:\n", pos_name))
      
      pos_results <- calculate_host_distances(tree, comprehensive_info, 
                                              sequences, pos_name)
      
      if (!is.null(pos_results)) {
        cat(sprintf("  Ratio: %.2f", pos_results$ratio))
        
        # Compare to all sequences ratio
        if (!is.null(control_results$all_sequences)) {
          diff <- pos_results$ratio - control_results$all_sequences$ratio
          cat(sprintf(" (%.1f%% %s than background)\n", 
                      abs(diff) / control_results$all_sequences$ratio * 100,
                      ifelse(diff > 0, "higher", "lower")))
        } else {
          cat("\n")
        }
      }
    }
  }
  
  cat("\n\n=== SUMMARY ===\n")
  cat("This control analysis reveals whether the within/between host clustering\n")
  cat("is specific to PTM sequences or a general feature of your phylogenetic tree.\n")
  cat("If PTM sequences show a significantly different pattern than random samples,\n")
  cat("it suggests the PTMs are specifically associated with host adaptation.\n")
  
  