# ==== Combined Protein Evolution Analysis Across All Segments ====
# This script combines the robust mapping and ancestral reconstruction methods 
# from the second script with the comprehensive downstream analysis from the first script.

# ==== Load Required Libraries ====
library(ape)
library(phangorn)
library(readr)
library(dplyr)
library(stringr)
library(ggplot2)
library(tidyr)
library(jsonlite)
library(fs)
library(MASS)  # For robust statistics

# Try to load visualization libraries, install if needed
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}
if (!requireNamespace("ggtree", quietly = TRUE)) {
  BiocManager::install("ggtree")
}
if (!requireNamespace("RColorBrewer", quietly = TRUE)) {
  install.packages("RColorBrewer")
}
if (!requireNamespace("aplot", quietly = TRUE)) {
  install.packages("aplot")
}

# ==== Global Settings ====
burnin_fraction <- 0.25
base_dir <- "processed_segments_fixed_final"
main_output_dir <- "site_specific_analysis_results"
debug_output_dir <- file.path(main_output_dir, "debug_output")
dir.create(main_output_dir, showWarnings = FALSE)
dir.create(debug_output_dir, showWarnings = FALSE)

# Debug log file
debug_log <- file.path(debug_output_dir, "processing_debug.log")
file.create(debug_log)
log_message <- function(...) {
  message <- paste0(...)
  cat(message, "\n")
  cat(message, "\n", file=debug_log, append=TRUE)
}

log_message("===== STARTING COMBINED ANALYSIS =====")
log_message("Base directory:", base_dir)
log_message("Output directory:", main_output_dir)

# ==== Prepare Overall Summary File ====
summary_all_file <- file.path(main_output_dir, "all_proteins_summary.txt")
file.create(summary_all_file)  # Create an empty summary file

# ==== Find All Segment Directories ====
log_message("Searching for segment directories...")
seg_dirs <- list.dirs(base_dir, recursive = FALSE)
seg_dirs <- seg_dirs[grepl("^seg", basename(seg_dirs))]
if (length(seg_dirs) == 0) {
  stop("No segment directories found. Check that base_dir is correct.")
}
log_message("Found", length(seg_dirs), "segment directories:", paste(basename(seg_dirs), collapse = ", "), "\n")

# ==== Flexible File Matching Function ====
find_matching_files <- function(tree_file) {
  dir_path <- dirname(tree_file)
  all_files <- list.files(dir_path, full.names = TRUE)
  
  # Extract protein name from the tree file.
  tree_basename <- basename(tree_file)
  raw_name <- sub("_tree\\.nex$", "", tree_basename, ignore.case = TRUE)
  protein_name <- sub("_protein$", "", raw_name, ignore.case = TRUE)
  log_message("Extracted protein name:", protein_name, "from file:", tree_basename)
  
  # Build regex patterns for related files
  branch_regex  <- paste0("^", protein_name, "(_protein)?_branch_rates\\.log$")
  site_regex    <- paste0("^", protein_name, "(_protein)?_site_rates\\.log$")
  mapping_regex <- paste0("^", protein_name, "(_protein)?_id_mapping\\.json$")
  
  branch_log_file <- NA
  if (length(grep(branch_regex, basename(all_files), ignore.case = TRUE)) > 0) {
    branch_candidates <- all_files[grep(branch_regex, basename(all_files), ignore.case = TRUE)]
    branch_log_file <- branch_candidates[1]
  }
  
  site_rate_file <- NA
  if (length(grep(site_regex, basename(all_files), ignore.case = TRUE)) > 0) {
    site_candidates <- all_files[grep(site_regex, basename(all_files), ignore.case = TRUE)]
    site_rate_file <- site_candidates[1]
  }
  
  mapping_file <- NA
  if (length(grep(mapping_regex, basename(all_files), ignore.case = TRUE)) > 0) {
    mapping_candidates <- all_files[grep(mapping_regex, basename(all_files), ignore.case = TRUE)]
    mapping_file <- mapping_candidates[1]
  }
  
  log_message("Found matching files for", protein_name, ":")
  log_message("  Tree file:         ", tree_file)
  log_message("  Branch rate file:  ", ifelse(is.na(branch_log_file), "Not found", branch_log_file))
  log_message("  Site rate file:    ", ifelse(is.na(site_rate_file), "Not found", site_rate_file))
  log_message("  Mapping file:      ", ifelse(is.na(mapping_file), "Not found", mapping_file))
  
  list(
    protein_name = protein_name,
    tree_file = tree_file,
    branch_log_file = branch_log_file,
    site_rate_file = site_rate_file,
    mapping_file = mapping_file
  )
}

# ==== Improved Label Matching Algorithm ====
match_strains_to_meta <- function(tip_labels, meta_data, mapping) {
  # Initialize result dataframe
  result <- data.frame(
    label = tip_labels,
    Final_Accession = NA_character_,
    Matched_Strain = NA_character_,
    match_method = NA_character_,  # Track which method made the match
    stringsAsFactors = FALSE
  )
  
  # Log the number of tip labels and metadata entries
  log_message("\n===== MATCHING DEBUG =====")
  log_message("Total tip labels to match: ", length(tip_labels))
  log_message("Total metadata entries: ", nrow(meta_data))
  
  # STEP 1: Direct mapping via JSON (including both key->value and value mapping)
  if (length(mapping) > 0) {
    log_message("\nStep 1: Direct mapping via JSON file")
    for (i in 1:nrow(result)) {
      curr_label <- result$label[i]
      
      # Case 1: Check if the tip label is a value in the mapping (reverse lookup)
      if (curr_label %in% mapping) {
        # Find all keys that map to this value
        matching_keys <- names(mapping)[mapping == curr_label]
        if (length(matching_keys) > 0) {
          # Extract accession from the key if possible
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
          
          # If we didn't find a match above, just use the tip label as the accession
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
      
      # Case 2: Check if the tip label is a key in the mapping
      if (is.na(result$Final_Accession[i]) && curr_label %in% names(mapping)) {
        mapped_value <- mapping[[curr_label]]
        
        # Try different patterns depending on the mapped value
        if (grepl("^H[1-9]N[1-9]__", mapped_value)) {
          # Handle case like "H1N1__AUO38178_"
          clean_value <- gsub("^(H[1-9]N[1-9]__)", "", mapped_value)
          clean_value <- gsub("_$", "", clean_value)
          idx <- match(clean_value, meta_data$Accession)
          if (!is.na(idx)) {
            result$Final_Accession[i] <- clean_value
            result$Matched_Strain[i] <- meta_data$Strain[idx]
            result$match_method[i] <- "json_key_h1n1_pattern"
          } else {
            # If we can't find it in metadata, still use it
            result$Final_Accession[i] <- clean_value
            result$match_method[i] <- "json_key_h1n1_pattern_no_meta"
          }
        } else {
          # Handle direct accession like "A01134331"
          idx <- match(mapped_value, meta_data$Accession)
          if (!is.na(idx)) {
            result$Final_Accession[i] <- mapped_value
            result$Matched_Strain[i] <- meta_data$Strain[idx]
            result$match_method[i] <- "json_key_direct"
          } else {
            # If we can't find it in metadata, still use it
            result$Final_Accession[i] <- mapped_value
            result$match_method[i] <- "json_key_direct_no_meta"
          }
        }
      }
    }
    
    json_matches <- sum(!is.na(result$match_method) & 
                          grepl("^json_", result$match_method), na.rm=TRUE)
    log_message("Total matches via JSON mapping: ", json_matches, " (", round(json_matches/nrow(result)*100, 1), "%)")
  }
  
  # Step 2: Direct accession matching for remaining unmatched
  rem <- which(is.na(result$Final_Accession))
  log_message("\nStep 2: Direct accession matching for ", length(rem), " remaining unmatched labels")
  
  if (length(rem) > 0) {
    acc_pat <- "^[A-Z]{2}[_]?[0-9]{6,}$"
    hit <- which(is.na(result$Final_Accession) & 
                   grepl(acc_pat, result$label) & 
                   result$label %in% meta_data$Accession)
    
    if (length(hit) > 0) {
      result$Final_Accession[hit] <- result$label[hit]
      result$Matched_Strain[hit] <- meta_data$Strain[match(result$label[hit], meta_data$Accession)]
      result$match_method[hit] <- "direct_accession"
      
      log_message("  Direct accession matches: ", length(hit))
    } else {
      log_message("  No direct accession matches found")
    }
  }
  
  # Step 3: HxNy__ACC_ pattern
  rem <- which(is.na(result$Final_Accession))
  log_message("\nStep 3: HxNy__ACC_ pattern matching for ", length(rem), " remaining unmatched labels")
  
  if (length(rem) > 0) {
    pattern_matches <- 0
    for (i in rem) {
      if (grepl("^H\\d+N\\d+__([A-Z]{2,}\\d{5,})_$", result$label[i])) {
        acc <- sub("^H\\d+N\\d+__([A-Z]{2,}\\d{5,})_$", "\\1", result$label[i])
        idx <- match(acc, meta_data$Accession)
        if (!is.na(idx)) {
          result$Final_Accession[i] <- acc
          result$Matched_Strain[i] <- meta_data$Strain[idx]
          result$match_method[i] <- "HxNy_pattern"
          pattern_matches <- pattern_matches + 1
        }
      }
    }
    log_message("  HxNy pattern matches: ", pattern_matches)
  }
  
  # Step 4: ARG via JSON mapping (for any remaining ARG-like patterns)
  rem <- which(is.na(result$Final_Accession))
  log_message("\nStep 4: ARG pattern matching for ", length(rem), " remaining unmatched labels")
  
  if (length(rem) > 0 && length(mapping) > 0) {
    arg_matches <- 0
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
            arg_matches <- arg_matches + 1
          }
        }
      }
    }
    log_message("  ARG pattern matches: ", arg_matches)
  }
  
  # Step 5: Clean up and try no underscore suffix
  rem <- which(is.na(result$Final_Accession))
  log_message("\nStep 5: No underscore suffix matching for ", length(rem), " remaining unmatched labels")
  
  if (length(rem) > 0) {
    no_underscore_matches <- 0
    # Prepare accessions without trailing underscores
    label_no_underscore <- gsub("_$", "", result$label[rem])
    meta_no_underscore <- gsub("_$", "", meta_data$Accession)
    
    for (i in seq_along(rem)) {
      idx <- match(label_no_underscore[i], meta_no_underscore)
      if (!is.na(idx)) {
        result$Final_Accession[rem[i]] <- meta_data$Accession[idx]
        result$Matched_Strain[rem[i]] <- meta_data$Strain[idx]
        result$match_method[rem[i]] <- "no_underscore"
        no_underscore_matches <- no_underscore_matches + 1
      }
    }
    log_message("  No underscore matches: ", no_underscore_matches)
  }
  
  # Matching results summary
  log_message("\nMatching results summary:")
  log_message("- Total tips:", nrow(result))
  log_message("- Successfully matched:", sum(!is.na(result$Final_Accession)))
  log_message("- Unmatched:", sum(is.na(result$Final_Accession)))
  
  # Summarize by match method
  for (method in unique(result$match_method[!is.na(result$match_method)])) {
    count <- sum(result$match_method == method, na.rm=TRUE)
    log_message("- Matched by ", method, ": ", count, " (", round(count/nrow(result)*100, 1), "%)")
  }
  
  return(result)
}

# ==== Function to fix the tree for ACE compatibility ====
fix_tree_for_ace <- function(tree) {
  # Check for problem branch lengths
  zero_branches <- sum(tree$edge.length <= 1e-6, na.rm=TRUE)
  na_branches <- sum(is.na(tree$edge.length))
  
  log_message("Tree has", zero_branches, "zero-length branches and", na_branches, "NA branches")
  
  if (zero_branches > 0 || na_branches > 0) {
    log_message("Setting small positive values for zero/NA branch lengths")
    # Replace zero/NA branch lengths with small positive values - this is critical for ACE
    tree$edge.length[tree$edge.length <= 1e-6 | is.na(tree$edge.length)] <- 0.01
  }
  
  # Check for negative branch lengths
  neg_branches <- sum(tree$edge.length < 0, na.rm=TRUE)
  if (neg_branches > 0) {
    log_message("Fixing", neg_branches, "negative branch lengths")
    tree$edge.length[tree$edge.length < 0] <- 0.01
  }
  
  # Make sure tree is binary
  if (!is.binary.tree(tree)) {
    log_message("Making tree binary")
    tree <- multi2di(tree)
  }
  
  # Make sure tree is rooted
  if (!is.rooted(tree)) {
    log_message("Rooting tree")
    tree <- root(tree, outgroup=1, resolve.root=TRUE)
  }
  
  # Check the tree
  log_message("Checking fixed tree:")
  log_message("  Zero branches after fix:", sum(tree$edge.length <= 1e-6))
  log_message("  Negative branches after fix:", sum(tree$edge.length < 0))
  log_message("  NA branches after fix:", sum(is.na(tree$edge.length)))
  log_message("  Min branch length:", min(tree$edge.length))
  log_message("  Max branch length:", max(tree$edge.length))
  
  return(tree)
}

# ==== Function that only ever uses ACE, with robust error handling ====
run_ace_with_fixes <- function(host_states, tree) {
  # First, fix any tree issues
  tree <- fix_tree_for_ace(tree)
  
  # Fix host state issues - make sure all tips have states
  missing_states <- setdiff(tree$tip.label, names(host_states))
  if (length(missing_states) > 0) {
    log_message("Found", length(missing_states), "tips without host states")
    for (tip in missing_states) {
      host_states[tip] <- "Unknown"
    }
  }
  
  # Make sure it's a factor with all levels present
  if (!is.factor(host_states)) {
    log_message("Converting host states to factor")
    host_states <- factor(host_states)
  }
  
  # Make sure "Unknown" is included in the levels
  if (!"Unknown" %in% levels(host_states)) {
    log_message("Adding 'Unknown' to host state levels")
    host_states <- factor(host_states, levels = c(levels(host_states), "Unknown"))
  }
  
  # Try different ACE variations, while always using ACE
  tryCatch({
    log_message("Attempting ACE with default parameters")
    ace_result <- ace(host_states, tree, type="discrete")
    log_message("ACE successful!")
    return(validate_ace_result(ace_result, tree, host_states))
  }, error = function(e) {
    log_message("First ACE attempt failed:", e$message)
    
    # Try with marginal reconstruction
    tryCatch({
      log_message("Trying ACE with marginal=TRUE")
      ace_result <- ace(host_states, tree, type="discrete", marginal=TRUE)
      log_message("ACE with marginal=TRUE successful!")
      return(validate_ace_result(ace_result, tree, host_states))
    }, error = function(e2) {
      log_message("Second ACE attempt failed:", e2$message)
      
      # Try with no confidence intervals
      tryCatch({
        log_message("Trying ACE with CI=FALSE")
        ace_result <- ace(host_states, tree, type="discrete", CI=FALSE)
        log_message("ACE with CI=FALSE successful!")
        return(validate_ace_result(ace_result, tree, host_states))
      }, error = function(e3) {
        log_message("Third ACE attempt failed:", e3$message)
        
        # Make more aggressive fixes to the tree
        log_message("Making more aggressive tree fixes")
        tree$edge.length <- pmax(tree$edge.length, 0.05)
        
        # Try with all optimizations
        tryCatch({
          log_message("Final ACE attempt with heavily optimized tree")
          ace_result <- ace(host_states, tree, type="discrete", CI=FALSE, marginal=TRUE)
          log_message("Final ACE attempt successful!")
          return(validate_ace_result(ace_result, tree, host_states))
        }, error = function(e4) {
          # Last-ditch effort - create a proper dummy result
          log_message("All ACE attempts failed:", e4$message)
          log_message("Creating dummy ACE result with most common state.")
          return(create_dummy_ace_result(tree, host_states))
        })
      })
    })
  })
}

# Helper function to validate and fix ACE results
validate_ace_result <- function(ace_result, tree, host_states) {
  n_internal <- tree$Nnode
  
  # Check if the result is properly structured
  if (is.null(ace_result) || !is.list(ace_result) || is.null(ace_result$lik.anc)) {
    log_message("Invalid ACE result: NULL or missing lik.anc")
    return(create_dummy_ace_result(tree, host_states))
  }
  
  if (!is.matrix(ace_result$lik.anc)) {
    log_message("Invalid ACE result: lik.anc is not a matrix")
    return(create_dummy_ace_result(tree, host_states))
  }
  
  # Check column names
  if (is.null(colnames(ace_result$lik.anc))) {
    log_message("ACE result missing column names, adding them")
    colnames(ace_result$lik.anc) <- levels(host_states)
  }
  
  # Check and fix dimensions
  if (nrow(ace_result$lik.anc) != n_internal) {
    log_message("ACE result has incorrect number of rows:", nrow(ace_result$lik.anc), 
                "expected:", n_internal)
    
    if (nrow(ace_result$lik.anc) > n_internal) {
      # Truncate
      log_message("Truncating ACE result rows")
      ace_result$lik.anc <- ace_result$lik.anc[1:n_internal, , drop = FALSE]
    } else {
      # Extend with neutral probability
      log_message("Extending ACE result rows")
      n_missing <- n_internal - nrow(ace_result$lik.anc)
      n_states <- ncol(ace_result$lik.anc)
      
      # Create extension with all states having equal probability
      extension <- matrix(1/n_states, nrow = n_missing, ncol = n_states,
                          dimnames = list(NULL, colnames(ace_result$lik.anc)))
      
      # If "Unknown" is a state, give it higher probability
      unknown_col <- match("Unknown", colnames(ace_result$lik.anc))
      if (!is.na(unknown_col)) {
        extension[, ] <- 0.1 / (n_states - 1)  # Small probability for other states
        extension[, unknown_col] <- 0.9        # High probability for Unknown
      }
      
      ace_result$lik.anc <- rbind(ace_result$lik.anc, extension)
    }
  }
  
  return(ace_result)
}

# Helper function to create a dummy ACE result when all else fails
create_dummy_ace_result <- function(tree, host_states) {
  log_message("Creating dummy ACE result")
  n_internal <- tree$Nnode
  state_levels <- levels(host_states)
  n_states <- length(state_levels)
  
  # Create a matrix with the right dimensions
  lik_anc <- matrix(0, nrow = n_internal, ncol = n_states,
                    dimnames = list(NULL, state_levels))
  
  # Find most common state among tips (excluding 'Unknown')
  tip_states <- as.character(host_states)
  tip_states <- tip_states[tip_states != "Unknown"]
  
  if (length(tip_states) > 0) {
    most_common <- names(sort(table(tip_states), decreasing = TRUE))[1]
    most_common_idx <- match(most_common, state_levels)
    if (!is.na(most_common_idx)) {
      lik_anc[, most_common_idx] <- 0.6  # Give most common state higher probability
    }
  }
  
  # Add some probability for Unknown
  unknown_idx <- match("Unknown", state_levels)
  if (!is.na(unknown_idx)) {
    lik_anc[, unknown_idx] <- 0.4
  }
  
  # Add small probabilities to other states to avoid all zeros
  for (i in 1:n_states) {
    if (lik_anc[1, i] == 0) {
      lik_anc[, i] <- 0.01
    }
  }
  
  # Normalize to ensure rows sum to 1
  row_sums <- rowSums(lik_anc)
  for (i in 1:n_internal) {
    if (row_sums[i] > 0) {
      lik_anc[i, ] <- lik_anc[i, ] / row_sums[i]
    } else {
      # If row sum is 0, distribute evenly with higher weight to Unknown
      if (!is.na(unknown_idx)) {
        lik_anc[i, ] <- 0.1 / (n_states - 1)
        lik_anc[i, unknown_idx] <- 0.9
      } else {
        lik_anc[i, ] <- 1 / n_states
      }
    }
  }
  
  # Create ACE result structure
  result <- list(
    lik.anc = lik_anc,
    loglik = 0,
    rates = rep(1, choose(n_states, 2)),
    se = NA
  )
  
  return(result)
}

# ==== Tree Visualization Function with Pie Charts ====
create_tree_visualization <- function(tree, all_states, ace_result, output_dir, protein_name) {
  # Check for required packages
  if (!requireNamespace("ggtree", quietly = TRUE) || 
      !requireNamespace("RColorBrewer", quietly = TRUE) ||
      !requireNamespace("aplot", quietly = TRUE)) {
    log_message("Cannot create tree visualization - packages missing")
    return(FALSE)
  }
  
  # Load the packages
  library(ggtree)
  library(ggplot2)
  library(RColorBrewer)
  
  log_message("\nCreating tree visualization with state information...")
  
  # Get the number of tips and internal nodes
  n_tips <- length(tree$tip.label)
  n_internal <- tree$Nnode
  
  # Create a data frame for node annotations
  node_data <- data.frame(
    node = 1:(n_tips + n_internal),
    state = all_states,
    is_tip = 1:length(all_states) <= n_tips,
    label = c(tree$tip.label, rep("", n_internal)),
    stringsAsFactors = FALSE
  )
  
  # Get unique states and assign colors
  unique_states <- unique(all_states)
  n_states <- length(unique_states)
  
  # Choose an appropriate color palette
  if (n_states <= 9) {
    state_colors <- brewer.pal(max(3, n_states), "Set1")
  } else {
    state_colors <- colorRampPalette(brewer.pal(9, "Set1"))(n_states)
  }
  names(state_colors) <- unique_states
  
  # Create a basic tree plot with tips colored by state
  p <- ggtree(tree) %<+% node_data +
    geom_tippoint(aes(color = state), size = 2) +
    geom_tiplab(offset = 0.001, size = 2, align = TRUE) +
    geom_nodepoint(aes(color = state, size = !is_tip), alpha = 0.8) +
    scale_color_manual(values = state_colors, name = "Host") +
    scale_size_manual(values = c("TRUE" = 1, "FALSE" = 3), guide = "none") +
    theme(legend.position = "right") +
    labs(title = paste0(protein_name, " Tree with Host States"))
  
  # Save the rectangular tree
  rect_file <- file.path(output_dir, paste0(protein_name, "_tree_rectangular.pdf"))
  ggsave(rect_file, p, width = 12, height = max(8, n_tips/30), limitsize = FALSE)
  log_message("Saved rectangular tree to", rect_file)
  
  # Create a circular tree for large trees
  p2 <- ggtree(tree, layout = "circular") %<+% node_data +
    geom_tippoint(aes(color = state), size = 2) +
    geom_nodepoint(aes(color = state, size = !is_tip), alpha = 0.8) +
    scale_color_manual(values = state_colors, name = "Host") +
    scale_size_manual(values = c("TRUE" = 1, "FALSE" = 3), guide = "none") +
    theme(legend.position = "right") +
    labs(title = paste0(protein_name, " Circular Tree"))
  
  # Save the circular tree
  circ_file <- file.path(output_dir, paste0(protein_name, "_tree_circular.pdf"))
  ggsave(circ_file, p2, width = 10, height = 10)
  log_message("Saved circular tree to", circ_file)
  
  # Try to create pie charts
  tryCatch({
    # Create a new plot with nodepie
    p_pie <- ggtree(tree) %<+% node_data
    
    # Add tip points and labels
    p_pie <- p_pie + 
      geom_tippoint(aes(color = state), size = 2) +
      geom_tiplab(offset = 0.001, size = 2, align = TRUE) +
      scale_color_manual(values = state_colors, name = "Host")
    
    # Add pie charts for internal nodes
    start_node <- n_tips + 1
    for (i in 1:nrow(ace_result$lik.anc)) {
      node_id <- start_node + i - 1
      probs <- ace_result$lik.anc[i,]
      
      # Only include non-zero probabilities
      if (sum(probs) > 0) {
        p_pie <- p_pie + nodepie(node_id, probs, colors = state_colors, size = 3)
      }
    }
    
    # Save the pie chart tree
    pie_file <- file.path(output_dir, paste0(protein_name, "_tree_with_pies.pdf"))
    ggsave(pie_file, p_pie, width = 12, height = max(8, n_tips/30), limitsize = FALSE)
    log_message("Saved tree with pie charts to", pie_file)
  }, error = function(e) {
    log_message("Error creating pie chart tree:", e$message)
  })
  
  # Create a host distribution plot
  host_counts <- table(node_data$state[node_data$is_tip])
  host_df <- data.frame(
    Host = names(host_counts),
    Count = as.numeric(host_counts),
    stringsAsFactors = FALSE
  )
  
  p_dist <- ggplot(host_df, aes(x = reorder(Host, -Count), y = Count, fill = Host)) +
    geom_bar(stat = "identity") +
    scale_fill_manual(values = state_colors) +
    labs(title = paste0(protein_name, " Host Distribution"),
         x = "Host", y = "Number of Tips") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  # Save the distribution plot
  dist_file <- file.path(output_dir, paste0(protein_name, "_host_distribution.pdf"))
  ggsave(dist_file, p_dist, width = 10, height = 6)
  log_message("Saved host distribution plot to", dist_file)
  
  return(TRUE)
}

# ==== Function to Process a Single Protein ====
process_protein <- function(tree_file) {
  cat("=============================================\n")
  cat("Processing tree file:", tree_file, "\n")
  file_matches <- find_matching_files(tree_file)
  protein_name <- file_matches$protein_name
  branch_log_file <- file_matches$branch_log_file
  site_rate_file <- file_matches$site_rate_file
  mapping_file <- file_matches$mapping_file
  segment_name <- basename(dirname(tree_file))
  
  if (is.na(site_rate_file)) {
    cat("Warning: Site rate file not found for", protein_name, "in", segment_name, ". Skipping analysis.\n")
    return(NULL)
  }
  if (is.na(branch_log_file)) {
    cat("Warning: Branch rate file not found for", protein_name, "in", segment_name, ". Continuing with default branch rates.\n")
  }
  
  protein_name_clean <- gsub("[^a-zA-Z0-9]", "_", protein_name)
  output_dir <- file.path(main_output_dir, paste0(segment_name, "_", protein_name_clean))
  dir.create(output_dir, showWarnings = FALSE)
  log_message("Processing protein", protein_name, "in segment", segment_name)
  
  # ==== Load Tree ====
  log_message("Loading tree from", tree_file, "...")
  trees <- tryCatch({
    if (grepl("\\.nex$", tree_file, ignore.case = TRUE)) {
      read.nexus(tree_file)
    } else {
      read.tree(tree_file)
    }
  }, error = function(e) {
    log_message("Error loading tree:", e$message)
    return(NULL)
  })
  if (is.null(trees)) {
    log_message("Failed to load tree. Skipping protein", protein_name)
    return(NULL)
  }
  if (inherits(trees, "multiPhylo")) {
    log_message("Loaded", length(trees), "trees")
    burnin_idx <- floor(length(trees) * burnin_fraction) + 1
    if (burnin_idx > length(trees)) burnin_idx <- 1
    trees_to_use <- trees[burnin_idx:length(trees)]
    cons_tree <- consensus(trees_to_use)
  } else {
    log_message("Loaded a single tree")
    cons_tree <- trees
  }
  if (is.null(cons_tree$edge.length) || any(is.na(cons_tree$edge.length))) {
    cons_tree$edge.length <- rep(1, nrow(cons_tree$edge))
  }
  
  # Fix the tree to ensure it's compatible with ACE
  cons_tree <- fix_tree_for_ace(cons_tree)
  
  cons_tree_rooted <- tryCatch({
    midpoint.root(cons_tree)
  }, error = function(e) {
    log_message("Midpoint rooting failed:", e$message, "\nFalling back to outgroup rooting...")
    root(cons_tree, outgroup = 1, resolve.root = TRUE)
  })
  cons_tree_binary <- multi2di(cons_tree_rooted)
  
  # ==== Load Mapping File (if available) ====
  name_map <- list()
  if (!is.na(mapping_file)) {
    log_message("Loading mapping file from", mapping_file, "...")
    name_map <- tryCatch({
      fromJSON(mapping_file)
    }, error = function(e) {
      warning("Error loading mapping file:", e$message)
      list()
    })
  } else {
    log_message("No mapping file found; proceeding without mapping.")
  }
  
  # ==== Load Site Rates ====
  log_message("Loading site rates from", site_rate_file, "...")
  site_rates_log <- tryCatch({
    site_log <- read.table(site_rate_file, header = TRUE)
    log_message("Site rate log has", nrow(site_log), "rows and", ncol(site_log), "columns")
    
    site_rate_cols <- grep("site_rates", tolower(colnames(site_log)), value = TRUE)
    if (length(site_rate_cols) == 0) {
      # If we can't find site_rates in column names, try other common patterns
      site_rate_cols <- grep("rate|omega|dN/dS|omega", tolower(colnames(site_log)), value = TRUE)
      if (length(site_rate_cols) == 0) {
        # As last resort, use all columns except first (usually iteration or sample number)
        site_rate_cols <- colnames(site_log)[-1]
        log_message("Using columns 2 onward as site rates")
      } else {
        log_message("Found rate columns using alternative patterns:", paste(site_rate_cols, collapse=", "))
      }
    } else {
      log_message("Found site rate columns:", paste(site_rate_cols, collapse=", "))
    }
    
    # Print first few rows of the rate columns to verify data
    log_message("First few values of first rate column:")
    print(head(site_log[, site_rate_cols[1], drop = FALSE]))
    
    site_log_burned <- site_log[(floor(nrow(site_log) * burnin_fraction) + 1):nrow(site_log), ]
    site_means <- colMeans(site_log_burned[, site_rate_cols, drop = FALSE])
    
    # Check if site_means has reasonable values
    log_message("Range of site means:", min(site_means), "to", max(site_means))
    
    site_indices <- 1:length(site_means)
    result <- data.frame(site = site_indices, rate = site_means, stringsAsFactors = FALSE)
    log_message("Created site rate data frame with", nrow(result), "rows")
    log_message("Site rate summary: min =", min(result$rate), ", max =", max(result$rate), ", mean =", mean(result$rate))
    
    # Add additional safety check for uniform rates
    if(max(result$rate) - min(result$rate) < 1e-6) {
      log_message("WARNING: All sites have nearly identical rates. Check input data format.")
    }
    
    result
  }, error = function(e) {
    log_message("Error loading site rates:", e$message)
    NULL
  })
  if (is.null(site_rates_log)) {
    log_message("Skipping protein", protein_name, "due to site rate loading error")
    return(NULL)
  }
  
  # ==== Load Branch Rates ====
  rate_vector <- rep(1, nrow(cons_tree_binary$edge))
  if (!is.na(branch_log_file)) {
    log_message("Loading branch rates from", branch_log_file, "...")
    branch_rates_log <- tryCatch({
      branch_log <- read.table(branch_log_file, header = TRUE)
      branch_log_burned <- branch_log[(floor(nrow(branch_log) * burnin_fraction) + 1):nrow(branch_log), ]
      branch_rate_cols <- grep("branch_rates", tolower(colnames(branch_log_burned)), value = TRUE)
      if (length(branch_rate_cols) == 0) {
        stop("No branch rate columns found in", branch_log_file)
      }
      branch_means <- colMeans(branch_log_burned[, branch_rate_cols, drop = FALSE])
      branch_means
    }, error = function(e) {
      log_message("Error loading branch rates:", e$message)
      NULL
    })
    if (!is.null(branch_rates_log)) {
      log_message("Branch rates loaded successfully.")
      rate_vector <- branch_rates_log
      if (length(rate_vector) < nrow(cons_tree_binary$edge)) {
        rate_vector <- c(rate_vector, rep(mean(rate_vector), nrow(cons_tree_binary$edge) - length(rate_vector)))
      } else if (length(rate_vector) > nrow(cons_tree_binary$edge)) {
        rate_vector <- rate_vector[1:nrow(cons_tree_binary$edge)]
      }
    } else {
      log_message("Using default branch rates (all 1.0).")
    }
  } else {
    log_message("No branch rate file provided; using default branch rates.")
  }
  zero_branches <- rate_vector <= 1e-10
  log_message("Found", sum(zero_branches), "branches with near-zero rates.")
  
  # ==== Load Metadata ====
  log_message("Loading metadata...")
  
  meta_h1n1 <- tryCatch(
    read_csv("C:/Users/Ashley/Documents/flu_project/H1N1/H1N1_seq_info.csv") %>% mutate(Strain = "H1N1"),
    error = function(e) { 
      log_message("Error loading H1N1 metadata from original path:", e$message)
      # Try an alternative path
      tryCatch(
        read_csv("H1N1_seq_info.csv") %>% mutate(Strain = "H1N1"),
        error = function(e2) {
          log_message("Also failed with local path. H1N1 metadata not loaded.")
          data.frame()
        }
      )
    }
  )
  
  meta_h5n1 <- tryCatch(
    read_csv("C:/Users/Ashley/Documents/flu_project/H5N1/H5N1_seq_info.csv") %>% mutate(Strain = "H5N1"),
    error = function(e) { 
      log_message("Error loading H5N1 metadata from original path:", e$message)
      # Try an alternative path
      tryCatch(
        read_csv("H5N1_seq_info.csv") %>% mutate(Strain = "H5N1"),
        error = function(e2) {
          log_message("Also failed with local path. H5N1 metadata not loaded.")
          data.frame()
        }
      )
    }
  )
  
  meta_h7n9 <- tryCatch(
    read_csv("C:/Users/Ashley/Documents/flu_project/H7N9/H7N9_seq_info.csv") %>% mutate(Strain = "H7N9"),
    error = function(e) { 
      log_message("Error loading H7N9 metadata from original path:", e$message)
      # Try an alternative path
      tryCatch(
        read_csv("H7N9_seq_info.csv") %>% mutate(Strain = "H7N9"),
        error = function(e2) {
          log_message("Also failed with local path. H7N9 metadata not loaded.")
          data.frame()
        }
      )
    }
  )
  
  all_meta <- bind_rows(meta_h1n1, meta_h5n1, meta_h7n9)
  if (nrow(all_meta) == 0) {
    log_message("Warning: No metadata loaded.")
  } else {
    all_meta$Accession <- str_trim(all_meta$Accession)
    all_meta$Accession <- gsub("\\.1$", "", all_meta$Accession)
  }
  
  # ==== Convert Tip Labels Using Mapping ====
  original_tip_labels <- cons_tree_binary$tip.label
  
  # Apply the enhanced matching function
  strain_matching <- match_strains_to_meta(original_tip_labels, all_meta, name_map)
  
  # Process the matched results and handle duplicates
  new_tip_labels <- sapply(original_tip_labels, function(tip) {
    idx <- which(strain_matching$label == tip)
    if(length(idx) > 0 && !is.na(strain_matching$Final_Accession[idx])) {
      return(strain_matching$Final_Accession[idx])
    } else if(length(name_map) > 0 && tip %in% names(name_map)) {
      return(gsub("^(H1N1__|H5N1__|H7N9__)", "", name_map[[tip]]))
    } else {
      return(gsub("^(H1N1__|H5N1__|H7N9__)", "", tip))
    }
  }, USE.NAMES = FALSE)
  
  # Check for duplicates in new labels
  if (length(unique(new_tip_labels)) < length(new_tip_labels)) {
    log_message("\nWARNING: Found duplicate tip labels after renaming!")
    dupe_table <- table(new_tip_labels)
    dupes <- names(dupe_table[dupe_table > 1])
    log_message("  Duplicate labels:", paste(dupes, collapse=", "))
    
    # For each duplicate, find the original labels
    for (dupe in dupes) {
      dupe_indices <- which(new_tip_labels == dupe)
      log_message("  '", dupe, "' appears", length(dupe_indices), "times. Original labels:")
      
      # Make labels unique by appending a suffix
      for (i in seq_along(dupe_indices)) {
        if (i > 1) {  # Keep the first instance as is
          new_tip_labels[dupe_indices[i]] <- paste0(dupe, "_dupe", i)
          log_message("    Renamed to: '", new_tip_labels[dupe_indices[i]], "'")
        }
      }
    }
  }
  
  cons_tree_binary$tip.label <- new_tip_labels
  
  # Create a data frame for tip metadata
  tip_df <- data.frame(
    label = new_tip_labels, 
    Matched_Strain = strain_matching$Matched_Strain[match(new_tip_labels, strain_matching$Final_Accession)],
    stringsAsFactors = FALSE
  )
  
  # Merge with metadata for additional information
  annotated_tips <- left_join(tip_df, all_meta, by = c("label" = "Accession"))
  
  # Fill in Strain if it was matched but not joined
  idx_missing <- which(is.na(annotated_tips$Strain) & !is.na(annotated_tips$Matched_Strain))
  if(length(idx_missing) > 0) {
    annotated_tips$Strain[idx_missing] <- annotated_tips$Matched_Strain[idx_missing]
  }
  
  no_match_count <- sum(is.na(annotated_tips$Strain) & is.na(annotated_tips$Host))
  log_message("Directly matched", nrow(tip_df) - no_match_count, "of", nrow(tip_df), "tips")
  if (no_match_count > 0) {
    log_message("Still", no_match_count, "unmatched tips after enhanced matching")
  }
  
  host_states <- annotated_tips$Host
  host_states[is.na(host_states)] <- "Unknown"
  names(host_states) <- new_tip_labels
  host_states <- factor(host_states)
  
  # ==== Ancestral State Reconstruction using robust ACE approach ====
  log_message("Performing ancestral state reconstruction with robust ACE approach...")
  ace_result <- run_ace_with_fixes(host_states, cons_tree_binary)
  
  # Check if ACE result is properly structured
  log_message("Checking ACE result structure...")
  valid_ace_result <- TRUE
  
  if (is.null(ace_result)) {
    log_message("ERROR: ACE result is NULL")
    valid_ace_result <- FALSE
  } else if (!is.list(ace_result)) {
    log_message("ERROR: ACE result is not a list")
    valid_ace_result <- FALSE
  } else if (is.null(ace_result$lik.anc)) {
    log_message("ERROR: ACE result has no lik.anc component")
    valid_ace_result <- FALSE
  } else if (!is.matrix(ace_result$lik.anc)) {
    log_message("ERROR: ACE result lik.anc is not a matrix")
    log_message("Type of lik.anc: ", class(ace_result$lik.anc))
    valid_ace_result <- FALSE
  } else {
    log_message("ACE result structure valid:")
    log_message("  lik.anc dimensions: ", paste(dim(ace_result$lik.anc), collapse=" x "))
    if (!is.null(colnames(ace_result$lik.anc))) {
      log_message("  lik.anc columns: ", paste(colnames(ace_result$lik.anc), collapse=", "))
    } else {
      log_message("  WARNING: lik.anc has no column names")
    }
  }
  
  # Initialize empty ancestral states
  n_internal_nodes <- cons_tree_binary$Nnode
  log_message("Tree has", n_internal_nodes, "internal nodes")
  ancestral_states <- rep("Unknown", n_internal_nodes)
  
  # Only extract ancestral states if the ACE result is valid
  if (valid_ace_result) {
    tryCatch({
      # Make sure we have the right number of rows in lik.anc
      if (nrow(ace_result$lik.anc) != n_internal_nodes) {
        log_message("WARNING: Number of rows in lik.anc (", nrow(ace_result$lik.anc), 
                    ") does not match number of internal nodes (", n_internal_nodes, ")")
        # Adjust the size if needed - either truncate or extend
        if (nrow(ace_result$lik.anc) > n_internal_nodes) {
          log_message("Truncating lik.anc to match the number of internal nodes")
          ace_result$lik.anc <- ace_result$lik.anc[1:n_internal_nodes, , drop = FALSE]
        } else {
          log_message("Extending lik.anc with Unknown states to match the number of internal nodes")
          # Create extension rows with all zeros except for "Unknown" column
          unknown_col <- match("Unknown", colnames(ace_result$lik.anc))
          if (is.na(unknown_col)) {
            # If no "Unknown" column, just add rows of NA
            extension <- matrix(NA, nrow = n_internal_nodes - nrow(ace_result$lik.anc), 
                                ncol = ncol(ace_result$lik.anc),
                                dimnames = list(NULL, colnames(ace_result$lik.anc)))
          } else {
            # If there is an "Unknown" column, set it to 1 for extension rows
            extension <- matrix(0, nrow = n_internal_nodes - nrow(ace_result$lik.anc), 
                                ncol = ncol(ace_result$lik.anc),
                                dimnames = list(NULL, colnames(ace_result$lik.anc)))
            extension[, unknown_col] <- 1
          }
          ace_result$lik.anc <- rbind(ace_result$lik.anc, extension)
        }
      }
      
      # Now extract the ancestral states with error handling
      ancestral_states <- sapply(1:nrow(ace_result$lik.anc), function(i) {
        p <- ace_result$lik.anc[i, ]
        if (all(is.na(p)) || sum(p, na.rm = TRUE) == 0) {
          return("Unknown")
        } else {
          max_idx <- which.max(p)
          if (length(max_idx) == 0) {
            return("Unknown")
          } else {
            col_name <- colnames(ace_result$lik.anc)[max_idx]
            if (is.null(col_name) || is.na(col_name)) {
              return("Unknown")
            } else {
              return(col_name)
            }
          }
        }
      })
      
      log_message("Successfully extracted ancestral states for", length(ancestral_states), "internal nodes")
      log_message("Distribution of ancestral states:")
      state_counts <- table(ancestral_states)
      for (state in names(state_counts)) {
        log_message("  ", state, ": ", state_counts[state])
      }
      
    }, error = function(e) {
      log_message("ERROR extracting ancestral states:", e$message)
      log_message("Using 'Unknown' for all internal nodes")
      ancestral_states <- rep("Unknown", n_internal_nodes)
    })
  } else {
    log_message("Using 'Unknown' for all internal nodes due to invalid ACE result")
  }
  
  # Combine tip states and internal node states
  all_node_nums <- c(1:length(cons_tree_binary$tip.label),
                     (length(cons_tree_binary$tip.label)+1):(length(cons_tree_binary$tip.label)+cons_tree_binary$Nnode))
  all_states <- rep("", length(all_node_nums))
  names(all_states) <- all_node_nums
  all_states[1:length(cons_tree_binary$tip.label)] <- as.character(host_states)
  all_states[(length(cons_tree_binary$tip.label)+1):length(all_states)] <- ancestral_states
  
  # ==== Identify Host Shifts ====
  edge_df <- as.data.frame(cons_tree_binary$edge)
  colnames(edge_df) <- c("parent", "child")
  edge_df$parent_state <- all_states[as.character(edge_df$parent)]
  edge_df$child_state <- all_states[as.character(edge_df$child)]
  edge_df$host_switch <- (edge_df$parent_state != edge_df$child_state) &
    (edge_df$parent_state != "Unknown") &
    (edge_df$child_state != "Unknown")
  edge_df$rate <- rate_vector
  edge_df$zero_branch <- zero_branches
  edge_df$group <- ifelse(edge_df$host_switch, "Host Shift", "No Shift")
  edge_df_nonzero <- edge_df[!edge_df$zero_branch, ]
  log_message("After excluding near-zero branches:", nrow(edge_df_nonzero), "of", nrow(edge_df), "branches remain.")
  n_shifts <- sum(edge_df_nonzero$host_switch)
  log_message("Identified", n_shifts, "host shifts among non-zero branches.")
  host_transitions <- edge_df_nonzero %>% filter(host_switch) %>%
    count(parent_state, child_state, sort = TRUE) %>%
    mutate(transition = paste(parent_state, "→", child_state))
  log_message("Host shifts by transition type:")
  print(host_transitions)
  transition_file <- file.path(output_dir, paste0(protein_name_clean, "_host_transitions.csv"))
  write_csv(host_transitions, transition_file)
  log_message("Saved host transitions to", transition_file)
  zero_file <- file.path(output_dir, paste0(protein_name_clean, "_zero_branches.csv"))
  write_csv(edge_df[edge_df$zero_branch, ], zero_file)
  log_message("Saved list of zero branches to", zero_file)
  
  # ==== Create Tree Visualization ====
  log_message("Creating tree visualization...")
  tryCatch({
    viz_result <- create_tree_visualization(cons_tree_binary, all_states, ace_result, output_dir, protein_name)
    if (viz_result) {
      log_message("Tree visualization created successfully.")
    }
  }, error = function(e) {
    log_message("Error creating tree visualization:", e$message)
    log_message("Continuing with analysis without visualization.")
  })
  
  # ==== Branch Rate Analysis ====
  log_message("Performing branch rate analysis...")
  branch_rate_summary <- tryCatch({
    if (n_shifts > 0 && n_shifts < nrow(edge_df_nonzero)) {
      edge_df_nonzero %>% group_by(group) %>%
        summarize(mean_rate = mean(rate),
                  median_rate = median(rate),
                  min_rate = min(rate),
                  max_rate = max(rate),
                  n_branches = n())
    } else {
      NULL
    }
  }, error = function(e) {
    log_message("Error in branch rate summary:", e$message)
    NULL
  })
  if (!is.null(branch_rate_summary) && n_shifts > 0) {
    print(branch_rate_summary)
    wilcox_result <- tryCatch({
      wilcox.test(edge_df_nonzero$rate[edge_df_nonzero$host_switch],
                  edge_df_nonzero$rate[!edge_df_nonzero$host_switch])
    }, error = function(e) {
      log_message("Error in Wilcoxon test:", e$message)
      list(p.value = NA)
    })
    log_message("Wilcoxon test p-value (non-zero branches):", wilcox_result$p.value)
  } else {
    log_message("Insufficient data for branch rate comparison.")
    wilcox_result <- list(p.value = NA)
  }
  
  # ==== Site Rate Analysis ====
  log_message("Performing site rate analysis...")
  n_sites <- nrow(site_rates_log)
  
  # Calculate both standard and robust statistics
  mean_rate_val <- mean(site_rates_log$rate)
  sd_rate_val <- sd(site_rates_log$rate)
  median_rate_val <- median(site_rates_log$rate)
  mad_rate_val <- mad(site_rates_log$rate)  # Median Absolute Deviation (more robust)
  
  log_message("Site rate statistics:")
  log_message("  Mean:", mean_rate_val)
  log_message("  Standard Deviation:", sd_rate_val)
  log_message("  Coefficient of Variation:", sd_rate_val/mean_rate_val * 100, "%")
  log_message("  Median:", median_rate_val)
  log_message("  Median Absolute Deviation:", mad_rate_val)
  
  # Add quartile-based categorization
  site_rates_ranked <- site_rates_log %>%
    arrange(rate) %>%
    mutate(rank = row_number(),
           percentile = rank / n() * 100,
           quartile = ntile(rate, 4),
           rate_category = factor(quartile, levels = 1:4, labels = c("Very Slow", "Slow", "Fast", "Very Fast")))
  
  site_rate_summary <- site_rates_ranked %>% group_by(rate_category) %>%
    summarize(n_sites = n(),
              min_rate = min(rate),
              max_rate = max(rate),
              mean_rate = mean(rate),
              median_rate = median(rate))
  
  log_message("Site rate statistics by category:")
  print(site_rate_summary)
  
  # Test for normality
  shapiro_test <- tryCatch({
    shapiro.test(site_rates_log$rate)
  }, error = function(e) {
    log_message("Normality test failed:", e$message)
    list(p.value = 0)  # If test fails, assume non-normality
  })
  
  log_message("Shapiro-Wilk normality test: W =", shapiro_test$statistic, 
              ", p-value =", shapiro_test$p.value)
  is_normal <- shapiro_test$p.value >= 0.05
  
  # Calculate z-scores
  site_rates_log$z_score <- (site_rates_log$rate - mean_rate_val) / sd_rate_val
  site_rates_log$robust_z <- (site_rates_log$rate - median_rate_val) / mad_rate_val
  
  # Identify sites under selection
  log_message("Identifying sites under selection...")
  
  # Determine which z-score to use based on normality
  z_score_col <- if(is_normal) "z_score" else "robust_z"
  reference_value <- if(is_normal) mean_rate_val else median_rate_val
  spread_value <- if(is_normal) sd_rate_val else mad_rate_val
  
  log_message("Using", if(is_normal) "standard" else "robust", "statistics for site selection")
  
  # Check for extreme cases
  if (spread_value <= 1e-6 || spread_value/reference_value < 0.001) {
    log_message("WARNING: Very low variation. All sites appear to evolve at similar rates.")
    fast_evolving_sites <- data.frame()
    very_fast_sites <- data.frame()
    slow_evolving_sites <- data.frame()
  } else if (any(is.nan(site_rates_log[[z_score_col]])) || any(is.infinite(site_rates_log[[z_score_col]]))) {
    log_message("WARNING: Invalid z-scores detected. Using percentile-based selection instead.")
    # Use percentile-based selection as fallback
    fast_evolving_sites <- site_rates_log %>% 
      arrange(desc(rate)) %>% 
      slice_head(prop = 0.025) %>%
      arrange(desc(rate))
    
    very_fast_sites <- site_rates_log %>% 
      arrange(desc(rate)) %>% 
      slice_head(prop = 0.01) %>%
      arrange(desc(rate))
    
    slow_evolving_sites <- site_rates_log %>% 
      arrange(rate) %>% 
      slice_head(prop = 0.025) %>%
      arrange(rate)
  } else {
    # Use z-score based selection with appropriate threshold
    fast_evolving_sites <- site_rates_log %>% filter(.data[[z_score_col]] > 1.96) %>% arrange(desc(rate))
    very_fast_sites <- site_rates_log %>% filter(.data[[z_score_col]] > 2.58) %>% arrange(desc(rate))
    slow_evolving_sites <- site_rates_log %>% filter(.data[[z_score_col]] < -1.96) %>% arrange(desc(rate))
    
    # Fallback to percentile-based if none detected despite high variation
    if ((nrow(fast_evolving_sites) == 0 && nrow(slow_evolving_sites) == 0) && 
        (spread_value/reference_value > 0.5)) {
      log_message("No sites detected with z-score method despite high variation. Using percentile approach.")
      fast_evolving_sites <- site_rates_log %>% 
        arrange(desc(rate)) %>% 
        slice_head(prop = 0.025) %>%
        arrange(desc(rate))
      
      very_fast_sites <- site_rates_log %>% 
        arrange(desc(rate)) %>% 
        slice_head(prop = 0.01) %>%
        arrange(desc(rate))
      
      slow_evolving_sites <- site_rates_log %>% 
        arrange(rate) %>% 
        slice_head(prop = 0.025) %>%
        arrange(rate)
    }
  }
  
  log_message("Fast-evolving sites:", nrow(fast_evolving_sites))
  log_message("Very fast sites:", nrow(very_fast_sites))
  log_message("Slow-evolving sites:", nrow(slow_evolving_sites))
  
  # Save site lists if they exist
  fast_sites_file <- if(nrow(fast_evolving_sites) > 0) 
    file.path(output_dir, paste0(protein_name_clean, "_fast_sites.csv")) else NA
  if (!is.na(fast_sites_file)) {
    log_message("Top 5 fast-evolving sites:")
    top_fast <- head(fast_evolving_sites, 5)
    for (i in 1:nrow(top_fast)) {
      comparison_val <- top_fast$rate[i] / reference_value
      z_val <- if("robust_z" %in% names(top_fast)) top_fast$robust_z[i] else top_fast$z_score[i]
      log_message(sprintf("  Site %d: rate = %.4f (%.2fx faster, z = %.2f)",
                          top_fast$site[i], top_fast$rate[i], comparison_val, z_val))
    }
    write_csv(fast_evolving_sites, fast_sites_file)
    log_message("Saved fast-evolving sites to", fast_sites_file)
  }
  
  slow_sites_file <- if(nrow(slow_evolving_sites) > 0) 
    file.path(output_dir, paste0(protein_name_clean, "_slow_sites.csv")) else NA
  if (!is.na(slow_sites_file)) {
    log_message("Top 5 slow-evolving sites:")
    top_slow <- head(slow_evolving_sites, 5)
    for (i in 1:nrow(top_slow)) {
      comparison_val <- reference_value / top_slow$rate[i]
      z_val <- if("robust_z" %in% names(top_slow)) top_slow$robust_z[i] else top_slow$z_score[i]
      log_message(sprintf("  Site %d: rate = %.4f (%.2fx slower, z = %.2f)",
                          top_slow$site[i], top_slow$rate[i], comparison_val, z_val))
    }
    write_csv(slow_evolving_sites, slow_sites_file)
    log_message("Saved slow-evolving sites to", slow_sites_file)
  }
  
  # ==== Create Site Rate Visualizations ====
  log_message("Creating site rate plot...")
  
  # Basic site rate plot
  site_rate_plot <- ggplot(site_rates_log, aes(x = site, y = rate)) +
    geom_point(alpha = 0.7) +
    geom_smooth(method = "loess", se = TRUE, color = "red") +
    labs(title = paste("Site-specific Rates for", protein_name, "Protein"),
         x = "Site Position", y = "Evolutionary Rate") +
    theme_minimal()
  site_plot_file <- file.path(output_dir, paste0(protein_name_clean, "_site_rates.pdf"))
  ggsave(site_plot_file, site_rate_plot, width = 10, height = 6)
  log_message("Saved site rate plot to", site_plot_file)
  
  # Quartile boxplot
  quartile_plot <- ggplot(site_rates_ranked, aes(x = rate_category, y = rate, fill = rate_category)) +
    geom_boxplot() +
    labs(title = paste("Site Rate Distribution for", protein_name, "Protein"),
         x = "Rate Category", y = "Evolutionary Rate") +
    theme_minimal() +
    theme(legend.position = "none")
  quartile_plot_file <- file.path(output_dir, paste0(protein_name_clean, "_rate_quartiles.pdf"))
  ggsave(quartile_plot_file, quartile_plot, width = 8, height = 6)
  log_message("Saved quartile plot to", quartile_plot_file)
  
  # Significant sites plot
  # Create dataframe for plotting with significant sites highlighted
  significant_sites <- c()
  if(nrow(fast_evolving_sites) > 0) significant_sites <- c(significant_sites, fast_evolving_sites$site)
  if(nrow(slow_evolving_sites) > 0) significant_sites <- c(significant_sites, slow_evolving_sites$site)
  
  site_rates_log$significant <- site_rates_log$site %in% significant_sites
  
  site_significance_plot <- ggplot(site_rates_log, aes(x = site, y = rate, color = significant)) +
    geom_point(alpha = 0.7, size = 1.5) +
    scale_color_manual(values = c("FALSE" = "gray50", "TRUE" = "red"),
                       labels = c("FALSE" = "Average", "TRUE" = "Significant"),
                       name = "Evolution Rate") +
    geom_hline(yintercept = reference_value, linetype = "dashed", color = "black") +
    geom_hline(yintercept = reference_value + 1.96 * spread_value, linetype = "dotted", color = "red") +
    geom_hline(yintercept = reference_value - 1.96 * spread_value, linetype = "dotted", color = "blue") +
    labs(title = paste("Significant Sites in", protein_name, "Protein"),
         subtitle = paste(nrow(fast_evolving_sites), "fast and", nrow(slow_evolving_sites), "slow sites"),
         x = "Site Position", y = "Evolutionary Rate") +
    theme_minimal()
  
  sig_plot_file <- file.path(output_dir, paste0(protein_name_clean, "_significant_sites.pdf"))
  ggsave(sig_plot_file, site_significance_plot, width = 10, height = 6)
  log_message("Saved site significance plot to", sig_plot_file)
  
  # ==== Branch Rate Visualizations ====
  branch_plot_file <- NA
  violin_plot_file <- NA
  log_message("Creating branch rate visualizations...")
  if (n_shifts > 0 && n_shifts < nrow(edge_df_nonzero)) {
    branch_rate_plot <- ggplot(edge_df_nonzero, aes(x = rate, fill = group)) +
      geom_density(alpha = 0.5) +
      labs(title = paste("Branch Rate Distribution for", protein_name, "Protein"),
           subtitle = paste(n_shifts, "host shifts among", nrow(edge_df_nonzero), "non-zero branches"),
           x = "Branch Rate", y = "Density", fill = "Branch Type") +
      theme_minimal()
    branch_plot_file <- file.path(output_dir, paste0(protein_name_clean, "_branch_rates.pdf"))
    ggsave(branch_plot_file, branch_rate_plot, width = 10, height = 6)
    log_message("Saved branch rate density plot to", branch_plot_file)
    
    branch_rate_plot_v <- ggplot(edge_df_nonzero, aes(x = rate, y = group, fill = group)) +
      geom_violin(alpha = 0.5) +
      labs(title = paste("Branch Rate Violin Plot for", protein_name, "Protein"),
           subtitle = paste(n_shifts, "host shifts among", nrow(edge_df_nonzero), "non-zero branches"),
           x = "Branch Rate", y = "Branch Type", fill = "Branch Type") +
      theme_minimal()
    violin_plot_file <- file.path(output_dir, paste0(protein_name_clean, "_branch_rates_violin.pdf"))
    ggsave(violin_plot_file, branch_rate_plot_v, width = 10, height = 6)
    log_message("Saved branch rate violin plot to", violin_plot_file)
  } else {
    log_message("Skipping branch rate visualizations (insufficient host shift data).")
  }
  
  # ==== Create Per-Protein Summary Report ====
  log_message("Generating per-protein summary report...")
  summary_file <- file.path(output_dir, paste0(protein_name_clean, "_summary.txt"))
  sink(summary_file)
  cat("=============================================\n")
  cat(segment_name, "-", protein_name, "PROTEIN EVOLUTIONARY ANALYSIS\n")
  cat("=============================================\n\n")
  cat("ANALYSIS OVERVIEW:\n")
  cat("- Sites analyzed:", n_sites, "\n")
  cat("- Non-zero branches:", nrow(edge_df_nonzero), "\n")
  cat("- Zero branches excluded:", sum(zero_branches), "\n")
  if (n_shifts > 0) {
    cat("- Host shifts identified:", n_shifts, " (", round(n_shifts / nrow(edge_df_nonzero) * 100, 1), "%)\n")
  } else {
    cat("- No host shifts identified\n")
  }
  cat("- Fast-evolving sites:", nrow(fast_evolving_sites), "\n")
  cat("- Slow-evolving sites:", nrow(slow_evolving_sites), "\n\n")
  if (n_shifts > 0) {
    cat("HOST SHIFT SUMMARY:\n")
    top_transitions <- head(host_transitions, 5)
    for (i in 1:nrow(top_transitions)) {
      cat("  * ", top_transitions$transition[i], ": ", top_transitions$n[i], " occurrences\n", sep = "")
    }
    cat("\n")
    if (!is.null(branch_rate_summary)) {
      cat("BRANCH RATE ANALYSIS:\n")
      cat("- Host shift branches (n =", branch_rate_summary$n_branches[branch_rate_summary$group == "Host Shift"], "): median rate =", 
          branch_rate_summary$median_rate[branch_rate_summary$group == "Host Shift"], "\n")
      cat("- Non-host shift branches (n =", branch_rate_summary$n_branches[branch_rate_summary$group == "No Shift"], "): median rate =", 
          branch_rate_summary$median_rate[branch_rate_summary$group == "No Shift"], "\n")
      cat("- Wilcoxon test p-value:", wilcox_result$p.value, "\n")
    }
  } else {
    cat("HOST SHIFT SUMMARY: No host shifts identified\n")
  }
  cat("\nSITE RATE ANALYSIS:\n")
  cat("- Mean site rate:", mean_rate_val, "\n")
  cat("- SD:", sd_rate_val, "\n")
  cat("- Median site rate:", median_rate_val, "\n")
  cat("- MAD:", mad_rate_val, "\n")
  cat("- Coefficient of variation:", round(sd_rate_val / mean_rate_val * 100, 1), "%\n")
  cat("- Range:", min(site_rates_log$rate), "to", max(site_rates_log$rate), "\n\n")
  cat("FILES GENERATED:\n")
  cat("- Site rate plot:", site_plot_file, "\n")
  cat("- Quartile plot:", quartile_plot_file, "\n")
  cat("- Site significance plot:", sig_plot_file, "\n")
  if (!is.na(branch_plot_file)) cat("- Branch rate density plot:", branch_plot_file, "\n")
  if (!is.na(violin_plot_file)) cat("- Branch rate violin plot:", violin_plot_file, "\n")
  if (n_shifts > 0) cat("- Host transitions table:", transition_file, "\n")
  cat("- Zero branch list:", zero_file, "\n")
  if (!is.na(fast_sites_file)) cat("- Fast-evolving sites table:", fast_sites_file, "\n")
  if (!is.na(slow_sites_file)) cat("- Slow-evolving sites table:", slow_sites_file, "\n")
  
  # Add interpretations if significant sites found
  if (nrow(fast_evolving_sites) > 0 || nrow(slow_evolving_sites) > 0) {
    cat("\nINTERPRETATION:\n")
    cat("This analysis has identified sites that evolve at significantly different rates\n")
    cat("compared to the protein average. Fast-evolving sites may be under positive\n")
    cat("selection or relaxed constraints, while slow-evolving sites typically represent\n")
    cat("functionally or structurally critical regions under purifying selection.\n")
  }
  
  sink()
  log_message("Per-protein summary report saved to", summary_file)
  
  # ==== Collect Results for Overall Summary ====
  # Now include significant_branch_difference computed from the wilcox p-value.
  results_row <- data.frame(
    segment = segment_name,
    protein = protein_name,
    n_sites = n_sites,
    n_branches = nrow(edge_df),
    n_nonzero_branches = nrow(edge_df_nonzero),
    n_zero_branches = sum(zero_branches),
    n_host_shifts = n_shifts,
    host_shift_percent = ifelse(nrow(edge_df_nonzero) > 0, round(n_shifts / nrow(edge_df_nonzero) * 100, 1), NA),
    wilcox_pvalue = wilcox_result$p.value,
    significant_branch_difference = ifelse(!is.null(wilcox_result) && !is.na(wilcox_result$p.value), wilcox_result$p.value < 0.05, NA),
    mean_site_rate = mean_rate_val,
    sd_site_rate = sd_rate_val,
    cv_site_rate = sd_rate_val / mean_rate_val * 100,
    min_site_rate = min(site_rates_log$rate),
    max_site_rate = max(site_rates_log$rate),
    n_fast_sites = nrow(fast_evolving_sites),
    n_very_fast_sites = nrow(very_fast_sites),
    n_slow_sites = nrow(slow_evolving_sites),
    top_fast_site = ifelse(nrow(fast_evolving_sites) > 0, fast_evolving_sites$site[1], NA),
    top_fast_rate = ifelse(nrow(fast_evolving_sites) > 0, fast_evolving_sites$rate[1], NA),
    top_slow_site = ifelse(nrow(slow_evolving_sites) > 0, slow_evolving_sites$site[1], NA),
    top_slow_rate = ifelse(nrow(slow_evolving_sites) > 0, slow_evolving_sites$rate[1], NA),
    stringsAsFactors = FALSE
  )
  if (n_shifts > 0 && nrow(host_transitions) > 0) {
    results_row$top_transition <- host_transitions$transition[1]
    results_row$top_transition_count <- host_transitions$n[1]
  } else {
    results_row$top_transition <- NA
    results_row$top_transition_count <- NA
  }
  
  log_message("Completed analysis for", protein_name, "in", segment_name)
  log_message("=============================================\n")
  
  return(results_row)
}

# ==== Main Processing Loop ====
log_message("==== Starting Batch Processing of All Segments ====\n")
all_results <- data.frame()
for (seg_dir in seg_dirs) {
  segment_name <- basename(seg_dir)
  log_message("Processing segment directory:", segment_name)
  # List only tree files ending with _tree.nex (case-insensitive)
  tree_files <- list.files(seg_dir, pattern = "_tree\\.nex$", full.names = TRUE, ignore.case = TRUE)
  if (length(tree_files) == 0) {
    log_message("No tree files found in", segment_name, ". Skipping this segment.\n")
    next
  }
  log_message("Found", length(tree_files), "tree file(s) in", segment_name)
  for (tree_file in tree_files) {
    results <- process_protein(tree_file)
    if (!is.null(results)) {
      all_results <- rbind(all_results, results)
    }
  }
}

# ==== Generate Overall Summary Report ====
log_message("\n==== Generating Overall Summary Report ====")
if (nrow(all_results) == 0) {
  log_message("No results were collected. Check if the file paths and file patterns are correct.")
} else {
  all_results <- all_results %>% arrange(segment, protein)
  results_csv <- file.path(main_output_dir, "all_proteins_results.csv")
  write_csv(all_results, results_csv)
  log_message("Saved complete results table to", results_csv)
  
  sink(summary_all_file, append = TRUE)
  cat("=======================================================\n")
  cat("SUMMARY OF EVOLUTIONARY ANALYSIS ACROSS ALL PROTEINS\n")
  cat("=======================================================\n\n")
  cat("Analysis includes", nrow(all_results), "proteins across", length(unique(all_results$segment)), "segments\n\n")
  
  cat("SITE RATE VARIATION:\n")
  cat("------------------\n")
  cv_sorted <- all_results %>% arrange(desc(cv_site_rate)) %>% 
    dplyr::select(segment, protein, n_sites, mean_site_rate, cv_site_rate, n_fast_sites, n_slow_sites)
  cat("Proteins with highest rate variation (CV):\n")
  for(i in 1:min(5, nrow(cv_sorted))) {
    cat(sprintf("%d. %s-%s: CV = %.1f%%, %d fast sites, %d slow sites\n",
                i, cv_sorted$segment[i], cv_sorted$protein[i],
                cv_sorted$cv_site_rate[i], cv_sorted$n_fast_sites[i], cv_sorted$n_slow_sites[i]))
  }
  cat("\n")
  
  selection_sorted <- all_results %>%
    mutate(total_selection_sites = n_fast_sites + n_slow_sites) %>%
    arrange(desc(total_selection_sites)) %>% 
    dplyr::select(segment, protein, n_sites, total_selection_sites, n_fast_sites, n_slow_sites)
  cat("Proteins with most sites under selection:\n")
  for(i in 1:min(5, nrow(selection_sorted))) {
    cat(sprintf("%d. %s-%s: %d sites (%.1f%%), %d fast, %d slow\n",
                i, selection_sorted$segment[i], selection_sorted$protein[i],
                selection_sorted$total_selection_sites[i],
                selection_sorted$total_selection_sites[i] / selection_sorted$n_sites[i] * 100,
                selection_sorted$n_fast_sites[i], selection_sorted$n_slow_sites[i]))
  }
  cat("\n")
  
  cat("HOST SHIFT ANALYSIS:\n")
  cat("------------------\n")
  host_shift_sorted <- all_results %>% filter(!is.na(n_host_shifts) & n_host_shifts > 0) %>%
    arrange(desc(host_shift_percent)) %>% 
    dplyr::select(segment, protein, n_nonzero_branches, n_host_shifts, host_shift_percent, significant_branch_difference, top_transition, top_transition_count)
  if (nrow(host_shift_sorted) > 0) {
    cat("Proteins with highest proportion of host shifts:\n")
    for(i in 1:min(5, nrow(host_shift_sorted))) {
      cat(sprintf("%d. %s-%s: %d shifts (%.1f%%), %s\n",
                  i, host_shift_sorted$segment[i], host_shift_sorted$protein[i],
                  host_shift_sorted$n_host_shifts[i], host_shift_sorted$host_shift_percent[i],
                  ifelse(host_shift_sorted$significant_branch_difference[i], "significant", "not significant")))
      if (!is.na(host_shift_sorted$top_transition[i])) {
        cat(sprintf("   Most common transition: %s (%d occurrences)\n",
                    host_shift_sorted$top_transition[i], host_shift_sorted$top_transition_count[i]))
      }
    }
  } else {
    cat("No proteins with host shifts were found.\n")
  }
  cat("\n")
  
  cat("SUMMARY BY SEGMENT:\n")
  cat("------------------\n")
  for (segment in unique(all_results$segment)) {
    segment_results <- all_results %>% filter(segment == segment)
    cat(segment, "- analyzed", nrow(segment_results), "proteins\n")
    fast_sites_summary <- segment_results %>% arrange(desc(n_fast_sites)) %>% 
      dplyr::select(protein, n_sites, n_fast_sites) %>% head(3)
    cat("  Proteins with most fast sites:\n")
    for (i in 1:nrow(fast_sites_summary)) {
      cat(sprintf("    * %s: %d sites (%.1f%%)\n",
                  fast_sites_summary$protein[i],
                  fast_sites_summary$n_fast_sites[i],
                  fast_sites_summary$n_fast_sites[i] / fast_sites_summary$n_sites[i] * 100))
    }
    host_shifts_summary <- segment_results %>% filter(!is.na(n_host_shifts) & n_host_shifts > 0) %>% 
      arrange(desc(n_host_shifts)) %>% dplyr::select(protein, n_host_shifts, host_shift_percent) %>% head(3)
    if (nrow(host_shifts_summary) > 0) {
      cat("  Proteins with most host shifts:\n")
      for (i in 1:nrow(host_shifts_summary)) {
        cat(sprintf("    * %s: %d shifts (%.1f%%)\n",
                    host_shifts_summary$protein[i],
                    host_shifts_summary$n_host_shifts[i],
                    host_shifts_summary$host_shift_percent[i]))
      }
    } else {
      cat("  No host shifts detected in this segment.\n")
    }
    cat("\n")
  }
  
  # ==== Report Interesting Sites Per Protein ====
  cat("INTERESTING SITES PER PROTEIN:\n")
  cat("---------------------------\n")
  for (i in 1:nrow(all_results)) {
    fast_site_info <- ifelse(!is.na(all_results$top_fast_site[i]), 
                             sprintf("site %s (rate = %.4f)", all_results$top_fast_site[i], all_results$top_fast_rate[i]),
                             "none identified")
    slow_site_info <- ifelse(!is.na(all_results$top_slow_site[i]), 
                             sprintf("site %s (rate = %.4f)", all_results$top_slow_site[i], all_results$top_slow_rate[i]),
                             "none identified")
    
    cat(sprintf("%s-%s: Top fast site: %s; Top slow site: %s\n",
                all_results$segment[i],
                all_results$protein[i],
                fast_site_info,
                slow_site_info))
  }
  
  cat("\nDETAILED RESULTS TABLE:\n")
  cat("---------------------\n")
  cat("For complete results, see:", results_csv, "\n\n")
  
  cat("KEY FINDINGS:\n")
  cat("------------\n")
  cat("Protein with highest rate variation:", cv_sorted$segment[1], "-", cv_sorted$protein[1], "\n")
  
  # Add interpretation if sites under selection were found
  if(sum(all_results$n_fast_sites) > 0 || sum(all_results$n_slow_sites) > 0) {
    cat("\nINTERPRETATION:\n")
    cat("This analysis has identified proteins with significant heterogeneity in evolutionary rates\n")
    cat("across sites. Fast-evolving sites may be under positive selection or relaxed constraints,\n")
    cat("while slow-evolving sites typically represent functionally or structurally critical\n")
    cat("regions under purifying selection. Host shifts appear to be associated with certain\n")
    cat("protein evolutionary patterns, which may provide insights into adaptation mechanisms.\n")
  }
  
  sink()
  
  log_message("\nOverall summary appended to", summary_all_file)
}

# ==== Final Success Message ====
log_message("\n==== ANALYSIS COMPLETE ====")
log_message("All output files saved to:", main_output_dir)
log_message("Debug log saved to:", debug_log)
log_message("Summary file:", summary_all_file)
log_message("Results CSV:", results_csv)
log_message("\nTo load and view the results in R:")
log_message('results <- read_csv("', results_csv, '")')
log_message('View(results)')
log_message("\nThank you for using the combined protein evolution analysis script!")