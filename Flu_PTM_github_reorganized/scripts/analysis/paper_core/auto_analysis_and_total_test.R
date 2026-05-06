# ==== Combined Protein Evolution Analysis Across All Segments with KS Test ====
# This script extends the original analysis to include Kolmogorov-Smirnov (KS) tests
# alongside Wilcoxon tests and adds an overall analysis to assess differences between groups
# independent of protein.

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

# Initialize data frame for overall branch rates
all_branch_rates <- data.frame()

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
  result <- data.frame(
    label = tip_labels,
    Final_Accession = NA_character_,
    Matched_Strain = NA_character_,
    match_method = NA_character_,
    stringsAsFactors = FALSE
  )
  
  log_message("\n===== MATCHING DEBUG =====")
  log_message("Total tip labels to match: ", length(tip_labels))
  log_message("Total metadata entries: ", nrow(meta_data))
  
  if (length(mapping) > 0) {
    log_message("\nStep 1: Direct mapping via JSON file")
    for (i in 1:nrow(result)) {
      curr_label <- result$label[i]
      
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
    
    json_matches <- sum(!is.na(result$match_method) & 
                          grepl("^json_", result$match_method), na.rm=TRUE)
    log_message("Total matches via JSON mapping: ", json_matches, " (", round(json_matches/nrow(result)*100, 1), "%)")
  }
  
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
  
  rem <- which(is.na(result$Final_Accession))
  log_message("\nStep 5: No underscore suffix matching for ", length(rem), " remaining unmatched labels")
  
  if (length(rem) > 0) {
    no_underscore_matches <- 0
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
  
  log_message("\nMatching results summary:")
  log_message("- Total tips:", nrow(result))
  log_message("- Successfully matched:", sum(!is.na(result$Final_Accession)))
  log_message("- Unmatched:", sum(is.na(result$Final_Accession)))
  
  for (method in unique(result$match_method[!is.na(result$match_method)])) {
    count <- sum(result$match_method == method, na.rm=TRUE)
    log_message("- Matched by ", method, ": ", count, " (", round(count/nrow(result)*100, 1), "%)")
  }
  
  return(result)
}

# ==== Function to fix the tree for ACE compatibility ====
fix_tree_for_ace <- function(tree) {
  zero_branches <- sum(tree$edge.length <= 1e-6, na.rm=TRUE)
  na_branches <- sum(is.na(tree$edge.length))
  
  log_message("Tree has", zero_branches, "zero-length branches and", na_branches, "NA branches")
  
  if (zero_branches > 0 || na_branches > 0) {
    log_message("Setting small positive values for zero/NA branch lengths")
    tree$edge.length[tree$edge.length <= 1e-6 | is.na(tree$edge.length)] <- 0.01
  }
  
  neg_branches <- sum(tree$edge.length < 0, na.rm=TRUE)
  if (neg_branches > 0) {
    log_message("Fixing", neg_branches, "negative branch lengths")
    tree$edge.length[tree$edge.length < 0] <- 0.01
  }
  
  if (!is.binary.tree(tree)) {
    log_message("Making tree binary")
    tree <- multi2di(tree)
  }
  
  if (!is.rooted(tree)) {
    log_message("Rooting tree")
    tree <- root(tree, outgroup=1, resolve.root=TRUE)
  }
  
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
  tree <- fix_tree_for_ace(tree)
  
  missing_states <- setdiff(tree$tip.label, names(host_states))
  if (length(missing_states) > 0) {
    log_message("Found", length(missing_states), "tips without host states")
    for (tip in missing_states) {
      host_states[tip] <- "Unknown"
    }
  }
  
  if (!is.factor(host_states)) {
    log_message("Converting host states to factor")
    host_states <- factor(host_states)
  }
  
  if (!"Unknown" %in% levels(host_states)) {
    log_message("Adding 'Unknown' to host state levels")
    host_states <- factor(host_states, levels = c(levels(host_states), "Unknown"))
  }
  
  tryCatch({
    log_message("Attempting ACE with default parameters")
    ace_result <- ace(host_states, tree, type="discrete")
    log_message("ACE successful!")
    return(validate_ace_result(ace_result, tree, host_states))
  }, error = function(e) {
    log_message("First ACE attempt failed:", e$message)
    
    tryCatch({
      log_message("Trying ACE with marginal=TRUE")
      ace_result <- ace(host_states, tree, type="discrete", marginal=TRUE)
      log_message("ACE with marginal=TRUE successful!")
      return(validate_ace_result(ace_result, tree, host_states))
    }, error = function(e2) {
      log_message("Second ACE attempt failed:", e2$message)
      
      tryCatch({
        log_message("Trying ACE with CI=FALSE")
        ace_result <- ace(host_states, tree, type="discrete", CI=FALSE)
        log_message("ACE with CI=FALSE successful!")
        return(validate_ace_result(ace_result, tree, host_states))
      }, error = function(e3) {
        log_message("Third ACE attempt failed:", e3$message)
        
        log_message("Making more aggressive tree fixes")
        tree$edge.length <- pmax(tree$edge.length, 0.05)
        
        tryCatch({
          log_message("Final ACE attempt with heavily optimized tree")
          ace_result <- ace(host_states, tree, type="discrete", CI=FALSE, marginal=TRUE)
          log_message("Final ACE attempt successful!")
          return(validate_ace_result(ace_result, tree, host_states))
        }, error = function(e4) {
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
  
  if (is.null(ace_result) || !is.list(ace_result) || is.null(ace_result$lik.anc)) {
    log_message("Invalid ACE result: NULL or missing lik.anc")
    return(create_dummy_ace_result(tree, host_states))
  }
  
  if (!is.matrix(ace_result$lik.anc)) {
    log_message("Invalid ACE result: lik.anc is not a matrix")
    return(create_dummy_ace_result(tree, host_states))
  }
  
  if (is.null(colnames(ace_result$lik.anc))) {
    log_message("ACE result missing column names, adding them")
    colnames(ace_result$lik.anc) <- levels(host_states)
  }
  
  if (nrow(ace_result$lik.anc) != n_internal) {
    log_message("ACE result has incorrect number of rows:", nrow(ace_result$lik.anc), 
                "expected:", n_internal)
    
    if (nrow(ace_result$lik.anc) > n_internal) {
      log_message("Truncating ACE result rows")
      ace_result$lik.anc <- ace_result$lik.anc[1:n_internal, , drop = FALSE]
    } else {
      n_missing <- n_internal - nrow(ace_result$lik.anc)
      n_states <- ncol(ace_result$lik.anc)
      
      extension <- matrix(1/n_states, nrow = n_missing, ncol = n_states,
                          dimnames = list(NULL, colnames(ace_result$lik.anc)))
      
      unknown_col <- match("Unknown", colnames(ace_result$lik.anc))
      if (!is.na(unknown_col)) {
        extension[, ] <- 0.1 / (n_states - 1)
        extension[, unknown_col] <- 0.9
      }
      
      ace_result$lik.anc <- rbind(ace_result$lik.anc, extension)
    }
  }
  
  return(ace_result)
}

# Helper function to create a dummy ACE result
create_dummy_ace_result <- function(tree, host_states) {
  log_message("Creating dummy ACE result")
  n_internal <- tree$Nnode
  state_levels <- levels(host_states)
  n_states <- length(state_levels)
  
  lik_anc <- matrix(0, nrow = n_internal, ncol = n_states,
                    dimnames = list(NULL, state_levels))
  
  tip_states <- as.character(host_states)
  tip_states <- tip_states[tip_states != "Unknown"]
  
  if (length(tip_states) > 0) {
    most_common <- names(sort(table(tip_states), decreasing = TRUE))[1]
    most_common_idx <- match(most_common, state_levels)
    if (!is.na(most_common_idx)) {
      lik_anc[, most_common_idx] <- 0.6
    }
  }
  
  unknown_idx <- match("Unknown", state_levels)
  if (!is.na(unknown_idx)) {
    lik_anc[, unknown_idx] <- 0.4
  }
  
  for (i in 1:n_states) {
    if (lik_anc[1, i] == 0) {
      lik_anc[, i] <- 0.01
    }
  }
  
  row_sums <- rowSums(lik_anc)
  for (i in 1:n_internal) {
    if (row_sums[i] > 0) {
      lik_anc[i, ] <- lik_anc[i, ] / row_sums[i]
    } else {
      if (!is.na(unknown_idx)) {
        lik_anc[i, ] <- 0.1 / (n_states - 1)
        lik_anc[i, unknown_idx] <- 0.9
      } else {
        lik_anc[i, ] <- 1 / n_states
      }
    }
  }
  
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
  if (!requireNamespace("ggtree", quietly = TRUE) || 
      !requireNamespace("RColorBrewer", quietly = TRUE) ||
      !requireNamespace("aplot", quietly = TRUE)) {
    log_message("Cannot create tree visualization - packages missing")
    return(FALSE)
  }
  
  library(ggtree)
  library(ggplot2)
  library(RColorBrewer)
  
  log_message("\nCreating tree visualization with state information...")
  
  n_tips <- length(tree$tip.label)
  n_internal <- tree$Nnode
  
  node_data <- data.frame(
    node = 1:(n_tips + n_internal),
    state = all_states,
    is_tip = 1:length(all_states) <= n_tips,
    label = c(tree$tip.label, rep("", n_internal)),
    stringsAsFactors = FALSE
  )
  
  unique_states <- unique(all_states)
  n_states <- length(unique_states)
  
  if (n_states <= 9) {
    state_colors <- brewer.pal(max(3, n_states), "Set1")
  } else {
    state_colors <- colorRampPalette(brewer.pal(9, "Set1"))(n_states)
  }
  names(state_colors) <- unique_states
  
  p <- ggtree(tree) %<+% node_data +
    geom_tippoint(aes(color = state), size = 2) +
    geom_tiplab(offset = 0.001, size = 2, align = TRUE) +
    geom_nodepoint(aes(color = state, size = !is_tip), alpha = 0.8) +
    scale_color_manual(values = state_colors, name = "Host") +
    scale_size_manual(values = c("TRUE" = 1, "FALSE" = 3), guide = "none") +
    theme(legend.position = "right") +
    labs(title = paste0(protein_name, " Tree with Host States"))
  
  rect_file <- file.path(output_dir, paste0(protein_name, "_tree_rectangular.pdf"))
  ggsave(rect_file, p, width = 12, height = max(8, n_tips/30), limitsize = FALSE)
  log_message("Saved rectangular tree to", rect_file)
  
  p2 <- ggtree(tree, layout = "circular") %<+% node_data +
    geom_tippoint(aes(color = state), size = 2) +
    geom_nodepoint(aes(color = state, size = !is_tip), alpha = 0.8) +
    scale_color_manual(values = state_colors, name = "Host") +
    scale_size_manual(values = c("TRUE" = 1, "FALSE" = 3), guide = "none") +
    theme(legend.position = "right") +
    labs(title = paste0(protein_name, " Circular Tree"))
  
  circ_file <- file.path(output_dir, paste0(protein_name, "_tree_circular.pdf"))
  ggsave(circ_file, p2, width = 10, height = 10)
  log_message("Saved circular tree to", circ_file)
  
  tryCatch({
    p_pie <- ggtree(tree) %<+% node_data
    
    p_pie <- p_pie + 
      geom_tippoint(aes(color = state), size = 2) +
      geom_tiplab(offset = 0.001, size = 2, align = TRUE) +
      scale_color_manual(values = state_colors, name = "Host")
    
    start_node <- n_tips + 1
    for (i in 1:nrow(ace_result$lik.anc)) {
      node_id <- start_node + i - 1
      probs <- ace_result$lik.anc[i,]
      
      if (sum(probs) > 0) {
        p_pie <- p_pie + nodepie(node_id, probs, colors = state_colors, size = 3)
      }
    }
    
    pie_file <- file.path(output_dir, paste0(protein_name, "_tree_with_pies.pdf"))
    ggsave(pie_file, p_pie, width = 12, height = max(8, n_tips/30), limitsize = FALSE)
    log_message("Saved tree with pie charts to", pie_file)
  }, error = function(e) {
    log_message("Error creating pie chart tree:", e$message)
  })
  
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
      site_rate_cols <- grep("rate|omega|dN/dS|omega", tolower(colnames(site_log)), value = TRUE)
      if (length(site_rate_cols) == 0) {
        site_rate_cols <- colnames(site_log)[-1]
        log_message("Using columns 2 onward as site rates")
      } else {
        log_message("Found rate columns using alternative patterns:", paste(site_rate_cols, collapse=", "))
      }
    } else {
      log_message("Found site rate columns:", paste(site_rate_cols, collapse=", "))
    }
    
    log_message("First few values of first rate column:")
    print(head(site_log[, site_rate_cols[1], drop = FALSE]))
    
    site_log_burned <- site_log[(floor(nrow(site_log) * burnin_fraction) + 1):nrow(site_log), ]
    site_means <- colMeans(site_log_burned[, site_rate_cols, drop = FALSE])
    
    log_message("Range of site means:", min(site_means), "to", max(site_means))
    
    site_indices <- 1:length(site_means)
    result <- data.frame(site = site_indices, rate = site_means, stringsAsFactors = FALSE)
    log_message("Created site rate data frame with", nrow(result), "rows")
    log_message("Site rate summary: min =", min(result$rate), ", max =", max(result$rate), ", mean =", mean(result$rate))
    
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
  
  strain_matching <- match_strains_to_meta(original_tip_labels, all_meta, name_map)
  
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
  
  if (length(unique(new_tip_labels)) < length(new_tip_labels)) {
    log_message("\nWARNING: Found duplicate tip labels after renaming!")
    dupe_table <- table(new_tip_labels)
    dupes <- names(dupe_table[dupe_table > 1])
    log_message("  Duplicate labels:", paste(dupes, collapse=", "))
    
    for (dupe in dupes) {
      dupe_indices <- which(new_tip_labels == dupe)
      log_message("  '", dupe, "' appears", length(dupe_indices), "times. Original labels:")
      
      for (i in seq_along(dupe_indices)) {
        if (i > 1) {
          new_tip_labels[dupe_indices[i]] <- paste0(dupe, "_dupe", i)
          log_message("    Renamed to: '", new_tip_labels[dupe_indices[i]], "'")
        }
      }
    }
  }
  
  cons_tree_binary$tip.label <- new_tip_labels
  
  tip_df <- data.frame(
    label = new_tip_labels, 
    Matched_Strain = strain_matching$Matched_Strain[match(new_tip_labels, strain_matching$Final_Accession)],
    stringsAsFactors = FALSE
  )
  
  annotated_tips <- left_join(tip_df, all_meta, by = c("label" = "Accession"))
  
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
  
  # ==== Ancestral State Reconstruction ====
  log_message("Performing ancestral state reconstruction with robust ACE approach...")
  ace_result <- run_ace_with_fixes(host_states, cons_tree_binary)
  
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
  
  n_internal_nodes <- cons_tree_binary$Nnode
  log_message("Tree has", n_internal_nodes, "internal nodes")
  ancestral_states <- rep("Unknown", n_internal_nodes)
  
  if (valid_ace_result) {
    tryCatch({
      if (nrow(ace_result$lik.anc) != n_internal_nodes) {
        log_message("WARNING: Number of rows in lik.anc (", nrow(ace_result$lik.anc), 
                    ") does not match number of internal nodes (", n_internal_nodes, ")")
        if (nrow(ace_result$lik.anc) > n_internal_nodes) {
          log_message("Truncating lik.anc to match the number of internal nodes")
          ace_result$lik.anc <- ace_result$lik.anc[1:n_internal_nodes, , drop = FALSE]
        } else {
          log_message("Extending lik.anc with Unknown states to match the number of internal nodes")
          unknown_col <- match("Unknown", colnames(ace_result$lik.anc))
          if (is.na(unknown_col)) {
            extension <- matrix(NA, nrow = n_internal_nodes - nrow(ace_result$lik.anc), 
                                ncol = ncol(ace_result$lik.anc),
                                dimnames = list(NULL, colnames(ace_result$lik.anc)))
          } else {
            extension <- matrix(0, nrow = n_internal_nodes - nrow(ace_result$lik.anc), 
                                ncol = ncol(ace_result$lik.anc),
                                dimnames = list(NULL, colnames(ace_result$lik.anc)))
            extension[, unknown_col] <- 1
          }
          ace_result$lik.anc <- rbind(ace_result$lik.anc, extension)
        }
      }
      
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
  edge_df$segment <- segment_name
  edge_df$protein <- protein_name
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
  
  # Collect branch rates for overall analysis
  all_branch_rates <<- rbind(all_branch_rates, edge_df_nonzero %>% 
                               dplyr::select(segment, protein, rate, group))
  
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
  
  wilcox_result <- list(p.value = NA)
  ks_result <- list(p.value = NA)
  if (!is.null(branch_rate_summary) && n_shifts > 0) {
    print(branch_rate_summary)
    wilcox_result <- tryCatch({
      wilcox.test(edge_df_nonzero$rate[edge_df_nonzero$host_switch],
                  edge_df_nonzero$rate[!edge_df_nonzero$host_switch])
    }, error = function(e) {
      log_message("Error in Wilcoxon test:", e$message)
      list(p.value = NA)
    })
    ks_result <- tryCatch({
      ks.test(edge_df_nonzero$rate[edge_df_nonzero$host_switch],
              edge_df_nonzero$rate[!edge_df_nonzero$host_switch])
    }, error = function(e) {
      log_message("Error in KS test:", e$message)
      list(p.value = NA)
    })
    log_message("Wilcoxon test p-value (non-zero branches):", wilcox_result$p.value)
    log_message("KS test p-value (non-zero branches):", ks_result$p.value)
  } else {
    log_message("Insufficient data for branch rate comparison.")
  }
  
  # ==== Site Rate Analysis ====
  log_message("Performing site rate analysis...")
  n_sites <- nrow(site_rates_log)
  
  mean_rate_val <- mean(site_rates_log$rate)
  sd_rate_val <- sd(site_rates_log$rate)
  median_rate_val <- median(site_rates_log$rate)
  mad_rate_val <- mad(site_rates_log$rate)
  
  log_message("Site rate statistics:")
  log_message("  Mean:", mean_rate_val)
  log_message("  Standard Deviation:", sd_rate_val)
  log_message("  Coefficient of Variation:", sd_rate_val/mean_rate_val * 100, "%")
  log_message("  Median:", median_rate_val)
  log_message("  Median Absolute Deviation:", mad_rate_val)
  
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
  
  shapiro_test <- tryCatch({
    shapiro.test(site_rates_log$rate)
  }, error = function(e) {
    log_message("Normality test failed:", e$message)
    list(p.value = 0)
  })
  
  log_message("Shapiro-Wilk normality test: W =", shapiro_test$statistic, 
              ", p-value =", shapiro_test$p.value)
  is_normal <- shapiro_test$p.value >= 0.05
  
  site_rates_log$z_score <- (site_rates_log$rate - mean_rate_val) / sd_rate_val
  site_rates_log$robust_z <- (site_rates_log$rate - median_rate_val) / mad_rate_val
  
  log_message("Identifying sites under selection...")
  
  z_score_col <- if(is_normal) "z_score" else "robust_z"
  reference_value <- if(is_normal) mean_rate_val else median_rate_val
  spread_value <- if(is_normal) sd_rate_val else mad_rate_val
  
  log_message("Using", if(is_normal) "standard" else "robust", "statistics for site selection")
  
  if (spread_value <= 1e-6 || spread_value/reference_value < 0.001) {
    log_message("WARNING: Very low variation. All sites appear to evolve at similar rates.")
    fast_evolving_sites <- data.frame()
    very_fast_sites <- data.frame()
    slow_evolving_sites <- data.frame()
  } else if (any(is.nan(site_rates_log[[z_score_col]])) || any(is.infinite(site_rates_log[[z_score_col]]))) {
    log_message("WARNING: Invalid z-scores detected. Using percentile-based selection instead.")
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
    fast_evolving_sites <- site_rates_log %>% filter(.data[[z_score_col]] > 1.96) %>% arrange(desc(rate))
    very_fast_sites <- site_rates_log %>% filter(.data[[z_score_col]] > 2.58) %>% arrange(desc(rate))
    slow_evolving_sites <- site_rates_log %>% filter(.data[[z_score_col]] < -1.96) %>% arrange(rate)
    
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
      comparison_val <- top_slow$rate[i] / reference_value
      z_val <- if("robust_z" %in% names(top_slow)) top_slow$robust_z[i] else top_slow$z_score[i]
      log_message(sprintf("  Site %d: rate = %.4f (%.2fx slower, z = %.2f)",
                          top_slow$site[i], top_slow$rate[i], comparison_val, z_val))
    }
    write_csv(slow_evolving_sites, slow_sites_file)
    log_message("Saved slow-evolving sites to", slow_sites_file)
  }
  
  very_fast_file <- if(nrow(very_fast_sites) > 0) 
    file.path(output_dir, paste0(protein_name_clean, "_very_fast_sites.csv")) else NA
  if (!is.na(very_fast_file)) {
    log_message("Top 5 very fast-evolving sites:")
    top_very_fast <- head(very_fast_sites, 5)
    for (i in 1:nrow(top_very_fast)) {
      comparison_val <- top_very_fast$rate[i] / reference_value
      z_val <- if("robust_z" %in% names(top_very_fast)) top_very_fast$robust_z[i] else top_very_fast$z_score[i]
      log_message(sprintf("  Site %d: rate = %.4f (%.2fx faster, z = %.2f)",
                          top_very_fast$site[i], top_very_fast$rate[i], comparison_val, z_val))
    }
    write_csv(very_fast_sites, very_fast_file)
    log_message("Saved very fast-evolving sites to", very_fast_file)
  }
  
  # ==== Save Summary ====
  summary_text <- paste0(
    "\nProtein: ", protein_name, " (Segment: ", segment_name, ")\n",
    "Number of sites: ", n_sites, "\n",
    "Site rate mean: ", round(mean_rate_val, 4), "\n",
    "Site rate SD: ", round(sd_rate_val, 4), "\n",
    "Site rate median: ", round(median_rate_val, 4), "\n",
    "Site rate MAD: ", round(mad_rate_val, 4), "\n",
    "Fast-evolving sites: ", nrow(fast_evolving_sites), "\n",
    "Very fast sites: ", nrow(very_fast_sites), "\n",
    "Slow-evolving sites: ", nrow(slow_evolving_sites), "\n",
    "Number of host shifts: ", n_shifts, "\n",
    "Wilcoxon test p-value: ", sprintf("%.4f", wilcox_result$p.value), "\n",
    "KS test p-value: ", sprintf("%.4f", ks_result$p.value), "\n"
  )
  
  cat(summary_text, file = summary_all_file, append = TRUE)
  log_message("Appended summary for", protein_name, "to", summary_all_file)
  
  # ==== Plot Site Rates ====
  p_site_rates <- ggplot(site_rates_log, aes(x = site, y = rate)) +
    geom_point(aes(color = .data[[z_score_col]]), size = 2) +
    scale_color_gradient2(low = "blue", mid = "grey", high = "red", midpoint = 0) +
    labs(title = paste0(protein_name, " Site Rates"),
         x = "Site", y = "Rate") +
    theme_minimal()
  
  site_plot_file <- file.path(output_dir, paste0(protein_name_clean, "_site_rates.pdf"))
  ggsave(site_plot_file, p_site_rates, width = 10, height = 6)
  log_message("Saved site rates plot to", site_plot_file)
  
  return(list(
    protein = protein_name,
    segment = segment_name,
    n_sites = n_sites,
    fast_sites = nrow(fast_evolving_sites),
    very_fast_sites = nrow(very_fast_sites),
    slow_sites = nrow(slow_evolving_sites),
    n_shifts = n_shifts,
    wilcox_p = wilcox_result$p.value,
    ks_p = ks_result$p.value
  ))
}

# ==== Process All Segments and Proteins ====
log_message("Processing all segment directories...")
results_list <- list()

for (seg_dir in seg_dirs) {
  log_message("\nProcessing segment directory:", seg_dir)
  tree_files <- list.files(seg_dir, pattern = "_tree\\.nex$", full.names = TRUE)
  
  if (length(tree_files) == 0) {
    log_message("No tree files found in", seg_dir)
    next
  }
  
  log_message("Found", length(tree_files), "tree files in", seg_dir)
  
  for (tree_file in tree_files) {
    result <- tryCatch({
      process_protein(tree_file)
    }, error = function(e) {
      log_message("Error processing", tree_file, ":", e$message)
      NULL
    })
    
    if (!is.null(result)) {
      results_list <- append(results_list, list(result))
    }
  }
}

# ==== Overall Analysis Across All Proteins ====
log_message("\n===== PERFORMING OVERALL ANALYSIS =====")

if (nrow(all_branch_rates) == 0) {
  log_message("No branch rates collected. Skipping overall analysis.")
} else {
  log_message("Total branch rates collected:", nrow(all_branch_rates))
  log_message("Number of host shifts:", sum(all_branch_rates$group == "Host Shift"))
  log_message("Number of non-shifts:", sum(all_branch_rates$group == "No Shift"))
  
  overall_summary <- all_branch_rates %>% 
    group_by(group) %>%
    summarize(mean_rate = mean(rate),
              median_rate = median(rate),
              min_rate = min(rate),
              max_rate = max(rate),
              n_branches = n())
  
  log_message("Overall branch rate summary:")
  print(overall_summary)
  
  overall_wilcox <- tryCatch({
    wilcox.test(rate ~ group, data = all_branch_rates)
  }, error = function(e) {
    log_message("Error in overall Wilcoxon test:", e$message)
    list(p.value = NA)
  })
  
  overall_ks <- tryCatch({
    ks.test(all_branch_rates$rate[all_branch_rates$group == "Host Shift"],
            all_branch_rates$rate[all_branch_rates$group == "No Shift"])
  }, error = function(e) {
    log_message("Error in overall KS test:", e$message)
    list(p.value = NA)
  })
  
  log_message("Overall Wilcoxon test p-value:", overall_wilcox$p.value)
  log_message("Overall KS test p-value:", overall_ks$p.value)
  
  overall_summary_text <- paste0(
    "\n===== OVERALL ANALYSIS =====\n",
    "Total branches analyzed: ", nrow(all_branch_rates), "\n",
    "Host shifts: ", sum(all_branch_rates$group == "Host Shift"), "\n",
    "Non-shifts: ", sum(all_branch_rates$group == "No Shift"), "\n",
    "Wilcoxon test p-value: ", sprintf("%.4f", overall_wilcox$p.value), "\n",
    "KS test p-value: ", sprintf("%.4f", overall_ks$p.value), "\n"
  )
  
  cat(overall_summary_text, file = summary_all_file, append = TRUE)
  log_message("Appended overall summary to", summary_all_file)
  
  # ==== Plot Overall Branch Rates ====
  p_overall_rates <- ggplot(all_branch_rates, aes(x = group, y = rate, fill = group)) +
    geom_violin(trim = FALSE) +
    geom_boxplot(width = 0.2, fill = "white") +
    labs(title = "Branch Rates Across All Proteins",
         x = "Group", y = "Rate") +
    theme_minimal()
  
  overall_plot_file <- file.path(main_output_dir, "overall_branch_rates.pdf")
  ggsave(overall_plot_file, p_overall_rates, width = 8, height = 6)
  log_message("Saved overall branch rates plot to", overall_plot_file)
}

# ==== Summarize Results ====
results_df <- bind_rows(results_list)
if (nrow(results_df) > 0) {
  log_message("\n===== FINAL SUMMARY =====")
  log_message("Processed", nrow(results_df), "proteins across", length(unique(results_df$segment)), "segments")
  log_message("Total fast-evolving sites:", sum(results_df$fast_sites))
  log_message("Total very fast sites:", sum(results_df$very_fast_sites))
  log_message("Total slow-evolving sites:", sum(results_df$slow_sites))
  log_message("Total host shifts:", sum(results_df$n_shifts))
  
  results_file <- file.path(main_output_dir, "protein_results_summary.csv")
  write_csv(results_df, results_file)
  log_message("Saved results summary to", results_file)
}

log_message("===== ANALYSIS COMPLETE =====")


