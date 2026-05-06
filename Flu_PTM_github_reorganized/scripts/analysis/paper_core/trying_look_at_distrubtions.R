# ==== Host Shift Branch Rate Distribution Analysis with Ridgeline Plot ====
# This script analyzes branch rate distributions for different types of host shifts
# across influenza segments and produces a ridgeline plot to visualize them.
# It uses ancestral state reconstruction to identify host shifts and categorizes them
# by transition type (e.g., Human to Swine).

# ==== Load Required Libraries ====
library(ape)
library(phangorn)
library(readr)
library(dplyr)
library(stringr)
library(ggplot2)
library(ggridges)  # For ridgeline plots
library(jsonlite)
library(fs)
library(stats)  # For KS test

# Install ggridges if not available
if (!requireNamespace("ggridges", quietly = TRUE)) {
  install.packages("ggridges")
}

# ==== Global Settings ====
burnin_fraction <- 0.25
base_dir <- "processed_segments_fixed_final"
main_output_dir <- "host_shift_analysis_results"
debug_output_dir <- file.path(main_output_dir, "debug_output")
dir.create(main_output_dir, showWarnings = FALSE)
dir.create(debug_output_dir, showWarnings = FALSE)

# Debug log file
debug_log <- file.path(debug_output_dir, "host_shift_debug.log")
file.create(debug_log)
log_message <- function(...) {
  message <- paste0(...)
  cat(message, "\n")
  cat(message, "\n", file = debug_log, append = TRUE)
}

log_message("===== STARTING HOST SHIFT BRANCH RATE ANALYSIS =====")
log_message("Base directory:", base_dir)
log_message("Output directory:", main_output_dir)

# Initialize data frame for branch rates by host shift type
all_host_shift_rates <- data.frame()

# ==== Find All Segment Directories ====
log_message("Searching for segment directories...")
seg_dirs <- list.dirs(base_dir, recursive = FALSE)
seg_dirs <- seg_dirs[grepl("^seg", basename(seg_dirs))]
if (length(seg_dirs) == 0) {
  stop("No segment directories found. Check that base_dir is correct.")
}
log_message("Found", length(seg_dirs), "segment directories:", paste(basename(seg_dirs), collapse = ", "), "\n")

# ==== File Matching Function ====
find_matching_files <- function(tree_file) {
  dir_path <- dirname(tree_file)
  all_files <- list.files(dir_path, full.names = TRUE)
  
  tree_basename <- basename(tree_file)
  raw_name <- sub("_tree\\.nex$", "", tree_basename, ignore.case = TRUE)
  protein_name <- sub("_protein$", "", raw_name, ignore.case = TRUE)
  log_message("Extracted protein name:", protein_name, "from file:", tree_basename)
  
  branch_regex <- paste0("^", protein_name, "(_protein)?_branch_rates\\.log$")
  mapping_regex <- paste0("^", protein_name, "(_protein)?_id_mapping\\.json$")
  
  branch_log_file <- NA
  if (length(grep(branch_regex, basename(all_files), ignore.case = TRUE)) > 0) {
    branch_candidates <- all_files[grep(branch_regex, basename(all_files), ignore.case = TRUE)]
    branch_log_file <- branch_candidates[1]
  }
  
  mapping_file <- NA
  if (length(grep(mapping_regex, basename(all_files), ignore.case = TRUE)) > 0) {
    mapping_candidates <- all_files[grep(mapping_regex, basename(all_files), ignore.case = TRUE)]
    mapping_file <- mapping_candidates[1]
  }
  
  log_message("Found matching files for", protein_name, ":")
  log_message("  Tree file:         ", tree_file)
  log_message("  Branch rate file:  ", ifelse(is.na(branch_log_file), "Not found", branch_log_file))
  log_message("  Mapping file:      ", ifelse(is.na(mapping_file), "Not found", mapping_file))
  
  list(
    protein_name = protein_name,
    tree_file = tree_file,
    branch_log_file = branch_log_file,
    mapping_file = mapping_file
  )
}

# ==== Label Matching Function ====
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
                          grepl("^json_", result$match_method), na.rm = TRUE)
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
  
  log_message("\nMatching results summary:")
  log_message("- Total tips:", nrow(result))
  log_message("- Successfully matched:", sum(!is.na(result$Final_Accession)))
  log_message("- Unmatched:", sum(is.na(result$Final_Accession)))
  
  return(result)
}

# ==== Fix Tree for ACE Compatibility ====
fix_tree_for_ace <- function(tree) {
  zero_branches <- sum(tree$edge.length <= 1e-6, na.rm = TRUE)
  na_branches <- sum(is.na(tree$edge.length))
  
  if (zero_branches > 0 || na_branches > 0) {
    log_message("Setting small positive values for zero/NA branch lengths")
    tree$edge.length[tree$edge.length <= 1e-6 | is.na(tree$edge.length)] <- 0.01
  }
  
  neg_branches <- sum(tree$edge.length < 0, na.rm = TRUE)
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
    tree <- root(tree, outgroup = 1, resolve.root = TRUE)
  }
  
  return(tree)
}

# ==== Robust ACE Function ====
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
    host_states <- factor(host_states)
  }
  
  if (!"Unknown" %in% levels(host_states)) {
    host_states <- factor(host_states, levels = c(levels(host_states), "Unknown"))
  }
  
  tryCatch({
    log_message("Attempting ACE with default parameters")
    ace_result <- ace(host_states, tree, type = "discrete")
    return(ace_result)
  }, error = function(e) {
    log_message("ACE failed:", e$message)
    log_message("Trying ACE with marginal=TRUE")
    tryCatch({
      ace_result <- ace(host_states, tree, type = "discrete", marginal = TRUE)
      return(ace_result)
    }, error = function(e2) {
      log_message("Second ACE attempt failed:", e2$message)
      log_message("Creating dummy ACE result")
      n_internal <- tree$Nnode
      state_levels <- levels(host_states)
      lik_anc <- matrix(1/length(state_levels), nrow = n_internal, ncol = length(state_levels),
                        dimnames = list(NULL, state_levels))
      return(list(lik.anc = lik_anc))
    })
  })
}

# ==== Process Single Protein ====
process_protein <- function(tree_file) {
  log_message("Processing tree file:", tree_file)
  file_matches <- find_matching_files(tree_file)
  protein_name <- file_matches$protein_name
  branch_log_file <- file_matches$branch_log_file
  mapping_file <- file_matches$mapping_file
  segment_name <- basename(dirname(tree_file))
  
  if (is.na(branch_log_file)) {
    log_message("Warning: Branch rate file not found for", protein_name, "in", segment_name, ". Skipping.")
    return(NULL)
  }
  
  protein_name_clean <- gsub("[^a-zA-Z0-9]", "_", protein_name)
  output_dir <- file.path(main_output_dir, paste0(segment_name, "_", protein_name_clean))
  dir.create(output_dir, showWarnings = FALSE)
  
  # Load Tree
  trees <- tryCatch({
    read.nexus(tree_file)
  }, error = function(e) {
    log_message("Error loading tree:", e$message)
    return(NULL)
  })
  if (is.null(trees)) {
    return(NULL)
  }
  
  if (inherits(trees, "multiPhylo")) {
    burnin_idx <- floor(length(trees) * burnin_fraction) + 1
    trees_to_use <- trees[burnin_idx:length(trees)]
    cons_tree <- consensus(trees_to_use)
  } else {
    cons_tree <- trees
  }
  
  if (is.null(cons_tree$edge.length)) {
    cons_tree$edge.length <- rep(1, nrow(cons_tree$edge))
  }
  
  cons_tree <- fix_tree_for_ace(cons_tree)
  cons_tree <- midpoint.root(cons_tree)
  cons_tree <- multi2di(cons_tree)
  
  # Load Mapping File
  name_map <- list()
  if (!is.na(mapping_file)) {
    name_map <- tryCatch({
      fromJSON(mapping_file)
    }, error = function(e) {
      log_message("Error loading mapping file:", e$message)
      list()
    })
  }
  
  # Load Branch Rates
  branch_rates_log <- tryCatch({
    branch_log <- read.table(branch_log_file, header = TRUE)
    branch_log_burned <- branch_log[(floor(nrow(branch_log) * burnin_fraction) + 1):nrow(branch_log), ]
    branch_rate_cols <- grep("branch_rates", tolower(colnames(branch_log_burned)), value = TRUE)
    if (length(branch_rate_cols) == 0) {
      stop("No branch rate columns found in", branch_log_file)
    }
    colMeans(branch_log_burned[, branch_rate_cols, drop = FALSE])
  }, error = function(e) {
    log_message("Error loading branch rates:", e$message)
    rep(1, nrow(cons_tree$edge))
  })
  
  rate_vector <- branch_rates_log
  if (length(rate_vector) != nrow(cons_tree$edge)) {
    rate_vector <- rep(mean(rate_vector), nrow(cons_tree$edge))
  }
  
  zero_branches <- rate_vector <= 1e-10
  
  # Load Metadata
  meta_h1n1 <- tryCatch(
    read_csv("C:/Users/Ashley/Documents/flu_project/H1N1/H1N1_seq_info.csv") %>% mutate(Strain = "H1N1"),
    error = function(e) { 
      log_message("Error loading H1N1 metadata:", e$message)
      data.frame()
    }
  )
  
  meta_h5n1 <- tryCatch(
    read_csv("C:/Users/Ashley/Documents/flu_project/H5N1/H5N1_seq_info.csv") %>% mutate(Strain = "H5N1"),
    error = function(e) { 
      log_message("Error loading H5N1 metadata:", e$message)
      data.frame()
    }
  )
  
  meta_h7n9 <- tryCatch(
    read_csv("C:/Users/Ashley/Documents/flu_project/H7N9/H7N9_seq_info.csv") %>% mutate(Strain = "H7N9"),
    error = function(e) { 
      log_message("Error loading H7N9 metadata:", e$message)
      data.frame()
    }
  )
  
  all_meta <- bind_rows(meta_h1n1, meta_h5n1, meta_h7n9)
  if (nrow(all_meta) > 0) {
    all_meta$Accession <- str_trim(all_meta$Accession)
    all_meta$Accession <- gsub("\\.1$", "", all_meta$Accession)
  }
  
  # Match Tip Labels
  original_tip_labels <- cons_tree$tip.label
  strain_matching <- match_strains_to_meta(original_tip_labels, all_meta, name_map)
  
  new_tip_labels <- sapply(original_tip_labels, function(tip) {
    idx <- which(strain_matching$label == tip)
    if (length(idx) > 0 && !is.na(strain_matching$Final_Accession[idx])) {
      return(strain_matching$Final_Accession[idx])
    } else if (length(name_map) > 0 && tip %in% names(name_map)) {
      return(gsub("^(H1N1__|H5N1__|H7N9__)", "", name_map[[tip]]))
    } else {
      return(gsub("^(H1N1__|H5N1__|H7N9__)", "", tip))
    }
  }, USE.NAMES = FALSE)
  
  if (length(unique(new_tip_labels)) < length(new_tip_labels)) {
    dupe_table <- table(new_tip_labels)
    dupes <- names(dupe_table[dupe_table > 1])
    for (dupe in dupes) {
      dupe_indices <- which(new_tip_labels == dupe)
      for (i in seq_along(dupe_indices)[-1]) {
        new_tip_labels[dupe_indices[i]] <- paste0(dupe, "_dupe", i)
      }
    }
  }
  
  cons_tree$tip.label <- new_tip_labels
  
  tip_df <- data.frame(
    label = new_tip_labels,
    Matched_Strain = strain_matching$Matched_Strain[match(new_tip_labels, strain_matching$Final_Accession)],
    stringsAsFactors = FALSE
  )
  
  annotated_tips <- left_join(tip_df, all_meta, by = c("label" = "Accession"))
  
  idx_missing <- which(is.na(annotated_tips$Strain) & !is.na(annotated_tips$Matched_Strain))
  if (length(idx_missing) > 0) {
    annotated_tips$Strain[idx_missing] <- annotated_tips$Matched_Strain[idx_missing]
  }
  
  host_states <- annotated_tips$Host
  host_states[is.na(host_states)] <- "Unknown"
  names(host_states) <- new_tip_labels
  host_states <- factor(host_states)
  
  # Ancestral State Reconstruction
  ace_result <- run_ace_with_fixes(host_states, cons_tree)
  
  n_internal_nodes <- cons_tree$Nnode
  ancestral_states <- rep("Unknown", n_internal_nodes)
  
  if (!is.null(ace_result$lik.anc) && is.matrix(ace_result$lik.anc)) {
    if (nrow(ace_result$lik.anc) != n_internal_nodes) {
      ace_result$lik.anc <- ace_result$lik.anc[1:n_internal_nodes, , drop = FALSE]
    }
    
    ancestral_states <- sapply(1:nrow(ace_result$lik.anc), function(i) {
      p <- ace_result$lik.anc[i, ]
      if (all(is.na(p)) || sum(p, na.rm = TRUE) == 0) {
        return("Unknown")
      } else {
        colnames(ace_result$lik.anc)[which.max(p)]
      }
    })
  }
  
  all_node_nums <- 1:(length(cons_tree$tip.label) + cons_tree$Nnode)
  all_states <- rep("", length(all_node_nums))
  names(all_states) <- all_node_nums
  all_states[1:length(cons_tree$tip.label)] <- as.character(host_states)
  all_states[(length(cons_tree$tip.label)+1):length(all_states)] <- ancestral_states
  
  # Identify Host Shifts
  edge_df <- data.frame(
    parent = cons_tree$edge[,1],
    child = cons_tree$edge[,2],
    parent_state = all_states[cons_tree$edge[,1]],
    child_state = all_states[cons_tree$edge[,2]],
    rate = rate_vector,
    zero_branch = zero_branches
  )
  
  edge_df$host_shift_type <- ifelse(
    edge_df$parent_state != edge_df$child_state &
      edge_df$parent_state != "Unknown" &
      edge_df$child_state != "Unknown",
    paste(edge_df$parent_state, "→", edge_df$child_state),
    "No Shift"
  )
  
  edge_df_nonzero <- edge_df[!edge_df$zero_branch, ]
  log_message("Identified", sum(edge_df_nonzero$host_shift_type != "No Shift"), "host shifts")
  
  # Collect rates for host shift types
  shift_rates <- edge_df_nonzero %>%
    filter(host_shift_type != "No Shift") %>%
    mutate(segment = segment_name, protein = protein_name) %>%
    dplyr::select(segment, protein, rate, host_shift_type)
  
  all_host_shift_rates <<- rbind(all_host_shift_rates, shift_rates)
  
  # Save host shift counts
  host_transitions <- edge_df_nonzero %>%
    filter(host_shift_type != "No Shift") %>%
    count(host_shift_type, sort = TRUE)
  
  transition_file <- file.path(output_dir, paste0(protein_name_clean, "_host_shift_types.csv"))
  write_csv(host_transitions, transition_file)
  log_message("Saved host shift types to", transition_file)
  
  return(list(
    protein = protein_name,
    segment = segment_name,
    n_shifts = sum(edge_df_nonzero$host_shift_type != "No Shift")
  ))
}

# ==== Process All Segments and Proteins ====
log_message("Processing all segment directories...")
results_list <- list()

for (seg_dir in seg_dirs) {
  tree_files <- list.files(seg_dir, pattern = "_tree\\.nex$", full.names = TRUE)
  if (length(tree_files) == 0) {
    log_message("No tree files found in", seg_dir)
    next
  }
  
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

# ==== Ridgeline Plot and Statistical Analysis ====
log_message("\n===== CREATING RIDGELINE PLOT AND STATISTICAL ANALYSIS =====")

if (nrow(all_host_shift_rates) > 0) {
  # Summarize host shift types
  shift_summary <- all_host_shift_rates %>%
    group_by(host_shift_type) %>%
    summarize(n_branches = n(),
              mean_rate = mean(rate),
              median_rate = median(rate))
  
  log_message("Host shift types summary:")
  print(shift_summary)
  
  # Perform KS tests between pairs of host shift types
  shift_types <- unique(all_host_shift_rates$host_shift_type)
  if (length(shift_types) > 1) {
    ks_results <- data.frame(
      comparison = character(),
      p_value = numeric(),
      stringsAsFactors = FALSE
    )
    
    for (i in 1:(length(shift_types)-1)) {
      for (j in (i+1):length(shift_types)) {
        type1 <- shift_types[i]
        type2 <- shift_types[j]
        rates1 <- all_host_shift_rates$rate[all_host_shift_rates$host_shift_type == type1]
        rates2 <- all_host_shift_rates$rate[all_host_shift_rates$host_shift_type == type2]
        
        if (length(rates1) > 10 && length(rates2) > 10) {
          ks_test <- tryCatch({
            ks.test(rates1, rates2)
          }, error = function(e) {
            log_message("KS test failed for", type1, "vs", type2, ":", e$message)
            list(p.value = NA)
          })
          
          ks_results <- rbind(ks_results, data.frame(
            comparison = paste(type1, "vs", type2),
            p_value = ks_test$p.value
          ))
        }
      }
    }
    
    log_message("KS test results for host shift type comparisons:")
    print(ks_results)
    
    ks_file <- file.path(main_output_dir, "host_shift_ks_tests.csv")
    write_csv(ks_results, ks_file)
    log_message("Saved KS test results to", ks_file)
  }
  
  # Create ridgeline plot
  p_ridgeline <- ggplot(all_host_shift_rates, aes(x = rate, y = host_shift_type, fill = host_shift_type)) +
    geom_density_ridges(scale = 3, alpha = 0.6) +
    labs(title = "Branch Rate Distributions by Host Shift Type",
         x = "Branch Rate", y = "Host Shift Type") +
    theme_minimal() +
    theme(legend.position = "none")
  
  ridge_plot_file <- file.path(main_output_dir, "host_shift_ridgeline.pdf")
  ggsave(ridge_plot_file, p_ridgeline, width = 10, height = max(6, length(unique(all_host_shift_rates$host_shift_type)) * 0.5))
  log_message("Saved ridgeline plot to", ridge_plot_file)
} else {
  log_message("No host shifts detected. Skipping ridgeline plot and statistical analysis.")
}

# ==== Final Summary ====
results_df <- bind_rows(results_list)
if (nrow(results_df) > 0) {
  log_message("\n===== FINAL SUMMARY =====")
  log_message("Processed", nrow(results_df), "proteins across", length(unique(results_df$segment)), "segments")
  log_message("Total host shifts:", sum(results_df$n_shifts))
  
  results_file <- file.path(main_output_dir, "host_shift_results_summary.csv")
  write_csv(results_df, results_file)
  log_message("Saved results summary to", results_file)
}

log_message("===== ANALYSIS COMPLETE =====")


# Replace the original ridgeline plot code with this:

# Filter for top 5 host shift types by number of shifts
top_shift_types <- all_host_shift_rates %>%
  group_by(host_shift_type) %>%
  summarize(n_shifts = n()) %>%
  arrange(desc(n_shifts)) %>%
  slice_head(n = 5) %>%
  pull(host_shift_type)

# [Previous script sections unchanged up to the analysis section]
# Replace the "Ridgeline Plot and Statistical Analysis" section with this:

# ==== Per-Protein Ridgeline Plots and Statistical Analysis ====
log_message("\n===== PER-PROTEIN ANALYSIS =====")

if (nrow(all_host_shift_rates) > 0) {
  # Group by protein and analyze
  proteins <- unique(all_host_shift_rates$protein)
  log_message("Analyzing", length(proteins), "proteins:", paste(proteins, collapse = ", "))
  
  for (prot in proteins) {
    log_message("\n=== Protein:", prot, "===")
    
    # Subset data for this protein
    prot_data <- all_host_shift_rates %>% filter(protein == prot)
    
    # Summarize host shift types for this protein
    shift_summary <- prot_data %>%
      group_by(host_shift_type) %>%
      summarize(
        n_branches = n(),
        mean_rate = mean(rate),
        median_rate = median(rate),
        sd_rate = sd(rate)
      ) %>%
      arrange(desc(n_branches))
    
    log_message("Host shift types summary for", prot, ":")
    print(shift_summary)
    
    # Filter for top 5 host shift types by number of shifts for this protein
    top_shift_types <- shift_summary %>%
      slice_head(n = 5) %>%
      pull(host_shift_type)
    
    # Subset data to include top 5 host shift types and add all shifts
    plot_data <- prot_data %>%
      filter(host_shift_type %in% top_shift_types) %>%
      bind_rows(
        prot_data %>%
          mutate(host_shift_type = "All Shifts")
      )
    
    # Create ridgeline plot for this protein
    p_ridgeline <- ggplot(plot_data, aes(x = rate, y = host_shift_type, fill = host_shift_type)) +
      geom_density_ridges(scale = 1, alpha = 0.6) +
      geom_vline(xintercept = median(prot_data$rate), linetype = "dashed", color = "black") +
      labs(title = paste("Branch Rate Distributions for", prot),
           x = "Branch Rate", y = "Host Shift Type") +
      theme_minimal() +
      theme(legend.position = "none")
    
    print(p_ridgeline)  # Display in RStudio plot window
  }
  
  # Overall KS tests (unchanged)
  log_message("\n===== OVERALL KS TEST RESULTS =====")
  shift_types <- unique(all_host_shift_rates$host_shift_type)
  if (length(shift_types) > 1) {
    ks_results <- data.frame(
      comparison = character(),
      p_value = numeric(),
      stringsAsFactors = FALSE
    )
    
    for (i in 1:(length(shift_types)-1)) {
      for (j in (i+1):length(shift_types)) {
        type1 <- shift_types[i]
        type2 <- shift_types[j]
        rates1 <- all_host_shift_rates$rate[all_host_shift_rates$host_shift_type == type1]
        rates2 <- all_host_shift_rates$rate[all_host_shift_rates$host_shift_type == type2]
        
        if (length(rates1) > 3 && length(rates2) > 3) {
          ks_test <- tryCatch({
            ks.test(rates1, rates2)
          }, error = function(e) {
            log_message("KS test failed for", type1, "vs", type2, ":", e$message)
            list(p.value = NA)
          })
          
          ks_results <- rbind(ks_results, data.frame(
            comparison = paste(type1, "vs", type2),
            p_value = ks_test$p.value
          ))
        }
      }
    }
    
    log_message("KS test results for host shift type comparisons:")
    print(ks_results)
  } else {
    log_message("Not enough host shift types for KS tests.")
  }
  
  # Overall summary
  overall_summary <- all_host_shift_rates %>%
    group_by(host_shift_type) %>%
    summarize(
      n_branches = n(),
      mean_rate = mean(rate),
      median_rate = median(rate),
      sd_rate = sd(rate)
    ) %>%
    arrange(desc(n_branches))
  
  log_message("\n===== OVERALL HOST SHIFT SUMMARY =====")
  print(overall_summary)
} else {
  log_message("No host shifts detected. No plots or summaries generated.")
}

# ==== Final Summary ====
results_df <- bind_rows(results_list)
if (nrow(results_df) > 0) {
  log_message("\n===== FINAL SUMMARY =====")
  log_message("Processed", nrow(results_df), "proteins across", length(unique(results_df$segment)), "segments")
  log_message("Total host shifts:", sum(results_df$n_shifts))
}

log_message("===== ANALYSIS COMPLETE =====")












# [Previous script sections unchanged up to the analysis section]
# Replace the "Ridgeline Plot and Statistical Analysis" section with this:


# [Previous script sections unchanged up to the analysis section]
# Replace the "Robust Linear Regression and Per-Protein Analysis" section with this:

# ==== Mixed-Effects Model and Outlier Detection Analysis ====
log_message("\n===== MIXED-EFFECTS MODEL AND OUTLIER DETECTION =====")

if (nrow(all_host_shift_rates) > 0) {
  # Load required libraries
  library(lme4)
  library(ape)
  library(emmeans)
  library(Biostrings)
  library(robust)
  
  # Filter host shift types with sufficient branches (>=10)
  log_message("Filtering host shift types with at least 10 branches...")
  valid_shifts <- all_host_shift_rates %>%
    group_by(host_shift_type) %>%
    summarize(n_branches = n()) %>%
    filter(n_branches >= 10) %>%
    pull(host_shift_type)
  
  filtered_data <- all_host_shift_rates %>%
    filter(host_shift_type %in% valid_shifts)
  
  log_message("Retained", length(valid_shifts), "host shift types:", paste(valid_shifts, collapse = ", "))
  log_message("Filtered data has", nrow(filtered_data), "branches across", length(unique(filtered_data$protein)), "proteins")
  
  if (nrow(filtered_data) == 0) {
    log_message("No host shift types with sufficient branches. Skipping analysis.")
  } else {
    # Add tree ID for random effects (assuming segment + protein uniquely identifies a tree)
    filtered_data$tree_id <- paste(filtered_data$segment, filtered_data$protein, sep = "_")
    
    # Fit mixed-effects model with phylogenetic correlation
    log_message("Fitting mixed-effects model...")
    # Placeholder for phylogenetic correlation matrix (requires tree)
    # For simplicity, we use random effects for tree_id
    lmm_model <- tryCatch({
      lmer(rate ~ protein + host_shift_type + (1 | tree_id), 
           data = filtered_data, 
           control = lmerControl(optimizer = "bobyqa"))
    }, error = function(e) {
      log_message("Error fitting mixed-effects model:", e$message)
      log_message("Trying without protein effect...")
      lmer(rate ~ host_shift_type + (1 | tree_id), 
           data = filtered_data, 
           control = lmerControl(optimizer = "bobyqa"))
    })
    
    # Print model summary
    log_message("Mixed-effects model results:")
    print(summary(lmm_model))
    
    # Post-hoc tests for significant host shift types
    log_message("Performing post-hoc tests for host_shift_type...")
    emm <- emmeans(lmm_model, ~ host_shift_type)
    emm_pairs <- pairs(emm, adjust = "tukey")
    sig_pairs <- summary(emm_pairs)$p.value < 0.05
    if (any(sig_pairs)) {
      log_message("Significant host shift type comparisons (p < 0.05):")
      print(summary(emm_pairs)[sig_pairs, ])
    } else {
      log_message("No significant host shift type comparisons (p < 0.05).")
    }
    
    # Outlier detection using standardized residuals
    log_message("Detecting outlier branches...")
    residuals <- resid(lmm_model)
    std_residuals <- residuals / sd(residuals)
    outlier_branches <- filtered_data[abs(std_residuals) > 2.5, ]
    
    if (nrow(outlier_branches) > 0) {
      log_message("Found", nrow(outlier_branches), "outlier branches (standardized residual > 2.5):")
      print(outlier_branches[, c("segment", "protein", "host_shift_type", "rate")])
    } else {
      log_message("No outlier branches detected.")
    }
    
    # Robust Mahalanobis distance for multivariate outliers (rate and predictors)
    log_message("Computing robust Mahalanobis distance for outliers...")
    numeric_data <- filtered_data[, "rate", drop = FALSE]
    robust_cov <- covRob(numeric_data, estim = "mcd")
    mah_dist <- mahalanobis(numeric_data, center = robust_cov$center, cov = robust_cov$cov)
    outlier_mah <- filtered_data[mah_dist > qchisq(0.975, df = ncol(numeric_data)), ]
    
    if (nrow(outlier_mah) > 0) {
      log_message("Found", nrow(outlier_mah), "outlier branches (Mahalanobis distance):")
      print(outlier_mah[, c("segment", "protein", "host_shift_type", "rate")])
    } else {
      log_message("No multivariate outliers detected.")
    }
    
    # Per-protein ridgeline plots and summaries
    log_message("\n===== PER-PROTEIN ANALYSIS =====")
    proteins <- unique(filtered_data$protein)
    log_message("Analyzing", length(proteins), "proteins:", paste(proteins, collapse = ", "))
    
    for (prot in proteins) {
      log_message("\n=== Protein:", prot, "===")
      
      # Subset data for this protein
      prot_data <- filtered_data %>% filter(protein == prot)
      
      # Summarize host shift types
      shift_summary <- prot_data %>%
        group_by(host_shift_type) %>%
        summarize(
          n_branches = n(),
          mean_rate = mean(rate),
          median_rate = median(rate),
          sd_rate = sd(rate)
        ) %>%
        arrange(desc(n_branches))
      
      log_message("Host shift types summary for", prot, ":")
      print(shift_summary)
      
      # Filter for top 5 host shift types
      top_shift_types <- shift_summary %>%
        slice_head(n = 5) %>%
        pull(host_shift_type)
      
      # Subset data for plot
      plot_data <- prot_data %>%
        filter(host_shift_type %in% top_shift_types) %>%
        bind_rows(
          prot_data %>%
            mutate(host_shift_type = "All Shifts")
        )
      
      # Create ridgeline plot
      p_ridgeline <- ggplot(plot_data, aes(x = rate, y = host_shift_type, fill = host_shift_type)) +
        geom_density_ridges(scale = 1, alpha = 0.6) +
        geom_vline(xintercept = median(prot_data$rate), linetype = "dashed", color = "black") +
        labs(title = paste("Branch Rate Distributions for", prot),
             x = "Branch Rate", y = "Host Shift Type") +
        theme_minimal() +
        theme(legend.position = "none")
      
      print(p_ridgeline)  # Display in RStudio plot window
    }
    
    # PTM site analysis for outliers and significant shifts
    log_message("\n===== PTM SITE ANALYSIS =====")
    log_message("Note: PTM analysis requires alignment files or MutSiteDeep output.")
    
    # Focus on chicken-to-duck and outlier branches
    for (prot in unique(filtered_data$protein)) {
      ptm_file <- file.path(base_dir, unique(filtered_data$segment[filtered_data$protein == prot]), 
                            paste0(prot, "_ptm_alignment.fasta"))
      
      if (file.exists(ptm_file)) {
        log_message("Found PTM alignment for", prot, ":", ptm_file)
        
        # Chicken-to-duck shifts
        prot_data <- filtered_data %>% 
          filter(protein == prot, host_shift_type == "Gallus gallus → Anas platyrhynchos")
        
        if (nrow(prot_data) > 0) {
          log_message("Summarizing PTM changes for", prot, "in Gallus gallus → Anas platyrhynchos:")
          log_message("Found", nrow(prot_data), "branches.")
          # Placeholder PTM analysis
          log_message("Parse", ptm_file, "for glycosylation/phosphorylation changes.")
        }
        
        # Outlier branches
        prot_outliers <- outlier_branches %>% filter(protein == prot)
        if (nrow(prot_outliers) > 0) {
          log_message("Summarizing PTM changes for", nrow(prot_outliers), "outlier branches in", prot, ":")
          log_message("Parse", ptm_file, "for PTM changes in branches with rates:", 
                      paste(prot_outliers$rate, collapse = ", "))
        }
      } else {
        log_message("PTM alignment not found for", prot, ". Expected:", ptm_file)
      }
    }
  }
  
} else {
  log_message("No host shifts detected. No analysis performed.")
}

# ==== Final Summary ====
results_df <- bind_rows(results_list)
if (nrow(results_df) > 0) {
  log_message("\n===== FINAL SUMMARY =====")
  log_message("Processed", nrow(results_df), "proteins across", length(unique(results_df$segment)), "segments")
  log_message("Total host shifts:", sum(results_df$n_shifts))
}

log_message("===== ANALYSIS COMPLETE =====")


