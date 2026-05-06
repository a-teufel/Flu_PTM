# install.packages("jsonlite")  # uncomment if needed
library(jsonlite)

# Path to your JSON file
json_file <- "NS1_protein_alignment_data.json"

# read the JSON
data     <- fromJSON(json_file)
mat_list <- data$alignment_matrix

# get sorted position indices
pos_indices <- sort(as.integer(names(mat_list)))

# prepare a vector to hold the consensus residues
consensus <- character(length(pos_indices))

for (i in seq_along(pos_indices)) {
  pos     <- as.character(pos_indices[i])
  entries <- mat_list[[pos]]
  
  # if there are no entries at this position, assign a gap (“-”)
  if (is.null(entries) || length(entries) == 0) {
    consensus[i] <- "-"
  } else {
    # extract residues whether entries is a data.frame/matrix or list
    if (is.data.frame(entries) || is.matrix(entries)) {
      residues <- entries[, 2]
    } else {
      residues <- vapply(entries, `[`, character(1), 2)
    }
    
    # tally frequencies and pick the most common residue
    freq_tab <- table(residues)
    max_n    <- max(freq_tab)
    tops     <- names(freq_tab)[freq_tab == max_n]
    
    # break ties alphabetically
    consensus[i] <- sort(tops)[1]
  }
}

consensus[183]
# collapse into one string and print
consensus_seq <- paste0(consensus, collapse = "")
cat("Consensus sequence (length", nchar(consensus_seq), "):\n")
cat(consensus_seq, "\n")


