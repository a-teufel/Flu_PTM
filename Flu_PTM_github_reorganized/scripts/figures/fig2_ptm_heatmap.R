########################################################################
#  LIBRARIES
########################################################################
rm(list = ls())

library(dplyr)
library(tidyr)
library(ggplot2)
library(jsonlite)
library(readr)
library(stringr)
library(purrr)
library(ggalluvial)
library(forcats)  # For factor manipulation
library(gridExtra) # For arranging multiple plots
library(viridis)   # For color palettes
library(scales)    # For formatting
library(DescTools) # For CramerV calculation
library(here)

########################################################################
#  CONSTANTS
########################################################################
SYMBOL_COLOURS <- c(
  A="#E69F00",   # Orange
  C="#1E90FF",   # Dodger Blue
  K="#56B4E9",   # Sky Blue
  M="#009E73",   # Emerald Green
  N="#F0E442",   # Butter Yellow
  O="#0072B2",   # Intense Blue
  Q="#D55E00",   # Vermilion
  R="#CC79A7",   # Rosacea
  S="#984EA3",   # Black
  T="#82B366",   # Muted Green
  Y="#8C564B",   # Brown
  U="#4DFFC3",   # Bright Turquoise
  H="#FF6B6B",   # Soft Red
  L="#6A5ACD",    # Slate Blue
  "-"="darkgray"  
)
STRAIN_ORDER <- c("H7N9","H1N1","H5N1")

########################################################################
#  PART 1 – helper functions
########################################################################
parse_nexus_file <- function(file){
  ln <- readLines(file)
  start <- grep("^MATRIX",ln)
  end   <- grep("^;",ln[(start+1):length(ln)])[1]+start
  mat   <- ln[(start+1):(end-1)];
  mat   <- mat[mat!=""]
  bind_rows(lapply(mat, function(line){
    pr <- strsplit(line,"\\s+")[[1]]; pr <- pr[pr!=""]
    if(length(pr)<=1) return(NULL)
    df <- data.frame(Strain=pr[1], stringsAsFactors=FALSE)
    for(i in seq_along(pr[-1])) df[[paste0("Site_",i)]] <- pr[-1][i]
    df
  }))
}

# Simplified function to map strains to CSV with no data manipulation
map_strains_to_csv <- function(cm, json, csv){
  jmap <- fromJSON(json)
  
  # Read CSV files - consistent column handling
  h1 <- read_csv(csv$h1n1, show_col_types=FALSE) %>% 
    mutate(StrainType="H1N1") 
  
  h5 <- read_csv(csv$h5n1, show_col_types=FALSE) %>% 
    mutate(StrainType="H5N1")
  
  h7 <- read_csv(csv$h7n9, show_col_types=FALSE) %>% 
    mutate(StrainType="H7N9")
  
  all <- bind_rows(h1, h5, h7)
  
  cm$Final_Accession <- cm$Matched_Strain <- NA_character_
  
  ## step 1 – direct accession
  acc_pat <- "^[A-Z]{2}[_]?[0-9]{6,}$"
  hit <- which(grepl(acc_pat, cm$Strain) & cm$Strain %in% all$Accession)
  cm$Final_Accession[hit] <- cm$Strain[hit]
  cm$Matched_Strain[hit]  <- all$StrainType[match(cm$Strain[hit], all$Accession)]
  
  ## step 2 – HxNy__ACC_
  rem <- which(is.na(cm$Final_Accession) &
                 grepl("^H\\d+N\\d+__([A-Z]{2,}\\d{5,})_$", cm$Strain))
  for(i in rem){
    acc <- sub("^H\\d+N\\d+__([A-Z]{2,}\\d{5,})_$", "\\1", cm$Strain[i])
    idx <- match(acc, all$Accession)
    if(!is.na(idx)){
      cm$Final_Accession[i] <- acc
      cm$Matched_Strain[i]  <- all$StrainType[idx]
    }
  }
  
  ## step 3 – ARG via JSON
  rem <- which(is.na(cm$Final_Accession) & grepl("^A[0-9]{8,}$", cm$Strain))
  for(i in rem){
    key <- names(jmap)[jmap == cm$Strain[i]]
    if(length(key) && grepl("__([A-Z]{2,}\\d{5,})_\\d+_", key)){
      acc <- sub(".*__([A-Z]{2,}\\d{5,})_\\d+_.*", "\\1", key)
      idx <- match(acc, all$Accession)
      if(!is.na(idx)){
        cm$Final_Accession[i] <- acc
        cm$Matched_Strain[i]  <- all$StrainType[idx]
      }
    }
  }
  
  # Join metadata
  join_meta <- function(d, m) {
    if(!nrow(d)) return(d) 
    left_join(d, m, by=c("Final_Accession"="Accession"))
  }
  
  list(H1N1=join_meta(filter(cm, Matched_Strain=="H1N1"), h1),
       H5N1=join_meta(filter(cm, Matched_Strain=="H5N1"), h5),
       H7N9=join_meta(filter(cm, Matched_Strain=="H7N9"), h7))
}

# Consistent handling of Year as character type
process_data <- function(df, strain){
  if(!nrow(df)) return(NULL)
  
  site_cols <- names(dplyr::select(df, starts_with("Site_")))
  
  df %>% 
    pivot_longer(cols = starts_with("Site_"),
                 names_to="Position", values_to="Symbol") %>%
    filter(Symbol %in% names(SYMBOL_COLOURS), Symbol != "X") %>%
    mutate(
      Strain = strain,
      Position = factor(Position, levels=site_cols),
      # Always ensure Year is a character type for consistency
      Year = as.character(ifelse(!is.na(Collection_Date), 
                                 str_extract(Collection_Date, "\\d{4}"), 
                                 "Unknown"))
    )
}

analyse_prot <- function(nex, json, csv, tag){
  mp <- map_strains_to_csv(parse_nexus_file(nex), json, csv)
  
  # Process each strain dataset separately
  h1_proc <- process_data(mp$H1N1, "H1N1")
  h5_proc <- process_data(mp$H5N1, "H5N1")
  h7_proc <- process_data(mp$H7N9, "H7N9")
  
  # Combine results
  all_freq <- bind_rows(h1_proc, h5_proc, h7_proc)
  
  # Add protein info
  all_freq$Protein <- sub(".*/\\s*", "", tag)
  
  return(all_freq)
}

########################################################################
#  PART 2 – Process data for all proteins
########################################################################
# Define file paths for your CSV files
csv_files <- list(
  h1n1=here("data/raw_sequences/H1N1/H1N1_seq_info.csv"),
  h5n1=here("data/raw_sequences/H5N1/H5N1_seq_info.csv"),
  h7n9=here("data/raw_sequences/H7N9/H7N9_seq_info.csv")
)

# Find all segment directories
processed_base <- file.path(normalizePath(file.path(here::here(), ".."), winslash = "/", mustWork = FALSE),
                            "cluster_output_final_converged", "processed_segments_fixed_final")
seg_dirs <- list.dirs(processed_base, recursive = FALSE, full.names = TRUE) %>%
  .[grepl("seg[1-8]$", basename(.))]
if(!length(seg_dirs)) stop("No seg1–seg8 folders located.")

dir.create(here("results/summary_figures"), showWarnings = FALSE, recursive = TRUE)

# Run analysis for each protein
all_results <- list()
for(seg in seg_dirs){
  for(nx in list.files(seg, "_multistate\\.txt$", full.names=TRUE)){
    js <- file.path(dirname(nx),
                    sub("_multistate\\.txt$", "_id_mapping.json", basename(nx)))
    if(!file.exists(js)){ warning("JSON missing for ", basename(nx)); next }
    tag <- paste(basename(seg),
                 sub("_multistate\\.txt$", "", basename(nx)), sep=" / ")
    message("Analysing ", tag)
    all_results[[tag]] <- analyse_prot(nx, js, csv_files, tag)
  }
}

# Check if we have results
if(length(all_results) == 0) stop("No protein analyses completed.")

# ── drop unwanted proteins ──
exclude_proteins <- c("HA1_domain", "HA2_domain")
all_results <- all_results[
  !sapply(all_results, function(x)
    any(grepl(paste(exclude_proteins, collapse="|"), names(x)))
  )
]

# Combine all data
all_data <- bind_rows(all_results)

# Create a time period category for more meaningful analysis
all_data <- all_data %>%
  mutate(
    TimePeriod = case_when(
      Year == "Unknown" ~ "Unknown",
      Year %in% c("2000", "2001", "2002", "2003", "2004") ~ "2000-2004",
      Year %in% c("2005", "2006", "2007", "2008", "2009") ~ "2005-2009",
      Year %in% c("2010", "2011", "2012", "2013", "2014") ~ "2010-2014",
      Year %in% c("2015", "2016", "2017", "2018", "2019") ~ "2015-2019",
      Year %in% c("2020", "2021", "2022", "2023", "2024") ~ "2020-2024",
      TRUE ~ "Other"
    )
  )

# Extract country from Geo_Location if needed
all_data <- all_data %>%
  mutate(Country = ifelse(is.na(Country) | Country == "", 
                          sub("^([^/]+).*", "\\1", Geo_Location), 
                          Country))
# Load required packages

# Load necessary libraries
library(dplyr)
library(tibble)
# Step 1: Create a contingency table of Symbol counts by Host
symbol_host_table <- all_data %>%
  filter(!is.na(Host), !is.na(Symbol)) %>%
  count(Host, Symbol) %>%
  tidyr::pivot_wider(names_from = Symbol, values_from = n, values_fill = 0) %>%
  column_to_rownames("Host") %>%
  as.matrix()

# Step 2: Perform the Chi-square test
chisq_result <- chisq.test(symbol_host_table)

# Step 3: View the result
print(chisq_result)

# Optional: Examine which combinations contribute most to the Chi-square statistic
chisq_contrib <- chisq_result$stdres^2  # squared standardized residuals
print(chisq_contrib)




# Create contingency table of Symbol counts by Host
symbol_host_table <- all_data %>%
  filter(!is.na(Host), !is.na(Symbol)) %>%
  count(Host, Symbol) %>%
  tidyr::pivot_wider(names_from = Symbol, values_from = n, values_fill = 0) %>%
  column_to_rownames("Host") %>%
  as.matrix()

# 1. Chi-square test of independence
chisq_result <- chisq.test(symbol_host_table)
print(chisq_result)

# 2. Inspect expected counts
print("Expected counts:")
print(chisq_result$expected)

# 3. If needed, simulate p-value
if(any(chisq_result$expected < 5)) {
  chisq_sim <- chisq.test(symbol_host_table, simulate.p.value = TRUE, B = 10000)
  print("Chi-square with simulated p-value:")
  print(chisq_sim)
}

# 4. Effect size (Cramer's V)
cramer_v <- DescTools::CramerV(symbol_host_table)
message("Cramer's V = ", round(cramer_v, 3))

# 5. Contribution of each cell to the chi-square statistic
chisq_contrib <- chisq_result$stdres^2
print("Squared standardized residuals (contributions):")
print(chisq_contrib)









########################################################################
#  LIBRARIES
########################################################################
library(dplyr)
library(tidyr)
library(tibble)
library(viridis)
library(pheatmap)

########################################################################
#  EMOJI MAPPING
########################################################################
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
  "Suricata suricatta"   = "\U1F43F",
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
  "Parus major"          = "\U1F426"

)
default_symbol <- "\U003F"  # question‐mark for anything missing

########################################################################
#  STEP 1: build Host × Symbol table, dropping the "-" symbol
########################################################################
symbol_host_table <- all_data %>%
  filter(!is.na(Host), !is.na(Symbol), Symbol != "-") %>%
  count(Host, Symbol) %>%
  pivot_wider(
    names_from   = Symbol,
    values_from  = n,
    values_fill  = list(n = 0)
  ) %>%
  column_to_rownames("Host") %>%
  as.matrix()

########################################################################
#  STEP 2: attach emojis to host names
########################################################################
host_names  <- rownames(symbol_host_table)
host_emojis <- sapply(host_names, function(h) {
  if(h %in% names(species_symbols)) species_symbols[[h]] else default_symbol
})
host_labels <- paste0(host_names, "  ", host_emojis)
rownames(symbol_host_table) <- host_labels

########################################################################
#  STEP 3: normalize rows to proportions
########################################################################
mat_norm <- prop.table(symbol_host_table, margin = 1)

########################################################################
#  STEP 4: draw the clustered heatmap
########################################################################
pheatmap(
  mat_norm,
  cluster_rows    = TRUE,
  cluster_cols    = TRUE,
  show_rownames   = TRUE,
  show_colnames   = TRUE,
  fontsize_row    = 6,
  color           = magma(100),
  main            = "Symbol – Host heatmap (no “-” symbols, with emojis)"
)





########################################################################
#  LIBRARIES
########################################################################
library(dplyr)
library(tidyr)
library(tibble)
library(viridis)
library(pheatmap)

########################################################################
#  EMOJI MAPPING
########################################################################
species_symbols <- list(
  "Sus scrofa"           = "\U1F416",
  "Homo sapiens"         = "\U1F9D1",
  "Mus musculus"         = "\U1F42D",
  "Gallus gallus"        = "\U1F414",
  "Phasianinae"          = "\U1F99A",
  "Phasianidae"          = "\U1F99A",
  "Panthera tigris"      = "\U1F405",
  "Mustela lutreola"     = "\U1F98A",
  "Anatidae"             = "\U1F986",
  "Passer montanus"      = "\U1F426",
  "Anas platyrhynchos"   = "\U1F986",
  "Arenaria interpres"   = "\U1F426",
  "Accipitriformes"      = "\U1F985",
  "Aves"                 = "\U1F426",
  "Psittacidae"          = "\U1F99C",
  "Columbidae"           = "\U1F54A",
  "Suricata suricatta"   = "\U1F43F",
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
  "Parus major"          = "\U1F426"
  
)
default_symbol <- "\U003F"

########################################################################
#  STEP 1: build Host × Symbol table, dropping the "-" symbol
########################################################################
al <- all_data

symbol_host_table <- all_data %>%
  filter(!is.na(Host), !is.na(Symbol), Symbol != "-") %>%
  count(Host, Symbol) %>%
  pivot_wider(
    names_from   = Symbol,
    values_from  = n,
    values_fill  = list(n = 0)
  ) %>%
  column_to_rownames("Host") %>%
  as.matrix()

########################################################################
#  STEP 2: add emojis to host names
########################################################################
host_names  <- rownames(symbol_host_table)
host_emojis <- sapply(host_names, function(h) {
  if(h %in% names(species_symbols)) species_symbols[[h]] else default_symbol
})
host_labels <- paste0(host_names, "  ", host_emojis)
rownames(symbol_host_table) <- host_labels

########################################################################
#  STEP 3: normalize rows to proportions
########################################################################
mat_norm <- prop.table(symbol_host_table, margin = 1)

########################################################################
#  STEP 4: compute per-host Strain proportions for annotation
########################################################################
host_strain_prop <- all_data %>%
  filter(!is.na(Host), !is.na(Strain), Host %in% host_names) %>%
  count(Host, Strain) %>%
  group_by(Host) %>%
  mutate(prop = n / sum(n)) %>%
  dplyr::select(-n) %>%
  pivot_wider(names_from = Strain, values_from = prop, values_fill = list(prop = 0)) %>%
  column_to_rownames("Host")

# match the emoji‐augmented row names
rownames(host_strain_prop) <- host_labels[match(rownames(host_strain_prop), host_names)]


########################################################################
#  STEP 5: draw the heatmap with row‐annotation
########################################################################
pheatmap(
  mat_norm,
  annotation_row   = host_strain_prop,
  annotation_colors  = list(
    H7N9 = c("white", "#A8E6CF"),
    H1N1 = c("white", "#66B2FF"),
    H5N1 = c("white", "#B266FF")
  ),
  cluster_rows     = TRUE,
  cluster_cols     = TRUE,
  show_rownames    = TRUE,
  show_colnames    = TRUE,
  fontsize_row     = 6,
  color            = viridis(100),
  main             = "Symbol – Host heatmap (with Strain proportions)"
)




# 1. Convert to presence/absence factors
host_strain_pa <- host_strain_prop > 0
host_strain_pa <- as.data.frame(host_strain_pa)
host_strain_pa[] <- lapply(host_strain_pa, function(col) {
  factor(ifelse(col, "present", "absent"), levels = c("absent", "present"))
})

# 2. Define discrete colors for absent/present
ann_colors <- list(
  H7N9 = c(absent = "white", present = "#A8E6CF"),
  H1N1 = c(absent = "white", present = "#66B2FF"),
  H5N1 = c(absent = "white", present = "#B266FF")
)

# 3. Draw the heatmap with binary strain‐annotation
pheatmap(
  mat_norm,
  annotation_row     = host_strain_pa,
  annotation_colors  = ann_colors,
  cluster_rows       = TRUE,
  cluster_cols       = TRUE,
  show_rownames      = TRUE,
  show_colnames      = TRUE,
  fontsize_row       = 12,
  color              = c("white",plasma(100000)),
  filename          = here("results/summary_figures/fig2_ptm_heatmap_all_proteins.png"),
  width             = 8,                 # inches
  height            = 10                 # inches
)





al<-all_data[all_data$Protein=="NA_neuraminidase",]

symbol_host_table <- al %>%
  filter(!is.na(Host), !is.na(Symbol), Symbol != "-") %>%
  count(Host, Symbol) %>%
  pivot_wider(
    names_from   = Symbol,
    values_from  = n,
    values_fill  = list(n = 0)
  ) %>%
  column_to_rownames("Host") %>%
  as.matrix()

########################################################################
#  STEP 2: add emojis to host names
########################################################################
host_names  <- rownames(symbol_host_table)
host_emojis <- sapply(host_names, function(h) {
  if(h %in% names(species_symbols)) species_symbols[[h]] else default_symbol
})
host_labels <- paste0(host_names, "  ", host_emojis)
rownames(symbol_host_table) <- host_labels

########################################################################
#  STEP 3: normalize rows to proportions
########################################################################
mat_norm <- prop.table(symbol_host_table, margin = 1)

########################################################################
#  STEP 4: compute per-host Strain proportions for annotation
########################################################################
host_strain_prop <- al %>%
  filter(!is.na(Host), !is.na(Strain), Host %in% host_names) %>%
  count(Host, Strain) %>%
  group_by(Host) %>%
  mutate(prop = n / sum(n)) %>%
  dplyr::select(-n) %>%
  pivot_wider(names_from = Strain, values_from = prop, values_fill = list(prop = 0)) %>%
  column_to_rownames("Host")

# match the emoji‐augmented row names
rownames(host_strain_prop) <- host_labels[match(rownames(host_strain_prop), host_names)]


########################################################################
#  STEP 5: draw the heatmap with row‐annotation
########################################################################
pheatmap(
  mat_norm,
  annotation_row   = host_strain_prop,
  annotation_colors  = list(
    H7N9 = c("white", "#A8E6CF"),
    H1N1 = c("white", "#66B2FF"),
    H5N1 = c("white", "#B266FF")
  ),
  cluster_rows     = TRUE,
  cluster_cols     = TRUE,
  show_rownames    = TRUE,
  show_colnames    = TRUE,
  fontsize_row     = 6,
  color            = viridis(100),
  main             = "Symbol – Host heatmap (with Strain proportions)"
)




# 1. Convert to presence/absence factors
host_strain_pa <- host_strain_prop > 0
host_strain_pa <- as.data.frame(host_strain_pa)
host_strain_pa[] <- lapply(host_strain_pa, function(col) {
  factor(ifelse(col, "present", "absent"), levels = c("absent", "present"))
})

# 2. Define discrete colors for absent/present
ann_colors <- list(
  H7N9 = c(absent = "white", present = "#A8E6CF"),
  H1N1 = c(absent = "white", present = "#66B2FF"),
  H5N1 = c(absent = "white", present = "#B266FF")
)

# 3. Draw the heatmap with binary strain‐annotation
pheatmap(
  mat_norm,
  annotation_row     = host_strain_pa,
  annotation_colors  = ann_colors,
  cluster_rows       = TRUE,
  cluster_cols       = TRUE,
  show_rownames      = TRUE,
  show_colnames      = TRUE,
  fontsize_row       = 12,
  color              = c("white",plasma(100000)),
  filename          = here("results/summary_figures/fig2_ptm_heatmap_NA_neuraminidase.png"),
  width             = 8,                 # inches
  height            = 10     
)


