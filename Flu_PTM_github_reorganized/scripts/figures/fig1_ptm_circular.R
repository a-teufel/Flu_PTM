########################################################################
#  LIBRARIES
########################################################################
library(dplyr)
library(tidyr)
library(ggplot2)
library(jsonlite)
library(readr)
library(stringr)
library(purrr)
library(geomtextpath)
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
  X= "white",   # Gray DDDDDD
  Y="#8C564B",   # Brown
  U="#4DFFC3",   # Bright Turquoise
  H="#FF6B6B",   # Soft Red
  L="#6A5ACD",    # Slate Blue
  "-"="darkgray"  
)
STRAIN_ORDER <- c("H7N9","H1N1","H5N1")       # inner → outer
RING_W       <- 0.6                           # Compact ring thickness
PROTEIN_GAP  <- 2.0                            # Gap between protein segments
DATA_RADIUS <- 2.6                           # Inner radius for data rings
PROTEIN_RADIUS <- 2.0                        # Radius for protein label ring
SEGMENT_RADIUS <- 1.2                        # Radius for segment ring

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

map_strains_to_csv <- function(cm,json,csv){
  jmap <- fromJSON(json)
  h1 <- read_csv(csv$h1n1,show_col_types=FALSE)%>%mutate(StrainType="H1N1")
  h5 <- read_csv(csv$h5n1,show_col_types=FALSE)%>%mutate(StrainType="H5N1")
  h7 <- read_csv(csv$h7n9,show_col_types=FALSE)%>%mutate(StrainType="H7N9")
  all <- bind_rows(h1,h5,h7)
  
  cm$Final_Accession <- cm$Matched_Strain <- NA_character_
  
  ## step 1 – direct accession
  acc_pat <- "^[A-Z]{2}[_]?[0-9]{6,}$"
  hit <- which(grepl(acc_pat,cm$Strain)&cm$Strain%in%all$Accession)
  cm$Final_Accession[hit] <- cm$Strain[hit]
  cm$Matched_Strain[hit]  <- all$StrainType[match(cm$Strain[hit],all$Accession)]
  
  ## step 2 – HxNy__ACC_
  rem <- which(is.na(cm$Final_Accession)&
                 grepl("^H\\d+N\\d+__([A-Z]{2,}\\d{5,})_$",cm$Strain))
  for(i in rem){
    acc <- sub("^H\\d+N\\d+__([A-Z]{2,}\\d{5,})_$","\\1",cm$Strain[i])
    idx <- match(acc,all$Accession)
    if(!is.na(idx)){
      cm$Final_Accession[i]<-acc
      cm$Matched_Strain[i] <-all$StrainType[idx]
    }
  }
  
  ## step 3 – ARG via JSON
  rem <- which(is.na(cm$Final_Accession)&grepl("^A[0-9]{8,}$",cm$Strain))
  for(i in rem){
    key <- names(jmap)[jmap==cm$Strain[i]]
    if(length(key)&&grepl("__([A-Z]{2,}\\d{5,})_\\d+_",key)){
      acc <- sub(".*__([A-Z]{2,}\\d{5,})_\\d+_.*","\\1",key)
      idx <- match(acc,all$Accession)
      if(!is.na(idx)){
        cm$Final_Accession[i]<-acc
        cm$Matched_Strain[i] <-all$StrainType[idx]
      }
    }
  }
  
  join_meta <- function(d,m) if(!nrow(d)) d else left_join(d,m,
                                                           by=c("Final_Accession"="Accession"))
  list(H1N1=join_meta(filter(cm,Matched_Strain=="H1N1"),h1),
       H5N1=join_meta(filter(cm,Matched_Strain=="H5N1"),h5),
       H7N9=join_meta(filter(cm,Matched_Strain=="H7N9"),h7))
}

symbol_freq <- function(df,strain){
  if(!nrow(df)) return(NULL)
  site_cols <- names(dplyr::select(df, starts_with("Site_")))
  df %>% pivot_longer(cols = starts_with("Site_"),
                      names_to="Position", values_to="Symbol") %>%
    filter(Symbol %in% names(SYMBOL_COLOURS)) %>%
    mutate(Strain=strain,
           Position=factor(Position,levels=site_cols)) %>%
    count(Position,Strain,Symbol,name="n") %>%
    group_by(Position,Strain) %>% mutate(Frac=n/sum(n)) %>% ungroup()
}

perform_tests <- function(h1,h5,h7){
  to_long <- function(d,s) d %>%
    pivot_longer(cols = starts_with("Site_"),
                 names_to="Position", values_to="Symbol") %>%
    filter(Symbol %in% names(SYMBOL_COLOURS)) %>% mutate(Strain=s)
  all <- bind_rows(to_long(h1,"H1N1"),to_long(h5,"H5N1"),to_long(h7,"H7N9"))
  map_dfr(unique(all$Position), function(pos){
    sub <- filter(all,Position==pos)
    if(length(unique(sub$Symbol))<=1) return(NULL)
    tab <- table(sub$Strain,sub$Symbol)
    if(nrow(tab)<2||ncol(tab)<2) return(NULL)
    p <- if(any(tab<5)&&nrow(tab)<=3&&ncol(tab)<=3)
      fisher.test(tab,workspace=5e7)$p.value
    else chisq.test(tab)$p.value
    tibble(Site=pos,Pvalue=p)
  }) %>% mutate(AdjP=p.adjust(Pvalue,"bonferroni"),
                Sig=AdjP<0.05)
}

analyse_prot <- function(nex,json,csv,tag){
  mp <- map_strains_to_csv(parse_nexus_file(nex),json,csv)
  list(
    id    = tag,
    freq  = bind_rows(symbol_freq(mp$H1N1,"H1N1"),
                      symbol_freq(mp$H5N1,"H5N1"),
                      symbol_freq(mp$H7N9,"H7N9")),
    stats = perform_tests(mp$H1N1,mp$H5N1,mp$H7N9)
  )
}

########################################################################
#  PART 2 – batch across seg1…seg8
########################################################################
# Define file paths for your CSV files
csv_files <- list(
  h1n1=here("data/raw_sequences/H1N1/H1N1_seq_info.csv"),
  h5n1=here("data/raw_sequences/H5N1/H5N1_seq_info.csv"),
  h7n9=here("data/raw_sequences/H7N9/H7N9_seq_info.csv")
)

# Find all segment directories in the paper-era processed folder
processed_base <- file.path(normalizePath(file.path(here::here(), ".."), winslash = "/", mustWork = FALSE),
                            "cluster_output_final_converged", "processed_segments_fixed_final")
seg_dirs <- list.dirs(processed_base, recursive = FALSE, full.names = TRUE) %>%
  .[grepl("seg[1-8]$", basename(.))]
if(!length(seg_dirs)) stop("No seg1–seg8 folders located.")

# Run analysis for each protein
results <- list()
for(seg in seg_dirs){
  for(nx in list.files(seg,"_multistate\\.txt$",full.names=TRUE)){
    js <- file.path(dirname(nx),
                    sub("_multistate\\.txt$","_id_mapping.json",basename(nx)))
    if(!file.exists(js)){ warning("JSON missing for ",basename(nx)); next }
    tag <- paste(basename(seg),
                 sub("_multistate\\.txt$","",basename(nx)), sep=" / ")
    message("Analysing ",tag)
    results[[tag]] <- analyse_prot(nx,js,csv_files,tag)
  }
}
if(!length(results)) stop("No protein analyses completed.")


# ── drop unwanted proteins from the results list ──
exclude_proteins <- c("HA1_domain","HA2_domain")
results <- results[
  !sapply(results, function(x)
    any(grepl(paste(exclude_proteins, collapse="|"), x$id))
  )
]


########################################################################
#  PART 3 – table of significant sites
########################################################################
sig_df <- bind_rows(lapply(results, function(x){
  if(!nrow(x$stats)) return(NULL)
  mutate(x$stats, Protein=x$id,
         SiteNum=as.integer(sub("Site_","",Site))) %>% filter(Sig)
}))
write_csv(sig_df,"significant_sites_table.csv")
message("✓ significant_sites_table.csv written (", nrow(sig_df), " rows)")

dir.create(here("results/summary_figures"), showWarnings = FALSE, recursive = TRUE)
write_csv(sig_df, here("results/summary_figures/significant_sites_table_fig1.csv"))

########################################################################
#  PART 4 – build circular stacked plot with protein‑ring labels
########################################################################
# ---------- composition data -------------
comp <- bind_rows(lapply(results, `[[`, "freq"), .id="Protein")

# Ensure all sites have equal representation
comp_full <- comp %>%
  complete(Protein, Position, Strain, Symbol = names(SYMBOL_COLOURS), 
           fill = list(n = 0, Frac = 0))

# Extract protein short names from the full protein names
comp_full <- comp_full %>%
  mutate(ShortName = sub(".*/\\s*", "", Protein))

# Extract segment numbers from protein names
comp_full <- comp_full %>%
  mutate(SegmentNumber = as.integer(sub("^seg(\\d+).*", "\\1", 
                                        sub("^(seg\\d+).*", "\\1", Protein))))

# Add PosNum column first
comp_full <- comp_full %>%
  mutate(PosNum = as.integer(sub("Site_", "", Position)))

# Count actual data points per protein to determine proper widths
protein_data_counts <- comp_full %>%
  filter(n > 0) %>%  # Only count actual data points, not empty positions
  group_by(ShortName) %>%
  summarize(data_count = n_distinct(Position)) %>%
  ungroup()

# Modified block geometry with spacing based on actual data count
prot_order <- unique(comp_full$ShortName)
blocks <- map_dfr(prot_order, function(p){
  # Get data count for this protein
  data_count <- protein_data_counts %>% 
    filter(ShortName == p) %>% 
    pull(data_count)
  
  # Default to 1 if no data found (avoid 0-width segments)
  if(length(data_count) == 0 || is.na(data_count) || data_count == 0) data_count <- 1
  
  tibble(
    Protein = p,
    n_sites = data_count,  # Use actual data count instead of all possible positions
    SegmentNumber = unique(filter(comp_full, ShortName == p)$SegmentNumber)[1]
  )
}) %>% mutate(
  start = cumsum(lag(n_sites, default = 0) + lag(rep(PROTEIN_GAP, n()), default = 0)) + 1,
  end = start + n_sites - 1,
  # Calculate midpoint for each protein segment
  mid = (start + end) / 2
)

# Calculate total plot width for proper circular layout
total_sites <- max(blocks$end) + PROTEIN_GAP

# Create a mapping from original Position to new Global position
position_mapping <- comp_full %>%
  dplyr::select(ShortName, Position, PosNum) %>%
  distinct() %>%
  group_by(ShortName) %>%
  arrange(PosNum) %>%
  # Create a sequential number within each protein
  mutate(LocalPos = row_number()) %>%
  # Join with blocks to get start position
  left_join(blocks, by = c("ShortName" = "Protein")) %>%
  # Calculate global position
  mutate(
    Global = start + LocalPos - 1
  ) %>%
  dplyr::select(ShortName, Position, Global)

# merge coords with new position mapping
comp_full <- comp_full %>% 
  # Use the position mapping instead of original blocks join
  left_join(position_mapping, by = c("ShortName", "Position")) %>%
  mutate(
    Ring = match(Strain, STRAIN_ORDER)
  ) %>%
  # Filter out positions with no data
  filter(!is.na(Global), n > 0)  # Only keep positions that have data

# stacked bar ymin/ymax - adjusted for new DATA_RADIUS
comp_full <- comp_full %>% 
  group_by(ShortName, Position, Strain) %>%
  arrange(Symbol) %>%
  mutate(
    ymin = DATA_RADIUS + (Ring-1)*RING_W + cumsum(lag(Frac, default=0)) * RING_W,
    ymax = DATA_RADIUS + (Ring-1)*RING_W + cumsum(    Frac         ) * RING_W
  ) %>% 
  ungroup() %>% 
  mutate(
    xmin = Global - 0.45,
    xmax = Global + 0.45
  )


# significance star coords - adjusted for new DATA_RADIUS and position mapping
sig_plot <- sig_df %>%
  # Extract short protein name and create Position column
  mutate(
    ShortName = sub(".*/\\s*", "", Protein),
    Position = paste0("Site_", SiteNum)  # Create Position column to match position_mapping
  ) %>%
  # Use position mapping with correct join columns
  left_join(position_mapping, by = c("ShortName", "Position")) %>%
  mutate(
    y = DATA_RADIUS + length(STRAIN_ORDER)*RING_W + 0.2
  ) %>%
  filter(!is.na(Global))  # Only keep positions that have a mapping

# Add dividers between protein segments - adjusted for new DATA_RADIUS
protein_dividers <- blocks %>%
  dplyr::select(Protein, start, end) %>%
  mutate(
    x_div_start = start - 0.5,
    x_div_end = end + 0.5
  ) %>%
  tidyr::pivot_longer(
    cols = c(x_div_start, x_div_end),
    names_to = "type",
    values_to = "x"
  ) %>%
  mutate(
    y_min = DATA_RADIUS,
    y_max = DATA_RADIUS + length(STRAIN_ORDER) * RING_W
  )

# Create segment label data in ring format
segment_labels <- blocks %>%
  group_by(SegmentNumber) %>%
  summarize(
    start_pos = min(start),
    end_pos = max(end),
    mid_pos = mean(c(min(start), max(end)))
  ) %>%
  mutate(
    # Calculate angles for the segment middle
    mid_angle = (mid_pos / total_sites) * 2 * pi,
    
    # Position text in the ring
    xmin = start_pos - 0.5,
    xmax = end_pos + 0.5,
    
    # Create the segment label ring - make thicker
    ymin = SEGMENT_RADIUS - 0.4,
    ymax = SEGMENT_RADIUS + 0.4,
    
    # Create label text
    label = paste("Segment", SegmentNumber),
    
    # Calculate coordinates for label text
    text_x = mid_pos,
    text_y = SEGMENT_RADIUS,
    
    # Calculate text angle for readability
    text_angle = (mid_angle * 180/pi) %% 360,
    text_angle = ifelse(text_angle > 90 & text_angle < 270, 
                        text_angle + 180, text_angle),
    
    # Adjust text alignment based on position
    hjust = 0.5
  )

# Black background bars for segment labels
segment_bg <- segment_labels %>%
  dplyr::select(xmin, xmax, ymin, ymax, SegmentNumber)

# Create protein label data in ring format
protein_label_bg <- blocks %>%
  mutate(
    # Position text in the ring
    xmin = start - 0.5,
    xmax = end + 0.5,
    
    # Create the protein label ring - make thicker
    ymin = PROTEIN_RADIUS - 0.4,
    ymax = PROTEIN_RADIUS + 0.4,
    
    # Calculate middle position for text
    text_x = mid,
    text_y = PROTEIN_RADIUS,
    
    # Calculate text angle based on position
    mid_angle = (mid / total_sites) * 2 * pi,
    text_angle = (mid_angle * 180/pi) %% 360,
    text_angle = ifelse(text_angle > 90 & text_angle < 270, 
                        text_angle + 180, text_angle),
    
    # Adjust text alignment
    hjust = 0.5
  )

# Strain labels - positioned within each ring - adjusted for new DATA_RADIUS
strain_labels <- tibble(
  Strain = STRAIN_ORDER,
  # Position in the middle of each ring
  y = DATA_RADIUS + seq(RING_W/2, by = RING_W, length.out = length(STRAIN_ORDER)),
  # Position at a fixed angle (right side of circle)
  angle_rad = 0,
  x = y * cos(angle_rad),
  # Use rectangular labels for better visibility
  background = TRUE
)

# Fix the crossing function issue by using expand_grid instead
segment_strain_labels <- segment_labels %>%
  dplyr::select(text_x, text_angle, hjust) %>%      # one row per segment
  expand_grid(Strain = STRAIN_ORDER) %>%  # expand to 3× as many rows
  mutate(
    ring = match(Strain, STRAIN_ORDER),
    # place labels just outside each ring
    y    = DATA_RADIUS + ring * RING_W + 0.1,
    x    = text_x,
    angle = text_angle
  )

label_x <- blocks$start - (PROTEIN_GAP/2)

# ─── build curved label paths - FIX FOR CROSSING ERROR ───
protein_strain_paths <- blocks %>%
  mutate(
    x_start = start - PROTEIN_GAP,
    x_end   = start
  ) %>%
  dplyr::select(Protein, x_start, x_end) %>%
  expand_grid(Strain = STRAIN_ORDER) %>%  # Use expand_grid instead of crossing
  mutate(
    ring  = match(Strain, STRAIN_ORDER),
    y     = DATA_RADIUS + ring * RING_W - 0.1,
    group = paste(Protein, Strain, sep = "_")
  )

# Now create the expanded paths using do() instead of reframe
protein_strain_paths <- protein_strain_paths %>%
  group_by(Protein, Strain, x_start, x_end, y, group) %>%
  do({
    tibble(
      x     = seq(.$x_start[1], .$x_end[1], length.out = 50),
      y     = .$y[1],
      label = .$Strain[1],
      group = .$group[1]
    )
  }) %>%
  ungroup()

# ——— define background colors for rings ———
ring_colors <- tibble(
  Ring   = seq_along(STRAIN_ORDER),
  Strain = STRAIN_ORDER,
  # three pastel hues, none of which match your SYMBOL_COLOURS
  color  = c("#A8E6CF",  # warm peach
             "#66B2FF",  # vibrant sky‑blue
             "#B266FF"), # medium lavender
  ymin   = DATA_RADIUS + (Ring - 1) * RING_W,
  ymax   = DATA_RADIUS +  Ring    * RING_W,
  xmin   = 0,
  xmax   = total_sites
)

########
# ——— build the plot ———
p2 <- ggplot() +
  # (1) segment labels background + text
  geom_rect(data = segment_bg,
            aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
            fill = "black", size = 0.5, inherit.aes = FALSE) +
  
  geom_text(data = segment_labels,
            aes(x = text_x, y = text_y, label = label),
            color = "white", fontface = "bold",
            size = 3, hjust = 0.5, inherit.aes = FALSE) +
  
  # (2) protein label bands + text
  geom_rect(data = protein_label_bg,
            aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
            fill = "#003366", inherit.aes = FALSE) +
  
  geom_textpath(data = protein_label_bg,
                aes(x = text_x, y = text_y, label = Protein),
                color = "white", fontface = "bold",
                size = 2, hjust = .5, angle=90) +
  
  # (3) three strain‐ring backgrounds in distinct greys
  geom_rect(data = ring_colors,
            aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
            fill = ring_colors$color,
            alpha = 0.9, inherit.aes = FALSE) +
  
  # (4) segment dividers
  geom_segment(data = protein_dividers,
               aes(x = x, xend = x, y = y_min, yend = y_max),
               color = "black", size = 1, alpha = 0.8,
               inherit.aes = FALSE) +
  
  # (5) your stacked PTM bars
  geom_rect(data = comp_full,
            aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = Symbol),
            colour = "transparent", inherit.aes = FALSE, size=0.8) +
  scale_fill_manual(values = SYMBOL_COLOURS, name = "Symbol") +
  
  # (6) significance stars
  geom_point(data = sig_plot,
             aes(x = Global, y = y),
             colour = "red", shape = 20, size = 2.5,
             inherit.aes = FALSE) +
  
  # (7) strain labels
  geom_label(data = strain_labels,
             aes(x = 0, y = y, label = Strain),
             hjust = 0, vjust=-0.1, size = 2, fontface = "bold",
             fill = "white", alpha = 0.8, label.size = 0.5,
             label.padding = unit(0.15, "lines"), angle = 90,
             inherit.aes = FALSE) +
  
  # (8) ring boundaries
  geom_hline(yintercept = DATA_RADIUS + (0:length(STRAIN_ORDER)) * RING_W-.01,
             color = "black", size = 1, alpha = 0.8,
             inherit.aes = FALSE) +
  
  # polar coords + clean up
  coord_curvedpolar(theta = "x", start = 0, direction = 1, clip = "off") +
  scale_x_continuous(limits = c(0, total_sites), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0, DATA_RADIUS + length(STRAIN_ORDER) * RING_W + 0.8),
                     expand = c(0, 0)) +
  theme_void() +
  theme(
    legend.position = c(0.9, 0.5),  # Position legend inside the plot, closer to circle
    plot.title      = element_text(hjust = 0.5, size = 14, face = "bold"),
    legend.title    = element_text(size = 10),
    legend.text     = element_text(size = 8),
    plot.margin     = margin(20, 20, 20, 20)
    #plot.margin = margin(-50, -20, -50, -20)  # Extremely negative margins to force crop
    
  ) 

# Save the figure
#ggsave(here("results/summary_figures/fig1_ptm_circular.pdf"), p2, width = 18, height = 16, dpi = 300)


