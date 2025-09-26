############################################################
# WARNINGS: SynGO FOLDER NAMES AND DATA
#
# DIA-NN results for each sample were manually run on the SynGO portal. 
# SynGO data can be found in MassIVE dataset MSV000099156, 
# ProteomeXchange dataset PXD068359, or at https://github.com/LarryThePharmacologist
# SynGO exports often create very long folder names 
# that may exceed Windows path length limits and break the script. 
# SynGO XLSX file names can be indistinguishable.  
# To avoid this, rename long SynGO output folders to shorter names 
# and update the metadata CSV accordingly or avoid 
# hierarchical/nested directories. 
############################################################
# Load libraries
library(readxl)
library(dplyr)
library(stringr)
library(tidyr)
library(lubridate)
library(circlize)
library(ComplexHeatmap)

# Set working directory
setwd("~/Scripps 2/GitHub Testing/For GitHub")

# Set paths
syngo_data_dir <- "SynGO data"
results_dir <- file.path(syngo_data_dir, "SynGO comparison results")
dir.create(results_dir, showWarnings = FALSE)

# Load metadata
metadata <- read.csv("Metadata-Naming-ephys-grouping-cap-and-protein-ids.csv", check.names = FALSE)
folder_names <- metadata$`SynGO reanalysis folder name`

# --- Preflight checks ---
stopifnot(file.exists("Metadata-Naming-ephys-grouping-cap-and-protein-ids.csv"))

# Confirm metadata columns exist
required_meta_cols <- c("SynGO reanalysis folder name", "Electrophysiology", "Names for Analysis")
missing_meta <- setdiff(required_meta_cols, names(metadata))
if (length(missing_meta) > 0) {
  stop("Metadata is missing required columns: ", paste(missing_meta, collapse = ", "),
       "\nCheck the CSV name and headers.")
}

# Check folders exist & files present
message("Checking SynGO folders under: ", normalizePath(syngo_data_dir, mustWork = FALSE))
missing_folders <- folder_names[!dir.exists(file.path(syngo_data_dir, folder_names))]
if (length(missing_folders) > 0) {
  warning("These folders listed in metadata don't exist under '", syngo_data_dir, "':\n- ",
          paste(missing_folders, collapse = "\n- "), 
          "\nFix syngo_data_dir or the 'SynGO reanalysis folder name' column.")
}

# Peek at one expected file
example_files <- file.path(syngo_data_dir, folder_names, "syngo_ontologies_with_annotations_matching_user_input.xlsx")
message("Found ", sum(file.exists(example_files)), " SynGO Excel files out of ", length(example_files), " listed.")

# Initialize data frames
compiled_CC <- data.frame()
compiled_BP <- data.frame()

# Loop over folders to compile CC and BP results
for (folder in folder_names) {
  file_path <- file.path(syngo_data_dir, folder, "syngo_ontologies_with_annotations_matching_user_input.xlsx")
  if (file.exists(file_path)) {
    syngo_data <- read_excel(file_path, sheet = "sheet1")
    syngo_data$Sample <- folder
    syngo_data$Electrophysiology <- metadata$Electrophysiology[metadata$`SynGO reanalysis folder name` == folder]
    
    syngo_data_filtered <- syngo_data %>% 
      filter(`GSEA 'gene cluster' FDR corrected p-value` < 0.05) %>%
      mutate(`-log10Q` = -log10(`GSEA 'gene cluster' FDR corrected p-value`))
    
    CC_data <- syngo_data_filtered %>% filter(`GO domain` == "CC")
    BP_data <- syngo_data_filtered %>% filter(`GO domain` == "BP")
    
    compiled_CC <- bind_rows(compiled_CC, CC_data)
    compiled_BP <- bind_rows(compiled_BP, BP_data)
  }
}

# Timestamp
timestamp <- format(Sys.Date(), "%m%d%y")

# Save compiled full tables
write.csv(compiled_CC, file.path(results_dir, paste0("Compiled_CC_SynGO_", timestamp, ".csv")), row.names = FALSE)
write.csv(compiled_BP, file.path(results_dir, paste0("Compiled_BP_SynGO_", timestamp, ".csv")), row.names = FALSE)

# Function to generate split heatmaps and export gene lists
generate_split_heatmaps <- function(compiled_df, go_type = "BP") {
  domain <- ifelse(go_type == "BP", "Biological Process", "Cellular Component")
  file_prefix <- ifelse(go_type == "BP", "BP", "CC")
  
  # Create wide matrix of -log10Q
  heatmap_data <- compiled_df %>%
    select(Sample, `GO term name`, `-log10Q`) %>%
    distinct() %>%
    pivot_wider(names_from = Sample, values_from = `-log10Q`, values_fill = 0)
  
  mat <- as.matrix(heatmap_data[, -1])
  rownames(mat) <- heatmap_data$`GO term name`
  
  # Map to simplified sample names
  simplified_names <- setNames(metadata$`Names for Analysis`, metadata$`SynGO reanalysis folder name`)
  colnames(mat) <- simplified_names[colnames(mat)]
  
  # Split into shared and partially shared
  present_counts <- rowSums(mat > 0)
  shared_mat <- mat[present_counts == ncol(mat), , drop = FALSE]
  partial_mat <- mat[present_counts < ncol(mat), , drop = FALSE]
  
  # Export gene lists per GO term; hierarchichal ok
  gene_info <- compiled_df %>%
    select(`GO term name`, `genes - your genelist input`, Sample) %>%
    distinct()
  
  write.csv(gene_info, file.path(results_dir, paste0(file_prefix, "_GO_term_gene_mappings_", timestamp, ".csv")), row.names = FALSE)
  
  # Set heatmap parameters
  col_fun <- colorRamp2(c(0, 2, 5, 9), c("white", "lightblue", "purple", "red"))
  cell_w_mm <- 8
  cell_h_mm <- 5
  
  # Create annotation data frame for current heatmap columns
  annotation_df <- metadata %>%
    select(`SynGO reanalysis folder name`, `Names for Analysis`, Electrophysiology) %>%
    rename(SampleFolder = `SynGO reanalysis folder name`,
           Sample = `Names for Analysis`) %>%
    filter(Sample %in% colnames(partial_mat)) %>%
    arrange(match(Sample, colnames(partial_mat)))
  
  # Define color mapping
  ephys_colors <- c(
    "Preserved" = "green",
    "None" = "red",
    "Lost" = "orange",
    "Neuron Torn" = "gray"
  )
  
  # Create column annotation
  ha_col <- HeatmapAnnotation(
    Electrophysiology = annotation_df$Electrophysiology,
    col = list(Electrophysiology = ephys_colors),
    annotation_legend_param = list(title = "Ephys Status"),
    gp = gpar(col = "black", lwd = 1)
  )
  
  plot_heatmap <- function(mat, title, suffix) {
    png(file.path(results_dir, paste0(file_prefix, "_", suffix, "_Heatmap_", timestamp, ".png")), 
        width = 700 + ncol(mat) * cell_w_mm, 
        height = 500 + nrow(mat) * cell_h_mm, 
        units = "mm", res = 300)
    
    draw(
      Heatmap(mat,
              name = "-log10(Q)",
              col = col_fun,
              cluster_rows = FALSE,
              cluster_columns = TRUE,
              show_column_names = TRUE,
              column_names_side = "top",   # <-- move sample names to top
              show_row_names = TRUE,
              rect_gp = gpar(col = "black", lwd = 1),
              row_names_gp = gpar(fontsize = 14),
              column_names_gp = gpar(fontsize = 16),
              column_names_rot = 45,        # keep horizontal for clean readout
              column_title = title,
              top_annotation = ha_col,
              heatmap_legend_param = list(
                title = "-log10(Q-value)",
                at = c(0, 2, 5, 9),
                labels_gp = gpar(fontsize = 10),
                title_gp = gpar(fontsize = 12),
                legend_width = unit(5, "cm")
              ),
              width = unit(ncol(mat) * cell_w_mm, "mm"),
              height = unit(nrow(mat) * cell_h_mm + 40, "mm")  # add some top buffer
      ),
      heatmap_legend_side = "left",
      annotation_legend_side = "left"
    )
    dev.off()
  }
  # ----  Limit to top 10 GO terms by average -log10(Q) across samples ----
  top_n <- 10  #  Can change this to your desired number of GO terms
  
  if (nrow(shared_mat) > top_n) {
    shared_avg <- rowMeans(shared_mat, na.rm = TRUE)
    shared_mat <- shared_mat[order(-shared_avg)[1:top_n], , drop = FALSE]
  }
  
  if (nrow(partial_mat) > top_n) {
    partial_avg <- rowMeans(partial_mat, na.rm = TRUE)
    partial_mat <- partial_mat[order(-partial_avg)[1:top_n], , drop = FALSE]
  }
  # ---------------------------------------------------------------------------------
  if (nrow(shared_mat) > 0) {
    plot_heatmap(shared_mat, paste("Shared", domain, "GO Terms"), "Shared")
  }
  if (nrow(partial_mat) > 0) {
    plot_heatmap(partial_mat, paste("Unique/Partial", domain, "GO Terms"), "Partial")
  }
}

# Run on both domains
generate_split_heatmaps(compiled_BP, go_type = "BP")
generate_split_heatmaps(compiled_CC, go_type = "CC")

cat("Split GO heatmaps and gene lists saved.\n")
