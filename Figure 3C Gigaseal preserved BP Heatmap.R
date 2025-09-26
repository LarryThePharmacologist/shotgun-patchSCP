############################################################
# NOTE:
# This script depends on the compiled SynGO results generated 
# in the Figure 5 script. 
# Run Figure 5 script first to produce the required inputs.
############################################################
# Load necessary libraries
library(readr)
library(dplyr)
library(readxl)
library(tidyr)
library(circlize)
library(ComplexHeatmap)

# Set working directory
setwd("~/Scripps 2/GitHub Testing/For GitHub")
timestamp <- format(Sys.Date(), "%m%d%y")

# Define paths (modify as needed)
compiled_BP_path <- "SynGO data/SynGO comparison results/Compiled_BP_SynGO_092225.csv"
metadata_path <- "Metadata-Naming-ephys-grouping-cap-and-protein-ids.csv"
syngo_data_dir <- "SynGO data"
results_dir <- file.path(syngo_data_dir, "SynGO comparison results/Heatmap_GP")
dir.create(results_dir, showWarnings = FALSE, recursive = TRUE)

# Load data
compiled_BP <- read_csv(compiled_BP_path)
metadata <- read_csv(metadata_path)

# Filter for Gigaseal Preserved samples
gigaseal_samples <- metadata %>%
  filter(Gigaseal_Status == "Gigaseal Preserved") %>%
  select(`SynGO reanalysis folder name`, `Names for Analysis`)

compiled_BP_filtered <- compiled_BP %>%
  filter(Sample %in% gigaseal_samples$`SynGO reanalysis folder name`,
         `GSEA 'gene cluster' FDR corrected p-value` < 0.05) %>%
  left_join(gigaseal_samples, by = c("Sample" = "SynGO reanalysis folder name")) %>%
  rename(ShortName = `Names for Analysis`)

# Count how many cells each GO term appears in
go_counts <- compiled_BP_filtered %>%
  group_by(`GO term name`) %>%
  summarize(cell_count = n_distinct(ShortName))

# Get top 15 GO terms by average significance (across all relevant cells)
top15_gos <- compiled_BP_filtered %>%
  group_by(`GO term name`) %>%
  summarize(mean_logQ = mean(-log10(`GSEA 'gene cluster' FDR corrected p-value`))) %>%
  arrange(desc(mean_logQ)) %>%
  slice_head(n = 15) %>%
  pull(`GO term name`)

# Prepare data matrix
heatmap_data_filtered <- compiled_BP_filtered %>%
  filter(`GO term name` %in% top15_gos) %>%
  mutate(logQ = -log10(`GSEA 'gene cluster' FDR corrected p-value`)) %>%
  select(ShortName, `GO term name`, logQ) %>%
  distinct() %>%
  pivot_wider(names_from = ShortName, values_from = logQ, values_fill = 0)

# Convert to matrix
BP_matrix <- as.matrix(heatmap_data_filtered[,-1])
rownames(BP_matrix) <- heatmap_data_filtered$`GO term name`

# Match order to metadata
column_order <- intersect(colnames(BP_matrix), metadata$`Names for Analysis`)
metadata_ordered <- metadata %>% filter(`Names for Analysis` %in% column_order)


# Color function for -log10(Q)
col_fun <- colorRamp2(c(0, 2, 5, 9), c("white", "lightblue", "purple", "red"))

# Determine number of rows and columns
n_cols <- ncol(BP_matrix)
n_rows <- nrow(BP_matrix)

# Define consistent cell size
cell_w_mm <- 8
cell_h_mm <- 5

# Generate and save Heatmap
png(file.path(results_dir, paste0("Top_BP_GO_Terms_Unique_Heatmap_Complex_GP", timestamp, ".png")), 
    width = 350 + n_cols * cell_w_mm, 
    height = 300 + n_rows * cell_h_mm, 
    units = "mm", res = 300)

ht_BP_U <- Heatmap(BP_matrix,
                   name = "-log10(Q)",
                   col = col_fun,
                   cluster_rows = TRUE,
                   cluster_columns = TRUE,
                   show_column_names = TRUE,
                   show_row_names = TRUE,
                   rect_gp = gpar(col = "black", lwd = 1),
                   row_names_gp = gpar(fontsize = 14),
                   column_names_rot = 45,
                   column_names_side = "top",
                   column_names_gp = gpar(fontsize = 16),
                   column_title = "Top 15 SynGO Terms (BP)",
                   column_title_gp = gpar(fontface = "bold", fontsize = 18),  
                   heatmap_legend_param = list(
                     title = "-log10(Q-value)",
                     at = c(0, 2, 5, 9),
                     labels_gp = gpar(fontsize = 14),
                     title_gp = gpar(fontsize = 16),
                     legend_width = unit(5, "cm")
                   ),
                   width = unit(n_cols * cell_w_mm, "mm"),
                   height = unit(n_rows * cell_h_mm + 40, "mm"))

draw(ht_BP_U, heatmap_legend_side = "left")
dev.off()

cat("ComplexHeatmap saved for partially shared BP terms.\n")
