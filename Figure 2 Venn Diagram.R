# Load necessary libraries
library(readr)
library(dplyr)
library(eulerr)
library(ggplot2)
library(RColorBrewer)

# Set working directory
setwd("~/Scripps 2/GitHub Testing/For GitHub")
output_dir <- "Within Group Comparison"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# Paths to files
diann_report_path <- "report.pg_matrix.tsv"
metadata_path <- "Metadata-Naming-ephys-grouping-cap-and-protein-ids.csv"

# Generate timestamp for output files
timestamp <- format(Sys.Date(), "%Y-%m-%d")

# Define DIA-NN samples of interest
samples_of_interest <- c(
  "/gpfs/group/home/titusj/1303_larry/20240709_Ast_Neo_Aurora_25um-75umID_FAIMS_CV-50_20mz_DIA_400-1000_50maxIT_800agc_25nce_36m_pilot-2_neuron_bcmk4.raw.dia",
  "/gpfs/group/home/titusj/1303_larry/20240709_Ast_Neo_Aurora_25um-75umID_FAIMS_CV-50_20mz_DIA_400-1000_50maxIT_800agc_25nce_36m_pilot-2_neuron_bcmk6.raw.dia",
  "/gpfs/group/home/titusj/1303_larry/20240709_Ast_Neo_Aurora_25um-75umID_FAIMS_CV-50_20mz_DIA_400-1000_50maxIT_800agc_25nce_36m_pilot-2_neuron_bcmk7.raw.dia"
)

# Load DIA-NN data
diann_data <- read_tsv(diann_report_path, col_types = cols())
metadata <- read_csv(metadata_path, col_types = cols())

# Map to friendly names
name_mapping <- setNames(metadata$`Names for Analysis`,
                         metadata$`Name in report.pg_matrix`)

# Filter DIA-NN data
diann_filtered <- diann_data %>%
  select(Genes, all_of(samples_of_interest))
colnames(diann_filtered)[-1] <- name_mapping[samples_of_interest]

# Extract lists of detected proteins per sample
protein_lists <- lapply(colnames(diann_filtered)[-1], function(col) {
  unique(diann_filtered$Genes[diann_filtered[[col]] > 0])
})
names(protein_lists) <- colnames(diann_filtered)[-1]

# Create named logical vectors for eulerr input
all_proteins <- unique(unlist(protein_lists))

eulerr_input <- lapply(protein_lists, function(gene_list) {
  all_proteins %in% gene_list
})
eulerr_df <- as.data.frame(eulerr_input)
colnames(eulerr_df) <- names(protein_lists)

# Make proportional Euler diagram
fit <- euler(eulerr_df)

# Define custom colors
custom_colors <- c("#1b9e77", "#d95f02", "#7570b3")  # Green, orange, purple

# Plot with custom colors
venn_plot <- plot(fit,
                  fills = list(fill = custom_colors, alpha = 0.6),
                  edges = list(col = "black"),
                  labels = list(font = 3, cex = 1.2),
                  quantities = list(font = 1, cex = 1.2),
                  main = list(label = "Total Protein Overlap", cex = 1.6, font = 3))

# Save plot
ggsave(filename = file.path(output_dir, paste0("Proportional_Venn_Diagram", timestamp, ".png")),
       plot = venn_plot, width = 8, height = 8, dpi = 300)