############################################################
# NOTE:
# This script depends on the compiled SynGO results generated 
# in the Figure 5 script. 
# Run Figure 5 script first to produce the required inputs.
############################################################

# Load necessary libraries
library(dplyr)
library(readr)
library(tidyr)
library(UpSetR)


# Set working directory
setwd("~/Scripps 2/GitHub Testing/For GitHub")

# Define paths (modify as needed)
compiled_BP_path <- "SynGO data/SynGO comparison results/Compiled_BP_SynGO_092225.csv"
metadata_path <- "Metadata-Naming-ephys-grouping-cap-and-protein-ids.csv"
output_upset_dir <- "SynGO data/SynGO comparison results/UpSet_Plots_GP"
dir.create(output_upset_dir, showWarnings = FALSE, recursive = TRUE)

# Load compiled BP data and metadata
compiled_BP <- read_csv(compiled_BP_path)
metadata <- read_csv(metadata_path)

# Generate timestamp for output files
timestamp <- format(Sys.Date(), "%Y-%m-%d")

# Map sample names to shorter labels
name_mapping <- setNames(metadata$`Names for Analysis`, metadata$`SynGO reanalysis folder name`)

# Filter metadata for "Gigaseal Preserved" cells only
gigaseal_preserved_samples <- metadata %>%
  filter(Gigaseal_Status == "Gigaseal Preserved") %>%
  pull(`SynGO reanalysis folder name`)

# Filter BP data for significant terms (FDR < 0.05) and Gigaseal Preserved samples only
significant_BP_filtered <- compiled_BP %>%
  filter(`GSEA 'gene cluster' FDR corrected p-value` < 0.05,
         Sample %in% gigaseal_preserved_samples) %>%
  select(Sample, `GO term name`) %>%
  distinct() %>%
  mutate(Sample = name_mapping[Sample])  # Rename for clarity

# Transform into wide format for UpSet plot
BP_upset_data <- significant_BP_filtered %>%
  mutate(present = 1) %>%
  pivot_wider(names_from = Sample, values_from = present, values_fill = 0)

# Convert to data frame for UpSet plot, set GO terms as rownames
BP_upset_matrix <- as.data.frame(BP_upset_data)
rownames(BP_upset_matrix) <- BP_upset_matrix$`GO term name`
BP_upset_matrix <- BP_upset_matrix[,-1]

# Generate and save UpSet plot
png(file.path(output_upset_dir, paste0("BP_SynGO_UpSet_GigasealPreserved_", timestamp, ".png")), 
    width = 9000, height = 3000, res = 300)

upset(BP_upset_matrix,
      nsets = ncol(BP_upset_matrix),
      order.by = "freq",
      sets.bar.color = "black",
      main.bar.color = "black",
      matrix.color = "black",
      mainbar.y.label = "Shared SynGO BPs",
      sets.x.label = "SynGO BPs per Neuron",
      point.size = 12,
      text.scale = c(6, 5, 5, 5, 6, 7.0),  # Y title, Y ticks, Set title, Set ticks, Set names, matrix labels
)

dev.off()

# Save intersection data to CSV with explicit GO term names
BP_upset_output <- BP_upset_matrix %>%
  mutate(`GO term name` = rownames(BP_upset_matrix)) %>% 
  relocate(`GO term name`)

write_csv(BP_upset_output, 
          file.path(output_upset_dir, paste0("BP_SynGO_UpSet_Matrix_GigasealPreserved_with_GO_names_", timestamp, ".csv")))

cat("UpSet plot and CSV file for BP domain (Gigaseal Preserved only) successfully generated!\n")
