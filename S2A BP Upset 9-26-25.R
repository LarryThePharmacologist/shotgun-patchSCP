############################################################
# NOTE:
# This script depends on the compiled SynGO results generated 
# in the Figure 5 script. 
# Run Figure 5 script first to produce the required inputs.
############################################################
# Load necessary libraries
library(dplyr)
library(readr)
library(UpSetR)
library(lubridate)
library(tidyr)

# Set working directories
setwd("~/Scripps 2/GitHub Testing/For GitHub")

# Define paths (modify as needed)
compiled_BP_path <- "SynGO data/SynGO comparison results/Compiled_BP_SynGO_092225.csv"
output_upset_dir <- "SynGO data/SynGO comparison results/UpSet_Plots_All_Neurons"
dir.create(output_upset_dir, showWarnings = FALSE, recursive = TRUE)

# Load previously compiled BP SynGO data
timestamp <- format(Sys.Date(), "%m%d%y")
compiled_BP <- read_csv(compiled_BP_path)

# Load Metadata to replace sample names
metadata_path <- "Metadata-Naming-ephys-grouping-cap-and-protein-ids.csv"
metadata <- read_csv(metadata_path)

# Map original SynGO sample folder names to shorter display names
name_mapping <- setNames(metadata$`Names for Analysis`, metadata$`SynGO reanalysis folder name`)

# Filter significant BP terms based on FDR < 0.05
significant_BP <- compiled_BP %>%
  filter(`GSEA 'gene cluster' FDR corrected p-value` < 0.05) %>%
  select(Sample, `GO term name`) %>%
  distinct() %>%
  mutate(Sample = name_mapping[Sample])  # Replace sample names here!

# Transform into wide format for UpSet plot
BP_upset_data <- significant_BP %>%
  mutate(present = 1) %>%
  pivot_wider(names_from = Sample, values_from = present, values_fill = 0)

# Set GO term names as rownames and remove them from columns
BP_upset_matrix <- as.data.frame(BP_upset_data)
rownames(BP_upset_matrix) <- BP_upset_matrix$`GO term name`
BP_upset_matrix <- BP_upset_matrix[,-1]

# Generate and save the UpSet plot with customized font sizes
png(file.path(output_upset_dir, paste0("BP_SynGO_UpSet_All_Neurons", timestamp, ".png")), 
    width = 9500, height = 7000, res = 300)

upset(BP_upset_matrix,
      nsets = ncol(BP_upset_matrix),
      nintersects = NA, 
      order.by = "freq",
      main.bar.color = "black",
      sets.bar.color = "black",
      matrix.color = "black",
      mainbar.y.label = "Shared SynGO BPs",
      sets.x.label = "SynGO BPs per Neuron",
      text.scale = c(6, 5, 5, 5, 6, 7.0),  # Y title, Y ticks, Set title, Set ticks, Set names, matrix labels
      point.size = 6) # dot size 
dev.off()


# Save the intersection matrix for future reference
write_csv(BP_upset_data, file.path(output_upset_dir, paste0("BP_SynGO_UpSet_Matrix_All_Neurons_", timestamp, ".csv")))

print("UpSet plot for all neurons generated successfully.")
