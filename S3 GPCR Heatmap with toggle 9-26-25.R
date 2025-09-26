# GPCR heatmap (IUPHAR + optional SynGO rescue)
# Expected files in the working directory:
#   - report.pg_matrix.tsv                          (DIA-NN report v1.8.1, wide)
#   - Metadata-Naming-ephys-grouping-cap-and-protein-ids.csv
#       (must contain columns: "Name in report.pg_matrix","Sample","Names for Analysis","Electrophysiology")
#   - G_to_P_db_Version_2024_4.csv               (Guide to Pharmacology export; downloaded in March 2025)
#   - SynGO-Derived_GPCR_Candidates__HGNC_Named_.csv   (optional rescue list)

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(stringr)
  library(pheatmap)
  library(RColorBrewer)
})

# ------------------------------
# User paths
# ------------------------------
setwd("~/Scripps 2/GitHub Testing/For GitHub")   # <-- adjust as needed

diann_path  <- "report.pg_matrix.tsv"
meta_path   <- "Metadata-Naming-ephys-grouping-cap-and-protein-ids.csv"
iuphar_path <- "G_to_P_db_Version_2024_4.csv"
rescue_path <- "SynGO-Derived_GPCR_Candidates__HGNC_Named_.csv"  # optional

out_dir   <- "TM_protein_detection"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
timestamp <- format(Sys.Date(), "%Y-%m-%d")

# ------------------------------
# Load DIA-NN report & metadata
# ------------------------------
stopifnot(file.exists(diann_path), file.exists(meta_path), file.exists(iuphar_path))

diann <- read_tsv(diann_path, col_types = cols(.default = "c"))
if (!"Genes" %in% names(diann)) stop("DIA-NN report must have a 'Genes' column.")

meta <- read_csv(meta_path, show_col_types = FALSE, progress = FALSE)
# ------------------------------
#  DIA-NN report sanity checks
# ------------------------------
required_cols <- c("Name in report.pg_matrix","Sample","Names for Analysis","Electrophysiology")
missing <- setdiff(required_cols, names(meta))
if (length(missing)) stop("Metadata missing required columns: ", paste(missing, collapse = ", "))

# Ensure 1:1 mapping
if (anyDuplicated(meta$`Name in report.pg_matrix`))
  stop("Duplicate 'Name in report.pg_matrix' in metadata; fix before proceeding.")
if (anyDuplicated(meta$Sample))
  stop("Duplicate 'Sample' in metadata; fix before proceeding.")

# The DIA-NN report must contain all mapped sample columns
missing_in_report <- setdiff(meta$`Name in report.pg_matrix`, names(diann))
if (length(missing_in_report))
  stop("These 'Name in report.pg_matrix' columns are not present in DIA-NN report: ",
       paste(missing_in_report, collapse = ", "))

# Pick a distinguishing protein-group/ID column for duplicate toggle 
pg_candidates <- c("Protein.Group")
pg_col <- intersect(pg_candidates, names(diann))

if (length(pg_col) == 0) {
  warning("No protein-group/ID column found in DIA-NN report; will synthesize a row ID.")
  diann$PG_ID <- seq_len(nrow(diann))
} else {
  diann$PG_ID <- diann[[pg_col[1]]]
}

# Keep Genes + sample columns
diann_sel <- diann %>% select(Genes, PG_ID, all_of(meta$`Name in report.pg_matrix`))

# Build rename map as NEW = OLD (names(new) = new, values = old)
rename_map <- setNames(meta$`Name in report.pg_matrix`, meta$Sample)

# Programmatic rename: !!! expects a named character vector c(new = "old")
diann_sel <- diann_sel %>% rename(!!!rename_map)

# ------------------------------
# Build binary detection matrix
# ------------------------------
sample_cols <- setdiff(names(diann_sel), c("Genes", "PG_ID"))

bin_mat <- diann_sel
bin_mat[sample_cols] <- lapply(bin_mat[sample_cols], function(x) {
  xnum <- suppressWarnings(as.numeric(x))
  ifelse(!is.na(xnum) & xnum > 0, 1L, 0L)
})

# If your DIA-NN 'Genes' are already HGNC symbols, uppercasing is enough.
# If they are mouse MGI symbols, insert an MGI→HGNC mapping step here.
bin_mat <- bin_mat %>% mutate(HGNC_Genes = toupper(Genes))

# ------------------------------
# IUPHAR: derive GPCR gene list (downloaded March 2025)
# ------------------------------
read_iuphar <- function(path) {
  # peek to look for headers
  x <- readr::read_csv(path, show_col_types = FALSE, progress = FALSE, n_max = 2)
  has_headers <- all(c("HGNC symbol","Family name") %in% names(x))
  if (!has_headers) {
    # try re-read with headers on row 2
    x <- readr::read_csv(path, show_col_types = FALSE, progress = FALSE, skip = 1)
  } else {
    # read full file if headers were already on row 1
    x <- readr::read_csv(path, show_col_types = FALSE, progress = FALSE)
  }
  x
}

# make IUPHAR r/reader friendly
iuphar <- read_iuphar(iuphar_path) %>%
  dplyr::mutate(
    # strip any HTML subscripts if present; harmless if absent
    `HGNC symbol` = gsub("<sub>|</sub>", "", .data$`HGNC symbol`),
    `Family name` = gsub("<sub>|</sub>", "", .data$`Family name`),
    hgnc   = toupper(.data$`HGNC symbol`),
    family = .data$`Family name`,
    # Normalize "Type" if present, else blank
    Type   = if ("Type" %in% names(.)) toupper(dplyr::coalesce(.data$Type, "")) else ""
  )

# keep only GPCR rows by IUPHARDB Type, if the column exists; otherwise skip this filter
gpcr_types <- ("GPCR")
if ("Type" %in% names(iuphar)) {
  iuphar <- dplyr::filter(iuphar, .data$Type %in% gpcr_types)
}

gpcr_targets <- iuphar %>%
  dplyr::filter(!is.na(hgnc), hgnc != "") %>%
  dplyr::distinct(hgnc, family, Type)


# ------------------------------
# Optional "rescue": SynGO GPCR annotations
# ------------------------------
if (file.exists(rescue_path)) {
  rescue <- read_csv(rescue_path, show_col_types = FALSE, progress = FALSE) %>%
    mutate(
      hgnc   = toupper(.data[[if ("hgnc_symbol" %in% names(.)) "hgnc_symbol" else names(.)[1]]]),
      family = if ("Family name" %in% names(.)) gsub("<sub>|</sub>", "", `Family name`) else NA_character_
    ) %>%
    filter(!is.na(hgnc), hgnc != "") %>%
    distinct(hgnc, family)
  gpcr_targets <- bind_rows(gpcr_targets, anti_join(rescue, gpcr_targets, by = "hgnc")) %>%
    distinct(hgnc, .keep_all = TRUE)
}

# ------------------------------
# Keep only GPCR rows from the DIA-NN matrix
# ------------------------------
bin_gpcr <- bin_mat %>% filter(HGNC_Genes %in% gpcr_targets$hgnc)

# Annotate family/type for rows
ann_rows <- gpcr_targets %>%
  transmute(HGNC_Genes = hgnc, Family = family, Type = Type)

bin_gpcr <- bin_mat %>%
  filter(HGNC_Genes %in% gpcr_targets$hgnc) %>%
  left_join(ann_rows, by = "HGNC_Genes")


# ---- Toggle: keep duplicates (per protein group) vs collapse to gene-level
# TRUE produces GRM3 & GRM3.1; FALSE produces one row per gene
keep_duplicates_as_unique <- FALSE  
if (keep_duplicates_as_unique) {
  # Keep multiple rows per gene (distinguished by PG_ID), add display names like GRM3, GRM3.1
  # Ensure we don't keep *literal* duplicate rows (same gene & same PG_ID)
  bin_gpcr <- bin_gpcr %>% distinct(HGNC_Genes, PG_ID, .keep_all = TRUE)
  
  # Create display name with suffix within each gene
  bin_gpcr <- bin_gpcr %>%
    group_by(HGNC_Genes) %>%
    arrange(HGNC_Genes, PG_ID, .by_group = TRUE) %>%
    mutate(HGNC_display = make.unique(HGNC_Genes, sep = ".")) %>%
    ungroup()
  
  # Build matrix using display names
  mat <- bin_gpcr %>%
    select(HGNC_display, all_of(meta$Sample)) %>%
    as.data.frame()
  rownames(mat) <- mat$HGNC_display
  mat$HGNC_display <- NULL
  
  # Row annotation keyed by display name
  row_ann <- bin_gpcr %>%
    select(HGNC_display, Family) %>%
    distinct(HGNC_display, .keep_all = TRUE) %>%
    tibble::column_to_rownames("HGNC_display")
  
} else {
  # Collapse to one row per gene (logical OR across sample columns)
  collapsed <- bin_gpcr %>%
    group_by(HGNC_Genes) %>%
    summarise(across(all_of(meta$Sample), ~ as.integer(any(.x > 0))), .groups = "drop")
  
  # Build matrix
  mat <- as.data.frame(collapsed)
  rownames(mat) <- mat$HGNC_Genes
  mat$HGNC_Genes <- NULL
  
  # Row annotation
  row_ann <- bin_gpcr %>%
    group_by(HGNC_Genes) %>%
    summarise(Family = first(na.omit(Family)), .groups = "drop") %>%
    tibble::column_to_rownames("HGNC_Genes")
}
# ------------------------------
# Order columns by electrophysiology; label with "Names for Analysis"
# ------------------------------
meta$Electrophysiology <- factor(meta$Electrophysiology,
                                 levels = c("Preserved","Lost","None","Neuron Torn"))
meta <- meta %>% arrange(Electrophysiology)

# Keep only columns that exist, and reorder to match meta$Sample
mat <- mat[, intersect(meta$Sample, colnames(mat)), drop = FALSE]

# Drop all-zero rows
mat <- mat[rowSums(mat) > 0, , drop = FALSE]

# Remap column names to friendly labels for plotting
name_map <- setNames(meta$`Names for Analysis`, meta$Sample)
colnames(mat) <- name_map[colnames(mat)]

# Align row annotation to mat
row_ann <- row_ann[rownames(mat), , drop = FALSE]

# ------------------------------
# Save CSV 
# ------------------------------
write_csv(
  tibble::tibble(Gene = rownames(mat)) %>% bind_cols(as.data.frame(mat)),
  file.path(out_dir, paste0("GPCR_detection_matrix_v2", timestamp, ".csv"))
)
# ------------------------------
# Heatmap
# ------------------------------
# Row annotation (Family)
row_ann <- bin_gpcr %>%
  select(HGNC_Genes, Family) %>%
  filter(HGNC_Genes %in% rownames(mat)) %>%
  distinct(HGNC_Genes, .keep_all = TRUE) %>%
  tibble::column_to_rownames("HGNC_Genes")

# Column annotation (Electrophysiology)
ann_col <- meta %>%
  select(Sample, `Names for Analysis`, Electrophysiology) %>%
  filter(`Names for Analysis` %in% colnames(mat)) %>%
  tibble::column_to_rownames("Names for Analysis")

annotation_colors <- list(
  Electrophysiology = c(Preserved="green", Lost="orange", None="red", `Neuron Torn`="grey")
)

# Distinct palette per family (stable across reruns)
families <- sort(unique(row_ann$Family))
if (length(families) == 0) families <- "Unknown"
family_pal <- setNames(
  colorRampPalette(brewer.pal(min(max(length(families), 3), 8), "Set3"))(length(families)),
  families
)

png(file.path(out_dir, paste0("heatmap_gpcrs_v2", timestamp, ".png")),
    width = 3600, height = 1800, res = 300)
pheatmap(as.matrix(mat),
         cluster_rows = FALSE, cluster_cols = TRUE,
         annotation_row = data.frame(Family = row_ann$Family,
                                     row.names = rownames(row_ann)),
         annotation_col = data.frame(Electrophysiology = ann_col$Electrophysiology,
                                     row.names = rownames(ann_col)),
         annotation_colors = c(list(Electrophysiology = annotation_colors$Electrophysiology),
                               list(Family = family_pal)),
         legend = FALSE, annotation_legend = TRUE,
         color = colorRampPalette(c("white","blue"))(50),
         main = "GPCR subunits detected across neurons",
         na_col = "white", fontsize_row = 10, fontsize_col = 12,
         angle_col = 45, annotation_names_row = FALSE, annotation_names_col = FALSE,
         border_color = "black")
dev.off()

cat("GPCR matrix and heatmap written to:", out_dir, "\n")
