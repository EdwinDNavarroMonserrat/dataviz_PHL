#!/usr/bin/env Rscript

# Script: Grandeur results visualization & QC tables
# Author: Edwin Daniel Navarro Monserrat
# Description:
#   Parses Grandeur workflow outputs (AMRFinderPlus, SeqSero2, Shigatyper, MLST, FastANI, QUAST, MASH, fastp)
#   and produces ggplot objects + reactable tables for taxonomy, serotyping, QC, fastp, and AMR/plasmids.
#   Thresholds were made in accordance to PulseNet guidelines.
#   Note: if you want to have pagination on the reactable (i.e. tabs) switch pagination = FALSE to TRUE
#   Currently made preferably for interactive use, but can be modified for non-interactive usage
#   Set base dir 

# Loading R packages: install if you haven't done so
library(tidyverse)
library(ggthemes)
library(reactable)
library(jsonlite)

# Helper to always show full table (no pagination)
reactable_full <- function(df, ...) {
  reactable(
    df,
    pagination = FALSE,
    defaultPageSize = max(1, nrow(df)),
    minRows = min(10, max(1, nrow(df))),
    ...
  )
}

# 1) Configuration
base_dir <- path.expand("insert path where grandeur results are located")

safe_read_tsv <- function(path, ...) {
  if (!file.exists(path)) {
    warning(sprintf("File not found: %s", path))
    return(tibble())
  }
  readr::read_tsv(path, show_col_types = FALSE, ...)
}

# 2) AMRFinderPlus
amr_tsv_df_raw <- safe_read_tsv(file.path(base_dir, "amrfinder", "amrfinderplus.txt"))

if ("Name" %in% names(amr_tsv_df_raw)) {
  amr_tsv_df <- amr_tsv_df_raw |> mutate(sample = as.character(Name))
} else {
  amr_tsv_df <- amr_tsv_df_raw |> mutate(sample = as.character(sample))
}

amr_tsv_df <- amr_tsv_df |> filter(!grepl("Undetermined", sample, ignore.case = TRUE))

amr_summary <- amr_tsv_df |>
  filter(Type %in% c("AMR","STRESS","VIRULENCE")) |>
  count(sample, Type, name = "Gene_Count")

AMR_VIR_STRESS <- ggplot(amr_summary, aes(x = reorder(sample, -Gene_Count),
                                          y = Gene_Count, fill = Type)) +
  geom_bar(stat = "identity", position = "stack", color = "black") +
  coord_flip() +
  scale_fill_manual(values = c("AMR" = "cyan2",
                               "STRESS" = "green4",
                               "VIRULENCE" = "red2")) +
  theme_clean() +
  labs(title = "AMR, Stress and Virulence Genes per Sample",
       x = "Sample", y = "Number of Genes", fill = "Gene Type")

amr_heatmap <- amr_tsv_df |>
  filter(Type == "AMR") |>
  select(sample, `Element symbol`, `% Identity to reference`) |>
  group_by(`Element symbol`) |>
  mutate(gene_freq = n()) |>
  ungroup() |>
  mutate(`Element symbol` = forcats::fct_reorder(`Element symbol`, gene_freq))

AMR_HM_plot <- ggplot(amr_heatmap,
                      aes(x = `Element symbol`, y = sample, fill = `% Identity to reference`)) +
  geom_tile(color = "black") +
  scale_fill_gradient2(low = "bisque1", mid = "gold", high = "red2", midpoint = 70, name = "% Identity") +
  theme_clean() +
  labs(title = "AMR Genes per Sample", x = "AMR Gene", y = "Sample")

amr_heatmap_all <- amr_tsv_df |> select(sample, `Element symbol`, Type, `% Identity to reference`)

AMR_STRESS_VIR_HM_plot <- ggplot(amr_heatmap_all,
                                 aes(x = `Element symbol`, y = sample,
                                     fill = `% Identity to reference`)) +
  geom_tile(color = "black", size = 0.5) +
  scale_fill_gradient2(low = "bisque1", mid = "gold", high = "red2", midpoint = 70) +
  facet_wrap(~Type, scales = "free", shrink = FALSE) +
  theme_clean() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, face = "bold")) +
  labs(title = "AMR, Stress, and Virulence Genes per Sample", x = "Gene", y = "Sample")

# 3) Grandeur summary
grandeur_results <- safe_read_tsv(file.path(base_dir, "grandeur_summary.tsv")) |> select(!any_of(c("file", "version")))

taxonomy_table <- grandeur_results |>
  select(sample, predicted_organism, mlst_matching_pubmlst_scheme,
         fastani_top_organism, fastani_top_ani_estimate,
         mash_organism, `mash_mash-distance`)

taxa_reactable <- reactable_full(
  taxonomy_table,
  searchable = TRUE, striped = TRUE, highlight = TRUE, bordered = TRUE, resizable = TRUE,
  defaultColDef = colDef(align = "center"),
  columns = list(
    fastani_top_ani_estimate = colDef(
      name = "FastANI % Identity", format = colFormat(digits = 2),
      style = function(value) {
        if (is.na(value)) return(NULL)
        if (value > 97)       list(background = "lightgreen", fontWeight = "bold")
        else if (value > 90)  list(background = "khaki")
        else                  list(background = "salmon")
      }
    ),
    `mash_mash-distance` = colDef(
      name = "Mash Distance", format = colFormat(digits = 4),
      style = function(value) {
        if (is.na(value)) return(NULL)
        if (value < 0.05) list(background = "lightgreen") else list(background = "salmon")
      }
    )
  )
)

serotyping_table <- grandeur_results |>
  select(sample,
         any_of(c("seqsero2_predicted_antigenic_profile",
                  "seqsero2_predicted_serotype",
                  "serotypefinder_Serotype_O",
                  "serotypefinder_Serotype_H",
                  "shigatyper_prediction",
                  "shigatyper_hit")))

serotype_table <- reactable_full(
  serotyping_table,
  searchable = TRUE, filterable = TRUE, striped = TRUE, highlight = TRUE, bordered = TRUE,
  defaultColDef = colDef(align = "center")
)

taxonomy_serotyping_table <- grandeur_results |>
  select(sample, predicted_organism, mlst_matching_pubmlst_scheme,
         fastani_top_organism, fastani_top_ani_estimate,
         mash_organism, `mash_mash-distance`,
         any_of(c("seqsero2_predicted_antigenic_profile",
                  "seqsero2_predicted_serotype",
                  "serotypefinder_Serotype_O",
                  "serotypefinder_Serotype_H",
                  "shigatyper_prediction",
                  "shigatyper_hit")))

taxonomy_serotype_reactable <- reactable_full(
  taxonomy_serotyping_table,
  searchable = TRUE, filterable = TRUE, striped = TRUE, highlight = TRUE, bordered = TRUE,
  groupBy = "predicted_organism"
)

# 4) QUAST QC
Quast <- safe_read_tsv(file.path(base_dir, "quast", "quast_report.tsv")) |>
  select(sample, `# contigs`, N50, `Total length`, `Coverage >= 1x (%)`, `Avg. coverage depth`) |>
  left_join(taxonomy_table |> select(sample, predicted_organism), by = "sample")

qc_rules <- tibble::tribble(
  ~group,          ~min_cov, ~min_n50, ~max_contigs, ~len_min_mb, ~len_max_mb,
  "Campylobacter",      20,     75000,        400,         1.4,        2.2,
  "Listeria",           20,    200000,        100,         2.8,        3.1,
  "Escherichia",        40,     75000,        500,         4.9,        5.9,
  "Shigella",           40,     75000,        500,         4.2,        4.9,
  "Salmonella",         30,    100000,        400,         4.4,        6.2,
  "Vibrio",             40,    100000,        300,         3.8,        5.3,
  "Cronobacter",        30,     75000,        400,           NA,         NA
)

organism_to_group <- function(x){
  dplyr::case_when(
    stringr::str_detect(x, regex("Campylobacter", TRUE)) ~ "Campylobacter",
    stringr::str_detect(x, regex("Listeria", TRUE))      ~ "Listeria",
    stringr::str_detect(x, regex("Escherichia", TRUE))   ~ "Escherichia",
    stringr::str_detect(x, regex("Shigella", TRUE))      ~ "Shigella",
    stringr::str_detect(x, regex("Salmonella", TRUE))    ~ "Salmonella",
    stringr::str_detect(x, regex("Vibrio", TRUE))        ~ "Vibrio",
    stringr::str_detect(x, regex("Cronobacter", TRUE))   ~ "Cronobacter",
    TRUE ~ NA_character_
  )
}

Quast <- Quast |> mutate(group = organism_to_group(predicted_organism)) |> left_join(qc_rules, by = "group")

qc_table_clean <- Quast |>
  mutate(
    cov_flag     = if_else(!is.na(min_cov)     & `Avg. coverage depth` >= min_cov, "PASS", "FAIL"),
    n50_flag     = if_else(!is.na(min_n50)     & N50 >= min_n50,                   "PASS", "FAIL"),
    contigs_flag = if_else(!is.na(max_contigs) & `# contigs` <= max_contigs,       "PASS", "FAIL"),
    length_flag  = if_else(!is.na(len_min_mb)  & !is.na(len_max_mb) & (`Total length`/1e6 >= len_min_mb & `Total length`/1e6 <= len_max_mb), "PASS", "FAIL")
  ) |>
  mutate(
    `Avg. coverage depth` = as.numeric(`Avg. coverage depth`),
    N50                    = as.numeric(N50),
    `# contigs`            = as.numeric(`# contigs`),
    `Total length`         = as.numeric(`Total length`)
  )

qc_table_display <- qc_table_clean |>
  transmute(sample, predicted_organism, `Avg. coverage depth`, N50, `# contigs`, `Total length`, cov_flag, n50_flag, contigs_flag, length_flag)

._qc <- qc_table_display

Quast_QC_reactable <- reactable_full(
  qc_table_display,
  searchable = TRUE, striped = TRUE, highlight = TRUE, bordered = TRUE,
  defaultColDef = colDef(align = "center"),
  columns = list(
    `Avg. coverage depth` = colDef(
      format = colFormat(digits = 1),
      style = \(value, index) {
        flag <- ._qc$cov_flag[index]
        if (is.na(flag)) return(list(backgroundColor = "lightgrey"))
        if (flag == "PASS") list(backgroundColor = "lightgreen")
        else list(backgroundColor = "salmon", fontWeight = "bold")
      }
    ),
    N50 = colDef(
      format = colFormat(separators = TRUE),
      style = \(value, index) {
        flag <- ._qc$n50_flag[index]
        if (is.na(flag)) return(list(backgroundColor = "lightgrey"))
        if (flag == "PASS") list(backgroundColor = "lightgreen")
        else list(backgroundColor = "salmon", fontWeight = "bold")
      }
    ),
    `# contigs` = colDef(
      style = \(value, index) {
        flag <- ._qc$contigs_flag[index]
        if (is.na(flag)) return(list(backgroundColor = "lightgrey"))
        if (flag == "PASS") list(backgroundColor = "lightgreen")
        else list(backgroundColor = "salmon", fontWeight = "bold")
      }
    ),
    `Total length` = colDef(
      format = colFormat(separators = TRUE),
      style = \(value, index) {
        flag <- ._qc$length_flag[index]
        if (is.na(flag)) return(list(backgroundColor = "lightgrey"))
        if (flag == "PASS") list(backgroundColor = "lightgreen")
        else list(backgroundColor = "salmon", fontWeight = "bold")
      }
    ),
    cov_flag     = colDef(show = FALSE),
    n50_flag     = colDef(show = FALSE),
    contigs_flag = colDef(show = FALSE),
    length_flag  = colDef(show = FALSE)
  )
)

# 5) fastp summary
fastp_dir <- file.path(base_dir, "fastp")
json_files <- if (dir.exists(fastp_dir)) list.files(fastp_dir, pattern = "\\.json$", full.names = TRUE) else character(0)

fastp_summary <- if (length(json_files) > 0) purrr::map_dfr(json_files, function(file) {
  dat <- jsonlite::fromJSON(file)
  tibble::tibble(
    sample                 = tools::file_path_sans_ext(basename(file)),
    total_reads_before     = dat$summary$before_filtering$total_reads,
    total_reads_after      = dat$summary$after_filtering$total_reads,
    percent_reads_passing  = (dat$summary$after_filtering$total_reads / dat$summary$before_filtering$total_reads) * 100,
    q30_rate_after         = dat$summary$after_filtering$q30_rate * 100,
    mean_read_length_after = dat$summary$after_filtering$read1_mean_length,
    reads_too_many_N       = dat$filtering_result$too_many_N_reads,
    reads_too_short        = dat$filtering_result$too_short_reads
  )
}) else tibble()

fastp_reactable <- reactable_full(
  fastp_summary,
  searchable = TRUE, striped = TRUE, highlight = TRUE, bordered = TRUE,
  defaultColDef = colDef(align = "center"),
  columns = list(
    percent_reads_passing = colDef(
      name = "% Reads Passing",
      format = colFormat(digits = 1),
      style = function(value) {
        if (is.na(value)) return(list(backgroundColor = "lightgrey"))
        if (value >= 90)   list(backgroundColor = "lightgreen")
        else               list(backgroundColor = "salmon", fontWeight = "bold")
      }
    ),
    q30_rate_after = colDef(
      name = "Q30 (%)",
      format = colFormat(digits = 1),
      style = function(value) {
        if (is.na(value)) return(list(backgroundColor = "lightgrey"))
        if (value >= 80)   list(backgroundColor = "lightgreen")
        else               list(backgroundColor = "khaki")
      }
    ),
    mean_read_length_after = colDef(
      name = "Mean Read Length (bp)",
      format = colFormat(digits = 0),
      style = function(value) {
        if (is.na(value)) return(list(backgroundColor = "lightgrey"))
        if (value >= 135)  list(backgroundColor = "lightgreen")
        else               list(backgroundColor = "salmon", fontWeight = "bold")
      }
    )
  )
)

# 6) AMR & plasmids summary
plasmidfinder_table <- grandeur_results |>
  select(sample, any_of(c("plasmidfinder_plasmid_(identity)"))) |>
  rename(plasmids = `plasmidfinder_plasmid_(identity)`) |>
  mutate(
    plasmids = gsub("\\[|\\]|'", "", plasmids),
    plasmids = gsub(",\\s*", ", ", plasmids)
  )

amr_summary_table <- amr_tsv_df |>
  filter(Type == "AMR") |>
  group_by(sample) |>
  summarise(
    AMR_gene_count = dplyr::n(),
    AMR_genes = paste(unique(`Element symbol`), collapse = ", "),
    .groups = "drop"
  )

amr_plasmid_table <- full_join(amr_summary_table, plasmidfinder_table, by = "sample") |> arrange(sample)

AMR_plasmid_reactable <- reactable_full(
  amr_plasmid_table,
  searchable = TRUE, striped = TRUE, highlight = TRUE, bordered = TRUE
)


