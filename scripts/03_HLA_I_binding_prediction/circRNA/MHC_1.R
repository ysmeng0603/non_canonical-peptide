# ============================================================
# circRNA HLA-I assignment and netMHCpan analysis
# ============================================================

suppressPackageStartupMessages({
  library(Biostrings)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(purrr)
  library(tibble)
  library(ggplot2)
  library(ComplexHeatmap)
  library(grid)
})


# ============================================================
# 1. Build the sample-circRNA-HLA-I table
# ============================================================

filter_ToTable <- read.delim(
  "./circRNA/part2/DEseq2/filter_ToTable.txt",
  header = TRUE
)

sample_circRNA <- filter_ToTable %>%
  dplyr::select(group, comb, matches("SRR478")) %>%
  tidyr::pivot_longer(cols = 3:ncol(.)) %>%
  dplyr::filter(value > 0, group != "unknown")

colnames(sample_circRNA)[2] <- "circRNA"

circbank_circbase_HostGene <- read.delim(
  "./circRNA/a_info/circbank_circbase_HostGene.txt"
)

sample_circRNA <- sample_circRNA %>%
  dplyr::left_join(
    circbank_circbase_HostGene,
    by = "circRNA",
    relationship = "many-to-many"
  ) %>%
  dplyr::select(circRNA, circbaseID, name) %>%
  unique()

unique_all_database_info <- read.delim(
  "./circRNA/a_info/unique_all_database_info.txt",
  header = TRUE
) %>%
  dplyr::filter(Database == "circAtlas3")

circbase <- sample_circRNA %>%
  dplyr::filter(!is.na(circbaseID))

colnames(circbase)[2] <- "DB_name"

circAtlas3 <- sample_circRNA %>%
  dplyr::filter(is.na(circbaseID)) %>%
  dplyr::left_join(unique_all_database_info, by = "circRNA") %>%
  dplyr::select(circRNA, DB_name, name)

sample_circRNA <- rbind(circbase, circAtlas3)


# ============================================================
# 2. Retain circRNAs with BSJ-spanning ORFs
# ============================================================

bsj_files <- c(
  "./circRNA_peptide3/part2_BSJ/seq_repeat_4/H_eRA_down/contain_BSJ_result/Files/select_H_eRA_down_list.txt",
  "./circRNA_peptide3/part2_BSJ/seq_repeat_4/H_eRA_up/contain_BSJ_result/Files/select_H_eRA_up_list.txt",
  "./circRNA_peptide3/part2_BSJ/seq_repeat_4/H_RA_down/contain_BSJ_result/Files/select_H_RA_down_list.txt",
  "./circRNA_peptide3/part2_BSJ/seq_repeat_4/H_RA_up/contain_BSJ_result/Files/select_H_RA_up_list.txt",
  "./circRNA_peptide3/part2_BSJ/seq_repeat_4/RA_eRA_down/contain_BSJ_result/Files/select_RA_eRA_down_list.txt"
)

contain_BSJ <- lapply(bsj_files, read.delim) %>%
  dplyr::bind_rows() %>%
  tidyr::separate(
    col = V2,
    into = c("col1", "col2"),
    sep = "\\|"
  ) %>%
  tidyr::separate(
    col = col2,
    into = c("col3", "col4"),
    sep = ":"
  ) %>%
  unique()

sample_circRNA <- sample_circRNA %>%
  dplyr::filter(DB_name %in% contain_BSJ$col3)


# ============================================================
# 3. Add HLA-I information
# ============================================================

HLA_1 <- read.delim2(
  "./circRNA_peptide2/part4_netMHCpan/MHC_1/HLA_1.txt",
  header = FALSE
)


colnames(HLA_1) <- c("name", "HLA_1")

sample_circ_HLA_1 <- merge(sample_circRNA, HLA_1, by = "name")

write.table(
  sample_circ_HLA_1,
  "./circRNA_peptide3/part4_netMHCpan/Info/sample_circ_HLA_1.txt",
  quote = FALSE,
  row.names = FALSE,
  sep = "\t"
)



# ============================================================
# 4. Generate sample-specific HLA-I protein FASTA files
# ============================================================

protein_source_dir <-
  "./circRNA_peptide3/part4_netMHCpan/all_select_p_result"

write_sample_fastas <- function(sample_info, output_dir) {

  sample_info <- sample_info %>%
    dplyr::mutate(
      circ_path = paste0(
        protein_source_dir,
        "/p_result_",
        DB_name,
        "_selected.fasta"
      )
    )

  samples <- unique(sample_info$name)

  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  for (sample in samples) {

    df <- sample_info[sample_info$name == sample, ]
    all_sequences <- AAStringSet()

    for (file_path in df$circ_path) {
      fasta_content <- readAAStringSet(filepath = file_path)
      all_sequences <- c(all_sequences, fasta_content)
    }

    writeXStringSet(
      all_sequences,
      file.path(output_dir, paste0(sample, ".fasta"))
    )
  }
}

write_sample_fastas(
  sample_circ_HLA_1,
  "./circRNA_peptide3/part4_netMHCpan/MHC_1/SRR_Files"
)



# ============================================================
# 5. Build the HLA-I short-ID mapping table
# ============================================================

collect_fasta_headers <- function(input_dir) {

  fasta_files <- list.files(
    input_dir,
    pattern = "\\.fasta$",
    full.names = TRUE
  )

  unlist(
    lapply(
      fasta_files,
      function(file) names(readAAStringSet(file))
    ),
    use.names = FALSE
  )
}

MHC1_headers <- collect_fasta_headers(
  "./circRNA_peptide3/part4_netMHCpan/MHC_1/SRR_Files"
)


write.table(
  data.frame(ori_ID = MHC1_headers),
  "./circRNA_peptide3/part4_netMHCpan/Info/short_ID/MHC1_short_ID.txt",
  quote = FALSE,
  row.names = FALSE,
  sep = "\t"
)


short_ID_map <- data.frame(
  ori_ID = unique(MHC1_headers),
  stringsAsFactors = FALSE
) %>%
  tidyr::separate(
    col = ori_ID,
    into = c("new1", "new2"),
    sep = ":"
  ) %>%
  dplyr::select(new1) %>%
  tidyr::separate(
    col = new1,
    into = c("new1", "new2"),
    sep = "\\|"
  ) %>%
  dplyr::select(new2) %>%
  tidyr::separate(
    col = new2,
    into = paste0("new", 1:4),
    sep = "_"
  ) %>%
  dplyr::mutate(
    Ori_ID = str_c(new1, new2, new3, new4, sep = "_"),
    short_ID = str_c(new1, new4, sep = "_")
  ) %>%
  dplyr::select(Ori_ID, short_ID) %>%
  unique() %>%
  na.omit()

write.table(
  short_ID_map,
  "./circRNA_peptide3/part4_netMHCpan/Info/short_ID/all_short_ID.txt",
  quote = FALSE,
  row.names = FALSE,
  sep = "\t"
)


# ============================================================
# 6. Generate HLA-I FASTA files with shortened sequence IDs
# ============================================================

clean_fasta <- function(input_file, output_file, name_map) {

  sequences <- readAAStringSet(input_file, format = "fasta")

  old_names <- names(sequences)
  new_names <- old_names

  for (i in seq_along(old_names)) {
    for (j in seq_len(nrow(name_map))) {
      if (grepl(name_map$Ori_ID[j], old_names[i])) {
        new_names[i] <- name_map$short_ID[j]
        break
      }
    }
  }

  names(sequences) <- new_names

  writeXStringSet(
    sequences,
    filepath = output_file,
    format = "fasta"
  )
}

write_short_id_fastas <- function(input_dir, output_dir, name_map) {

  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  fasta_files <- list.files(
    input_dir,
    pattern = "\\.fasta$",
    full.names = TRUE
  )

  for (fasta_file in fasta_files) {
    clean_fasta(
      input_file = fasta_file,
      output_file = file.path(output_dir, basename(fasta_file)),
      name_map = name_map
    )
  }
}

write_short_id_fastas(
  "./circRNA_peptide3/part4_netMHCpan/MHC_1/SRR_Files",
  "./circRNA_peptide3/part4_netMHCpan/MHC_1/SRR_files_shortID",
  short_ID_map
)



# ============================================================
# 7. Build the netMHCpan MHC-I task table
# ============================================================

MHC1_file_paths <- sample_circ_HLA_1 %>%
  dplyr::select(name, HLA_1) %>%
  unique() %>%
  dplyr::mutate(
    index = seq_len(n()),
    file_path = paste0(
      "./circRNA_peptide3/part4_netMHCpan/MHC_1/SRR_files_shortID/",
      name,
      ".fasta"
    )
  ) %>%
  dplyr::select(index, file_path, name, HLA_1)

colnames(MHC1_file_paths)[3] <- "file_name"

write.table(
  MHC1_file_paths,
  "./circRNA_peptide3/part4_netMHCpan/MHC_1/bash_need/MHC1_file_paths.txt",
  quote = FALSE,
  row.names = FALSE,
  sep = "\t"
)


# ============================================================
# 8. Parse and filter netMHCpan MHC-I results
# ============================================================

netmhcpan_result_dir <-
  "./circRNA_peptide3/part4_netMHCpan/MHC_1/netMHCpan_1_result"

filter_result_dir <-
  "./circRNA_peptide3/part4_netMHCpan/MHC_1/filter_result"

dir.create(filter_result_dir, recursive = TRUE, showWarnings = FALSE)

parse_netmhcpan_result <- function(file) {

  data_read <- read.table(
    file,
    sep = "\t",
    header = FALSE,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )

  ID <- data_read[3:nrow(data_read), 1:3]
  colnames(ID) <- c("Pos", "Peptide", "ID")

  col <- data_read[1:2, 4:ncol(data_read)]
  rownames(col) <- c("HLA_type", "field")

  col <- col %>%
    t() %>%
    as.data.frame(stringsAsFactors = FALSE)

  col$HLA_type <- ifelse(
    str_length(col$HLA_type) == 0,
    NA,
    col$HLA_type
  )

  col <- col %>%
    tidyr::fill(HLA_type, .direction = "down") %>%
    dplyr::mutate(
      name = str_c(HLA_type, field, sep = "_"),
      name = ifelse(field %in% c("Ave", "NB"), field, name)
    )

  df <- data_read[3:nrow(data_read), 4:ncol(data_read)]
  colnames(df) <- col$name
  df <- cbind(ID, df)

  HLA_num <- unique(na.omit(col$HLA_type))

  df <- purrr::map_dfr(
    HLA_num,
    function(hla_type) {
      df %>%
        dplyr::select(Pos, Peptide, ID, matches(hla_type)) %>%
        dplyr::rename_with(~ str_replace(.x, hla_type, "")) %>%
        dplyr::mutate(type = hla_type, .after = 3)
    }
  )

  df$`_EL_Rank` <- as.numeric(df$`_EL_Rank`)

  df %>%
    dplyr::filter(`_EL_Rank` < 2)
}

input_files <- list.files(
  netmhcpan_result_dir,
  pattern = "\\.txt$",
  full.names = TRUE
)

for (file in input_files) {

  filtered_result <- parse_netmhcpan_result(file)

  original_file_name <- tools::file_path_sans_ext(
    basename(file)
  )

  write.table(
    filtered_result,
    file.path(
      filter_result_dir,
      paste0("f_", original_file_name, ".txt")
    ),
    row.names = FALSE,
    quote = FALSE,
    sep = "\t"
  )
}


# ============================================================
# 9. Combine filtered netMHCpan results
# ============================================================

filtered_files <- list.files(
  filter_result_dir,
  pattern = "\\.txt$",
  full.names = TRUE
)

comb_f_HLA_1 <- dplyr::bind_rows(
  lapply(filtered_files, read.delim, header = TRUE)
)

write.table(
  comb_f_HLA_1,
  "./circRNA_peptide3/part4_netMHCpan/MHC_1/combined_f.txt",
  quote = FALSE,
  row.names = FALSE,
  sep = "\t"
)

utest <- unique(comb_f_HLA_1)

u2test <- utest[, -1, drop = FALSE] %>%
  unique()

u2test <- u2test[, c(1, 2, 3, 7), drop = FALSE]

u2test <- na.omit(u2test)
u2test$X_EL_Rank <- as.numeric(as.character(u2test$X_EL_Rank))
u2test <- na.omit(u2test)

write.table(
  u2test,
  file = "./circRNA_peptide3/part4_netMHCpan/MHC_1/u2test.txt",
  quote = FALSE,
  row.names = FALSE
)


# ============================================================
# 10. Select top netMHCpan results for the heatmap
# ============================================================

heatmap_data <- u2test %>%
  dplyr::arrange(X_EL_Rank) %>%
  dplyr::slice_head(n = 150)

write.table(
  heatmap_data,
  "./circRNA_peptide3/part4_netMHCpan/MHC_1/heatmap_data/heatmap_data.txt"
)


# ============================================================
# 11. Plot the MHC-I binding heatmap
# ============================================================

heatmap_data$Score <- -log10(
  heatmap_data$X_EL_Rank
)

heatmap_matrix <- heatmap_data %>%
  dplyr::select(Peptide, type, Score) %>%
  unique() %>%
  tidyr::pivot_wider(
    names_from = type,
    values_from = Score
  ) %>%
  tibble::column_to_rownames(var = "Peptide")

heatmap_matrix[is.na(heatmap_matrix)] <- 0
heatmap_matrix <- as.matrix(heatmap_matrix)

# Sort HLA alleles by locus and allele number
colnames_vec <- colnames(heatmap_matrix)

hla_type <- gsub(
  "^HLA-([A-Z])[0-9]+:.*",
  "\\1",
  colnames_vec
)

num1 <- as.numeric(
  gsub(
    "^HLA-[A-Z]([0-9]+):.*",
    "\\1",
    colnames_vec
  )
)

num2 <- as.numeric(
  gsub(
    "^HLA-[A-Z][0-9]+:([0-9]+)$",
    "\\1",
    colnames_vec
  )
)

num2[is.na(num2)] <- 0

order_idx <- order(
  factor(hla_type, levels = c("A", "B", "C")),
  num1,
  num2
)

heatmap_matrix <- heatmap_matrix[, order_idx, drop = FALSE]
manual_order <- colnames(heatmap_matrix)

group <- factor(
  substr(colnames(heatmap_matrix), 1, 5),
  levels = c("HLA-A", "HLA-B", "HLA-C")
)

top_annotation <- HeatmapAnnotation(
  cluster = anno_block(
    gp = gpar(
      fill = c("#598F91", "#c0e2c7", "#eedeab"),
      col = NA
    ),
    labels = c("HLA-A", "HLA-B", "HLA-C"),
    labels_gp = gpar(
      col = "black",
      fontsize = 12
    )
  )
)

pdf(
  "./circRNA_peptide3/part4_netMHCpan/picture/MHC1.pdf",
  width = 18,
  height = 15
)

p2 <- Heatmap(
  heatmap_matrix,
  col = colorRampPalette(
    c("white", "#ede2cc", "#edae93", "#dd6670")
  )(100),
  show_row_names = TRUE,
  row_names_side = "left",
  show_column_names = TRUE,
  column_order = manual_order,
  column_names_rot = 45,
  column_title = NULL,
  top_annotation = top_annotation,
  column_split = group,
  rect_gp = gpar(col = "black", lwd = 0.5),
  width = ncol(heatmap_matrix) * unit(10, "mm"),
  name = "Protein Group\n-log10(EL_Rank)\nScore"
)

print(p2)

dev.off()
