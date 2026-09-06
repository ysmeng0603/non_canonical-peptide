# ============================================================
# lncRNA-derived peptide analysis for HLA class I
# ============================================================

suppressPackageStartupMessages({
  library(Biostrings)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(purrr)
  library(tibble)
  library(ggplot2)
})


# ============================================================
# 1. Prepare lncRNA ORF sequences and short-ID mapping
# ============================================================

orf_dir <- "./lncRNA_known/ORFfinder/fasta/p_result"

HC_eRA_result <- readAAStringSet(
  file.path(orf_dir, "HC_eRA_result.fasta")
)

HC_RA_result <- readAAStringSet(
  file.path(orf_dir, "HC_RA_result.fasta")
)

RA_eRA_result <- readAAStringSet(
  file.path(orf_dir, "RA_eRA_result.fasta")
)

# Merge ORFs and deduplicate by sequence header
all_seqs <- c(
  HC_eRA_result,
  HC_RA_result,
  RA_eRA_result
)

all_seqs <- all_seqs[!duplicated(names(all_seqs))]

# Build ORF ID table
all_need <- data.frame(
  ORF_ID = c(
    names(HC_eRA_result),
    names(HC_RA_result),
    names(RA_eRA_result)
  ),
  stringsAsFactors = FALSE
) %>%
  distinct(ORF_ID, .keep_all = TRUE) %>%
  mutate(
    new_ID = paste0("LncID_", row_number()),
    ENST_id = str_extract(ORF_ID, "(?<=_)[^:]+")
  )

write.table(
  all_need %>% select(ORF_ID, new_ID),
  "./lncRNA_known/netMHCpan/ori_shortID_map.txt",
  quote = FALSE,
  row.names = FALSE,
  sep = "\t"
)


# ============================================================
# 2. Identify expressed upregulated lncRNAs in each sample
# ============================================================

lncRNA_transcript_counts <- read.delim(
  "./lncRNA_known/netMHCpan/lncRNA_transcript_counts.txt",
  comment.char = "#"
)

all_lncRNA_group <- read.delim(
  "./lncRNA_known/ORFfinder/all_lncRNA_group.txt"
)

need_up <- all_lncRNA_group %>%
  filter(significance == "Up")

expressed_lncRNA <- lncRNA_transcript_counts %>%
  select(Geneid, matches("X.mnt.")) %>%
  pivot_longer(
    cols = -Geneid,
    names_to = "name",
    values_to = "value"
  ) %>%
  filter(value > 0) %>%
  filter(
    Geneid %in% all_need$ENST_id,
    Geneid %in% need_up$lncRNA_id
  ) %>%
  mutate(
    SRR_id = str_extract(name, "SRR\\d+")
  )


# ============================================================
# 3. Generate sample-specific lncRNA ORF FASTA files
# ============================================================

orf2new <- all_need %>%
  distinct(ORF_ID, new_ID)

idx <- match(names(all_seqs), orf2new$ORF_ID)
keep <- !is.na(idx)

all_seqs <- all_seqs[keep]
names(all_seqs) <- orf2new$new_ID[idx[keep]]

expr_by_sample <- expressed_lncRNA %>%
  select(SRR_id, ENST_id = Geneid) %>%
  distinct()

orf_by_sample <- expr_by_sample %>%
  inner_join(
    all_need %>% select(ENST_id, new_ID),
    by = "ENST_id"
  ) %>%
  distinct(SRR_id, new_ID)

sample_fasta_dir <- "./lncRNA_known/netMHCpan/HLA_1/SRR_Files_only_up"

dir.create(
  sample_fasta_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

for (srr in unique(orf_by_sample$SRR_id)) {

  ids <- orf_by_sample %>%
    filter(SRR_id == srr) %>%
    pull(new_ID) %>%
    unique()

  ids <- intersect(ids, names(all_seqs))

  fa_sub <- all_seqs[ids]

  writeXStringSet(
    fa_sub,
    filepath = file.path(sample_fasta_dir, paste0(srr, ".fasta")),
    format = "fasta"
  )
}


# ============================================================
# 4. Prepare HLA-I information for netMHCpan
# ============================================================

HLA_1 <- read.delim2(
  "./circRNA_peptide2/part4_netMHCpan/MHC_1/HLA_1.txt",
  header = FALSE
)

colnames(HLA_1) <- c("file_name", "HLA_1")

file_info <- data.frame(
  file_path = list.files(
    sample_fasta_dir,
    pattern = "\\.fasta$",
    full.names = TRUE
  ),
  stringsAsFactors = FALSE
) %>%
  mutate(
    file_name = str_extract(file_path, "SRR\\d+")
  )

final_for_netMHCpan_1 <- file_info %>%
  left_join(HLA_1, by = "file_name") %>%
  tidyr::drop_na() %>%
  mutate(index = row_number()) %>%
  select(index, file_path, file_name, HLA_1)

write.table(
  final_for_netMHCpan_1,
  "./lncRNA_known/netMHCpan/HLA_1/bash_need/MHC1_file_paths_only_up.txt",
  quote = FALSE,
  row.names = FALSE,
  sep = "\t"
)


# ============================================================
# 5. Parse and filter netMHCpan HLA-I results
# ============================================================

netmhcpan_result_dir <- "./lncRNA_known/netMHCpan/HLA_1/netMHCpan_1_result"
filter_result_dir <- "./lncRNA_known/netMHCpan/HLA_1/filter_result"

dir.create(
  filter_result_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

input_files <- list.files(
  netmhcpan_result_dir,
  pattern = "\\.txt$",
  full.names = TRUE
)

parse_netmhcpan_mhc1 <- function(file) {

  data_read <- read.table(
    file,
    sep = "\t",
    header = FALSE,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )

  ID <- data_read[3:nrow(data_read), 1:3, drop = FALSE]
  colnames(ID) <- c("Pos", "Peptide", "ID")

  col <- data_read[1:2, 4:ncol(data_read), drop = FALSE]
  rownames(col) <- c("HLA_type", "field")

  col <- as.data.frame(t(col), stringsAsFactors = FALSE)

  col$HLA_type <- ifelse(
    str_length(col$HLA_type) == 0,
    NA,
    col$HLA_type
  )

  col <- col %>%
    fill(HLA_type, .direction = "down") %>%
    mutate(
      name = str_c(HLA_type, field, sep = "_")
    ) %>%
    mutate(
      name = ifelse(field %in% c("Ave", "NB"), field, name)
    )

  df <- data_read[3:nrow(data_read), 4:ncol(data_read), drop = FALSE]
  colnames(df) <- col$name
  df <- cbind(ID, df)

  HLA_num <- unique(stats::na.omit(col$HLA_type))

  df <- map_dfr(
    HLA_num,
    function(hla_type) {
      df %>%
        select(Pos, Peptide, ID, matches(hla_type)) %>%
        rename_with(
          ~ str_replace(.x, fixed(hla_type), "")
        ) %>%
        mutate(type = hla_type, .after = 3)
    }
  )

  df$`_EL_Rank` <- as.numeric(df$`_EL_Rank`)

  df %>%
    filter(`_EL_Rank` < 2)
}

for (file in input_files) {

  filtered <- parse_netmhcpan_mhc1(file)

  original_file_name <- tools::file_path_sans_ext(
    basename(file)
  )

  new_file_name <- paste0(
    "f_",
    original_file_name,
    ".txt"
  )

  write.table(
    filtered,
    file.path(filter_result_dir, new_file_name),
    row.names = FALSE,
    quote = FALSE,
    sep = "\t"
  )
}


# ============================================================
# 6. Combine filtered HLA-I results
# ============================================================

filtered_files <- list.files(
  filter_result_dir,
  pattern = "\\.txt$",
  full.names = TRUE
)

combined_f_1 <- bind_rows(
  lapply(
    filtered_files,
    read.delim,
    header = TRUE,
    stringsAsFactors = FALSE
  )
)

write.table(
  combined_f_1,
  "./lncRNA_known/netMHCpan/HLA_1/combined_f.txt",
  quote = FALSE,
  row.names = FALSE,
  sep = "\t"
)


# ============================================================
# 7. Summarize HLA-I binding ranks
# ============================================================

df_all <- combined_f_1 %>%
  transmute(
    Peptide,
    ID,
    type,
    Rank = as.numeric(X_EL_Rank),
    source = "HLA_1"
  ) %>%
  mutate(
    locus = case_when(
      grepl("^HLA-A", type) ~ "HLA-A",
      grepl("^HLA-B", type) ~ "HLA-B",
      grepl("^HLA-C", type) ~ "HLA-C",
      TRUE ~ "Other"
    )
  ) %>%
  filter(locus != "Other")

df_med <- df_all %>%
  group_by(source, locus, Peptide) %>%
  summarise(
    Rank_median = stats::median(Rank, na.rm = TRUE),
    n = n(),
    .groups = "drop"
  ) %>%
  mutate(
    locus = factor(
      locus,
      levels = c("HLA-A", "HLA-B", "HLA-C")
    )
  )

lab_df <- df_med %>%
  group_by(locus) %>%
  summarise(
    n_peptide = n(),
    med = stats::median(Rank_median, na.rm = TRUE),
    .groups = "drop"
  )

write.table(
  df_med,
  "./lncRNA_known/netMHCpan/boxplot.txt",
  quote = FALSE,
  row.names = FALSE,
  sep = "\t"
)


# ============================================================
# 8. Plot HLA-I binding-rank distribution
# ============================================================

my_cols <- c(
  "HLA-A" = "#e3bec6",
  "HLA-B" = "#538085",
  "HLA-C" = "#faf1e2"
)

pdf(
  "./lncRNA_known/netMHCpan/boxplot.pdf",
  width = 10,
  height = 8
)

p_box <- ggplot(
  df_med,
  aes(
    x = locus,
    y = Rank_median,
    fill = locus
  )
) +
  geom_boxplot(
    alpha = 0.85,
    outlier.shape = NA,
    width = 0.65,
    linewidth = 0.6,
    colour = "black"
  ) +
  geom_text(
    data = lab_df,
    aes(
      x = locus,
      y = med,
      label = sprintf("%.2f", med)
    ),
    inherit.aes = FALSE,
    vjust = -3.5,
    size = 8
  ) +
  geom_text(
    data = lab_df,
    aes(
      x = locus,
      y = Inf,
      label = paste0("n=", n_peptide)
    ),
    inherit.aes = FALSE,
    hjust = 0.7,
    vjust = 0.65,
    size = 7
  ) +
  scale_fill_manual(values = my_cols) +
  coord_flip(clip = "off") +
  scale_y_continuous(
    expand = expansion(mult = c(0.02, 0.12))
  ) +
  theme_classic(base_size = 14) +
  theme(
    legend.position = "none",
    axis.title.y = element_blank(),
    axis.title.x = element_text(size = 24),
    axis.text.y = element_text(size = 22, colour = "black"),
    axis.text.x = element_text(size = 22, colour = "black"),
    plot.margin = margin(10, 50, 10, 10)
  ) +
  ylab("Median Rank(%)")

print(p_box)
dev.off()


# ============================================================
# 9. Calculate peptide prevalence across HLA-I samples
# ============================================================

File_paths <- list.files(
  filter_result_dir,
  pattern = "\\.txt$",
  full.names = TRUE
)

df_list <- lapply(
  File_paths,
  function(file_path) {

    df <- read.delim(
      file_path,
      header = TRUE,
      stringsAsFactors = FALSE
    )

    file_name <- tools::file_path_sans_ext(
      basename(file_path)
    )

    df %>%
      select(Peptide, type) %>%
      mutate(
        type = substr(type, 1, 5),
        SRR = file_name
      ) %>%
      unique()
  }
)

comb_MHC1 <- bind_rows(df_list)
total_SRR <- length(df_list)

MHC1_count <- comb_MHC1 %>%
  group_by(type, Peptide) %>%
  summarise(
    count = n(),
    .groups = "drop"
  ) %>%
  mutate(
    percentage = count / total_SRR * 100,
    percntg_range = cut(
      percentage,
      breaks = c(0, 10, 20, 30, 100),
      labels = c("0-10%", "10-20%", "20-30%", ">30%")
    )
  )

write.table(
  MHC1_count,
  "./lncRNA_known/netMHCpan/plot/sample_barplot/need_file/MHC1_for_pic_1.txt",
  quote = FALSE,
  row.names = FALSE,
  sep = "\t"
)

final_MHC1 <- MHC1_count %>%
  group_by(type, percntg_range) %>%
  summarise(
    counts = n(),
    .groups = "drop_last"
  ) %>%
  mutate(
    range_percntg = counts / sum(counts) * 100
  ) %>%
  ungroup()

write.table(
  final_MHC1,
  "./lncRNA_known/netMHCpan/plot/sample_barplot/need_file/final_MHC1.txt",
  quote = FALSE,
  row.names = FALSE,
  sep = "\t"
)


# ============================================================
# 10. Plot HLA-I peptide prevalence
# ============================================================

final_MHC1 <- final_MHC1 %>%
  mutate(HLA_type = "HLA-I")

colnames(final_MHC1)[1] <- "HLA"
colnames(final_MHC1)[2] <- "Range"
colnames(final_MHC1)[4] <- "Percentage"

Colors <- c(
  "0-10%" = "#F3BDA5",
  "10-20%" = "#B7E4C7",
  "20-30%" = "#FDFD96",
  ">30%" = "#D8B4E2"
)

draw_barplot <- function(need_data, need_group) {

  ggplot(
    need_data,
    aes(
      x = HLA,
      y = counts,
      fill = Range
    )
  ) +
    geom_bar(
      stat = "identity",
      position = "stack",
      width = 0.6
    ) +
    scale_fill_manual(values = Colors) +
    labs(
      title = need_group,
      x = ""
    ) +
    scale_y_continuous(
      labels = scales::comma
    ) +
    theme_classic() +
    theme(
      panel.border = element_rect(
        color = "black",
        fill = NA,
        linewidth = 1
      ),
      strip.background = element_blank(),
      plot.title = element_text(
        hjust = 0.5,
        size = 22
      ),
      legend.position = "right",
      legend.title = element_text(size = 16),
      legend.text = element_text(size = 13),
      axis.text = element_text(size = 15),
      axis.title.y = element_text(size = 18),
      axis.line = element_blank(),
      axis.ticks = element_line(color = "black")
    )
}

HLA_1_plot <- draw_barplot(
  final_MHC1,
  "HLA-I"
)

pdf(
  "./lncRNA_known/netMHCpan/plot/sample_barplot/plot/new_MHC_1_barplot.pdf",
  width = 6,
  height = 6
)

print(HLA_1_plot)
dev.off()


# ============================================================
# 11. Identify MS-supported HLA-I peptides
# ============================================================

containHuman <- read.delim(
  "./lncRNA_known/MaxQuant/result/z_filter_result/containHuman.txt"
)

combined_f <- read.delim(
  "./lncRNA_known/netMHCpan/HLA_1/combined_f.txt"
)

ori_shortID_map <- read.delim(
  "./lncRNA_known/netMHCpan/ori_shortID_map.txt"
)

map2 <- ori_shortID_map %>%
  mutate(
    orf_key = str_extract(
      ORF_ID,
      "ORF\\d+_ENST\\d+\\.\\d+"
    )
  ) %>%
  select(orf_key, new_ID) %>%
  distinct()

all_MS <- containHuman %>%
  mutate(
    new = str_extract(
      Proteins,
      "(?<=lcl\\|)[^:]+"
    ),
    new = str_c(
      new,
      Sequence,
      sep = "&"
    )
  ) %>%
  select(new, everything())

all_MS <- all_MS %>%
  mutate(
    orf_key = sub("&.*$", "", new),
    peptide = sub("^.*&", "", new)
  ) %>%
  left_join(
    map2,
    by = "orf_key"
  ) %>%
  mutate(
    new = ifelse(
      !is.na(new_ID),
      paste0(new_ID, "&", peptide),
      new
    )
  ) %>%
  select(-orf_key, -peptide)

for_plot <- combined_f %>%
  mutate(
    new = str_c(
      ID,
      Peptide,
      sep = "&"
    )
  ) %>%
  select(new, everything()) %>%
  filter(new %in% all_MS$new)

write.table(
  for_plot,
  "./lncRNA_known/netMHCpan/plot/MS_heatmap/need_file/for_plot_1.txt",
  quote = FALSE,
  row.names = FALSE,
  sep = "\t"
)


# ============================================================
# 12. Plot MS-supported HLA-I peptide binding heatmap
# ============================================================

suppressPackageStartupMessages({
  library(ComplexHeatmap)
  library(circlize)
  library(grid)
})

z_u_MS_MHC <- for_plot %>%
  select(
    Peptide,
    type,
    X_EL_Rank
  ) %>%
  mutate(
    HLA = str_sub(type, 1, 5)
  ) %>%
  unique()

plot_file <- z_u_MS_MHC %>%
  select(
    Peptide,
    type,
    X_EL_Rank
  ) %>%
  mutate(
    log = -log10(X_EL_Rank + 1e-4)
  )

heatmap_df <- plot_file %>%
  select(Peptide, type, log) %>%
  unique() %>%
  arrange(desc(log)) %>%
  slice_head(n = 300) %>%
  pivot_wider(
    names_from = type,
    values_from = log
  ) %>%
  column_to_rownames("Peptide")

heatmap_df[is.na(heatmap_df)] <- 0
heatmap_matrix <- as.matrix(heatmap_df)

# Sort HLA-I alleles as HLA-A, HLA-B, and HLA-C
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
  factor(
    hla_type,
    levels = c("A", "B", "C")
  ),
  num1,
  num2
)

heatmap_matrix <- heatmap_matrix[, order_idx, drop = FALSE]

hla_group <- str_extract(
  colnames(heatmap_matrix),
  "^HLA-[ABC]"
)

hla_group <- factor(
  hla_group,
  levels = c("HLA-A", "HLA-B", "HLA-C")
)

pdf(
  "./lncRNA_known/netMHCpan/plot/MS_heatmap/new_MHC1.pdf",
  width = 8,
  height = 6
)

p_heatmap <- Heatmap(
  heatmap_matrix,
  col = colorRampPalette(
    c(
      "white",
      "#ede2cc",
      "#edae93",
      "#dd6670"
    )
  )(100),
  show_row_names = TRUE,
  row_names_side = "right",
  show_column_names = TRUE,
  column_names_rot = 45,
  column_title = NULL,
  column_split = hla_group,
  rect_gp = gpar(
    col = "black",
    lwd = 0.5
  ),
  width = ncol(heatmap_matrix) * unit(8, "mm"),
  name = "HLA-I\nMS-supported peptides"
)

print(p_heatmap)
dev.off()


# ============================================================
# 13. Prepare HLA-I peptide-frequency table for word cloud
# ============================================================

peptide_Cloud <- combined_f %>%
  mutate(
    new = str_c(
      ID,
      Peptide,
      sep = "&"
    )
  ) %>%
  select(new, everything()) %>%
  filter(new %in% all_MS$new)

table_count <- peptide_Cloud %>%
  group_by(Peptide) %>%
  summarise(
    Num = n(),
    .groups = "drop"
  )

write.table(
  peptide_Cloud,
  "./lncRNA_known/netMHCpan/plot/wordCloud/need_file/peptide_Cloud_HLA_1.txt",
  quote = FALSE,
  row.names = FALSE,
  sep = "\t"
)

write.table(
  table_count,
  "./lncRNA_known/netMHCpan/plot/wordCloud/need_file/table_count_HLA_1.txt",
  quote = FALSE,
  row.names = FALSE,
  sep = "\t"
)


# ============================================================
# 14. Prepare MS-supported ORFs for NetChop
# ============================================================

peptide_Cloud_HLA_1 <- read.delim(
  "./lncRNA_known/netMHCpan/plot/wordCloud/need_file/peptide_Cloud_HLA_1.txt"
)

fa1 <- readAAStringSet(
  "./lncRNA_known/MaxQuant/need_file/only_lncRNA/HC_eRA_result.fasta"
)

fa2 <- readAAStringSet(
  "./lncRNA_known/MaxQuant/need_file/only_lncRNA/HC_RA_result.fasta"
)

fa3 <- readAAStringSet(
  "./lncRNA_known/MaxQuant/need_file/only_lncRNA/RA_eRA_result.fasta"
)

all_fa <- c(
  fa1,
  fa2,
  fa3
)

all_fa_dedup <- all_fa[!duplicated(names(all_fa))]

names(all_fa_dedup) <- sub(
  " .*",
  "",
  names(all_fa_dedup)
)

names(all_fa_dedup) <- str_extract(
  names(all_fa_dedup),
  "(?<=lcl\\|).+?(?=:\\d+)"
)

need_ORF_fasta <- peptide_Cloud_HLA_1 %>%
  select(ID) %>%
  distinct()

map2_netchop <- ori_shortID_map %>%
  mutate(
    orf_key = str_extract(
      ORF_ID,
      "ORF\\d+_ENST\\d+\\.\\d+"
    )
  ) %>%
  select(orf_key, new_ID) %>%
  distinct() %>%
  transmute(
    orf_key,
    ID = new_ID
  )

need_ORF_fasta <- need_ORF_fasta %>%
  left_join(
    map2_netchop,
    by = "ID"
  )

need_ORF_fasta <- all_fa_dedup[
  need_ORF_fasta$orf_key
]

# Deduplicate by amino-acid sequence
seqs <- as.character(need_ORF_fasta)
keep <- !duplicated(seqs)
need_ORF_fasta_uniq <- need_ORF_fasta[keep]

old_names <- names(need_ORF_fasta_uniq)
new_names <- paste0(
  "nCp_",
  seq_along(old_names)
)

names(need_ORF_fasta_uniq) <- new_names

id_map <- data.frame(
  short_id = new_names,
  original_id = old_names,
  stringsAsFactors = FALSE
)

write.csv(
  id_map,
  "./lncRNA_known/netChop/need_file/NetChop_ID_mapping.csv",
  row.names = FALSE
)

writeXStringSet(
  need_ORF_fasta_uniq,
  "./lncRNA_known/netChop/need_file/netChop_input.fasta"
)


# ============================================================
# 15. Plot HLA-I binding-rank ridge plot
# ============================================================

suppressPackageStartupMessages({
  library(ggridges)
})

u2test <- combined_f_1 %>%
  tidyr::drop_na() %>%
  mutate(
    X_EL_Rank = as.numeric(X_EL_Rank),
    Immunogenicity_Score = X_EL_Rank,
    HLA = str_sub(type, 1, 5)
  )

u2test_summary <- u2test %>%
  group_by(HLA, Peptide) %>%
  summarise(
    Upper_Quartile_Score = quantile(
      Immunogenicity_Score,
      0.75
    ),
    .groups = "drop"
  )

plot_file <- u2test_summary

sampleColor <- c(
  '#F6EB5E', '#CD8280',
  '#F0F0F0', '#54D0B4',
  '#D2D1D1', '#78C3ED',
  '#E69F00', '#CBDFF3',
  '#3BFFB8', '#ECAD67',
  '#B3B2B2', '#959595',
  '#E98FBD', '#2672B2'
)

peak_data <- plot_file %>%
  group_by(HLA) %>%
  reframe(
    density = density(Upper_Quartile_Score)$x,
    density_y = density(Upper_Quartile_Score)$y
  ) %>%
  group_by(HLA) %>%
  slice(which.max(density_y)) %>%
  rename(peak = density) %>%
  select(HLA, peak)

pdf(
  "./lncRNA_known/netMHCpan/ridge_1.pdf",
  width = 6.88,
  height = 4.86
)

RidgePlot <- ggplot(
  plot_file,
  aes(
    x = Upper_Quartile_Score,
    y = HLA
  )
) +
  ggridges::geom_density_ridges(
    aes(fill = HLA),
    alpha = 0.8,
    show.legend = FALSE,
    scale = 1
  ) +
  scale_fill_manual(
    values = sampleColor[c(1, 6, 10)]
  ) +
  geom_point(
    data = peak_data,
    aes(
      x = peak,
      y = HLA
    ),
    shape = 21,
    size = 3,
    fill = "white",
    color = "black"
  ) +
  geom_text(
    data = peak_data,
    aes(
      x = peak,
      y = HLA,
      label = round(peak, 2)
    ),
    vjust = -1,
    size = 4,
    color = "black",
    fontface = "bold"
  ) +
  geom_vline(
    xintercept = 0.5,
    linetype = "dashed"
  ) +
  labs(
    x = "",
    y = ""
  ) +
  ggridges::theme_ridges(grid = FALSE) +
  theme(
    legend.position = "none",
    axis.text.y = element_text(size = 14),
    axis.ticks.x = element_blank(),
    axis.text.x = element_blank(),
    plot.margin = unit(c(0.5, 0, 0, 0), "cm")
  )

print(RidgePlot)
dev.off()


# ============================================================
# 16. Compare MS-supported HLA-I peptides with IEDB
# ============================================================

suppressPackageStartupMessages({
  library(readxl)
})

IEDB_select <- read.csv(
  "./circRNA/a_info/IEDB_select.csv"
)

need_peptide <- data.frame(
  Peptide = unique(for_plot$Peptide),
  stringsAsFactors = FALSE
) %>%
  mutate(
    IS_IEDB = Peptide %in% IEDB_select$Epitopes...Epitope
  )

need_2 <- need_peptide %>%
  filter(!IS_IEDB)

peptide <- paste(
  need_2$Peptide,
  collapse = ","
)

writeLines(
  peptide,
  "./A_paper/lncRNA_not_IEDB.txt"
)

lncRNA_Not_IEDB <- read_excel(
  "./A_paper/lncRNA_Not_IEDB.xlsx"
)

table_copy <- lncRNA_Not_IEDB %>%
  rename(
    Entry = `Entry`,
    Match = `Match`,
    Entry_Name = `Entry Name`,
    Protein_Names = `Protein Names`,
    Gene_Names = `Gene Names`,
    Organism = `Organism`,
    Length = `Length`
  ) %>%
  fill(
    Entry,
    Entry_Name,
    Protein_Names,
    Gene_Names,
    Organism,
    Length,
    .direction = "down"
  ) %>%
  filter(
    !is.na(Match),
    Match != "NA",
    str_detect(Match, ":")
  ) %>%
  mutate(
    Peptide = str_trim(
      str_replace(
        Match,
        ".*:\\s*",
        ""
      )
    ),
    Match_Type = if_else(
      str_detect(Entry, "-"),
      "isoform",
      "canonical"
    )
  )

peptide_detail <- table_copy %>%
  select(
    Peptide,
    Entry,
    Entry_Name,
    Protein_Names,
    Gene_Names,
    Organism,
    Match_Type
  ) %>%
  distinct()

peptide_class <- peptide_detail %>%
  group_by(Peptide) %>%
  summarise(
    has_canonical = any(Match_Type == "canonical"),
    has_isoform = any(Match_Type == "isoform"),
    Final_Class = case_when(
      has_canonical & !has_isoform ~ "canonical",
      !has_canonical & has_isoform ~ "isoform",
      has_canonical & has_isoform ~ "both",
      TRUE ~ "unknown"
    ),
    Matched_Entries = paste(
      unique(Entry),
      collapse = "; "
    ),
    Matched_Genes = paste(
      unique(Gene_Names),
      collapse = "; "
    ),
    n_matches = n(),
    .groups = "drop"
  ) %>%
  select(1:6)

write.table(
  peptide_class,
  "./A_paper/lncRNA_not_IEDB.csv",
  quote = FALSE,
  row.names = FALSE,
  sep = ","
)
