# ============================================================
# Alternative splicing peptide workflow
# rMATS -> JCAST -> sample-specific FASTA -> netMHCpan
# ============================================================

suppressPackageStartupMessages({
  library(Biostrings)
  library(data.table)
  library(dplyr)
  library(ggplot2)
  library(stringr)
  library(tibble)
  library(tidyr)
})


# ============================================================
# 1. Build sample-level support tables from rMATS results
# ============================================================

read_rmats_core <- function(path, type) {
  x <- data.table::fread(path)

  selected_cols <- c(
    "ID", "GeneID", "geneSymbol", "chr", "strand",
    "IJC_SAMPLE_1", "SJC_SAMPLE_1", "IJC_SAMPLE_2", "SJC_SAMPLE_2",
    "FDR", "IncLevel1", "IncLevel2", "IncLevelDifference"
  )

  selected_cols <- intersect(selected_cols, names(x))
  x <- x[, ..selected_cols]

  if ("ID" %in% names(x)) {
    data.table::setnames(x, "ID", "event_id")
  } else {
    data.table::setnames(x, names(x)[1], "event_id")
  }

  x[, event_type := type]
  x %>% dplyr::select(event_type, dplyr::everything())
}


build_long_support_one <- function(
    rmats_dir,
    group1_name,
    group2_name,
    bam_SRR_group,
    out_dir,
    fdr_cut = 0.05,
    dpsi_cut = -0.10,
    present_cut = 0
) {
  stopifnot(dir.exists(rmats_dir))
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  files <- c(
    SE   = file.path(rmats_dir, "SE.MATS.JC.txt"),
    MXE  = file.path(rmats_dir, "MXE.MATS.JC.txt"),
    A3SS = file.path(rmats_dir, "A3SS.MATS.JC.txt"),
    A5SS = file.path(rmats_dir, "A5SS.MATS.JC.txt"),
    RI   = file.path(rmats_dir, "RI.MATS.JC.txt")
  )

  if (!all(file.exists(files))) {
    miss <- names(files)[!file.exists(files)]
    stop("Missing files in ", rmats_dir, ": ", paste(miss, collapse = ", "))
  }

  rmats_list <- mapply(
    read_rmats_core,
    path = files,
    type = names(files),
    SIMPLIFY = FALSE
  )

  rmats_all <- data.table::rbindlist(
    rmats_list,
    use.names = TRUE,
    fill = TRUE
  )

  filter_data <- rmats_all %>%
    dplyr::mutate(
      FDR = as.numeric(FDR),
      IncLevelDifference = as.numeric(IncLevelDifference)
    ) %>%
    dplyr::filter(
      !is.na(FDR),
      FDR < fdr_cut,
      !is.na(IncLevelDifference),
      IncLevelDifference < dpsi_cut
    )

  if (nrow(filter_data) == 0) {
    warning("No events left after filtering in: ", rmats_dir)
    return(invisible(NULL))
  }

  # Sample order must match the rMATS --b1/--b2 BAM lists.
  all_SRR <- list(
    healthy = bam_SRR_group$Sample[bam_SRR_group$Group == "healthy"],
    eRA     = bam_SRR_group$Sample[bam_SRR_group$Group == "eRA"],
    RA      = bam_SRR_group$Sample[bam_SRR_group$Group == "RA"]
  )

  if (!group1_name %in% names(all_SRR)) {
    stop("group1_name not found in all_SRR: ", group1_name)
  }

  if (!group2_name %in% names(all_SRR)) {
    stop("group2_name not found in all_SRR: ", group2_name)
  }

  suppose_num_1 <- length(all_SRR[[group1_name]])
  suppose_num_2 <- length(all_SRR[[group2_name]])

  max_n1 <- max(
    stringr::str_count(filter_data$IJC_SAMPLE_1, ",") + 1,
    na.rm = TRUE
  )

  max_n2 <- max(
    stringr::str_count(filter_data$IJC_SAMPLE_2, ",") + 1,
    na.rm = TRUE
  )

  if (max_n1 != suppose_num_1 || max_n2 != suppose_num_2) {
    stop(
      "Sample number mismatch: ", basename(rmats_dir),
      " rMATS(group1=", max_n1, ", group2=", max_n2,
      ") vs bam_SRR_group(", group1_name, "=", suppose_num_1,
      ", ", group2_name, "=", suppose_num_2, ")"
    )
  }

  if (any(
    stringr::str_count(filter_data$IJC_SAMPLE_1, ",") + 1 != suppose_num_1,
    na.rm = TRUE
  )) {
    warning(
      basename(rmats_dir),
      ": replicate number in IJC_SAMPLE_1 differs from ", suppose_num_1
    )
  }

  if (any(
    stringr::str_count(filter_data$IJC_SAMPLE_2, ",") + 1 != suppose_num_2,
    na.rm = TRUE
  )) {
    warning(
      basename(rmats_dir),
      ": replicate number in IJC_SAMPLE_2 differs from ", suppose_num_2
    )
  }

  filter_data <- filter_data %>%
    tidyr::separate(
      col = IJC_SAMPLE_1,
      into = paste0("IJC_1_", all_SRR[[group1_name]]),
      sep = ",",
      convert = TRUE,
      fill = "right",
      remove = FALSE
    ) %>%
    tidyr::separate(
      col = SJC_SAMPLE_1,
      into = paste0("SJC_1_", all_SRR[[group1_name]]),
      sep = ",",
      convert = TRUE,
      fill = "right",
      remove = FALSE
    ) %>%
    tidyr::separate(
      col = IJC_SAMPLE_2,
      into = paste0("IJC_2_", all_SRR[[group2_name]]),
      sep = ",",
      convert = TRUE,
      fill = "right",
      remove = FALSE
    ) %>%
    tidyr::separate(
      col = SJC_SAMPLE_2,
      into = paste0("SJC_2_", all_SRR[[group2_name]]),
      sep = ",",
      convert = TRUE,
      fill = "right",
      remove = FALSE
    )

  filter_data2 <- filter_data %>%
    dplyr::mutate(
      dplyr::across(
        dplyr::matches("^(IJC|SJC)_[12]_SRR"),
        ~ dplyr::if_else(is.na(.x), 0L, as.integer(.x))
      )
    )

  long_support <- filter_data2 %>%
    dplyr::select(
      event_id,
      event_type,
      geneSymbol,
      GeneID,
      chr,
      strand,
      FDR,
      IncLevelDifference,
      dplyr::matches("^IJC_[12]_SRR"),
      dplyr::matches("^SJC_[12]_SRR")
    ) %>%
    tidyr::pivot_longer(
      cols = dplyr::matches("^(IJC|SJC)_[12]_SRR"),
      names_to = "key",
      values_to = "count"
    ) %>%
    tidyr::separate(
      key,
      into = c("kind", "group", "SRR"),
      sep = "_"
    ) %>%
    tidyr::pivot_wider(
      names_from = kind,
      values_from = count,
      values_fill = 0
    ) %>%
    dplyr::mutate(
      support = IJC + SJC,
      present = support > present_cut
    ) %>%
    dplyr::filter(present) %>%
    dplyr::mutate(
      group_label = dplyr::if_else(
        group == "1",
        group1_name,
        group2_name
      )
    )

  out_file <- file.path(
    out_dir,
    paste0(group1_name, "_", group2_name, "_count.txt")
  )

  write.table(
    long_support,
    out_file,
    quote = FALSE,
    row.names = FALSE,
    sep = "\t"
  )

  message(
    "Done: ", basename(rmats_dir), " -> ", out_file,
    " (events=", dplyr::n_distinct(long_support$event_id),
    ", rows=", nrow(long_support), ")"
  )

  invisible(long_support)
}


bam_SRR_group <- read.delim("~/A_info/bam_SRR_group.txt")

rmats_count_dir <- "./AS_rMATS_jcast/rMATS/Count"

contrasts <- list(
  list(
    dir = "./AS_rMATS_jcast/rMATS/result/hc_vs_eRA",
    g1 = "healthy",
    g2 = "eRA"
  ),
  list(
    dir = "./AS_rMATS_jcast/rMATS/result/hc_vs_RA",
    g1 = "healthy",
    g2 = "RA"
  ),
  list(
    dir = "./AS_rMATS_jcast/rMATS/result/RA_vs_eRA",
    g1 = "RA",
    g2 = "eRA"
  )
)

rmats_support <- lapply(
  contrasts,
  function(x) {
    build_long_support_one(
      rmats_dir = x$dir,
      group1_name = x$g1,
      group2_name = x$g2,
      bam_SRR_group = bam_SRR_group,
      out_dir = rmats_count_dir,
      fdr_cut = 0.05,
      dpsi_cut = -0.10,
      present_cut = 0
    )
  }
)


# ============================================================
# 2. Parse JCAST T1/T2 protein FASTA files
# ============================================================

read_jcast_fasta_biostrings <- function(fa) {
  aa <- Biostrings::readAAStringSet(fa)
  hdr <- names(aa)
  seq <- as.character(aa)

  keypart <- sub("^>", "", hdr)
  keypart <- vapply(
    strsplit(keypart, " ", fixed = TRUE),
    `[`,
    character(1),
    1
  )

  parts <- strsplit(keypart, "\\|")

  get_part <- function(p, i) {
    if (length(p) >= i) p[[i]] else NA_character_
  }

  tibble::tibble(
    header = hdr,
    seq = seq,
    kb = vapply(parts, get_part, character(1), i = 1),
    uniprot = vapply(parts, get_part, character(1), i = 2),
    uniprot_name = vapply(parts, get_part, character(1), i = 3),
    gene_ensg = vapply(parts, get_part, character(1), i = 4),
    type_order = vapply(parts, get_part, character(1), i = 5),
    event_id = as.integer(vapply(parts, get_part, character(1), i = 6)),
    chr = vapply(parts, get_part, character(1), i = 7),
    sjc_raw = vapply(parts, get_part, character(1), i = 11),
    sjc = as.integer(sub(
      "^r",
      "",
      vapply(parts, get_part, character(1), i = 11)
    )),
    tier = vapply(parts, get_part, character(1), i = 12),
    event_type = sub(
      "[0-9]+$",
      "",
      vapply(parts, get_part, character(1), i = 5)
    )
  )
}


jcast_base <- "./AS_rMATS_jcast/jcast/result/T1_and_T2_info"

jcast_groups <- c("HC_eRA", "HC_RA", "RA_eRA")

jcast_info <- lapply(
  jcast_groups,
  function(group) {
    x <- read_jcast_fasta_biostrings(
      file.path(
        jcast_base,
        "fasta",
        paste0(group, "_T1_T2.fasta")
      )
    )

    x %>%
      dplyr::mutate(
        short_ID = paste0(group, "_", dplyr::row_number())
      ) %>%
      dplyr::select(short_ID, dplyr::everything())
  }
)

names(jcast_info) <- jcast_groups

for (group in names(jcast_info)) {
  write.table(
    jcast_info[[group]],
    file.path(jcast_base, paste0(group, "_T1_T2_info.txt")),
    quote = FALSE,
    row.names = FALSE,
    sep = "\t"
  )
}

all_mapping <- dplyr::bind_rows(jcast_info)

write.table(
  all_mapping,
  file.path(jcast_base, "all_shortID_map.txt"),
  quote = FALSE,
  row.names = FALSE,
  sep = "\t"
)


# ============================================================
# 3. Build sample-specific JCAST FASTA files
# ============================================================

write_srr_fastas_one <- function(
    count_path,
    info_path,
    out_dir,
    sjc_cut = 116
) {
  count_df <- read.delim(
    count_path,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )

  info_df <- read.delim(
    info_path,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )

  count_df <- count_df %>%
    dplyr::mutate(
      event_id = as.character(event_id),
      event_key = paste(event_type, event_id, sep = ":")
    )

  info_df <- info_df %>%
    dplyr::mutate(
      event_id = as.character(event_id),
      event_key = paste(event_type, event_id, sep = ":")
    )

  # Retain sample-level events with SJC > 116.
  count_df <- count_df %>%
    dplyr::filter(SJC > sjc_cut)

  selected <- count_df %>%
    dplyr::semi_join(
      info_df %>% dplyr::select(event_key),
      by = "event_key"
    )

  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  all_SRR <- unique(selected$SRR)

  for (now_SRR in all_SRR) {
    need_keys <- selected %>%
      dplyr::filter(SRR == now_SRR) %>%
      dplyr::pull(event_key) %>%
      unique()

    selected_fasta <- info_df %>%
      dplyr::filter(event_key %in% need_keys)

    if (nrow(selected_fasta) == 0) next

    fa_lines <- as.vector(rbind(
      paste0(">", selected_fasta$short_ID),
      selected_fasta$seq
    ))

    writeLines(
      fa_lines,
      file.path(out_dir, paste0(now_SRR, ".fasta"))
    )
  }

  message(
    "Done: ", basename(out_dir),
    " | sjc_cut=", sjc_cut,
    " | SRR files=", length(all_SRR),
    " | kept rows=", nrow(selected)
  )

  invisible(selected)
}


jobs <- list(
  list(
    name = "HC_eRA",
    count = "./AS_rMATS_jcast/rMATS/Count/healthy_eRA_count.txt",
    info = "./AS_rMATS_jcast/jcast/result/T1_and_T2_info/HC_eRA_T1_T2_info.txt",
    out = "./AS_rMATS_jcast/netMHCpan/SRR_files_shortID/HC_eRA",
    sjc_cut = 116
  ),
  list(
    name = "HC_RA",
    count = "./AS_rMATS_jcast/rMATS/Count/healthy_RA_count.txt",
    info = "./AS_rMATS_jcast/jcast/result/T1_and_T2_info/HC_RA_T1_T2_info.txt",
    out = "./AS_rMATS_jcast/netMHCpan/SRR_files_shortID/HC_RA",
    sjc_cut = 116
  ),
  list(
    name = "RA_eRA",
    count = "./AS_rMATS_jcast/rMATS/Count/RA_eRA_count.txt",
    info = "./AS_rMATS_jcast/jcast/result/T1_and_T2_info/RA_eRA_T1_T2_info.txt",
    out = "./AS_rMATS_jcast/netMHCpan/SRR_files_shortID/RA_eRA",
    sjc_cut = 116
  )
)

sample_fasta_results <- lapply(
  jobs,
  function(j) {
    message("Running: ", j$name)

    write_srr_fastas_one(
      count_path = j$count,
      info_path = j$info,
      out_dir = j$out,
      sjc_cut = j$sjc_cut
    )
  }
)


# ============================================================
# 4. Merge sample FASTA files and deduplicate by sequence
# ============================================================

merge_srr_fastas_bySeq_keepID_biostrings <- function(in_dirs, out_dir) {
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  fa_all <- unlist(
    lapply(
      in_dirs,
      function(d) {
        list.files(d, pattern = "\\.fasta$", full.names = TRUE)
      }
    ),
    use.names = FALSE
  )

  if (length(fa_all) == 0) {
    stop("No fasta files found in input dirs.")
  }

  srr <- sub("\\.fasta$", "", basename(fa_all))

  df <- data.frame(
    path = fa_all,
    SRR = srr,
    stringsAsFactors = FALSE
  )

  srr_list <- split(df$path, df$SRR)

  for (now_SRR in names(srr_list)) {
    paths <- srr_list[[now_SRR]]

    # Input-directory order determines which ID is retained for duplicates.
    dir_rank <- match(dirname(paths), in_dirs)
    paths <- paths[order(dir_rank, paths)]

    aa_list <- lapply(paths, Biostrings::readAAStringSet)
    aa <- do.call(c, aa_list)

    if (length(aa) == 0) next

    seq_chr <- as.character(aa)
    keep_idx <- !duplicated(seq_chr)
    aa_unique <- aa[keep_idx]

    Biostrings::writeXStringSet(
      aa_unique,
      filepath = file.path(out_dir, paste0(now_SRR, ".fasta"))
    )
  }

  message("Merged fasta written to: ", out_dir)
}


in_dirs <- c(
  "./AS_rMATS_jcast/netMHCpan/SRR_files_shortID/HC_eRA",
  "./AS_rMATS_jcast/netMHCpan/SRR_files_shortID/HC_RA",
  "./AS_rMATS_jcast/netMHCpan/SRR_files_shortID/RA_eRA"
)

merged_srr_dir <- "./AS_rMATS_jcast/netMHCpan/SRR_Files_only_up"

merge_srr_fastas_bySeq_keepID_biostrings(
  in_dirs,
  merged_srr_dir
)


# ============================================================
# 5. Match sample FASTA files with HLA types
# ============================================================

HLA_1 <- read.delim2(
  "/public/workspace/ysmeng/A_info/HLA_1.txt",
  header = FALSE
)

HLA_2 <- read.delim2(
  "/public/workspace/ysmeng/A_info/HLA_2.txt",
  header = FALSE
)

colnames(HLA_1) <- c("file_name", "HLA_1")
colnames(HLA_2) <- c("file_name", "HLA_2")

file_info <- list.files(
  merged_srr_dir,
  full.names = TRUE
) %>%
  as.data.frame(stringsAsFactors = FALSE) %>%
  dplyr::mutate(
    file_name = stringr::str_extract(.[[1]], "SRR\\d+")
  )

colnames(file_info)[1] <- "file_path"

final_for_netMHCpan_1 <- file_info %>%
  dplyr::left_join(HLA_1, by = "file_name") %>%
  tidyr::drop_na() %>%
  dplyr::mutate(index = dplyr::row_number()) %>%
  dplyr::select(index, file_path, file_name, HLA_1)

final_for_netMHCpan_2 <- file_info %>%
  dplyr::left_join(HLA_2, by = "file_name") %>%
  tidyr::drop_na() %>%
  dplyr::mutate(index = dplyr::row_number()) %>%
  dplyr::select(index, file_path, file_name, HLA_2)


# ============================================================
# 6. Split netMHCpan input tables into three files
# ============================================================

write_three_parts <- function(df, out_dir) {
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  n <- nrow(df)
  split_size <- ceiling(n / 3)

  if (n == 0) {
    parts <- list(df, df, df)
  } else {
    idx <- split(
      seq_len(n),
      pmin(ceiling(seq_len(n) / split_size), 3)
    )

    parts <- lapply(1:3, function(i) {
      if (as.character(i) %in% names(idx)) {
        df[idx[[as.character(i)]], , drop = FALSE]
      } else {
        df[0, , drop = FALSE]
      }
    })
  }

  out_files <- file.path(
    out_dir,
    c("file_105.txt", "file_106.txt", "file_107.txt")
  )

  for (i in seq_along(parts)) {
    write.table(
      parts[[i]],
      out_files[i],
      quote = FALSE,
      row.names = FALSE,
      col.names = TRUE,
      sep = "\t"
    )
  }

  invisible(parts)
}


write_three_parts(
  final_for_netMHCpan_1,
  "./AS_rMATS_jcast/netMHCpan/HLA_1/bash_need"
)

write_three_parts(
  final_for_netMHCpan_2,
  "./AS_rMATS_jcast/netMHCpan/HLA_2/bash_need"
)


# ============================================================
# 7. Parse and filter netMHCpan HLA-I results
# ============================================================

netmhcpan_result_dir <- "./AS_rMATS_jcast/netMHCpan/HLA_1/netMHCpan_1_result"
netmhcpan_filter_dir <- "./AS_rMATS_jcast/netMHCpan/HLA_1/filter_result"

dir.create(netmhcpan_filter_dir, recursive = TRUE, showWarnings = FALSE)

input_files <- list.files(
  netmhcpan_result_dir,
  pattern = "\\.xls$",
  full.names = TRUE
)

for (file in input_files) {
  data_read <- read.delim(file, header = FALSE)

  ID <- data_read[3:nrow(data_read), 1:3]
  colnames(ID) <- c("Pos", "Peptide", "ID")

  col <- data_read[1:2, 4:ncol(data_read)]
  rownames(col) <- c("HLA_type", "field")
  col <- as.data.frame(t(col), stringsAsFactors = FALSE)

  col$HLA_type <- ifelse(
    stringr::str_length(col$HLA_type) == 0,
    NA,
    col$HLA_type
  )

  col <- col %>%
    tidyr::fill(HLA_type, .direction = "down") %>%
    dplyr::mutate(
      name = stringr::str_c(HLA_type, field, sep = "_")
    ) %>%
    dplyr::mutate(
      name = ifelse(field %in% c("Ave", "NB"), field, name)
    )

  df <- data_read[3:nrow(data_read), 4:ncol(data_read)]
  colnames(df) <- col$name
  df <- cbind(ID, df)

  HLA_num <- unique(stats::na.omit(col$HLA_type))

  df_long <- dplyr::bind_rows(
    lapply(
      HLA_num,
      function(hla_type) {
        df %>%
          dplyr::select(Pos, Peptide, ID, dplyr::matches(hla_type)) %>%
          dplyr::rename_with(~ stringr::str_replace(.x, hla_type, "")) %>%
          dplyr::mutate(type = hla_type, .after = 3)
      }
    )
  )

  df_long$`_EL_Rank` <- as.numeric(df_long$`_EL_Rank`)
  filtered <- df_long[df_long$`_EL_Rank` < 2, ]

  original_file_name <- tools::file_path_sans_ext(basename(file))
  new_file_name <- paste0("f_", original_file_name, ".txt")

  write.table(
    filtered,
    file.path(netmhcpan_filter_dir, new_file_name),
    row.names = FALSE,
    quote = FALSE,
    sep = "\t"
  )
}


# ============================================================
# 8. Summarize HLA-I binding ranks
# ============================================================

combined_f_1 <- read.delim(
  "./AS_rMATS_jcast/netMHCpan/HLA_1/combined_f.txt"
)

df_all <- combined_f_1 %>%
  dplyr::transmute(
    Peptide,
    ID,
    type,
    Rank = as.numeric(X_EL_Rank),
    source = "HLA_1"
  ) %>%
  dplyr::mutate(
    locus = dplyr::case_when(
      grepl("^HLA-A", type) ~ "HLA-A",
      grepl("^HLA-B", type) ~ "HLA-B",
      grepl("^HLA-C", type) ~ "HLA-C",
      TRUE ~ "Other"
    )
  ) %>%
  dplyr::filter(locus != "Other")


df_med <- df_all %>%
  dplyr::group_by(source, locus, Peptide) %>%
  dplyr::summarise(
    Rank_median = stats::median(Rank, na.rm = TRUE),
    n = dplyr::n(),
    .groups = "drop"
  ) %>%
  dplyr::mutate(
    locus = factor(locus, levels = c("HLA-A", "HLA-B", "HLA-C"))
  )

lab_df <- df_med %>%
  dplyr::group_by(locus) %>%
  dplyr::summarise(
    n_peptide = dplyr::n(),
    med = stats::median(Rank_median, na.rm = TRUE),
    .groups = "drop"
  )


# ============================================================
# 9. Plot HLA-I binding-rank distributions
# ============================================================

my_cols <- c(
  "HLA-A" = "#e3bec6",
  "HLA-B" = "#538085",
  "HLA-C" = "#faf1e2"
)

pdf(
  "./AS_rMATS_jcast/netMHCpan/boxplot.pdf",
  width = 10,
  height = 8
)

p_rank <- ggplot(
  df_med,
  aes(x = locus, y = Rank_median, fill = locus)
) +
  geom_violin(
    alpha = 0.85,
    width = 0.65,
    linewidth = 0.6,
    colour = "black"
  ) +
  geom_text(
    data = lab_df,
    aes(x = locus, y = med, label = sprintf("%.2f", med)),
    inherit.aes = FALSE,
    vjust = -3.5,
    size = 8
  ) +
  geom_text(
    data = lab_df,
    aes(x = locus, y = Inf, label = paste0("n=", n_peptide)),
    inherit.aes = FALSE,
    hjust = 0.6,
    vjust = 0.55,
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

print(p_rank)
dev.off()


ridge <- df_all %>%
  dplyr::group_by(locus, Peptide) %>%
  dplyr::summarise(
    rank_median = median(Rank),
    .groups = "drop"
  )

write.table(
  ridge,
  "./AS_rMATS_jcast/netMHCpan/ridge_1.txt",
  quote = FALSE,
  row.names = FALSE,
  sep = "\t"
)


# ============================================================
# 10. Summarize peptide prevalence across samples
# ============================================================

filter_result_dir <- "./AS_rMATS_jcast/netMHCpan/HLA_1/filter_result"

filter_files <- list.files(
  filter_result_dir,
  full.names = TRUE
)

sample_peptides <- lapply(
  filter_files,
  function(file_path) {
    df <- read.delim(file_path, header = TRUE)
    file_name <- tools::file_path_sans_ext(basename(file_path))

    df %>%
      dplyr::select(Peptide, type) %>%
      dplyr::mutate(
        type = substr(type, 1, 5),
        SRR = file_name
      ) %>%
      unique()
  }
)

comb_MHC1 <- dplyr::bind_rows(sample_peptides)
total_SRR <- length(sample_peptides)

MHC1_count <- comb_MHC1 %>%
  dplyr::group_by(type, Peptide) %>%
  dplyr::summarise(
    count = dplyr::n(),
    .groups = "drop"
  ) %>%
  dplyr::mutate(
    percentage = count / total_SRR * 100,
    percntg_range = cut(
      percentage,
      breaks = c(0, 10, 20, 30, 100),
      labels = c("0-10%", "10-20%", "20-30%", ">30%")
    )
  )

write.table(
  MHC1_count,
  "./AS_rMATS_jcast/netMHCpan/plot/sample_barplot/need_file/MHC1_for_pic_1.txt",
  quote = FALSE,
  row.names = FALSE,
  sep = "\t"
)

count2 <- MHC1_count %>%
  dplyr::group_by(type, percntg_range) %>%
  dplyr::summarise(counts = dplyr::n()) %>%
  dplyr::mutate(
    range_percntg = counts / sum(counts) * 100
  )

write.table(
  count2,
  "./AS_rMATS_jcast/netMHCpan/plot/sample_barplot/need_file/final_MHC1.txt",
  quote = FALSE,
  row.names = FALSE,
  sep = "\t"
)


# ============================================================
# 11. Plot peptide prevalence across HLA-I loci
# ============================================================

final_MHC1 <- read.delim(
  "./AS_rMATS_jcast/netMHCpan/plot/sample_barplot/need_file/final_MHC1.txt"
) %>%
  dplyr::mutate(HLA_type = "HLA−I")

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
    aes(x = HLA, y = counts, fill = Range)
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
    scale_y_continuous(labels = scales::comma) +
    theme_classic() +
    theme(
      panel.border = element_rect(
        color = "black",
        fill = NA,
        linewidth = 1
      ),
      strip.background = element_blank(),
      plot.title = element_text(hjust = 0.5, size = 22),
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
  "./AS_rMATS_jcast/netMHCpan/plot/sample_barplot/plot/paper_MHC_1_barplot.pdf",
  width = 6,
  height = 6
)

print(HLA_1_plot)
dev.off()


# ============================================================
# 12. Retain netMHCpan peptides derived from significant AS events
# ============================================================

combined_f <- read.delim(
  "~/circpeptide4/AS_rMATS_jcast/netMHCpan/HLA_1/combined_f.txt"
)

all_events_sig <- read.delim(
  "~/circpeptide4/AS_rMATS_jcast/rMATS/result/all_events_sig.txt"
)

all_shortID_map <- read.delim(
  "~/circpeptide4/AS_rMATS_jcast/jcast/result/T1_and_T2_info/all_shortID_map.txt"
)

all_events_sig <- all_events_sig %>%
  dplyr::mutate(
    event_id = as.character(event_id),
    event_type = as.character(event_type),
    contrast = as.character(contrast)
  )

all_shortID_map <- all_shortID_map %>%
  dplyr::mutate(
    event_id = as.character(event_id),
    event_type = as.character(event_type),
    short_ID = as.character(short_ID),
    contrast = dplyr::case_when(
      stringr::str_detect(short_ID, "^HC_eRA_") ~ "hc_vs_eRA",
      stringr::str_detect(short_ID, "^HC_RA_") ~ "hc_vs_RA",
      stringr::str_detect(short_ID, "^RA_eRA_") ~ "RA_vs_eRA",
      TRUE ~ NA_character_
    )
  )

combined_f <- combined_f %>%
  dplyr::mutate(
    ID = as.character(ID)
  )

combined_anno <- combined_f %>%
  dplyr::left_join(
    all_shortID_map,
    by = c("ID" = "short_ID")
  )

sig_key <- all_events_sig %>%
  dplyr::distinct(
    contrast,
    event_type,
    event_id,
    FDR,
    delta_psi,
    abs_dpsi,
    group
  )

combined_f_sig <- combined_anno %>%
  dplyr::inner_join(
    sig_key,
    by = c("contrast", "event_type", "event_id")
  )

combined_f_sig_new <- combined_f_sig[
  ,
  names(combined_f),
  drop = FALSE
]

write.table(
  combined_f_sig_new,
  "~/circpeptide4/AS_rMATS_jcast/netMHCpan/HLA_1/combined_f_sigAS.txt",
  quote = FALSE,
  row.names = FALSE,
  sep = "\t"
)
