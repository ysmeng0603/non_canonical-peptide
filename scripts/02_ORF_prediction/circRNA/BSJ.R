# ============================================================
# Identification of circRNA ORFs spanning back-splice junctions
# ============================================================

suppressPackageStartupMessages({
  library(Biostrings)
  library(dplyr)
  library(stringr)
  library(ggplot2)
})

# -----------------------------------------------------------------------------
# Configuration
# -----------------------------------------------------------------------------

base_dir <- "./circRNA_peptide3"

groups <- c(
  "H_eRA_down",
  "H_eRA_up",
  "H_RA_down",
  "H_RA_up",
  "RA_eRA_down"
)

# -----------------------------------------------------------------------------
# Step 1: Generate FASTA path lists for circRNAs in the Venn intersections
# -----------------------------------------------------------------------------

prepare_fasta_paths <- function(group) {
  input_file <- file.path(
    base_dir,
    "part1_coding_prob/A_venn/data",
    paste0(group, ".txt")
  )

  fasta_dir <- file.path(
    base_dir,
    "get_sequence/sequence_file",
    group,
    "fasta"
  )

  output_file <- file.path(
    base_dir,
    "part2_BSJ/seq_repeat_4",
    group,
    paste0(group, "_FilePath.txt")
  )

  dat <- read.table(
    input_file,
    quote = "\"",
    comment.char = "",
    header = FALSE
  )

  dat <- dat %>%
    dplyr::mutate(
      path1 = paste0(fasta_dir, "/"),
      path3 = ".fasta"
    ) %>%
    dplyr::select(path1, V1, path3) %>%
    dplyr::mutate(full_path = stringr::str_c(path1, V1, path3))

  dat <- dat[4]

  write.table(
    dat,
    output_file,
    quote = FALSE,
    row.names = FALSE,
    col.names = FALSE,
    sep = "\t"
  )
}

invisible(lapply(groups, prepare_fasta_paths))

# -----------------------------------------------------------------------------
# Step 2: Repeat each circRNA sequence four times
# -----------------------------------------------------------------------------

repeat_fasta_four_times <- function(group) {
  input_folder <- file.path(
    base_dir,
    "part2_BSJ/seq_repeat_4",
    group,
    "data_seq4/select_raw"
  )

  output_folder <- file.path(
    base_dir,
    "part2_BSJ/seq_repeat_4",
    group,
    "data_seq4/repeat4"
  )

  fasta_files <- list.files(
    input_folder,
    pattern = "\\.fasta$",
    full.names = TRUE
  )

  for (file in fasta_files) {
    file_name <- basename(file)
    fasta <- readDNAStringSet(file)

    repeated_sequences <- sapply(
      as.character(fasta),
      function(seq) paste(rep(seq, 4), collapse = "")
    )

    names(repeated_sequences) <- names(fasta)
    output_fasta_obj <- DNAStringSet(repeated_sequences)
    output_file <- file.path(output_folder, file_name)

    writeXStringSet(output_fasta_obj, filepath = output_file)
    cat("Processed file:", file_name, "->", output_file, "\n")
  }

  cat("All files processed. Output directory:", output_folder, "\n")
}

invisible(lapply(groups, repeat_fasta_four_times))

# -----------------------------------------------------------------------------
# Step 3: Identify ORFs containing the BSJ sequence
# -----------------------------------------------------------------------------

identify_bsj_orfs <- function(group) {
  result_file_path <- file.path(
    base_dir,
    "part2_BSJ/seq_repeat_4",
    group,
    paste0("result_", group)
  )

  bsj_file_path <- file.path(
    base_dir,
    "part2_BSJ/seq_repeat_4",
    group,
    "BSJ_files"
  )

  summary_file <- file.path(
    base_dir,
    "part2_BSJ/seq_repeat_4",
    group,
    "contain_BSJ_result",
    paste0(group, "_summary_result.txt")
  )

  result_files <- list.files(result_file_path, full.names = TRUE)
  bsj_files <- list.files(bsj_file_path, full.names = TRUE)

  output_file <- file(summary_file, "w")
  on.exit(close(output_file), add = TRUE)

  for (i in seq_along(result_files)) {
    result_file <- result_files[i]
    bsj_file <- bsj_files[i]

    if (!file.exists(result_file) || !file.exists(bsj_file)) {
      cat(
        "Result file or BSJ file does not exist for:",
        basename(result_file),
        "\n"
      )
      next
    }

    result_data <- readDNAStringSet(result_file)
    bsj_sequence <- readLines(bsj_file)

    cat("Processing result file:", basename(result_file), "\n")

    for (j in seq_along(result_data)) {
      seq_name <- names(result_data)[j]
      seq <- result_data[[j]]
      match_positions <- matchPattern(bsj_sequence, seq)

      if (length(match_positions) > 0) {
        cat("Sequence", seq_name, "contains the BSJ sequence\n")
        cat(
          "Sequence", seq_name, "contains the BSJ sequence\n",
          file = output_file
        )
      } else {
        cat("BSJ sequence not found in", seq_name, "\n")
        cat(
          "BSJ sequence not found in", seq_name, "\n",
          file = output_file
        )
      }
    }
  }
}

invisible(lapply(groups, identify_bsj_orfs))

# -----------------------------------------------------------------------------
# Step 4: Retain ORFs spanning the BSJ
# -----------------------------------------------------------------------------

extract_bsj_positive_orfs <- function(group) {
  summary_file <- file.path(
    base_dir,
    "part2_BSJ/seq_repeat_4",
    group,
    "contain_BSJ_result",
    paste0(group, "_summary_result.txt")
  )

  select_file <- file.path(
    base_dir,
    "part2_BSJ/seq_repeat_4",
    group,
    "contain_BSJ_result/Files",
    paste0("select_", group, "_list.txt")
  )

  summary_result <- read.table(
    summary_file,
    quote = "\"",
    comment.char = ""
  )

  select_result <- summary_result %>%
    dplyr::filter(V4 == "contains")

  write.table(
    select_result,
    select_file,
    quote = FALSE,
    row.names = FALSE,
    sep = "\t"
  )

  select_result
}

selected_orf_lists <- setNames(
  lapply(groups, extract_bsj_positive_orfs),
  groups
)

# -----------------------------------------------------------------------------
# Step 5: Filter nucleotide and amino-acid ORF FASTA files using BSJ-positive IDs
# -----------------------------------------------------------------------------

filter_fasta_by_conditions <- function(
    input_directory,
    output_directory,
    filter_conditions,
    sequence_type = c("DNA", "AA")) {

  sequence_type <- match.arg(sequence_type)
  input_files <- list.files(input_directory, full.names = TRUE)

  for (input_file in input_files) {
    sequences <- if (sequence_type == "DNA") {
      readDNAStringSet(input_file)
    } else {
      readAAStringSet(input_file)
    }

    selected_indices <- logical(length(sequences))

    for (j in seq_along(sequences)) {
      sequence_name <- names(sequences[j])

      for (condition in filter_conditions) {
        if (grepl(condition, sequence_name)) {
          selected_indices[j] <- TRUE
          break
        }
      }
    }

    selected_sequences <- sequences[selected_indices]
    file_name <- tools::file_path_sans_ext(basename(input_file))

    output_file <- file.path(
      output_directory,
      paste0(file_name, "_selected.fasta")
    )

    writeXStringSet(selected_sequences, output_file)
  }
}

for (group in groups) {
  select_file <- file.path(
    base_dir,
    "part2_BSJ/seq_repeat_4",
    group,
    "contain_BSJ_result/Files",
    paste0("select_", group, "_list.txt")
  )

  select_list <- read.delim(select_file, header = TRUE)
  filter_conditions <- select_list$V3

  result_path <- file.path(
    base_dir,
    "part2_BSJ/seq_repeat_4",
    group,
    paste0("result_", group)
  )

  p_result_path <- file.path(
    base_dir,
    "part2_BSJ/seq_repeat_4",
    group,
    paste0("p_result_", group)
  )

  result_output_directory <- file.path(
    base_dir,
    "part2_BSJ/seq_repeat_4",
    group,
    "contain_BSJ_result/Files/select_result"
  )

  p_result_output_directory <- file.path(
    base_dir,
    "part2_BSJ/seq_repeat_4",
    group,
    "contain_BSJ_result/Files/select_p_result"
  )

  filter_fasta_by_conditions(
    input_directory = result_path,
    output_directory = result_output_directory,
    filter_conditions = filter_conditions,
    sequence_type = "DNA"
  )

  filter_fasta_by_conditions(
    input_directory = p_result_path,
    output_directory = p_result_output_directory,
    filter_conditions = filter_conditions,
    sequence_type = "AA"
  )
}

# -----------------------------------------------------------------------------
# Step 6: Calculate and plot the proportion of BSJ-spanning ORFs
# -----------------------------------------------------------------------------

setwd("./circRNA_peptide3/part2_BSJ")

summary_results <- setNames(
  lapply(groups, function(group) {
    read.table(
      file.path(
        base_dir,
        "part2_BSJ/seq_repeat_4",
        group,
        "contain_BSJ_result",
        paste0(group, "_summary_result.txt")
      ),
      quote = "\"",
      comment.char = ""
    )
  }),
  groups
)

df_list <- vector("list", length(groups))

for (i in seq_along(groups)) {
  group <- groups[i]
  df_list[[i]] <- as.data.frame(table(summary_results[[group]]$V1))

  df_list[[i]]$Group <- if (group == "RA_eRA_down") {
    "RA_eRA"
  } else {
    group
  }
}

comb <- do.call(rbind, df_list)

comb <- comb %>%
  dplyr::mutate(
    Var1 = dplyr::case_when(
      Var1 == "BSJ" ~ "Negative",
      Var1 == "Sequence" ~ "Positive"
    )
  ) %>%
  dplyr::group_by(Group) %>%
  dplyr::mutate(
    Sum = sum(Freq),
    Percentage = round(Freq / Sum * 100, 2),
    label = paste0(Percentage, "%")
  )

df_new <- comb %>%
  dplyr::mutate(
    new_group = dplyr::case_when(
      Group %in% c("H_eRA_down", "H_eRA_up") ~ "HC_eRA",
      Group %in% c("H_RA_down", "H_RA_up") ~ "HC_RA",
      Group %in% c("RA_eRA", "RA_eRA_down") ~ "RA_eRA",
      TRUE ~ Group
    )
  ) %>%
  dplyr::group_by(new_group, Var1) %>%
  dplyr::summarise(
    Freq = sum(Freq),
    .groups = "drop"
  ) %>%
  dplyr::group_by(new_group) %>%
  dplyr::mutate(
    Sum = sum(Freq),
    Percentage = round(Freq / Sum * 100, 2),
    label = paste0(Percentage, "%")
  ) %>%
  dplyr::ungroup() %>%
  dplyr::rename(Group = new_group) %>%
  dplyr::select(Var1, Freq, Group, Sum, Percentage, label)

write.table(
  df_new,
  "./circRNA_peptide3/part2_BSJ/plot/new_plot_file.txt",
  quote = FALSE,
  row.names = FALSE,
  sep = "\t"
)

pdf(
  "./circRNA_peptide3/part2_BSJ/paper/new_BSJ_pie_4.pdf",
  width = 8,
  height = 8
)

p <- ggplot(df_new, aes(x = "", y = Percentage, fill = Var1)) +
  geom_bar(stat = "identity", width = 1) +
  coord_polar(theta = "y") +
  labs(fill = "Result") +
  scale_fill_manual(
    values = c(
      "Positive" = "#eeb8c3",
      "Negative" = "#8abcd1"
    )
  ) +
  facet_wrap(~Group, strip.position = "top") +
  theme_minimal() +
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.grid = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks = element_blank(),
    strip.text = element_text(size = 18, face = "bold"),
    legend.title = element_text(size = 13),
    legend.text = element_text(size = 11)
  ) +
  geom_text(
    aes(label = label),
    position = position_stack(vjust = 0.5),
    size = 5
  )

print(p)
dev.off()
