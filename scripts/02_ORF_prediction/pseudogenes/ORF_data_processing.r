# ============================================================
# CPAT and CPC2 coding-potential analysis for pseudogenes
# ============================================================

library(Biostrings)
library(dplyr)
library(stringr)
library(ggplot2)


# ============================================================
# Directory settings
# ============================================================

step5_dir <- "./pseudogenes/result/step5_coding/result"

CPAT_dir <- file.path(step5_dir, "CPAT")
CPC2_dir <- file.path(step5_dir, "CPC2")
filter_dir <- file.path(step5_dir, "filter")
fasta_outdir <- file.path(filter_dir, "fasta")

pseudogene_fasta_file <-
  "./pseudogenes/result/step4_seq/DE_pseudogene_fasta/all_DE_pseudogene.transcripts.fasta"

protein_fasta_dir <-
  "./pseudogenes/result/step6_ORFfinder/result/fasta/p_result"


# ============================================================
# 1. Filter CPAT results
# ============================================================

CPAT_files <- c(
  HC_eRA = file.path(
    CPAT_dir,
    "HC_eRA/CPAT_H_eRA_output.txt.ORF_prob.best.tsv"
  ),
  HC_RA = file.path(
    CPAT_dir,
    "HC_RA/CPAT_H_RA_output.txt.ORF_prob.best.tsv"
  ),
  RA_eRA = file.path(
    CPAT_dir,
    "RA_eRA/CPAT_RA_eRA_output.txt.ORF_prob.best.tsv"
  )
)


CPAT_all <- lapply(
  CPAT_files,
  read.delim
)


# Extract the first two pipe-delimited fields from sequence IDs
extract_first2 <- function(x) {
  sub("^([^|]+\\|[^|]+).*", "\\1", x)
}


CPAT_coding <- lapply(
  CPAT_all,
  function(x) {

    x %>%
      dplyr::filter(Coding_prob >= 0.364) %>%
      dplyr::mutate(
        new_col = extract_first2(seq_ID)
      )

  }
)


# ============================================================
# 2. Filter CPC2 results
# ============================================================

CPC2_files <- c(
  HC_eRA = file.path(
    CPC2_dir,
    "H_eRA_output.txt.txt"
  ),
  HC_RA = file.path(
    CPC2_dir,
    "H_RA_output.txt.txt"
  ),
  RA_eRA = file.path(
    CPC2_dir,
    "RA_eRA_output.txt.txt"
  )
)


CPC2_all <- lapply(
  CPC2_files,
  read.delim,
  header = FALSE,
  comment.char = "#"
)


CPC2_coding <- lapply(
  CPC2_all,
  function(x) {

    x %>%
      dplyr::filter(V8 == "coding") %>%
      dplyr::rename(seq_ID = V1) %>%
      dplyr::mutate(
        new_col = extract_first2(seq_ID)
      )

  }
)


# ============================================================
# 3. Identify the intersection of CPAT and CPC2
# ============================================================

coding_intersect <- lapply(
  names(CPAT_coding),
  function(group) {

    inner_join(
      CPAT_coding[[group]],
      CPC2_coding[[group]],
      by = "new_col",
      suffix = c(".CPAT", ".CPC2")
    )

  }
)

names(coding_intersect) <- names(CPAT_coding)


# Export intersection results
dir.create(
  filter_dir,
  showWarnings = FALSE
)

for (group in names(coding_intersect)) {

  outfile <- file.path(
    filter_dir,
    paste0(
      group,
      "_CPAT_CPC2_intersect.tsv"
    )
  )

  write.table(
    coding_intersect[[group]],
    file = outfile,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )

}


# ============================================================
# 4. Summarize CPAT and CPC2 overlap
# ============================================================

venn_list <- list()
venn_stat <- data.frame()


for (group in names(CPAT_coding)) {

  # Unique IDs predicted by each tool
  cpat_id <- unique(
    CPAT_coding[[group]]$new_col
  )

  cpc2_id <- unique(
    CPC2_coding[[group]]$new_col
  )

  # Calculate Venn regions
  both_id <- intersect(
    cpat_id,
    cpc2_id
  )

  cpat_only_id <- setdiff(
    cpat_id,
    cpc2_id
  )

  cpc2_only_id <- setdiff(
    cpc2_id,
    cpat_id
  )

  # Input for Venn diagrams
  venn_list[[group]] <- list(
    CPAT = cpat_id,
    CPC2 = cpc2_id
  )

  # Venn statistics
  tmp_stat <- data.frame(
    Group = group,
    CPAT_only = length(cpat_only_id),
    CPC2_only = length(cpc2_only_id),
    Intersection = length(both_id),
    CPAT_total = length(cpat_id),
    CPC2_total = length(cpc2_id),
    Union_total = length(
      union(cpat_id, cpc2_id)
    )
  )

  venn_stat <- dplyr::bind_rows(
    venn_stat,
    tmp_stat
  )

}

print(venn_stat)


# ============================================================
# 5. Extract FASTA sequences
# ============================================================

pseudogenes_fasta <- readDNAStringSet(
  pseudogene_fasta_file
)


# Extract ENST IDs from sequence identifiers
extract_enst <- function(x) {
  str_extract(
    x,
    "ENST\\d+\\.\\d+"
  )
}


# Export FASTA sequences and ID mapping tables
export_orf_fasta <- function(df, fasta, prefix) {

  ids_full <- df$new_col
  ids_enst <- extract_enst(ids_full)

  # Build an index using the first two fields of FASTA headers
  fasta_key <- extract_first2(
    names(fasta)
  )

  names(fasta) <- fasta_key

  # Identify IDs that cannot be matched
  missing_ids <- setdiff(
    ids_full,
    names(fasta)
  )

  if (length(missing_ids) > 0) {

    warning(
      length(missing_ids),
      " IDs not found in fasta."
    )

  }

  # Retain matched IDs
  keep <- ids_full %in% names(fasta)

  ids_full_keep <- ids_full[keep]
  ids_enst_keep <- ids_enst[keep]

  # Build ID mapping table
  map <- data.frame(
    ENST = ids_enst_keep,
    full_ID = ids_full_keep,
    stringsAsFactors = FALSE
  )

  write.table(
    map,
    file = file.path(
      fasta_outdir,
      paste0(
        prefix,
        "_ID_mapping.tsv"
      )
    ),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )

  # Extract corresponding sequences
  fa_sub <- fasta[
    ids_full_keep
  ]

  # Rename FASTA headers using ENST IDs
  names(fa_sub) <- ids_enst_keep

  writeXStringSet(
    fa_sub,
    filepath = file.path(
      fasta_outdir,
      paste0(
        prefix,
        "_ENST.fasta"
      )
    ),
    format = "fasta"
  )

}


dir.create(
  fasta_outdir,
  showWarnings = FALSE
)


for (group in names(coding_intersect)) {

  export_orf_fasta(
    coding_intersect[[group]],
    pseudogenes_fasta,
    group
  )

}


# ============================================================
# 6. Summarize predicted ORF protein lengths
# ============================================================

groups <- c(
  "HC_eRA",
  "HC_RA",
  "RA_eRA"
)


protein_data <- lapply(
  groups,
  function(group) {

    fasta <- readAAStringSet(
      file.path(
        protein_fasta_dir,
        paste0(
          group,
          "_result.fasta"
        )
      )
    )

    data.frame(
      header = names(fasta),
      length = width(fasta),
      stringsAsFactors = FALSE
    )

  }
)

names(protein_data) <- groups


# Merge predicted ORFs
all_orfs <- do.call(
  rbind,
  protein_data
) %>%
  unique() %>%
  dplyr::mutate(
    ENST_id = str_extract(
      header,
      "ENST\\d+\\.\\d+"
    )
  )


# ============================================================
# 7. Calculate peptide-length distribution
# ============================================================

Length_plot_file <- all_orfs %>%
  dplyr::select(
    header,
    length
  ) %>%
  unique()


# Length intervals
breaks <- c(
  0,
  50,
  100,
  200,
  500,
  800,
  Inf
)

labels <- c(
  "1-50",
  "51-100",
  "101-200",
  "201-500",
  "501-800",
  ">800"
)


length_summary <- Length_plot_file %>%
  dplyr::mutate(
    length_bin = cut(
      length,
      breaks = breaks,
      labels = labels,
      include.lowest = TRUE,
      right = TRUE
    )
  ) %>%
  dplyr::count(
    length_bin,
    name = "count"
  ) %>%
  dplyr::mutate(
    proportion = count / sum(count),
    percentage = scales::percent(
      proportion,
      accuracy = 0.1
    )
  ) %>%
  dplyr::arrange(
    length_bin
  )


print(length_summary)


write.table(
  length_summary,
  "./pseudogenes/result/pesduo_length.csv",
  sep = ",",
  row.names = FALSE
)


# ============================================================
# 8. Plot peptide-length distribution
# ============================================================

pdf(
  "./pseudogenes/result/length.pdf",
  width = 12,
  height = 8
)


p <- ggplot(
  length_summary,
  aes(
    x = length_bin,
    y = count
  )
) +
  geom_col(
    fill = "#7593af",
    width = 0.7
  ) +
  geom_text(
    aes(
      label = paste(
        count,
        "\n(",
        percentage,
        ")"
      )
    ),
    vjust = -0.5,
    size = 5,
    fontface = "bold",
    colour = "black"
  ) +
  theme_bw() +
  labs(
    title = "Distribution of Predicted pseudogenes Peptide Lengths",
    x = "Length ",
    y = "Number of Peptides"
  ) +
  theme(
    plot.title = element_text(
      hjust = 0.5,
      size = 16,
      face = "bold"
    ),
    plot.subtitle = element_text(
      hjust = 0.5,
      size = 13
    ),
    axis.title = element_text(
      size = 14,
      face = "bold"
    ),
    axis.text = element_text(
      size = 12
    )
  )


print(p)

dev.off()