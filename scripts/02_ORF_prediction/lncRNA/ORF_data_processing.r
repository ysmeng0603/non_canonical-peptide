# ============================================================
# CPAT and CPC2 coding-potential analysis for lncRNAs
# ============================================================

library(Biostrings)
library(dplyr)
library(stringr)
library(ggplot2)


# ============================================================
# 1. Filter CPAT results
# ============================================================

CPAT_files <- c(
  HC_eRA = "~/lncRNA_known/coding/CPAT_CPC2/result/CPAT/HC_eRA/CPAT_H_eRA_output.txt.ORF_prob.best.tsv",
  HC_RA  = "~/lncRNA_known/coding/CPAT_CPC2/result/CPAT/HC_RA/CPAT_H_RA_output.txt.ORF_prob.best.tsv",
  RA_eRA = "~/lncRNA_known/coding/CPAT_CPC2/result/CPAT/RA_eRA/CPAT_RA_eRA_output.txt.ORF_prob.best.tsv"
)

CPAT_all <- lapply(CPAT_files, read.delim)

CPAT_coding <- lapply(
  CPAT_all,
  function(x) {
    x %>%
      dplyr::filter(Coding_prob >= 0.364)
  }
)


# ============================================================
# 2. Filter CPC2 results
# ============================================================

CPC2_files <- c(
  HC_eRA = "~/lncRNA_known/coding/CPAT_CPC2/result/CPC2/H_eRA_output.txt.txt",
  HC_RA  = "~/lncRNA_known/coding/CPAT_CPC2/result/CPC2/H_RA_output.txt.txt",
  RA_eRA = "~/lncRNA_known/coding/CPAT_CPC2/result/CPC2/RA_eRA_output.txt.txt"
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
      dplyr::rename(seq_ID = V1)
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
      by = "seq_ID",
      suffix = c(".CPAT", ".CPC2")
    )
  }
)

names(coding_intersect) <- names(CPAT_coding)


# Export intersection results
outdir <- "./lncRNA_known/coding/CPAT_CPC2/result/filter"

dir.create(
  outdir,
  showWarnings = FALSE
)

for (group in names(coding_intersect)) {

  outfile <- file.path(
    outdir,
    paste0(group, "_CPAT_CPC2_intersect.tsv")
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
# 4. Extract FASTA sequences of coding lncRNAs
# ============================================================

lncRNA_fasta <- readDNAStringSet(
  "./reference/gencode.v49.lncRNA_transcripts.fa"
)

fasta_outdir <- "./lncRNA_known/coding/CPAT_CPC2/result/filter/fasta"


# Extract ENST IDs from full FASTA headers
extract_enst <- function(x) {
  str_extract(x, "ENST\\d+\\.\\d+")
}


# Export FASTA sequences and ID mapping tables
export_orf_fasta <- function(df, fasta, prefix) {

  ids_full <- df$seq_ID
  ids_enst <- extract_enst(ids_full)

  # Build ID mapping table
  map <- data.frame(
    ENST = ids_enst,
    full_ID = ids_full,
    stringsAsFactors = FALSE
  )

  write.table(
    map,
    file = file.path(
      fasta_outdir,
      paste0(prefix, "_ID_mapping.tsv")
    ),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )

  # Extract corresponding transcript sequences
  fa_sub <- fasta[ids_full]

  # Rename FASTA headers using ENST IDs
  names(fa_sub) <- ids_enst

  writeXStringSet(
    fa_sub,
    filepath = file.path(
      fasta_outdir,
      paste0(prefix, "_ENST.fasta")
    ),
    format = "fasta"
  )
}


for (group in names(coding_intersect)) {

  export_orf_fasta(
    coding_intersect[[group]],
    lncRNA_fasta,
    group
  )
}


# ============================================================
# 5. Summarize predicted ORF protein lengths
# ============================================================

protein_fasta_dir <- "./lncRNA_known/coding/ORFfinder/result/fasta/p_result"

groups <- c(
  "HC_eRA",
  "HC_RA",
  "RA_eRA"
)


# Read protein FASTA files
protein_data <- lapply(
  groups,
  function(group) {

    fasta <- readAAStringSet(
      file.path(
        protein_fasta_dir,
        paste0(group, "_result.fasta")
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
all_orfs <- do.call(rbind, protein_data) %>%
  unique() %>%
  dplyr::mutate(
    ENST_id = str_extract(header, "ENST\\d+\\.\\d+")
  )


# ============================================================
# 6. Calculate peptide-length distribution
# ============================================================

Length_plot_file <- all_orfs %>%
  dplyr::select(header, length) %>%
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
  dplyr::arrange(length_bin)


print(length_summary)


write.table(
  length_summary,
  "./lncRNA_known/lnc_length.csv",
  sep = ",",
  row.names = FALSE
)


# ============================================================
# 7. Plot the peptide-length distribution
# ============================================================

pdf(
  "./lncRNA_known/length.pdf",
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
    title = "Distribution of Predicted LncRNA Peptide Lengths",
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