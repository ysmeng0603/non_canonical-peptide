# ============================================================
# Pseudogene-derived peptide preparation for HLA-I prediction
# ============================================================

library(Biostrings)
library(dplyr)
library(stringr)
library(tidyr)
library(purrr)


# ============================================================
# 1. Prepare ORF sequences and short-ID mapping
# ============================================================

protein_dir <- "./pseudogenes/result/step6_ORFfinder/result/fasta/p_result"

HC_eRA_result <- readAAStringSet(
  file.path(protein_dir, "HC_eRA_result.fasta")
)

HC_RA_result <- readAAStringSet(
  file.path(protein_dir, "HC_RA_result.fasta")
)

RA_eRA_result <- readAAStringSet(
  file.path(protein_dir, "RA_eRA_result.fasta")
)


# Merge ORF sequences and deduplicate by sequence name
all_seqs <- c(
  HC_eRA_result,
  HC_RA_result,
  RA_eRA_result
)

all_seqs <- all_seqs[
  !duplicated(names(all_seqs))
]


# Collect unique ORF IDs
all_need <- data.frame(
  ORF_ID = unique(c(
    names(HC_eRA_result),
    names(HC_RA_result),
    names(RA_eRA_result)
  )),
  stringsAsFactors = FALSE
)


# Build short-ID mapping
all_need <- all_need %>%
  dplyr::mutate(
    new_ID = paste0(
      "LncID_",
      seq_len(nrow(all_need))
    )
  )


write.table(
  all_need,
  "./pseudogenes/result/step8_netMHCpan/ori_shortID_map.txt",
  quote = FALSE,
  row.names = FALSE,
  sep = "\t"
)


# Extract transcript IDs from ORF headers
all_need <- all_need %>%
  dplyr::mutate(
    ENST_id = str_extract(
      ORF_ID,
      "(?<=_)[^:]+"
    )
  )


# ============================================================
# 2. Identify expressed upregulated pseudogenes in each sample
# ============================================================

pseudogenes_transcript_counts <- read.delim(
  "./pseudogenes/result/step2_feature_count/pseudogene.featureCounts.strict.txt",
  comment.char = "#"
)

all_pseudogenes_group <- read.delim(
  "./pseudogenes/result/step4_seq/DE_pseudogene_gene_transcript_group_map.tsv"
)


# Retain upregulated pseudogenes
need_up <- all_pseudogenes_group %>%
  dplyr::filter(
    significance == "Up"
  )


# Convert expression matrix to long format
expression_data <- pseudogenes_transcript_counts %>%
  dplyr::select(
    Geneid,
    matches("X.mnt.")
  ) %>%
  tidyr::pivot_longer(
    cols = 2:ncol(.)
  ) %>%
  dplyr::filter(
    value > 0
  ) %>%
  dplyr::filter(
    Geneid %in% need_up$pseudogene_id
  ) %>%
  dplyr::mutate(
    SRR_id = str_extract(
      name,
      "SRR\\d+"
    )
  )


# ============================================================
# 3. Map expressed pseudogenes to predicted ORFs
# ============================================================

orf2new <- all_need %>%
  dplyr::distinct(
    ORF_ID,
    new_ID
  )


# Retain ORFs with valid short-ID mappings
idx <- match(
  names(all_seqs),
  orf2new$ORF_ID
)

keep <- !is.na(idx)

all_seqs <- all_seqs[keep]

names(all_seqs) <- orf2new$new_ID[
  idx[keep]
]


# Build sample-level pseudogene expression table
expr_by_sample <- expression_data %>%
  dplyr::select(
    SRR_id,
    pseudogene_id = Geneid
  ) %>%
  dplyr::distinct()


# Map pseudogene IDs to transcript IDs
ENSG_ENST <- need_up %>%
  dplyr::select(
    pseudogene_id,
    ENST_id
  ) %>%
  unique()


expr_by_sample <- left_join(
  expr_by_sample,
  ENSG_ENST,
  by = "pseudogene_id"
)


# Map expressed transcripts to predicted ORFs
orf_by_sample <- expr_by_sample %>%
  inner_join(
    all_need %>%
      dplyr::select(
        ENST_id,
        new_ID
      ),
    by = "ENST_id"
  ) %>%
  distinct(
    SRR_id,
    new_ID
  )


# ============================================================
# 4. Generate sample-specific FASTA files
# ============================================================

sample_fasta_dir <- "./pseudogenes/result/step8_netMHCpan/result/HLA_1/SRR_Files_only_up"

dir.create(
  sample_fasta_dir,
  showWarnings = FALSE
)


for (srr in unique(orf_by_sample$SRR_id)) {

  ids <- orf_by_sample %>%
    dplyr::filter(
      SRR_id == srr
    ) %>%
    dplyr::pull(
      new_ID
    ) %>%
    unique()

  # Retain IDs present in the ORF sequence set
  ids <- intersect(
    ids,
    names(all_seqs)
  )

  fa_sub <- all_seqs[
    ids
  ]

  outfile <- file.path(
    sample_fasta_dir,
    paste0(
      srr,
      ".fasta"
    )
  )

  writeXStringSet(
    fa_sub,
    filepath = outfile,
    format = "fasta"
  )

}


# ============================================================
# 5. Match sample FASTA files with HLA-I genotypes
# ============================================================

HLA_1 <- read.delim2(
  "./pseudogenes/result/step8_netMHCpan/HLA_1.txt",
  header = FALSE
)

colnames(HLA_1) <- c(
  "file_name",
  "HLA_1"
)


file_info <- data.frame(
  file_path = list.files(
    sample_fasta_dir,
    full.names = TRUE
  ),
  stringsAsFactors = FALSE
) %>%
  dplyr::mutate(
    file_name = str_extract(
      file_path,
      "SRR\\d+"
    )
  )


final_for_netMHCpan_1 <- left_join(
  file_info,
  HLA_1,
  by = "file_name"
) %>%
  na.omit() %>%
  dplyr::mutate(
    index = seq_len(n())
  ) %>%
  dplyr::select(
    index,
    file_path,
    file_name,
    HLA_1
  )


write.table(
  final_for_netMHCpan_1,
  "./pseudogenes/result/step8_netMHCpan/result/HLA_1/bash_need/MHC1_file_paths_only_up.txt",
  quote = FALSE,
  row.names = FALSE,
  sep = "\t"
)


# ============================================================
# 6. Parse and filter netMHCpan HLA-I results
# ============================================================

netmhcpan_result_dir <- "./pseudogenes/result/step8_netMHCpan/result/HLA_1/netMHCpan_1_result"

filter_result_dir <- "./pseudogenes/result/step8_netMHCpan/result/HLA_1/filter_result"

dir.create(
  filter_result_dir,
  recursive = TRUE,
  showWarnings = FALSE
)


input_files <- list.files(
  netmhcpan_result_dir,
  pattern = "\\.xls$",
  full.names = TRUE
)


for (file in input_files) {

  data_read <- read.delim(
    file,
    header = FALSE
  )

  # Sequence information
  ID <- data_read[
    3:nrow(data_read),
    1:3
  ]

  colnames(ID) <- c(
    "Pos",
    "Peptide",
    "ID"
  )


  # Parse netMHCpan column headers
  col <- data_read[
    1:2,
    4:ncol(data_read)
  ]

  rownames(col) <- c(
    "HLA_type",
    "field"
  )

  col <- col %>%
    t() %>%
    as.data.frame()

  col$HLA_type <- ifelse(
    str_length(col$HLA_type) == 0,
    NA,
    col$HLA_type
  )

  col <- col %>%
    tidyr::fill(
      HLA_type,
      .direction = "down"
    ) %>%
    dplyr::mutate(
      name = str_c(
        HLA_type,
        field,
        sep = "_"
      )
    ) %>%
    dplyr::mutate(
      name = ifelse(
        field %in% c("Ave", "NB"),
        field,
        name
      )
    )


  # Parse prediction values
  df <- data_read[
    3:nrow(data_read),
    4:ncol(data_read)
  ]

  colnames(df) <- col$name

  df <- cbind(
    ID,
    df
  )


  HLA_num <- unique(
    col$HLA_type
  )

  HLA_num <- HLA_num[
    !is.na(HLA_num)
  ]


  df <- purrr::map_dfr(
    HLA_num,
    function(hla_type) {

      df %>%
        dplyr::select(
          Pos,
          Peptide,
          ID,
          matches(hla_type)
        ) %>%
        dplyr::rename_with(
          ~ str_replace(
            .,
            fixed(hla_type),
            ""
          )
        ) %>%
        dplyr::mutate(
          type = hla_type,
          .after = 3
        )

    }
  )


  # Retain netMHCpan binders with EL_Rank < 2
  df$`_EL_Rank` <- as.numeric(
    df$`_EL_Rank`
  )

  filtered <- df %>%
    dplyr::filter(
      `_EL_Rank` < 2
    )


  original_file_name <- tools::file_path_sans_ext(
    basename(file)
  )

  output_file <- file.path(
    filter_result_dir,
    paste0(
      "f_",
      original_file_name,
      ".txt"
    )
  )


  write.table(
    filtered,
    output_file,
    row.names = FALSE,
    quote = FALSE,
    sep = "\t"
  )

}


# ============================================================
# 7. Identify MS-supported HLA-I-binding peptides
# ============================================================

combined_f <- read.delim(
  "~/circpeptide4/pseudogenes/result/step8_netMHCpan/result/HLA_1/combined_f.txt"
)

containHuman <- read.delim(
  "./pseudogenes/result/step7_Maxquant/result/contain_pseudogenes/z_filter_result/containHuman.txt"
)

ori_shortID_map <- read.delim(
  "./pseudogenes/result/step8_netMHCpan/ori_shortID_map.txt"
)


# Build ORF short-ID mapping
map2 <- ori_shortID_map %>%
  dplyr::mutate(
    orf_key = stringr::str_extract(
      ORF_ID,
      "ORF\\d+_ENST\\d+\\.\\d+"
    )
  ) %>%
  dplyr::select(
    orf_key,
    new_ID
  ) %>%
  dplyr::distinct()


# Convert MS-supported peptide IDs to short-ID format
all_MS <- containHuman %>%
  dplyr::mutate(
    new = str_extract(
      Proteins,
      "(?<=lcl\\|)[^:]+"
    )
  ) %>%
  dplyr::mutate(
    new = str_c(
      new,
      Sequence,
      sep = "&"
    )
  ) %>%
  dplyr::select(
    new,
    everything()
  )


all_MS <- all_MS %>%
  dplyr::mutate(
    orf_key = sub(
      "&.*$",
      "",
      new
    ),
    peptide = sub(
      "^.*&",
      "",
      new
    )
  ) %>%
  dplyr::left_join(
    map2,
    by = "orf_key"
  ) %>%
  dplyr::mutate(
    new = ifelse(
      !is.na(new_ID),
      paste0(
        new_ID,
        "&",
        peptide
      ),
      new
    )
  ) %>%
  dplyr::select(
    -orf_key,
    -peptide
  )


# Match netMHCpan predictions with MS-supported peptides
MS_HLA1 <- combined_f %>%
  dplyr::mutate(
    new = str_c(
      ID,
      Peptide,
      sep = "&"
    )
  ) %>%
  dplyr::select(
    new,
    everything()
  ) %>%
  dplyr::filter(
    new %in% all_MS$new
  )


# ============================================================
# 8. Add UniProt sequence-attribution information
# ============================================================

is_contain_sp <- read.delim(
  "~/circpeptide4/pseudogenes/result/step7_Maxquant/result/contain_pseudogenes/z_filter_result/is_contain_sp.txt"
)


is_contain_sp <- is_contain_sp %>%
  dplyr::select(
    Sequence,
    hit_any_sp
  ) %>%
  unique()

colnames(is_contain_sp)[1] <- "Peptide"


MQ_netMHCpan <- left_join(
  MS_HLA1,
  is_contain_sp,
  by = "Peptide"
)


write.table(
  MQ_netMHCpan,
  "./pseudogenes/result/MQ_netMHCpan.txt",
  quote = FALSE,
  row.names = FALSE,
  sep = "\t"
)
