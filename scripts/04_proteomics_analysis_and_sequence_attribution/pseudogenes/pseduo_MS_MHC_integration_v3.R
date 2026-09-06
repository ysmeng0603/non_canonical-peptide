# 1. Build the pseudogene ORF dictionary

library(data.table)
library(Biostrings)
library(dplyr)
library(tidyr)
library(stringr)
library(readr)


output_dir <- paste0(
  "./MS_unspecific/pseudogenes/",
  "result/contain_pseudogenes/pseudogene_integrated_analysis"
)

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

fasta_dir <- paste0(
  "./MS_unspecific/pseudogenes/",
  "need_file/contain_human"
)

id_map_file <- paste0(
  "./circpeptide4/pseudogenes/",
  "result/step8_netMHCpan/ID_map.txt"
)


# Read the three comparison FASTA files

fasta_files <- list.files(
  fasta_dir,
  pattern = "^contain_human_(HC_eRA|HC_RA|RA_eRA)\\.fasta$",
  full.names = TRUE
)


if (length(fasta_files) != 3) {
  stop("Expected three comparison FASTA files but did not find all three")
}


# Read the NetMHCpan ID mapping

id_map <- fread(id_map_file)

setnames(
  id_map,
  old = "new_ID",
  new = "netmhc_ID"
)

# Extract pseudogene ORFs

pseudo_orf_dict <- rbindlist(
  lapply(fasta_files, function(f) {

    aa <- readAAStringSet(f)

    header <- names(aa)
    fasta_ID <- sub("\\s.*$", "", header)

    keep <- grepl("^lcl\\|ORF[0-9]+_ENST", fasta_ID)

    data.table(
      comparison_group = sub(
        "^contain_human_|\\.fasta$",
        "",
        basename(f)
      ),
      fasta_ID = fasta_ID[keep],
      orf_key = sub(
        ":[0-9]+:[0-9]+$",
        "",
        sub("^lcl\\|", "", fasta_ID[keep])
      ),
      protein_sequence = as.character(aa[keep])
    )
  })
)

pseudo_orf_dict[, protein_length := nchar(protein_sequence)]

pseudo_orf_dict <- merge(
  pseudo_orf_dict,
  id_map,
  by = "orf_key",
  all.x = TRUE,
  sort = FALSE
)


# Save outputs

saveRDS(
  pseudo_orf_dict,
  file.path(output_dir, "01_pseudogene_ORF_dictionary.rds")
)

fwrite(
  pseudo_orf_dict,
  file.path(output_dir, "01_pseudogene_ORF_dictionary.tsv"),
  sep = "\t"
)


# 2. Extract high-confidence pseudogene-associated PSMs


# 2.1 Standardize comparison-group names

pseudo_orf_dict[
  ,
  comparison_group := sub("\\.fasta$", "", comparison_group)
]

saveRDS(
  pseudo_orf_dict,
  file.path(output_dir, "01_pseudogene_ORF_dictionary.rds")
)


# 2.2 Define the six MaxQuant result paths


mq_root <- paste0(
  "./MS_unspecific/pseudogenes/",
  "result/contain_pseudogenes"
)

group_map <- data.table(
  dataset = rep(c("PXD037581", "PXD044963"), each = 3),
  comparison_group = rep(
    c("HC_eRA", "HC_RA", "RA_eRA"),
    times = 2
  ),
  quantification = c(
    rep("TMT10plex", 3),
    rep("LabelFree", 3)
  )
)

group_map[
  ,
  msms_file := file.path(
    mq_root,
    dataset,
    paste0(comparison_group, "_result"),
    "combined/txt/msms.txt"
  )
]


# 2.3 Filter high-confidence PSMs containing pseudogene ORFs

read_pseudo_msms <- function(dataset_name,
                             group_name,
                             quant_name,
                             msms_path) {

  x <- fread(
    msms_path,
    check.names = TRUE,
    showProgress = FALSE
  )

  x <- x[
    (is.na(Reverse) | Reverse != "+") &
      (is.na(Contaminant) | Contaminant != "+") &
      Length >= 7 &
      PEP < 0.01 &
      Score > 40 &
      Delta.score >= 15 &
      !is.na(Proteins) &
      grepl(
        "(^|;)(lcl\\|)?ORF[0-9]+_ENST",
        Proteins
      )
  ]

  x[, `:=`(
    dataset = dataset_name,
    comparison_group = group_name,
    quantification = quant_name
  )]

  x[, PSM_uid := paste(
    dataset,
    comparison_group,
    Raw.file,
    Scan.number,
    id,
    sep = "|"
  )]

  x[, .(
    PSM_uid,
    dataset,
    comparison_group,
    quantification,
    Sequence,
    Length,
    Proteins,
    PEP,
    Score,
    Delta_score = Delta.score,
    PIF,
    Raw_file = Raw.file,
    Scan_number = Scan.number,
    msms_id = id,
    Evidence_ID = Evidence.ID
  )]
}


# Read and combine the six MaxQuant groups

pseudo_PSM <- rbindlist(
  lapply(seq_len(nrow(group_map)), function(i) {

    read_pseudo_msms(
      dataset_name = group_map$dataset[i],
      group_name = group_map$comparison_group[i],
      quant_name = group_map$quantification[i],
      msms_path = group_map$msms_file[i]
    )
  }),
  use.names = TRUE
)

pseudo_PSM[, psm_row := .I]


# 2.4 Map PSMs to pseudogene ORFs

# Split the Proteins field and extract pseudogene ORF entries
pseudo_PSM_source <- pseudo_PSM[
  ,
  {
    protein_ids <- trimws(
      unlist(strsplit(Proteins, ";", fixed = TRUE))
    )

    protein_ids <- protein_ids[
      !is.na(protein_ids) & protein_ids != ""
    ]

    is_pseudo <- grepl(
      "^(lcl\\|)?ORF[0-9]+_ENST",
      protein_ids
    )

    pseudo_ids <- protein_ids[is_pseudo]
    other_ids  <- protein_ids[!is_pseudo]

    .(
      fasta_ID = pseudo_ids,
      has_non_pseudogene_protein = length(other_ids) > 0,
      non_pseudogene_proteins = if (length(other_ids) > 0) {
        paste(other_ids, collapse = ";")
      } else {
        NA_character_
      }
    )
  },
  by = psm_row
]


# Join the source mapping back to PSM information

pseudo_PSM_source <- merge(
  pseudo_PSM_source,
  pseudo_PSM,
  by = "psm_row",
  all.x = TRUE,
  sort = FALSE
)


# Standardize ORF identifiers

pseudo_PSM_source[
  ,
  orf_key := sub(
    ":[0-9]+:[0-9]+$",
    "",
    sub("^lcl\\|", "", fasta_ID)
  )
]


# Join the ORF dictionary
pseudo_PSM_source <- merge(
  pseudo_PSM_source,
  pseudo_orf_dict[
    ,
    .(
      comparison_group,
      orf_key,
      netmhc_ID,
      protein_sequence,
      protein_length
    )
  ],
  by = c("comparison_group", "orf_key"),
  all.x = TRUE,
  sort = FALSE
)


# Confirm that each peptide occurs in the corresponding ORF protein sequence

pseudo_PSM_source[
  ,
  peptide_in_ORF :=
    !is.na(protein_sequence) &
    mapply(
      grepl,
      pattern = Sequence,
      x = protein_sequence,
      MoreArgs = list(fixed = TRUE)
    )
]


# Retain valid PSM-ORF mappings

pseudo_PSM_source_valid <- pseudo_PSM_source[
  !is.na(protein_sequence) &
    peptide_in_ORF == TRUE
]


# Build QC summary

step2_qc <- data.table(
  metric = c(
    "high-confidence PSMs",
    "unique peptides",
    "PSM-ORF mapping rows",
    "unique pseudogene ORFs",
    "unmapped ORF rows",
    "peptide absent from ORF",
    "PSMs also matching non-pseudogene proteins"
  ),
  value = c(
    uniqueN(pseudo_PSM$PSM_uid),
    uniqueN(pseudo_PSM$Sequence),
    nrow(pseudo_PSM_source),
    uniqueN(pseudo_PSM_source_valid$orf_key),
    sum(is.na(pseudo_PSM_source$protein_sequence)),
    sum(!pseudo_PSM_source$peptide_in_ORF),
    uniqueN(
      pseudo_PSM_source[
        has_non_pseudogene_protein == TRUE,
        PSM_uid
      ]
    )
  )
)


# Save outputs

saveRDS(
  pseudo_PSM,
  file.path(
    output_dir,
    "02a_pseudogene_high_confidence_PSM.rds"
  )
)

saveRDS(
  pseudo_PSM_source_valid,
  file.path(
    output_dir,
    "02b_pseudogene_valid_PSM_ORF_mapping.rds"
  )
)

fwrite(
  step2_qc,
  file.path(
    output_dir,
    "02c_pseudogene_PSM_ORF_mapping_QC.tsv"
  ),
  sep = "\t"
)


# 3. Integrate NetMHCpan HLA-I predictions


# 3.1 Summarize MS evidence by ORF-peptide pair


safe_min <- function(x) {
  if (all(is.na(x))) NA_real_ else min(x, na.rm = TRUE)
}

safe_max <- function(x) {
  if (all(is.na(x))) NA_real_ else max(x, na.rm = TRUE)
}

pseudo_source_peptide <- pseudo_PSM_source_valid[
  ,
  .(
    peptide_length = unique(Length)[1],
    protein_length = unique(protein_length)[1],

    PSM_count = uniqueN(PSM_uid),

    raw_file_count = uniqueN(
      paste(dataset, comparison_group, Raw_file, sep = "|")
    ),

    datasets = paste(
      sort(unique(dataset)),
      collapse = ";"
    ),

    comparison_groups = paste(
      sort(unique(comparison_group)),
      collapse = ";"
    ),

    quantifications = paste(
      sort(unique(quantification)),
      collapse = ";"
    ),

    min_PEP = safe_min(PEP),
    max_Score = safe_max(Score),
    max_Delta_score = safe_max(Delta_score),
    max_PIF = safe_max(PIF),

    has_non_pseudogene_protein =
      any(has_non_pseudogene_protein)
  ),
  by = .(
    orf_key,
    netmhc_ID,
    peptide = Sequence
  )
]


# 3.2 Read and standardize the HLA-I prediction table
netmhc_file <- paste0(
  "./circpeptide4/pseudogenes/",
  "result/step8_netMHCpan/result/HLA_1/combined_f.txt"
)

hla_I <- fread(
  netmhc_file,
  select = c(
    "Peptide",
    "ID",
    "type",
    "_EL-score",
    "_EL_Rank"
  ),
  showProgress = TRUE
)

setnames(
  hla_I,
  old = c(
    "Peptide",
    "ID",
    "type",
    "_EL-score",
    "_EL_Rank"
  ),
  new = c(
    "peptide",
    "netmhc_ID",
    "HLA_allele",
    "EL_score",
    "EL_Rank"
  )
)

hla_I[, peptide := toupper(trimws(peptide))]
hla_I[, netmhc_ID := trimws(netmhc_ID)]
hla_I[, HLA_allele := trimws(HLA_allele)]

hla_I[, EL_score := suppressWarnings(as.numeric(EL_score))]
hla_I[, EL_Rank := suppressWarnings(as.numeric(EL_Rank))]


# Retain valid HLA-I binders
hla_I <- hla_I[
  !is.na(EL_Rank) &
    EL_Rank < 2 &
    grepl("^[A-Z]+$", peptide) &
    netmhc_ID != "ID" &
    HLA_allele != "type"
]

hla_I[
  ,
  binder_class := fifelse(
    EL_Rank < 0.5,
    "SB",
    "WB"
  )
]


# Retain the best prediction for each ORF-peptide-HLA combination

setorder(
  hla_I,
  netmhc_ID,
  peptide,
  HLA_allele,
  EL_Rank,
  -EL_score
)

hla_I_best <- hla_I[
  ,
  .SD[1],
  by = .(
    netmhc_ID,
    peptide,
    HLA_allele
  )
]


# 3.3 Strictly integrate MS and HLA-I evidence
pseudo_MS_HLA <- merge(
  pseudo_source_peptide,
  hla_I_best,
  by = c("netmhc_ID", "peptide"),
  all = FALSE,
  allow.cartesian = TRUE,
  sort = FALSE
)


# Summarize HLA-I evidence for each ORF-peptide pair

hla_summary <- pseudo_MS_HLA[
  ,
  .(
    HLA_allele_count = uniqueN(HLA_allele),

    HLA_alleles = paste(
      sort(unique(HLA_allele)),
      collapse = ";"
    )
  ),
  by = .(
    orf_key,
    netmhc_ID,
    peptide
  )
]


# Select the best HLA-I prediction

setorder(
  pseudo_MS_HLA,
  orf_key,
  netmhc_ID,
  peptide,
  EL_Rank,
  -EL_score
)

best_HLA <- pseudo_MS_HLA[
  ,
  .SD[1],
  by = .(
    orf_key,
    netmhc_ID,
    peptide
  )
][
  ,
  .(
    orf_key,
    netmhc_ID,
    peptide,
    best_HLA_allele = HLA_allele,
    best_EL_score = EL_score,
    best_EL_Rank = EL_Rank,
    best_binder_class = binder_class
  )
]

# Generate one row per ORF-peptide pair

pseudo_source_HLA_summary <- merge(
  pseudo_source_peptide,
  hla_summary,
  by = c("orf_key", "netmhc_ID", "peptide"),
  all = FALSE,
  sort = FALSE
)

pseudo_source_HLA_summary <- merge(
  pseudo_source_HLA_summary,
  best_HLA,
  by = c("orf_key", "netmhc_ID", "peptide"),
  all.x = TRUE,
  sort = FALSE
)


# Build QC summary

step3_qc <- data.table(
  metric = c(
    "MS ORF-peptide pairs",
    "HLA-I supported ORF-peptide pairs",
    "unique HLA-I supported peptides",
    "unique HLA-I supported ORFs",
    "best SB pairs",
    "best WB pairs",
    "pairs also matching non-pseudogene proteins",
    "duplicated final ORF-peptide keys"
  ),
  value = c(
    nrow(pseudo_source_peptide),
    nrow(pseudo_source_HLA_summary),
    uniqueN(pseudo_source_HLA_summary$peptide),
    uniqueN(pseudo_source_HLA_summary$orf_key),
    sum(
      pseudo_source_HLA_summary$best_binder_class == "SB"
    ),
    sum(
      pseudo_source_HLA_summary$best_binder_class == "WB"
    ),
    sum(
      pseudo_source_HLA_summary$
        has_non_pseudogene_protein
    ),
    pseudo_source_HLA_summary[
      ,
      .N,
      by = .(orf_key, peptide)
    ][N > 1, .N]
  )
)


# Save outputs

saveRDS(
  pseudo_source_HLA_summary,
  file.path(
    output_dir,
    "03a_pseudogene_MS_HLA_source_peptide_summary.rds"
  )
)

fwrite(
  pseudo_source_HLA_summary,
  file.path(
    output_dir,
    "03b_pseudogene_MS_HLA_source_peptide_summary.tsv"
  ),
  sep = "\t"
)

fwrite(
  step3_qc,
  file.path(
    output_dir,
    "03c_pseudogene_MS_HLA_QC.tsv"
  ),
  sep = "\t"
)


# 4. Classify candidates by exact matching against the complete search FASTA

candidate_peptides <- sort(
  unique(pseudo_source_HLA_summary$peptide)
)


# 4.1 Function for exact peptide matching in a FASTA file

search_peptides_in_fasta <- function(fasta_file, peptides) {

  aa <- readAAStringSet(fasta_file)

  fasta_header <- names(aa)
  fasta_id <- sub("\\s.*$", "", fasta_header)

  group_name <- sub(
    "\\.fasta$",
    "",
    sub(
      "^contain_human_",
      "",
      basename(fasta_file)
    )
  )

  rbindlist(
    lapply(peptides, function(pep) {

      hit_index <- which(
        vcountPattern(
          AAString(pep),
          aa,
          fixed = TRUE
        ) > 0
      )

      if (length(hit_index) == 0) {
        return(NULL)
      }

      hit_id <- fasta_id[hit_index]

      protein_type <- fcase(
        grepl(
          "^lcl\\|ORF[0-9]+_ENST",
          hit_id
        ),
        "pseudogene_ORF",

        grepl(
          "^sp\\|[^|]+-[0-9]+\\|",
          hit_id
        ),
        "SwissProt_isoform",

        grepl("^sp\\|", hit_id),
        "SwissProt_canonical",

        grepl("^tr\\|", hit_id),
        "TrEMBL",

        default = "other_search_fasta"
      )

      accession <- ifelse(
        grepl("^(sp|tr)\\|", hit_id),
        sub(
          "^[^|]+\\|([^|]+)\\|.*$",
          "\\1",
          hit_id
        ),
        sub(
          ":[0-9]+:[0-9]+$",
          "",
          sub("^lcl\\|", "", hit_id)
        )
      )

      data.table(
        peptide = pep,
        fasta_group = group_name,
        protein_id = hit_id,
        accession = accession,
        protein_type = protein_type,
        fasta_header = fasta_header[hit_index]
      )
    }),
    use.names = TRUE
  )
}


# 4.2 Search the three FASTA files used for MaxQuant
fasta_exact_hits <- rbindlist(
  lapply(
    fasta_files,
    search_peptides_in_fasta,
    peptides = candidate_peptides
  ),
  use.names = TRUE
)

fasta_exact_hits <- unique(
  fasta_exact_hits,
  by = c(
    "peptide",
    "fasta_group",
    "protein_id"
  )
)


# Collapse duplicate protein hits occurring in multiple comparison FASTA files

fasta_exact_hits_unique <- fasta_exact_hits[
  ,
  .(
    fasta_groups = paste(
      sort(unique(fasta_group)),
      collapse = ";"
    ),
    fasta_header = fasta_header[1]
  ),
  by = .(
    peptide,
    protein_id,
    accession,
    protein_type
  )
]


# 4.3 Generate peptide-level FASTA matching classes
peptide_fasta_classification <- fasta_exact_hits_unique[
  ,
  {
    non_pseudo <- protein_type != "pseudogene_ORF"

    .(
      pseudogene_ORF_match_count =
        uniqueN(accession[protein_type == "pseudogene_ORF"]),

      non_pseudogene_match_count =
        uniqueN(protein_id[non_pseudo]),

      exact_SwissProt_canonical =
        any(protein_type == "SwissProt_canonical"),

      exact_SwissProt_isoform =
        any(protein_type == "SwissProt_isoform"),

      exact_TrEMBL =
        any(protein_type == "TrEMBL"),

      exact_other_search_fasta =
        any(protein_type == "other_search_fasta"),

      matched_non_pseudogene_accessions =
        if (any(non_pseudo)) {
          paste(
            sort(unique(accession[non_pseudo])),
            collapse = ";"
          )
        } else {
          NA_character_
        }
    )
  },
  by = peptide
]


# Assign the final sequence-matching class

peptide_fasta_classification[
  ,
  full_fasta_class := fcase(
    exact_SwissProt_canonical,
    "shared_with_SwissProt_canonical",

    exact_SwissProt_isoform,
    "shared_with_SwissProt_isoform",

    exact_TrEMBL,
    "shared_with_TrEMBL",

    exact_other_search_fasta,
    "shared_with_other_search_fasta",

    default = "pseudogene_only_in_search_fasta"
  )
]


# 4.4 Merge sequence classification back to ORF-peptide pairs


pseudo_source_HLA_classified <- merge(
  pseudo_source_HLA_summary,
  peptide_fasta_classification,
  by = "peptide",
  all.x = TRUE,
  sort = FALSE
)

pseudo_source_HLA_classified[
  ,
  full_fasta_has_non_pseudogene_match :=
    full_fasta_class != "pseudogene_only_in_search_fasta"
]


# Build QC summary

step4_qc <- data.table(
  metric = c(
    "HLA-I supported ORF-peptide pairs",
    "unique peptides",
    "pairs with non-pseudogene FASTA match",
    "pairs pseudogene-only in search FASTA",
    "missing FASTA classification",
    "MaxQuant and full-FASTA disagreement pairs"
  ),
  value = c(
    nrow(pseudo_source_HLA_classified),
    uniqueN(pseudo_source_HLA_classified$peptide),

    sum(
      pseudo_source_HLA_classified$
        full_fasta_has_non_pseudogene_match,
      na.rm = TRUE
    ),

    sum(
      pseudo_source_HLA_classified$
        full_fasta_class ==
        "pseudogene_only_in_search_fasta",
      na.rm = TRUE
    ),

    sum(
      is.na(
        pseudo_source_HLA_classified$
          full_fasta_class
      )
    ),

    sum(
      pseudo_source_HLA_classified$
        has_non_pseudogene_protein !=
        pseudo_source_HLA_classified$
        full_fasta_has_non_pseudogene_match,
      na.rm = TRUE
    )
  )
)


# Save outputs


saveRDS(
  fasta_exact_hits_unique,
  file.path(
    output_dir,
    "04a_pseudogene_full_FASTA_exact_hits.rds"
  )
)

saveRDS(
  pseudo_source_HLA_classified,
  file.path(
    output_dir,
    "04b_pseudogene_MS_HLA_FASTA_classified.rds"
  )
)

fwrite(
  peptide_fasta_classification,
  file.path(
    output_dir,
    "04c_pseudogene_peptide_FASTA_classification.tsv"
  ),
  sep = "\t"
)

fwrite(
  step4_qc,
  file.path(
    output_dir,
    "04d_pseudogene_full_FASTA_QC.tsv"
  ),
  sep = "\t"
)


# 5. Generate the publication-oriented pseudogene candidate table


collapse_unique <- function(x) {
  x <- sort(unique(x[!is.na(x) & x != ""]))
  if (length(x) == 0) NA_character_ else paste(x, collapse = ";")
}


# 5.1 Add mass-error information from evidence.txt


evidence_all <- rbindlist(
  lapply(seq_len(nrow(group_map)), function(i) {

    evidence_file <- file.path(
      dirname(group_map$msms_file[i]),
      "evidence.txt"
    )

    x <- fread(
      evidence_file,
      check.names = TRUE,
      showProgress = FALSE
    )

    mass_error <- if ("Mass.error..ppm." %in% names(x)) {
      suppressWarnings(as.numeric(x$Mass.error..ppm.))
    } else {
      rep(NA_real_, nrow(x))
    }

    data.table(
      dataset = group_map$dataset[i],
      comparison_group = group_map$comparison_group[i],
      Evidence_ID = x$id,
      evidence_mass_error_ppm = mass_error
    )
  })
)

pseudo_PSM_evidence <- merge(
  pseudo_PSM,
  evidence_all,
  by = c(
    "dataset",
    "comparison_group",
    "Evidence_ID"
  ),
  all.x = TRUE,
  sort = FALSE
)


# 5.2 Summarize MS evidence at the peptide level

peptide_MS_summary <- pseudo_PSM_evidence[
  ,
  {
    lfq_error <- abs(
      evidence_mass_error_ppm[
        quantification == "LabelFree" &
          !is.na(evidence_mass_error_ppm)
      ]
    )

    has_lfq <- any(quantification == "LabelFree")

    .(
      n_PSM = uniqueN(PSM_uid),

      n_raw_files = uniqueN(
        paste(
          dataset,
          comparison_group,
          Raw_file,
          sep = "|"
        )
      ),

      datasets = collapse_unique(dataset),
      comparisons = collapse_unique(comparison_group),

      min_PEP = safe_min(PEP),
      max_Score = safe_max(Score),
      max_Delta_score = safe_max(Delta_score),
      max_PIF = safe_max(PIF),

      PIF_availability = ifelse(
        any(!is.na(PIF)),
        "available",
        "unavailable"
      ),

      median_LFQ_abs_mass_error_ppm =
        if (length(lfq_error) > 0) {
          median(lfq_error)
        } else {
          NA_real_
        },

      mass_error_reporting_status =
        if (!has_lfq) {
          "not_reported_TMT_only"
        } else if (length(lfq_error) > 0) {
          "reported_from_LFQ_evidence"
        } else {
          "LFQ_evidence_mass_error_unavailable"
        }
    )
  },
  by = .(
    peptide = Sequence
  )
]

# 5.3 Summarize pseudogene ORF sources

peptide_source <- unique(
  pseudo_source_HLA_classified[
    ,
    .(
      peptide,
      orf_key,
      full_fasta_class
    )
  ]
)

peptide_source[
  ,
  pseudogene_transcript_id :=
    sub("^ORF[0-9]+_", "", orf_key)
]

peptide_source[
  ,
  pseudogene_ORF_id :=
    sub("_ENST.*$", "", orf_key)
]

peptide_source_summary <- peptide_source[
  ,
  .(
    pseudogene_transcript_ids =
      collapse_unique(pseudogene_transcript_id),

    pseudogene_ORF_ids =
      collapse_unique(pseudogene_ORF_id),

    pseudogene_source_ids =
      collapse_unique(orf_key),

    n_pseudogene_ORFs = uniqueN(orf_key),

    pseudogene_mapping_class = ifelse(
      uniqueN(orf_key) == 1,
      "single_pseudogene_ORF",
      "multiple_pseudogene_ORFs"
    ),

    source_uniqueness_class = fcase(
      full_fasta_class[1] ==
        "pseudogene_only_in_search_fasta",
      "pseudogene_only_in_search_fasta",

      full_fasta_class[1] ==
        "shared_with_SwissProt_canonical",
      "shared_with_SwissProt",

      default = full_fasta_class[1]
    )
  ),
  by = peptide
]


# 5.4 Summarize HLA-I predictions at the peptide level

peptide_HLA_summary <- pseudo_MS_HLA[
  ,
  .(
    HLA_alleles = collapse_unique(HLA_allele),
    n_HLA_alleles = uniqueN(HLA_allele)
  ),
  by = peptide
]

best_peptide_HLA <- pseudo_MS_HLA[
  order(EL_Rank, -EL_score),
  .SD[1],
  by = peptide
][
  ,
  .(
    peptide,
    best_EL_allele = HLA_allele,
    best_EL_rank = EL_Rank,
    best_EL_score = EL_score,
    best_EL_binder_class = binder_class
  )
]


# 5.5 Build the publication table

pseudogene_publication_table <- merge(
  peptide_source_summary,
  peptide_MS_summary,
  by = "peptide",
  all.x = TRUE
)

pseudogene_publication_table <- merge(
  pseudogene_publication_table,
  peptide_HLA_summary,
  by = "peptide",
  all.x = TRUE
)

pseudogene_publication_table <- merge(
  pseudogene_publication_table,
  best_peptide_HLA,
  by = "peptide",
  all.x = TRUE
)


# Add interpretation fields for publication

pseudogene_publication_table[
  ,
  strict_MS_MHC_support := TRUE
]

pseudogene_publication_table[
  ,
  source_conclusion := fcase(
    source_uniqueness_class ==
      "pseudogene_only_in_search_fasta",

    paste0(
      "peptide sequence detected in a pseudogene ORF ",
      "with no exact non-pseudogene match in the ",
      "MaxQuant search FASTA"
    ),

    source_uniqueness_class ==
      "shared_with_SwissProt",

    paste0(
      "peptide sequence shared with Swiss-Prot; ",
      "shotgun MS cannot establish ",
      "pseudogene-specific origin"
    ),

    default = "source assignment requires cautious interpretation"
  )
]

pseudogene_publication_table[
  ,
  supports_natural_HLA_presentation := FALSE
]

pseudogene_publication_table[
  ,
  recommended_claim := fcase(
    source_uniqueness_class ==
      "pseudogene_only_in_search_fasta",

    paste0(
      "MS-supported and NetMHCpan-predicted HLA-I-binding ",
      "pseudogene ORF candidate peptide"
    ),

    default = paste0(
      "MS-supported and NetMHCpan-predicted HLA-I-binding ",
      "peptide mapped to a pseudogene ORF"
    )
  )
]

pseudogene_publication_table[
  ,
  manual_spectrum_review_priority := fcase(
    source_uniqueness_class ==
      "pseudogene_only_in_search_fasta" &
      best_EL_binder_class == "SB",
    "P1_review_pseudogene_candidate_SB",

    source_uniqueness_class ==
      "pseudogene_only_in_search_fasta" &
      best_EL_binder_class == "WB",
    "P2_review_pseudogene_candidate_WB",

    source_uniqueness_class ==
      "shared_with_SwissProt" &
      best_EL_binder_class == "SB",
    "P3_review_shared_source_SB",

    default = "P4_review_shared_source_WB"
  )
]


# Arrange publication-table columns
setcolorder(
  pseudogene_publication_table,
  c(
    "peptide",
    "pseudogene_transcript_ids",
    "pseudogene_ORF_ids",
    "pseudogene_source_ids",
    "pseudogene_mapping_class",
    "n_pseudogene_ORFs",

    "n_PSM",
    "n_raw_files",
    "datasets",
    "comparisons",
    "min_PEP",
    "max_Score",
    "max_Delta_score",
    "max_PIF",
    "PIF_availability",
    "median_LFQ_abs_mass_error_ppm",
    "mass_error_reporting_status",

    "HLA_alleles",
    "n_HLA_alleles",
    "best_EL_allele",
    "best_EL_rank",
    "best_EL_score",
    "best_EL_binder_class",
    "strict_MS_MHC_support",

    "source_uniqueness_class",
    "source_conclusion",
    "supports_natural_HLA_presentation",
    "recommended_claim",
    "manual_spectrum_review_priority"
  )
)

setorder(
  pseudogene_publication_table,
  source_uniqueness_class,
  best_EL_rank,
  -max_Score
)


# 5.6 Save final results
stopifnot(
  nrow(pseudogene_publication_table) == 293,
  uniqueN(pseudogene_publication_table$peptide) == 293,
  sum(
    pseudogene_publication_table$
      source_uniqueness_class ==
      "pseudogene_only_in_search_fasta"
  ) == 1
)

final_file <- file.path(
  output_dir,
  "05a_293_pseudogene_publication_candidate_summary_HLA_fixed.tsv"
)

core_file <- file.path(
  output_dir,
  "05b_pseudogene_only_candidate.tsv"
)

fwrite(
  pseudogene_publication_table,
  final_file,
  sep = "\t",
  na = "NA"
)

fwrite(
  pseudogene_publication_table[
    source_uniqueness_class ==
      "pseudogene_only_in_search_fasta"
  ],
  core_file,
  sep = "\t",
  na = "NA"
)


# 6. Split candidate peptides into UniProt query batches

pseduo_293 <- read_delim("./MS_unspecific/pseudogenes/result/contain_pseudogenes/pseudogene_integrated_analysis/05a_293_pseudogene_publication_candidate_summary_HLA_fixed.tsv",
                         delim = "\t", escape_double = FALSE,
                         trim_ws = TRUE)


# Extract, clean, and deduplicate peptide sequences
peptides <- unique(
  toupper(trimws(pseduo_293$peptide))
)

peptides <- peptides[
  !is.na(peptides) &
    peptides != ""
]


# Split peptides into batches of 100

peptide_batches <- split(
  peptides,
  ceiling(seq_along(peptides) / 100)
)


# Save each batch as a comma-separated single-line text file

output_dir <- "./MS_unspecific/pseudogenes/result/contain_pseudogenes/IEDB_uniprot_search"

dir.create(
  output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

for (i in seq_along(peptide_batches)) {

  peptide_text <- paste(
    peptide_batches[[i]],
    collapse = ","
  )

  writeLines(
    peptide_text,
    file.path(
      output_dir,
      sprintf("pseduo_peptide_batch_%02d.txt", i)
    )
  )
}


# 7. Process UniProt peptide-search results
# Parse UniProt peptide-search output
uniprot_result <- function(input_file, output_dir) {

  # Read the current UniProt Excel result file
  lncRNA_uniprot <- readxl::read_excel(input_file)

  table_copy <- lncRNA_uniprot %>%
    dplyr::rename(
      Entry_Name   = `Entry Name`,
      Protein_Names = `Protein Names`,
      Gene_Names   = `Gene Names`
    ) %>%
    tidyr::fill(
      Entry,
      Entry_Name,
      Protein_Names,
      Gene_Names,
      Organism,
      Length,
      .direction = "down"
    ) %>%
    dplyr::filter(
      !is.na(Match),
      Match != "NA",
      stringr::str_detect(Match, ":")
    ) %>%
    dplyr::mutate(
      Peptide = stringr::str_trim(
        stringr::str_replace(Match, ".*:\\s*", "")
      ),
      Match_Type = dplyr::if_else(
        stringr::str_detect(Entry, "-"),
        "isoform",
        "canonical"
      )
    )

  peptide_detail <- table_copy %>%
    dplyr::select(
      Peptide,
      Entry,
      Entry_Name,
      Protein_Names,
      Gene_Names,
      Organism,
      Match_Type
    ) %>%
    dplyr::distinct()

  peptide_class <- peptide_detail %>%
    dplyr::group_by(Peptide) %>%
    dplyr::summarise(
      has_canonical = any(Match_Type == "canonical"),
      has_isoform   = any(Match_Type == "isoform"),

      Final_Class = dplyr::case_when(
        has_canonical & !has_isoform ~ "canonical",
        !has_canonical & has_isoform ~ "isoform",
        has_canonical & has_isoform  ~ "both",
        TRUE                         ~ "unknown"
      ),

      Matched_Entries = paste(unique(Entry), collapse = "; "),
      Matched_Genes   = paste(unique(Gene_Names), collapse = "; "),
      n_matches       = dplyr::n(),
      .groups = "drop"
    )

  # Derive the sample name from the input filename
  sample_name <- tools::file_path_sans_ext(
    basename(input_file)
  )

  output_file <- file.path(
    output_dir,
    paste0(sample_name, "_output.csv")
  )

  readr::write_csv(
    peptide_class,
    output_file,
    na = ""
  )

  message(
    "Completed: ", basename(input_file),
    " -> ", output_file,
    "; peptides: ", dplyr::n_distinct(peptide_detail$Peptide)
  )

  return(peptide_class)
}
# Run UniProt result parsing
input_dir <- "./MS_unspecific/pseudogenes/result/contain_pseudogenes/IEDB_uniprot_search"


output_dir <- "./MS_unspecific/pseudogenes/result/contain_pseudogenes/IEDB_uniprot_search"

sample_names <- c("pseduo_peptide_batch_01","pseduo_peptide_batch_02","pseduo_peptide_batch_03")

full_file_paths <- file.path(
  input_dir,
  paste0(sample_names, ".xlsx")
)

results <- lapply(
  full_file_paths,
  uniprot_result,
  output_dir = output_dir
)


pseduo_batch_03_output <- read.csv("./MS_unspecific/pseudogenes/result/contain_pseudogenes/IEDB_uniprot_search/pseduo_peptide_batch_03_output.csv")
pseduo_batch_02_output <- read.csv("./MS_unspecific/pseudogenes/result/contain_pseudogenes/IEDB_uniprot_search/pseduo_peptide_batch_02_output.csv")
pseduo_batch_01_output <- read.csv("./MS_unspecific/pseudogenes/result/contain_pseudogenes/IEDB_uniprot_search/pseduo_peptide_batch_01_output.csv")
pseduo_uniprot_output <- rbind(pseduo_batch_03_output,pseduo_batch_02_output,pseduo_batch_01_output)
setDT(pseduo_uniprot_output)

sp_tr_map <- read.delim("./MS_unspecific/sp_tr_map.txt")

setDT(sp_tr_map)

# Add row identifiers
pseduo_uniprot_output[, row_id := .I]

# 7.1 Expand Matched_Entries
circRNA_entries_long <- pseduo_uniprot_output[
  ,
  .(
    accession = trimws(
      unlist(
        strsplit(Matched_Entries, ";", fixed = TRUE)
      )
    )
  ),
  by = row_id
]

circRNA_entries_long[
  ,
  accession := sub(
    "^(sp|tr)\\|",
    "",
    accession
  )
]

circRNA_entries_long[
  ,
  parent_accession := sub(
    "-[0-9]+$",
    "",
    accession
  )
]

# 7.2 Collect accessions requiring database annotation
needed_entries <- unique(
  c(
    circRNA_entries_long$accession,
    circRNA_entries_long$parent_accession
  )
)


# 7.3 Subset the large accession map to relevant entries

uniprot_map_small <- sp_tr_map[
  Entry %chin% needed_entries,
  .(
    Entry,
    database
  )
]


uniprot_map_small <- unique(
  uniprot_map_small,
  by = "Entry"
)

# 7.4 Match exact accessions
circRNA_entries_long[
  uniprot_map_small,
  on = .(accession = Entry),
  database_exact := i.database
]

# 7.5 Annotate isoforms using parent accessions when needed

circRNA_entries_long[
  uniprot_map_small,
  on = .(parent_accession = Entry),
  database_parent := i.database
]

circRNA_entries_long[
  ,
  database_final := fifelse(
    !is.na(database_exact),
    database_exact,
    database_parent
  )
]
# Add sp| or tr| prefixes
circRNA_entries_long[
  ,
  accession_with_db := fifelse(
    !is.na(database_final),
    paste0(database_final, "|", accession),
    accession
  )
]

# 7.6 Classify accession matching type
circRNA_entries_long[
  ,
  match_type := fifelse(
    !is.na(database_exact),
    "exact_accession",
    fifelse(
      is.na(database_exact) &
        !is.na(database_parent) &
        grepl("-[0-9]+$", accession),
      "matched_by_parent",
      "not_found"
    )
  )
]
# 7.7 Collapse accession annotations back to row level

circRNA_entries_summary <- circRNA_entries_long[
  ,
  .(
    Matched_Entries_with_db = paste(
      accession_with_db,
      collapse = "; "
    ),

    has_sp_match = any(
      database_final == "sp",
      na.rm = TRUE
    ),

    has_tr_match = any(
      database_final == "tr",
      na.rm = TRUE
    ),

    all_entries_annotated = all(
      match_type != "not_found"
    ),

    unmatched_entries = {
      x <- unique(
        accession[match_type == "not_found"]
      )

      if (length(x) == 0) {
        NA_character_
      } else {
        paste(x, collapse = "; ")
      }
    }
  ),
  by = row_id
]

# 7.8 Merge annotations back to the UniProt result table
all_circrNA_search_with_db <- merge(
  pseduo_uniprot_output,
  circRNA_entries_summary,
  by = "row_id",
  all.x = TRUE,
  sort = FALSE
)

setorder(
  all_circrNA_search_with_db,
  row_id
)

all_circrNA_search_with_db[
  ,
  row_id := NULL
]


# Replace the original Matched_Entries column
all_circrNA_search_with_db[
  ,
  Matched_Entries_original := Matched_Entries
]

all_circrNA_search_with_db[
  ,
  Matched_Entries := Matched_Entries_with_db
]

all_circrNA_search_with_db[
  ,
  Matched_Entries_with_db := NULL
]


write.table(all_circrNA_search_with_db,"./MS_unspecific/pseudogenes/result/contain_pseudogenes/IEDB_uniprot_search/all_pseduo_search_sp_tr.txt",quote = FALSE,row.names = FALSE,sep = "\t")


# 8. Generate the publication-oriented UniProt annotation table


# Summarize UniProt matches by peptide

tidy_uniprot <- function(df) {

  df %>%
    dplyr::select(
      Peptide,
      Matched_Entries,
      Matched_Genes
    ) %>%

    # Expand to one UniProt entry per row
    tidyr::separate_rows(
      Matched_Entries,
      sep = ";\\s*"
    ) %>%

    dplyr::mutate(
      Matched_Entries = stringr::str_trim(
        Matched_Entries
      ),

      entry_class = dplyr::case_when(
        stringr::str_detect(
          Matched_Entries,
          "^sp\\|.+-[0-9]+$"
        ) ~ "sp_isoform",

        stringr::str_detect(
          Matched_Entries,
          "^sp\\|"
        ) ~ "sp_canonical",

        stringr::str_detect(
          Matched_Entries,
          "^tr\\|"
        ) ~ "tr",

        TRUE ~ "unknown"
      )
    ) %>%

    dplyr::distinct(
      Peptide,
      Matched_Entries,
      .keep_all = TRUE
    ) %>%

    dplyr::group_by(Peptide) %>%

    dplyr::summarise(
      has_sp_canonical = any(
        entry_class == "sp_canonical"
      ),

      has_sp_isoform = any(
        entry_class == "sp_isoform"
      ),

      has_tr_match = any(
        entry_class == "tr"
      ),

      n_sp_canonical_matches = sum(
        entry_class == "sp_canonical"
      ),

      n_sp_isoform_matches = sum(
        entry_class == "sp_isoform"
      ),

      n_tr_matches = sum(
        entry_class == "tr"
      ),

      sp_canonical_entries = paste(
        Matched_Entries[
          entry_class == "sp_canonical"
        ],
        collapse = "; "
      ),

      sp_isoform_entries = paste(
        Matched_Entries[
          entry_class == "sp_isoform"
        ],
        collapse = "; "
      ),

      tr_entries = paste(
        Matched_Entries[
          entry_class == "tr"
        ],
        collapse = "; "
      ),

      matched_genes = paste(
        unique(
          Matched_Genes[
            !is.na(Matched_Genes) &
              Matched_Genes != ""
          ]
        ),
        collapse = "; "
      ),

      .groups = "drop"
    ) %>%

    dplyr::mutate(
      n_total_UniProt_matches =
        n_sp_canonical_matches +
        n_sp_isoform_matches +
        n_tr_matches,

      UniProt_match_class = dplyr::case_when(
        has_sp_canonical &
          has_sp_isoform &
          has_tr_match ~
          "SP_canonical+SP_isoform+TrEMBL",

        has_sp_canonical &
          has_sp_isoform ~
          "SP_canonical+SP_isoform",

        has_sp_canonical &
          has_tr_match ~
          "SP_canonical+TrEMBL",

        has_sp_isoform &
          has_tr_match ~
          "SP_isoform+TrEMBL",

        has_sp_canonical ~
          "SP_canonical",

        has_sp_isoform ~
          "SP_isoform",

        has_tr_match ~
          "TrEMBL_only",

        TRUE ~
          "no_UniProt_match"
      )
    ) %>%

    dplyr::select(
      peptide = Peptide,
      UniProt_match_class,
      has_sp_canonical,
      has_sp_isoform,
      has_tr_match,
      n_sp_canonical_matches,
      n_sp_isoform_matches,
      n_tr_matches,
      n_total_UniProt_matches,
      sp_canonical_entries,
      sp_isoform_entries,
      tr_entries,
      matched_genes
    )
}


# Run final UniProt summarization
all_pseduo_search_sp_tr <- read.delim2("~/MS_unspecific/pseudogenes/result/contain_pseudogenes/IEDB_uniprot_search/all_pseduo_search_sp_tr.txt")

pseduo_uniprot_publication  <- tidy_uniprot(
  all_pseduo_search_sp_tr
)


write.table(pseduo_uniprot_publication ,"./MS_unspecific/pseudogenes/result/contain_pseudogenes/IEDB_uniprot_search/new_all_pseduo_search_sp_tr.txt",quote = FALSE,row.names = FALSE,sep = "\t")

write.csv(pseduo_uniprot_publication,"./MS_unspecific/pseudogenes/result/contain_pseudogenes/IEDB_uniprot_search/new_all_pseduo_search_sp_tr.csv",quote = FALSE,row.names = FALSE,sep = ",")
