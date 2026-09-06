# ============================================================
# Alternative-splicing peptide MS and HLA-I integration
# ============================================================

library(data.table)
library(dplyr)
library(tidyr)
library(stringr)
library(readxl)
library(readr)

# ============================================================
# 1. Input paths and sample metadata
# ============================================================

MQ_ROOT <- paste0(
  "./MS_unspecific/",
  "AS/result/contain_jcast"
)

MAP_FILE <- paste0(
  "./circpeptide4/AS_rMATS_jcast/",
  "jcast/result/T1_and_T2_info/all_shortID_map.txt"
)

NETMHC_FILE <- paste0(
  "./circpeptide4/AS_rMATS_jcast/",
  "netMHCpan/HLA_1/combined_f.txt"
)

OUTDIR <- file.path(
  MQ_ROOT,
  "AS_integrated_analysis"
)

dir.create(
  OUTDIR,
  recursive = TRUE,
  showWarnings = FALSE
)

group_map <- data.frame(
  Dataset = c(
    "PXD037581",
    "PXD037581",
    "PXD037581",
    "PXD044963",
    "PXD044963",
    "PXD044963"
  ),

  AS_group = c(
    "HC_eRA",
    "HC_RA",
    "RA_eRA",
    "HC_eRA",
    "HC_RA",
    "RA_eRA"
  ),

  Comparison = c(
    "HC_vs_eRA",
    "HC_vs_estRA",
    "estRA_vs_eRA",
    "HC_vs_eRA",
    "HC_vs_estRA",
    "estRA_vs_eRA"
  ),

  Acquisition_type = c(
    "TMT10plex",
    "TMT10plex",
    "TMT10plex",
    "LabelFree",
    "LabelFree",
    "LabelFree"
  ),

  stringsAsFactors = FALSE
)

group_map$txt_dir <- file.path(
  MQ_ROOT,
  group_map$Dataset,
  paste0(group_map$AS_group, "_result"),
  "combined",
  "txt"
)

group_map$msms_file <- file.path(
  group_map$txt_dir,
  "msms.txt"
)

group_map$evidence_file <- file.path(
  group_map$txt_dir,
  "evidence.txt"
)

# ============================================================
# 2. JCAST source dictionary
# ============================================================

AS_map <- data.table::fread(
  MAP_FILE,
  sep = "\t",
  header = TRUE,
  check.names = FALSE,
  data.table = FALSE,
  colClasses = "character",
  na.strings = c("", "NA", "NaN")
)

required_map_cols <- c(
  "short_ID",
  "header",
  "seq",
  "kb",
  "uniprot",
  "uniprot_name",
  "gene_ensg",
  "type_order",
  "event_id",
  "chr",
  "sjc_raw",
  "sjc",
  "tier",
  "event_type"
)

missing_map_cols <- setdiff(
  required_map_cols,
  names(AS_map)
)

if (length(missing_map_cols) == 0) {
} else {
  stop("The source mapping table is missing required columns.")
}

AS_map$MaxQuant_source_id <- sub(
  "\\s.*$",
  "",
  AS_map$header
)

AS_map$AS_group <- sub(
  "_[0-9]+$",
  "",
  AS_map$short_ID
)

AS_map$event_isoform_number <- ifelse(
  grepl("[12]$", AS_map$type_order),
  sub("^.*([12])$", "\\1", AS_map$type_order),
  NA_character_
)

# ============================================================
# 3. High-confidence AS-related PSMs
# ============================================================

msms_keep_cols <- c(
  "Raw file",
  "Scan number",
  "Scan index",
  "Sequence",
  "Length",
  "Missed cleavages",
  "Modifications",
  "Modified sequence",
  "Proteins",
  "Charge",
  "Type",
  "Mass error [ppm]",
  "Mass error [Da]",
  "Simple mass error [ppm]",
  "Retention time",
  "PEP",
  "Score",
  "Delta score",
  "PIF",
  "Reverse",
  "Contaminant",
  "Number of matches",
  "Intensity coverage",
  "Peak coverage",
  "id",
  "Evidence ID"
)

read_one_msms <- function(i) {

  x <- data.table::fread(
    group_map$msms_file[i],
    select = msms_keep_cols,
    check.names = FALSE,
    data.table = TRUE,
    na.strings = c(
      "",
      "NA",
      "NaN"
    ),
    showProgress = TRUE
  )

  x[, Dataset := group_map$Dataset[i]]
  x[, AS_group := group_map$AS_group[i]]
  x[, Comparison := group_map$Comparison[i]]
  x[, Acquisition_type := group_map$Acquisition_type[i]]

  x[, PSM_uid := paste(
    Dataset,
    AS_group,
    id,
    sep = "||"
  )]

  return(x)
}

all_msms <- data.table::rbindlist(
  lapply(
    seq_len(nrow(group_map)),
    read_one_msms
  ),
  use.names = TRUE,
  fill = TRUE
)

numeric_cols <- c(
  "Scan number",
  "Scan index",
  "Length",
  "Missed cleavages",
  "Charge",
  "Mass error [ppm]",
  "Mass error [Da]",
  "Simple mass error [ppm]",
  "Retention time",
  "PEP",
  "Score",
  "Delta score",
  "PIF",
  "Number of matches",
  "Intensity coverage",
  "Peak coverage",
  "id",
  "Evidence ID"
)

numeric_cols <- intersect(
  numeric_cols,
  names(all_msms)
)

all_msms[
  ,
  (numeric_cols) := lapply(
    .SD,
    function(z) suppressWarnings(as.numeric(z))
  ),
  .SDcols = numeric_cols
]

character_cols <- c(
  "Raw file",
  "Sequence",
  "Modifications",
  "Modified sequence",
  "Proteins",
  "Type",
  "Reverse",
  "Contaminant",
  "Dataset",
  "AS_group",
  "Comparison",
  "Acquisition_type",
  "PSM_uid"
)

character_cols <- intersect(
  character_cols,
  names(all_msms)
)

all_msms[
  ,
  (character_cols) := lapply(
    .SD,
    as.character
  ),
  .SDcols = character_cols
]

is_positive_flag <- function(x) {

  y <- toupper(
    trimws(
      as.character(x)
    )
  )

  !is.na(y) &
    y %in% c(
      "+",
      "TRUE",
      "T",
      "1"
    )
}

all_msms[
  ,
  reverse_flag := is_positive_flag(Reverse)
]

all_msms[
  ,
  contaminant_flag := is_positive_flag(Contaminant)
]

AS_ID_PATTERN <- paste0(
  "ENSG[0-9]+\\.[0-9]+\\|",
  "(A3SS|A5SS|MXE|RI|SE)[12]\\|"
)

all_msms[
  ,
  has_JCAST_AS_source :=
    !is.na(Proteins) &
    grepl(
      AS_ID_PATTERN,
      Proteins,
      perl = TRUE
    )
]

all_msms[
  ,
  sequence_length_mismatch :=
    !is.na(Sequence) &
    !is.na(Length) &
    nchar(Sequence) != Length
]

all_msms[
  ,
  pass_01_clean :=
    !reverse_flag &
    !contaminant_flag
]

all_msms[
  ,
  pass_02_length :=
    pass_01_clean &
    !is.na(Length) &
    Length >= 7
]

all_msms[
  ,
  pass_03_PEP :=
    pass_02_length &
    !is.na(PEP) &
    PEP < 0.01
]

all_msms[
  ,
  pass_04_Score :=
    pass_03_PEP &
    !is.na(Score) &
    Score > 40
]

all_msms[
  ,
  pass_05_Delta :=
    pass_04_Score &
    !is.na(`Delta score`) &
    `Delta score` >= 15
]

all_msms[
  ,
  pass_06_JCAST_source :=
    pass_05_Delta &
    has_JCAST_AS_source
]

filter_attrition <- all_msms[
  ,
  .(
    n_total_PSM = .N,

    n_non_reverse_non_contaminant =
      sum(
        pass_01_clean,
        na.rm = TRUE
      ),

    n_length_ge_7 =
      sum(
        pass_02_length,
        na.rm = TRUE
      ),

    n_PEP_lt_0.01 =
      sum(
        pass_03_PEP,
        na.rm = TRUE
      ),

    n_Score_gt_40 =
      sum(
        pass_04_Score,
        na.rm = TRUE
      ),

    n_Delta_ge_15 =
      sum(
        pass_05_Delta,
        na.rm = TRUE
      ),

    n_high_conf_with_JCAST_AS_source =
      sum(
        pass_06_JCAST_source,
        na.rm = TRUE
      ),

    n_unique_high_conf_peptides =
      uniqueN(
        Sequence[
          pass_06_JCAST_source
        ]
      ),

    n_sequence_length_mismatch =
      sum(
        sequence_length_mismatch,
        na.rm = TRUE
      )
  ),
  by = .(
    Dataset,
    AS_group,
    Comparison,
    Acquisition_type
  )
]
AS_PSM_high_conf <- all_msms[
  pass_06_JCAST_source == TRUE
]

setorder(
  AS_PSM_high_conf,
  Dataset,
  AS_group,
  `Raw file`,
  `Scan number`,
  Sequence
)

data.table::fwrite(
  filter_attrition,
  file.path(OUTDIR, "04a_AS_PSM_filter_attrition.tsv"),
  sep = "\t", quote = FALSE, na = "NA"
)

data.table::fwrite(
  AS_PSM_high_conf,
  file.path(OUTDIR, "04b_AS_PSM_high_conf_quality_and_source.tsv"),
  sep = "\t", quote = FALSE, na = "NA"
)

saveRDS(
  AS_PSM_high_conf,
  file.path(OUTDIR, "04b_AS_PSM_high_conf_quality_and_source.rds"),
  compress = TRUE
)

# ============================================================
# 4. Strict JCAST source mapping
# ============================================================

psm_core_cols <- c(
  "PSM_uid",
  "Dataset",
  "AS_group",
  "Comparison",
  "Acquisition_type",
  "Raw file",
  "Scan number",
  "Scan index",
  "Sequence",
  "Length",
  "Missed cleavages",
  "Modifications",
  "Modified sequence",
  "Proteins",
  "Charge",
  "Type",
  "Mass error [ppm]",
  "Mass error [Da]",
  "Simple mass error [ppm]",
  "Retention time",
  "PEP",
  "Score",
  "Delta score",
  "PIF",
  "Number of matches",
  "Intensity coverage",
  "Peak coverage",
  "id",
  "Evidence ID"
)

psm_core_cols <- intersect(
  psm_core_cols,
  names(AS_PSM_high_conf)
)

AS_PSM_core <- copy(
  AS_PSM_high_conf[
    ,
    ..psm_core_cols
  ]
)

protein_list <- strsplit(
  AS_PSM_core$Proteins,
  split = ";",
  fixed = TRUE
)

n_sources_per_PSM <- lengths(
  protein_list
)

if (any(n_sources_per_PSM == 0)) {
  stop("PSMs without protein-source assignments were detected.")
}

AS_PSM_source_all <- AS_PSM_core[
  rep(
    seq_len(nrow(AS_PSM_core)),
    n_sources_per_PSM
  )
]

AS_PSM_source_all[
  ,
  MaxQuant_source_id :=
    trimws(
      unlist(
        protein_list,
        use.names = FALSE
      )
    )
]

AS_PSM_source_all <- AS_PSM_source_all[
  !is.na(MaxQuant_source_id) &
    MaxQuant_source_id != ""
]

AS_PSM_source_all <- unique(
  AS_PSM_source_all,
  by = c(
    "PSM_uid",
    "MaxQuant_source_id"
  )
)

AS_PSM_source_all[
  ,
  is_JCAST_AS_source :=
    grepl(
      AS_ID_PATTERN,
      MaxQuant_source_id,
      perl = TRUE
    )
]
PSM_source_context <- AS_PSM_source_all[
  ,
  .(
    n_all_protein_sources =
      uniqueN(MaxQuant_source_id),

    n_JCAST_AS_sources =
      uniqueN(
        MaxQuant_source_id[
          is_JCAST_AS_source
        ]
      ),

    n_non_JCAST_sources =
      uniqueN(
        MaxQuant_source_id[
          !is_JCAST_AS_source
        ]
      ),

    has_non_JCAST_source =
      any(
        !is_JCAST_AS_source
      ),

    has_multiple_JCAST_sources =
      uniqueN(
        MaxQuant_source_id[
          is_JCAST_AS_source
        ]
      ) > 1
  ),
  by = .(
    PSM_uid,
    Dataset,
    AS_group
  )
]
AS_PSM_JCAST <- AS_PSM_source_all[
  is_JCAST_AS_source == TRUE
]

AS_PSM_JCAST <- merge(
  AS_PSM_JCAST,
  PSM_source_context,
  by = c(
    "PSM_uid",
    "Dataset",
    "AS_group"
  ),
  all.x = TRUE,
  sort = FALSE
)

AS_map_key <- as.data.table(
  AS_map[
    ,
    c(
      "AS_group",
      "MaxQuant_source_id",
      "short_ID",
      "header",
      "seq",
      "kb",
      "uniprot",
      "uniprot_name",
      "gene_ensg",
      "type_order",
      "event_isoform_number",
      "event_id",
      "chr",
      "sjc_raw",
      "sjc",
      "tier",
      "event_type"
    )
  ]
)

setnames(
  AS_map_key,
  old = c(
    "header",
    "seq"
  ),
  new = c(
    "JCAST_header",
    "source_protein_sequence"
  )
)

map_key_duplicate <- AS_map_key[
  ,
  .N,
  by = .(
    AS_group,
    MaxQuant_source_id
  )
][
  N > 1
]

if (nrow(map_key_duplicate) > 0) {
  stop(
    "AS_group + MaxQuant_source_id is not unique in the source mapping table."
  )
}

AS_map_global <- unique(
  AS_map_key[
    ,
    .(
      MaxQuant_source_id,
      source_group_in_map = AS_group
    )
  ],
  by = "MaxQuant_source_id"
)

AS_PSM_JCAST[
  ,
  source_group_in_map :=
    AS_map_global$source_group_in_map[
      match(
        MaxQuant_source_id,
        AS_map_global$MaxQuant_source_id
      )
    ]
]

AS_PSM_source_mapped <- merge(
  AS_PSM_JCAST,
  AS_map_key,
  by = c(
    "AS_group",
    "MaxQuant_source_id"
  ),
  all.x = TRUE,
  sort = FALSE
)
AS_PSM_source_mapped[
  ,
  source_mapping_status :=
    fifelse(
      is.na(source_group_in_map),
      "source_not_in_dictionary",

      fifelse(
        source_group_in_map != AS_group,
        "AS_group_mismatch",

        fifelse(
          is.na(short_ID),
          "exact_join_failed",
          "mapped_exact"
        )
      )
    )
]
AS_PSM_source_mapped[
  ,
  peptide_sequence :=
    toupper(
      trimws(
        Sequence
      )
    )
]

AS_PSM_source_mapped[
  ,
  source_protein_sequence :=
    toupper(
      gsub(
        "\\s+",
        "",
        source_protein_sequence
      )
    )
]

AS_PSM_source_mapped[
  ,
  peptide_in_source_sequence :=
    mapply(
      function(peptide, protein_sequence) {

        if (
          is.na(peptide) ||
          peptide == "" ||
          is.na(protein_sequence) ||
          protein_sequence == ""
        ) {
          return(FALSE)
        }

        grepl(
          peptide,
          protein_sequence,
          fixed = TRUE
        )
      },

      peptide_sequence,
      source_protein_sequence
    )
]

AS_PSM_source_mapped[
  ,
  valid_JCAST_source_row :=
    source_mapping_status == "mapped_exact" &
    peptide_in_source_sequence
]

AS_PSM_source_valid <- AS_PSM_source_mapped[
  valid_JCAST_source_row == TRUE
]

AS_PSM_source_invalid <- AS_PSM_source_mapped[
  valid_JCAST_source_row == FALSE
]

setorder(
  AS_PSM_source_valid,
  Dataset,
  AS_group,
  `Raw file`,
  `Scan number`,
  Sequence,
  short_ID
)

data.table::fwrite(
  AS_PSM_source_valid,
  file.path(OUTDIR, "05d_AS_PSM_JCAST_source_valid.tsv"),
  sep = "\t", quote = FALSE, na = "NA"
)

saveRDS(
  AS_PSM_source_valid,
  file.path(OUTDIR, "05d_AS_PSM_JCAST_source_valid.rds"),
  compress = TRUE
)

# ============================================================
# 5. AS event-pair dictionary
# ============================================================

AS_dictionary <- as.data.table(
  AS_map
)

header_parts <- strsplit(
  AS_dictionary$MaxQuant_source_id,
  split = "|",
  fixed = TRUE
)

token_count <- lengths(
  header_parts
)

if (any(token_count != 12)) {
  stop(
    "Unexpected JCAST accession structure: expected 12 pipe-delimited fields."
  )
}

header_matrix <- do.call(
  rbind,
  header_parts
)

header_info <- as.data.table(
  header_matrix
)

setnames(
  header_info,
  c(
    "db_prefix",
    "header_uniprot",
    "header_uniprot_name",
    "header_gene_ensg",
    "header_type_order",
    "header_event_id",
    "header_chr",
    "anchor_region",
    "alternative_region",
    "strand_phase",
    "header_sjc_raw",
    "header_tier"
  )
)

AS_dictionary <- cbind(
  AS_dictionary,
  header_info
)

if (any(token_count != 12)) {
  stop("Unexpected JCAST accession structure: expected 12 pipe-delimited fields.")
}
AS_dictionary[
  ,
  AS_event_key := paste(
    AS_group,
    gene_ensg,
    event_type,
    event_id,
    chr,
    anchor_region,
    alternative_region,
    strand_phase,
    sjc_raw,
    tier,
    sep = "||"
  )
]

event_isoform_duplicate <- AS_dictionary[
  ,
  .N,
  by = .(
    AS_event_key,
    event_isoform_number
  )
][
  N > 1
]

if (nrow(event_isoform_duplicate) > 0) {

  stop(
    "AS_event_key + event_isoform_number is not unique."
  )
}

AS_dictionary[
  ,
  target_partner_isoform_number :=
    fifelse(
      event_isoform_number == "1",
      "2",
      "1"
    )
]

partner_lookup <- AS_dictionary[
  ,
  .(
    AS_event_key,

    partner_isoform_number =
      event_isoform_number,

    partner_short_ID =
      short_ID,

    partner_type_order =
      type_order,

    partner_protein_sequence =
      seq
  )
]

AS_source_dictionary <- merge(
  AS_dictionary,
  partner_lookup,

  by.x = c(
    "AS_event_key",
    "target_partner_isoform_number"
  ),

  by.y = c(
    "AS_event_key",
    "partner_isoform_number"
  ),

  all.x = TRUE,
  sort = FALSE
)

AS_source_dictionary[
  ,
  source_event_pair_status :=
    fifelse(
      !is.na(partner_short_ID),
      "complete_1_2_pair",
      "partner_isoform_unavailable"
    )
]

data.table::fwrite(
  AS_source_dictionary,
  file.path(OUTDIR, "06d_AS_source_dictionary_with_partner.tsv"),
  sep = "\t", quote = FALSE, na = "NA"
)

saveRDS(
  AS_source_dictionary,
  file.path(OUTDIR, "06d_AS_source_dictionary_with_partner.rds"),
  compress = TRUE
)

# ============================================================
# 6. AS isoform sequence specificity and MS support
# ============================================================

AS_source_peptide_base <- unique(
  AS_PSM_source_valid[
    ,
    .(
      short_ID,
      AS_group,
      MaxQuant_source_id,
      peptide = toupper(
        trimws(Sequence)
      ),
      source_protein_sequence =
        toupper(
          gsub(
            "\\s+",
            "",
            source_protein_sequence
          )
        )
    )
  ],
  by = c(
    "short_ID",
    "peptide"
  )
)

partner_info <- AS_source_dictionary[
  ,
  .(
    short_ID,

    AS_event_key,

    gene_ensg,

    event_type,

    event_id,

    type_order,

    event_isoform_number,

    tier,

    source_event_pair_status,

    partner_short_ID,

    partner_type_order,

    partner_protein_sequence =
      toupper(
        gsub(
          "\\s+",
          "",
          partner_protein_sequence
        )
      )
  )
]

partner_duplicate <- partner_info[
  ,
  .N,
  by = short_ID
][
  N > 1
]

if (nrow(partner_duplicate) > 0) {
  stop("short_ID is not unique in partner_info.")
}

AS_source_peptide_specificity <- merge(
  AS_source_peptide_base,
  partner_info,
  by = "short_ID",
  all.x = TRUE,
  sort = FALSE
)

if (
  nrow(AS_source_peptide_specificity) !=
  nrow(AS_source_peptide_base)
) {
  stop("Row count changed after joining the event-pair dictionary.")
}

AS_source_peptide_specificity[
  ,
  peptide_in_source_sequence :=
    mapply(
      function(peptide, protein_sequence) {

        if (
          is.na(peptide) ||
          peptide == "" ||
          is.na(protein_sequence) ||
          protein_sequence == ""
        ) {
          return(FALSE)
        }

        grepl(
          peptide,
          protein_sequence,
          fixed = TRUE
        )
      },
      peptide,
      source_protein_sequence
    )
]

AS_source_peptide_specificity[
  ,
  peptide_in_partner_sequence :=
    mapply(
      function(
    peptide,
    partner_sequence,
    pair_status
      ) {

        if (
          is.na(pair_status) ||
          pair_status != "complete_1_2_pair"
        ) {
          return(NA)
        }

        if (
          is.na(partner_sequence) ||
          partner_sequence == ""
        ) {
          return(NA)
        }

        grepl(
          peptide,
          partner_sequence,
          fixed = TRUE
        )
      },
    peptide,
    partner_protein_sequence,
    source_event_pair_status
    )
]

AS_source_peptide_specificity[
  ,
  AS_isoform_specificity_class :=
    fcase(

      !peptide_in_source_sequence,
      "sequence_mapping_unresolved",

      source_event_pair_status ==
        "partner_isoform_unavailable",
      "partner_isoform_unavailable",

      source_event_pair_status ==
        "complete_1_2_pair" &
        peptide_in_partner_sequence == FALSE,
      "AS_isoform_specific_to_source_model",

      source_event_pair_status ==
        "complete_1_2_pair" &
        peptide_in_partner_sequence == TRUE,
      "shared_between_event_isoforms",

      default =
        "sequence_mapping_unresolved"
    )
]

collapse_unique <- function(x) {

  x <- as.character(x)

  x <- x[
    !is.na(x) &
      x != ""
  ]

  if (length(x) == 0) {
    return(NA_character_)
  }

  paste(
    sort(
      unique(x)
    ),
    collapse = ";"
  )
}

AS_source_peptide_PSM_summary <-
  AS_PSM_source_valid[
    ,
    .(
      n_PSM =
        uniqueN(PSM_uid),

      n_raw_files =
        uniqueN(`Raw file`),

      raw_files =
        collapse_unique(`Raw file`),

      datasets =
        collapse_unique(Dataset),

      comparisons =
        collapse_unique(Comparison),

      acquisition_types =
        collapse_unique(Acquisition_type),

      min_PEP =
        if (
          all(is.na(PEP))
        ) {
          NA_real_
        } else {
          min(
            PEP,
            na.rm = TRUE
          )
        },

      max_Score =
        if (
          all(is.na(Score))
        ) {
          NA_real_
        } else {
          max(
            Score,
            na.rm = TRUE
          )
        },

      max_Delta_score =
        if (
          all(is.na(`Delta score`))
        ) {
          NA_real_
        } else {
          max(
            `Delta score`,
            na.rm = TRUE
          )
        },

      max_PIF =
        if (
          all(is.na(PIF))
        ) {
          NA_real_
        } else {
          max(
            PIF,
            na.rm = TRUE
          )
        },

      n_PIF_available =
        sum(
          !is.na(PIF)
        ),

      has_non_JCAST_source =
        any(
          has_non_JCAST_source,
          na.rm = TRUE
        ),

      max_n_non_JCAST_sources =
        if (
          all(is.na(n_non_JCAST_sources))
        ) {
          NA_integer_
        } else {
          max(
            n_non_JCAST_sources,
            na.rm = TRUE
          )
        },

      max_n_JCAST_AS_sources =
        if (
          all(is.na(n_JCAST_AS_sources))
        ) {
          NA_integer_
        } else {
          max(
            n_JCAST_AS_sources,
            na.rm = TRUE
          )
        }
    ),
    by = .(
      short_ID,
      peptide = Sequence
    )
  ]

AS_source_peptide_PSM_summary[
  ,
  peptide :=
    toupper(
      trimws(peptide)
    )
]

AS_MS_source_peptide <- merge(
  AS_source_peptide_specificity,
  AS_source_peptide_PSM_summary,
  by = c(
    "short_ID",
    "peptide"
  ),
  all.x = TRUE,
  sort = FALSE
)

data.table::fwrite(
  AS_MS_source_peptide,
  file.path(OUTDIR, "07d_AS_MS_source_peptide.tsv"),
  sep = "\t", quote = FALSE, na = "NA"
)

saveRDS(
  AS_MS_source_peptide,
  file.path(OUTDIR, "07d_AS_MS_source_peptide.rds"),
  compress = TRUE
)

# ============================================================
# 7. HLA-I prediction table
# ============================================================

netmhc_header <- data.table::fread(
  NETMHC_FILE,
  nrows = 0,
  check.names = FALSE,
  data.table = FALSE
)

netmhc_names <- trimws(
  names(netmhc_header)
)

find_first_column <- function(candidates, available) {

  hit <- candidates[
    candidates %in% available
  ]

  if (length(hit) == 0) {
    return(NA_character_)
  }

  hit[1]
}

netmhc_column_map <- c(
  Pos = find_first_column(
    c("Pos", "pos"),
    netmhc_names
  ),

  peptide = find_first_column(
    c("Peptide", "peptide"),
    netmhc_names
  ),

  short_ID = find_first_column(
    c("ID", "id"),
    netmhc_names
  ),

  HLA_allele = find_first_column(
    c("type", "Type", "Allele", "allele"),
    netmhc_names
  ),

  core = find_first_column(
    c("_core", "core"),
    netmhc_names
  ),

  icore = find_first_column(
    c("_icore", "icore"),
    netmhc_names
  ),

  EL_score = find_first_column(
    c(
      "_EL-score",
      "EL-score",
      "_EL_score",
      "EL_score"
    ),
    netmhc_names
  ),

  EL_Rank = find_first_column(
    c(
      "_EL_Rank",
      "EL_Rank",
      "_EL-rank",
      "EL-rank"
    ),
    netmhc_names
  ),

  BA_score = find_first_column(
    c(
      "_BA-score",
      "BA-score",
      "_BA_score",
      "BA_score"
    ),
    netmhc_names
  ),

  BA_Rank = find_first_column(
    c(
      "_BA_Rank",
      "BA_Rank",
      "_BA-rank",
      "BA-rank"
    ),
    netmhc_names
  )
)

if (any(is.na(netmhc_column_map))) {

  stop(
    "The following NetMHCpan fields were not recognized: ",
    paste(
      names(netmhc_column_map)[
        is.na(netmhc_column_map)
      ],
      collapse = ", "
    )
  )
}

netmhc_raw <- data.table::fread(
  NETMHC_FILE,
  select = unname(netmhc_column_map),
  check.names = FALSE,
  data.table = TRUE,
  na.strings = c(
    "",
    "NA",
    "NaN"
  ),
  showProgress = TRUE
)

data.table::setnames(
  netmhc_raw,
  old = unname(netmhc_column_map),
  new = names(netmhc_column_map)
)

netmhc_raw[
  ,
  peptide :=
    toupper(
      trimws(
        as.character(peptide)
      )
    )
]

netmhc_raw[
  ,
  short_ID :=
    trimws(
      as.character(short_ID)
    )
]

netmhc_raw[
  ,
  HLA_allele :=
    trimws(
      as.character(HLA_allele)
    )
]

netmhc_raw[
  ,
  core :=
    trimws(
      as.character(core)
    )
]

netmhc_raw[
  ,
  icore :=
    trimws(
      as.character(icore)
    )
]

numeric_netmhc_cols <- c(
  "Pos",
  "EL_score",
  "EL_Rank",
  "BA_score",
  "BA_Rank"
)

netmhc_raw[
  ,
  (numeric_netmhc_cols) := lapply(
    .SD,
    function(x) {
      suppressWarnings(
        as.numeric(x)
      )
    }
  ),
  .SDcols = numeric_netmhc_cols
]

netmhc_raw[
  ,
  binder_class :=
    data.table::fcase(

      is.na(EL_Rank),
      "UNRESOLVED",

      EL_Rank < 0.5,
      "SB",

      EL_Rank < 2,
      "WB",

      default = "NB"
    )
]
AS_sequence_lookup <- unique(
  as.data.table(AS_map)[
    ,
    .(
      short_ID,

      AS_group,

      gene_ensg,

      type_order,

      event_id,

      event_type,

      tier,

      source_protein_sequence =
        toupper(
          gsub(
            "\\s+",
            "",
            seq
          )
        )
    )
  ],
  by = "short_ID"
)

netmhc_checked <- merge(
  netmhc_raw,
  AS_sequence_lookup,
  by = "short_ID",
  all.x = TRUE,
  sort = FALSE
)

netmhc_checked[
  ,
  peptide_in_source_sequence :=
    mapply(
      function(peptide, protein_sequence) {

        if (
          is.na(peptide) ||
          peptide == "" ||
          is.na(protein_sequence) ||
          protein_sequence == ""
        ) {
          return(FALSE)
        }

        grepl(
          peptide,
          protein_sequence,
          fixed = TRUE
        )
      },

      peptide,
      source_protein_sequence
    )
]

netmhc_I_exact_unique <- unique(
  netmhc_checked[
    short_ID != "ID" &
      !is.na(short_ID) &
      short_ID != "" &
      !is.na(peptide) &
      peptide != "" &
      !is.na(HLA_allele) &
      HLA_allele != "" &
      !is.na(EL_Rank) &
      !is.na(source_protein_sequence) &
      peptide_in_source_sequence == TRUE,

    .(
      short_ID,
      peptide,
      HLA_allele,
      Pos,
      core,
      icore,
      EL_score,
      EL_Rank,
      BA_score,
      BA_Rank,
      binder_class
    )
  ]
)
data.table::setorder(
  netmhc_I_exact_unique,
  short_ID,
  peptide,
  HLA_allele,
  EL_Rank,
  -EL_score,
  BA_Rank,
  Pos
)

netmhc_I_best <- unique(
  netmhc_I_exact_unique,
  by = c(
    "short_ID",
    "peptide",
    "HLA_allele"
  )
)

saveRDS(
  netmhc_I_best,
  file.path(OUTDIR, "08g_HLA_classI_best_prediction_by_source_peptide_allele.rds"),
  compress = FALSE
)

# ============================================================
# 8. Strict MS-HLA-I integration
# ============================================================

if (!exists("AS_MS_source_peptide")) {
  AS_MS_source_peptide <- readRDS(
    file.path(
      OUTDIR,
      "07d_AS_MS_source_peptide.rds"
    )
  )
}

if (!exists("netmhc_I_best")) {
  netmhc_I_best <- readRDS(
    file.path(
      OUTDIR,
      "08g_HLA_classI_best_prediction_by_source_peptide_allele.rds"
    )
  )
}

setDT(AS_MS_source_peptide)
setDT(netmhc_I_best)

heavy_sequence_cols <- intersect(
  c(
    "source_protein_sequence",
    "partner_protein_sequence"
  ),
  names(AS_MS_source_peptide)
)

AS_MS_for_join <- copy(
  AS_MS_source_peptide[
    ,
    setdiff(
      names(AS_MS_source_peptide),
      heavy_sequence_cols
    ),
    with = FALSE
  ]
)

HLA_I_for_join <- copy(
  netmhc_I_best[
    ,
    .(
      short_ID,
      peptide,
      HLA_allele,
      Pos,
      core,
      icore,
      EL_score,
      EL_Rank,
      BA_score,
      BA_Rank,
      binder_class
    )
  ]
)

AS_MS_for_join[
  ,
  short_ID :=
    trimws(
      as.character(short_ID)
    )
]

AS_MS_for_join[
  ,
  peptide :=
    toupper(
      trimws(
        as.character(peptide)
      )
    )
]

HLA_I_for_join[
  ,
  short_ID :=
    trimws(
      as.character(short_ID)
    )
]

HLA_I_for_join[
  ,
  peptide :=
    toupper(
      trimws(
        as.character(peptide)
      )
    )
]

HLA_I_for_join[
  ,
  HLA_allele :=
    trimws(
      as.character(HLA_allele)
    )
]

non_class_I_alleles <- unique(
  HLA_I_for_join[
    !grepl(
      "^HLA-[ABC]",
      HLA_allele
    ),
    HLA_allele
  ]
)

if (length(non_class_I_alleles) > 0) {
  stop("Non-HLA-I alleles were detected in the HLA-I prediction table.")
}
MS_duplicate_key <- AS_MS_for_join[
  ,
  .N,
  by = .(
    short_ID,
    peptide
  )
][
  N > 1
]

HLA_duplicate_key <- HLA_I_for_join[
  ,
  .N,
  by = .(
    short_ID,
    peptide,
    HLA_allele
  )
][
  N > 1
]

if (
  nrow(MS_duplicate_key) > 0 ||
  nrow(HLA_duplicate_key) > 0
) {
  stop(
    "Duplicate join keys remain before MS-HLA-I integration."
  )
}
AS_strict_MS_HLA_I <- merge(
  AS_MS_for_join,
  HLA_I_for_join,
  by = c(
    "short_ID",
    "peptide"
  ),
  all = FALSE,
  sort = FALSE,
  allow.cartesian = TRUE
)

AS_strict_MS_HLA_I[
  ,
  strict_MS_MHC_support := TRUE
]

AS_strict_MS_HLA_I[
  ,
  supports_natural_HLA_presentation := FALSE
]
invalid_binder_rows <- AS_strict_MS_HLA_I[
  !binder_class %in%
    c(
      "SB",
      "WB"
    ) |
    is.na(EL_Rank) |
    EL_Rank >= 2
]

if (nrow(invalid_binder_rows) > 0) {
  stop(
    "Non-SB/WB records were detected after strict MS-HLA-I integration."
  )
}
setorder(
  AS_strict_MS_HLA_I,
  short_ID,
  peptide,
  EL_Rank,
  -EL_score,
  BA_Rank,
  HLA_allele,
  Pos,
  na.last = TRUE
)

AS_strict_MS_HLA_I_best <- unique(
  AS_strict_MS_HLA_I,
  by = c(
    "short_ID",
    "peptide"
  )
)

data.table::fwrite(
  AS_strict_MS_HLA_I_best,
  file.path(OUTDIR, "09e_AS_strict_MS_HLA_classI_best_per_source_peptide.tsv"),
  sep = "\t", quote = FALSE, na = "NA"
)

saveRDS(
  AS_strict_MS_HLA_I,
  file.path(OUTDIR, "09f_AS_strict_MS_HLA_classI_by_allele.rds"),
  compress = TRUE
)

saveRDS(
  AS_strict_MS_HLA_I_best,
  file.path(OUTDIR, "09e_AS_strict_MS_HLA_classI_best_per_source_peptide.rds"),
  compress = TRUE
)

# ============================================================
# 9. Source-peptide HLA-I summary
# ============================================================

if (!exists("AS_strict_MS_HLA_I")) {
  AS_strict_MS_HLA_I <- readRDS(
    file.path(
      OUTDIR,
      "09f_AS_strict_MS_HLA_classI_by_allele.rds"
    )
  )
}

if (!exists("AS_MS_for_join")) {

  AS_MS_source_peptide <- readRDS(
    file.path(
      OUTDIR,
      "07d_AS_MS_source_peptide.rds"
    )
  )

  setDT(AS_MS_source_peptide)

  heavy_sequence_cols <- intersect(
    c(
      "source_protein_sequence",
      "partner_protein_sequence"
    ),
    names(AS_MS_source_peptide)
  )

  AS_MS_for_join <- copy(
    AS_MS_source_peptide[
      ,
      setdiff(
        names(AS_MS_source_peptide),
        heavy_sequence_cols
      ),
      with = FALSE
    ]
  )
}

setDT(AS_strict_MS_HLA_I)
setDT(AS_MS_for_join)

strict_hla_sorted <- copy(
  AS_strict_MS_HLA_I
)

setorder(
  strict_hla_sorted,
  short_ID,
  peptide,
  EL_Rank,
  -EL_score,
  BA_Rank,
  HLA_allele,
  Pos,
  na.last = TRUE
)

collapse_ordered_unique <- function(x) {

  x <- as.character(x)

  x <- x[
    !is.na(x) &
      x != ""
  ]

  if (length(x) == 0) {
    return(NA_character_)
  }

  paste(
    unique(x),
    collapse = ";"
  )
}

HLA_I_source_peptide_summary <-
  strict_hla_sorted[
    ,
    .(
      n_HLA_I_alleles =
        uniqueN(HLA_allele),

      HLA_I_alleles =
        collapse_ordered_unique(
          HLA_allele
        ),

      n_SB_alleles =
        uniqueN(
          HLA_allele[
            binder_class == "SB"
          ]
        ),

      SB_alleles =
        collapse_ordered_unique(
          HLA_allele[
            binder_class == "SB"
          ]
        ),

      n_WB_alleles =
        uniqueN(
          HLA_allele[
            binder_class == "WB"
          ]
        ),

      WB_alleles =
        collapse_ordered_unique(
          HLA_allele[
            binder_class == "WB"
          ]
        ),

      best_EL_allele =
        HLA_allele[1],

      best_EL_rank =
        EL_Rank[1],

      best_EL_score =
        EL_score[1],

      best_EL_binder_class =
        binder_class[1],

      best_BA_rank =
        BA_Rank[1],

      best_BA_score =
        BA_score[1],

      best_core =
        core[1],

      best_icore =
        icore[1],

      best_prediction_position =
        Pos[1],

      strict_MS_MHC_support =
        TRUE,

      supports_natural_HLA_presentation =
        FALSE
    ),
    by = .(
      short_ID,
      peptide
    )
  ]

AS_MS_HLA_I_source_summary <- merge(
  AS_MS_for_join,
  HLA_I_source_peptide_summary,
  by = c(
    "short_ID",
    "peptide"
  ),
  all = FALSE,
  sort = FALSE
)

AS_MS_HLA_I_source_summary[
  ,
  primary_translation_tier :=
    tier == "T1"
]

AS_MS_HLA_I_source_summary[
  ,
  MaxQuant_source_background :=
    fifelse(
      has_non_JCAST_source,
      "JCAST_plus_non_JCAST_sources",
      "JCAST_only_in_search_database"
    )
]

data.table::fwrite(
  AS_MS_HLA_I_source_summary,
  file.path(OUTDIR, "10d_AS_MS_HLA_classI_source_peptide_summary.tsv"),
  sep = "\t", quote = FALSE, na = "NA"
)

saveRDS(
  AS_MS_HLA_I_source_summary,
  file.path(OUTDIR, "10d_AS_MS_HLA_classI_source_peptide_summary.rds"),
  compress = TRUE
)

# ============================================================
# 10. Peptide-level AS-MS-HLA-I master table
# ============================================================

if (!exists("AS_MS_HLA_I_source_summary")) {
  AS_MS_HLA_I_source_summary <- readRDS(
    file.path(
      OUTDIR,
      "10d_AS_MS_HLA_classI_source_peptide_summary.rds"
    )
  )
}

if (!exists("AS_strict_MS_HLA_I")) {
  AS_strict_MS_HLA_I <- readRDS(
    file.path(
      OUTDIR,
      "09f_AS_strict_MS_HLA_classI_by_allele.rds"
    )
  )
}

if (!exists("AS_PSM_source_valid")) {
  AS_PSM_source_valid <- readRDS(
    file.path(
      OUTDIR,
      "05d_AS_PSM_JCAST_source_valid.rds"
    )
  )
}

if (!exists("AS_source_dictionary")) {
  AS_source_dictionary <- readRDS(
    file.path(
      OUTDIR,
      "06d_AS_source_dictionary_with_partner.rds"
    )
  )
}

setDT(AS_MS_HLA_I_source_summary)
setDT(AS_strict_MS_HLA_I)
setDT(AS_PSM_source_valid)
setDT(AS_source_dictionary)
collapse_sorted_unique <- function(x) {

  x <- as.character(x)

  x <- x[
    !is.na(x) &
      trimws(x) != ""
  ]

  if (length(x) == 0) {
    return(NA_character_)
  }

  paste(
    sort(unique(x)),
    collapse = ";"
  )
}

collapse_ordered_unique <- function(x) {

  x <- as.character(x)

  x <- x[
    !is.na(x) &
      trimws(x) != ""
  ]

  if (length(x) == 0) {
    return(NA_character_)
  }

  paste(
    unique(x),
    collapse = ";"
  )
}

source_metadata <- unique(
  AS_source_dictionary[
    ,
    .(
      short_ID,
      source_kb = kb,
      source_uniprot = uniprot,
      source_uniprot_name = uniprot_name
    )
  ],
  by = "short_ID"
)

AS_source_level_for_peptide <- merge(
  AS_MS_HLA_I_source_summary,
  source_metadata,
  by = "short_ID",
  all.x = TRUE,
  sort = FALSE
)

AS_peptide_source_summary <-
  AS_source_level_for_peptide[
    ,
    .(
      n_AS_source_peptide_pairs =
        .N,

      n_AS_short_IDs =
        uniqueN(short_ID),

      AS_short_IDs =
        collapse_sorted_unique(short_ID),

      AS_groups =
        collapse_sorted_unique(AS_group),

      AS_event_keys =
        collapse_sorted_unique(AS_event_key),

      n_AS_events =
        uniqueN(AS_event_key),

      gene_ensg =
        collapse_sorted_unique(gene_ensg),

      source_uniprot_accessions =
        collapse_sorted_unique(source_uniprot),

      source_uniprot_names =
        collapse_sorted_unique(source_uniprot_name),

      event_types =
        collapse_sorted_unique(event_type),

      type_orders =
        collapse_sorted_unique(type_order),

      tiers =
        collapse_sorted_unique(tier),

      n_T1_source_pairs =
        sum(
          tier == "T1",
          na.rm = TRUE
        ),

      n_T2_source_pairs =
        sum(
          tier == "T2",
          na.rm = TRUE
        ),

      n_isoform_specific_source_pairs =
        sum(
          AS_isoform_specificity_class ==
            "AS_isoform_specific_to_source_model"
        ),

      n_shared_between_isoform_source_pairs =
        sum(
          AS_isoform_specificity_class ==
            "shared_between_event_isoforms"
        ),

      n_partner_unavailable_source_pairs =
        sum(
          AS_isoform_specificity_class ==
            "partner_isoform_unavailable"
        ),

      AS_isoform_specificity_classes =
        collapse_sorted_unique(
          AS_isoform_specificity_class
        ),

      has_AS_isoform_specific_source_model =
        any(
          AS_isoform_specificity_class ==
            "AS_isoform_specific_to_source_model"
        ),

      any_JCAST_only_in_search_database =
        any(
          !has_non_JCAST_source,
          na.rm = TRUE
        ),

      all_source_pairs_have_non_JCAST_source =
        all(
          has_non_JCAST_source
        )
    ),
    by = peptide
  ]

AS_peptide_source_summary[
  ,
  peptide_AS_specificity_summary :=
    fcase(

      has_AS_isoform_specific_source_model,
      "has_AS_isoform_specific_source_model",

      n_shared_between_isoform_source_pairs ==
        n_AS_source_peptide_pairs,
      "shared_between_event_isoforms_only",

      n_partner_unavailable_source_pairs ==
        n_AS_source_peptide_pairs,
      "partner_isoform_unavailable_only",

      default =
        "mixed_non_specific_source_context"
    )
]

AS_peptide_source_summary[
  ,
  MaxQuant_source_background :=
    fcase(

      all_source_pairs_have_non_JCAST_source,
      "all_AS_source_pairs_shared_with_non_JCAST",

      any_JCAST_only_in_search_database,
      "at_least_one_JCAST_only_source_pair_in_search_database",

      default =
        "source_background_unresolved"
    )
]

strict_source_peptide_keys <- unique(
  AS_MS_HLA_I_source_summary[
    ,
    .(
      short_ID,
      peptide
    )
  ],
  by = c(
    "short_ID",
    "peptide"
  )
)

AS_PSM_for_peptide <- copy(
  AS_PSM_source_valid
)

AS_PSM_for_peptide[
  ,
  peptide :=
    toupper(
      trimws(Sequence)
    )
]
strict_supported_PSM_sources <- merge(
  AS_PSM_for_peptide,
  strict_source_peptide_keys,
  by = c(
    "short_ID",
    "peptide"
  ),
  all = FALSE,
  sort = FALSE,
  allow.cartesian = TRUE
)

strict_supported_PSM_unique <- unique(
  strict_supported_PSM_sources,
  by = c(
    "PSM_uid",
    "peptide"
  )
)

AS_peptide_PSM_summary <-
  strict_supported_PSM_unique[
    ,
    .(
      n_PSM =
        uniqueN(PSM_uid),

      n_raw_files =
        uniqueN(`Raw file`),

      raw_files =
        collapse_sorted_unique(`Raw file`),

      datasets =
        collapse_sorted_unique(Dataset),

      comparisons =
        collapse_sorted_unique(Comparison),

      acquisition_types =
        collapse_sorted_unique(Acquisition_type),

      min_PEP =
        if (all(is.na(PEP))) {
          NA_real_
        } else {
          min(
            PEP,
            na.rm = TRUE
          )
        },

      max_Score =
        if (all(is.na(Score))) {
          NA_real_
        } else {
          max(
            Score,
            na.rm = TRUE
          )
        },

      max_Delta_score =
        if (all(is.na(`Delta score`))) {
          NA_real_
        } else {
          max(
            `Delta score`,
            na.rm = TRUE
          )
        },

      max_PIF =
        if (all(is.na(PIF))) {
          NA_real_
        } else {
          max(
            PIF,
            na.rm = TRUE
          )
        },

      n_PIF_available =
        sum(
          !is.na(PIF)
        ),

      PIF_availability =
        fifelse(
          sum(!is.na(PIF)) > 0,
          "available",
          "not_available"
        )
    ),
    by = peptide
  ]

HLA_peptide_sorted <- copy(
  AS_strict_MS_HLA_I
)

setorder(
  HLA_peptide_sorted,
  peptide,
  EL_Rank,
  -EL_score,
  BA_Rank,
  HLA_allele,
  short_ID,
  Pos,
  na.last = TRUE
)
AS_peptide_HLA_summary <-
  HLA_peptide_sorted[
    ,
    .(
      n_source_HLA_prediction_pairs =
        .N,

      n_HLA_I_alleles =
        uniqueN(HLA_allele),

      HLA_I_alleles =
        collapse_ordered_unique(HLA_allele),

      n_SB_prediction_pairs =
        sum(
          binder_class == "SB"
        ),

      n_WB_prediction_pairs =
        sum(
          binder_class == "WB"
        ),

      n_SB_alleles =
        uniqueN(
          HLA_allele[
            binder_class == "SB"
          ]
        ),

      SB_alleles =
        collapse_ordered_unique(
          HLA_allele[
            binder_class == "SB"
          ]
        ),

      n_WB_alleles =
        uniqueN(
          HLA_allele[
            binder_class == "WB"
          ]
        ),

      WB_alleles =
        collapse_ordered_unique(
          HLA_allele[
            binder_class == "WB"
          ]
        ),

      best_EL_source_short_ID =
        short_ID[1],

      best_EL_allele =
        HLA_allele[1],

      best_EL_rank =
        EL_Rank[1],

      best_EL_score =
        EL_score[1],

      best_EL_binder_class =
        binder_class[1],

      best_BA_rank =
        BA_Rank[1],

      best_BA_score =
        BA_score[1],

      best_core =
        core[1],

      best_icore =
        icore[1],

      strict_MS_MHC_support =
        TRUE,

      supports_natural_HLA_presentation =
        FALSE
    ),
    by = peptide
  ]

AS_peptide_master_base <- merge(
  AS_peptide_source_summary,
  AS_peptide_PSM_summary,
  by = "peptide",
  all.x = TRUE,
  sort = FALSE
)

AS_peptide_master_base <- merge(
  AS_peptide_master_base,
  AS_peptide_HLA_summary,
  by = "peptide",
  all.x = TRUE,
  sort = FALSE
)

AS_peptide_master_base[
  ,
  translation_tier_support :=
    fcase(

      n_T1_source_pairs > 0 &
        n_T2_source_pairs == 0,
      "T1_only",

      n_T1_source_pairs == 0 &
        n_T2_source_pairs > 0,
      "T2_only",

      n_T1_source_pairs > 0 &
        n_T2_source_pairs > 0,
      "T1_and_T2",

      default =
        "tier_unresolved"
    )
]

AS_peptide_master_base[
  ,
  recommended_claim_pre_UniProt :=
    fifelse(
      has_AS_isoform_specific_source_model,

      paste0(
        "MS-supported, source-consistent predicted HLA-I-binding ",
        "peptide with sequence specificity to at least one modeled ",
        "JCAST AS isoform"
      ),

      paste0(
        "MS-supported, source-consistent predicted HLA-I-binding ",
        "peptide mapped to a JCAST AS-derived protein sequence"
      )
    )
]

AS_peptide_master_base[
  ,
  source_conclusion_pre_UniProt :=
    fifelse(
      all_source_pairs_have_non_JCAST_source,

      paste0(
        "All strict JCAST source mappings also had non-JCAST ",
        "protein matches in the MaxQuant search database"
      ),

      paste0(
        "At least one strict JCAST source mapping had no non-JCAST ",
        "co-match in the MaxQuant Proteins field; independent ",
        "UniProt sequence verification is still required"
      )
    )
]

data.table::fwrite(
  AS_peptide_master_base,
  file.path(OUTDIR, "11c_AS_peptide_master_base.tsv"),
  sep = "\t", quote = FALSE, na = "NA"
)

saveRDS(
  AS_peptide_master_base,
  file.path(OUTDIR, "11c_AS_peptide_master_base.rds"),
  compress = TRUE
)

# ============================================================
# 11. Evidence-level annotation
# ============================================================

evidence_keep_cols <- c(
  "Sequence",
  "Length",
  "Proteins",
  "Leading proteins",
  "Leading razor protein",
  "Raw file",
  "Charge",
  "Mass error [ppm]",
  "Mass error [Da]",
  "Uncalibrated mass error [ppm]",
  "PIF",
  "PEP",
  "MS/MS count",
  "MS/MS scan number",
  "Score",
  "Delta score",
  "Reverse",
  "Potential contaminant",
  "id"
)

read_one_evidence <- function(i) {

  x <- fread(
    group_map$evidence_file[i],
    select = evidence_keep_cols,
    check.names = FALSE,
    data.table = TRUE,
    na.strings = c("", "NA", "NaN"),
    showProgress = TRUE
  )

  x[, Dataset := group_map$Dataset[i]]
  x[, AS_group := group_map$AS_group[i]]
  x[, Comparison := group_map$Comparison[i]]
  x[, Acquisition_type := group_map$Acquisition_type[i]]

  x[
    ,
    Evidence_uid := fifelse(
      is.na(id),
      NA_character_,
      paste(
        Dataset,
        AS_group,
        as.integer(id),
        sep = "||"
      )
    )
  ]

  x
}

all_evidence <- rbindlist(
  lapply(
    seq_len(nrow(group_map)),
    read_one_evidence
  ),
  use.names = TRUE,
  fill = TRUE
)

evidence_duplicate_uid <- all_evidence[
  !is.na(Evidence_uid),
  .N,
  by = Evidence_uid
][
  N > 1
]

if (nrow(evidence_duplicate_uid) > 0) {

  stop("Evidence_uid is not unique.")
}

numeric_evidence_cols <- c(
  "Length",
  "Charge",
  "Mass error [ppm]",
  "Mass error [Da]",
  "Uncalibrated mass error [ppm]",
  "PIF",
  "PEP",
  "MS/MS count",
  "MS/MS scan number",
  "Score",
  "Delta score",
  "id"
)

numeric_evidence_cols <- intersect(
  numeric_evidence_cols,
  names(all_evidence)
)

all_evidence[
  ,
  (numeric_evidence_cols) := lapply(
    .SD,
    function(x) {
      suppressWarnings(
        as.numeric(x)
      )
    }
  ),
  .SDcols = numeric_evidence_cols
]

evidence_annotation <- all_evidence[
  ,
  .(
    Evidence_uid,

    evidence_sequence =
      toupper(
        trimws(Sequence)
      ),

    evidence_raw_file =
      `Raw file`,

    evidence_proteins =
      Proteins,

    leading_proteins =
      `Leading proteins`,

    leading_razor_protein =
      `Leading razor protein`,

    evidence_mass_error_ppm =
      `Mass error [ppm]`,

    evidence_mass_error_Da =
      `Mass error [Da]`,

    evidence_uncalibrated_mass_error_ppm =
      `Uncalibrated mass error [ppm]`,

    evidence_PIF =
      PIF,

    evidence_PEP =
      PEP,

    evidence_Score =
      Score,

    evidence_Delta_score =
      `Delta score`
  )
]

if (!exists("AS_MS_HLA_I_source_summary")) {

  AS_MS_HLA_I_source_summary <- readRDS(
    file.path(
      OUTDIR,
      "10d_AS_MS_HLA_classI_source_peptide_summary.rds"
    )
  )
}

if (!exists("AS_PSM_source_valid")) {

  AS_PSM_source_valid <- readRDS(
    file.path(
      OUTDIR,
      "05d_AS_PSM_JCAST_source_valid.rds"
    )
  )
}

if (!exists("AS_peptide_master_base")) {

  AS_peptide_master_base <- readRDS(
    file.path(
      OUTDIR,
      "11c_AS_peptide_master_base.rds"
    )
  )
}

setDT(AS_MS_HLA_I_source_summary)
setDT(AS_PSM_source_valid)
setDT(AS_peptide_master_base)

strict_source_peptide_keys <- unique(
  AS_MS_HLA_I_source_summary[
    ,
    .(
      short_ID,
      peptide
    )
  ],
  by = c(
    "short_ID",
    "peptide"
  )
)

strict_PSM_source_base <- copy(
  AS_PSM_source_valid
)

strict_PSM_source_base[
  ,
  peptide :=
    toupper(
      trimws(Sequence)
    )
]

strict_PSM_source_base[
  ,
  Evidence_uid := fifelse(
    is.na(`Evidence ID`),
    NA_character_,
    paste(
      Dataset,
      AS_group,
      as.integer(`Evidence ID`),
      sep = "||"
    )
  )
]

strict_PSM_sources <- merge(
  strict_PSM_source_base,
  strict_source_peptide_keys,
  by = c(
    "short_ID",
    "peptide"
  ),
  all = FALSE,
  sort = FALSE,
  allow.cartesian = TRUE
)

strict_PSM_evidence_source <- merge(
  strict_PSM_sources,
  evidence_annotation,
  by = "Evidence_uid",
  all.x = TRUE,
  sort = FALSE
)

strict_PSM_evidence_source[
  ,
  evidence_record_found :=
    !is.na(evidence_sequence)
]

strict_PSM_evidence_source[
  ,
  evidence_sequence_match :=
    evidence_record_found &
    peptide == evidence_sequence
]

strict_PSM_evidence_source[
  ,
  leading_razor_available :=
    !is.na(leading_razor_protein) &
    trimws(leading_razor_protein) != ""
]

strict_PSM_evidence_source[
  ,
  leading_razor_is_JCAST_AS :=
    leading_razor_available &
    grepl(
      AS_ID_PATTERN,
      leading_razor_protein,
      perl = TRUE
    )
]

strict_PSM_evidence_source[
  ,
  leading_razor_matches_this_JCAST_source :=
    leading_razor_available &
    leading_razor_protein ==
    MaxQuant_source_id
]

first_non_missing_numeric <- function(x) {

  x <- x[
    !is.na(x)
  ]

  if (length(x) == 0) {
    return(NA_real_)
  }

  as.numeric(x[1])
}

first_non_missing_character <- function(x) {

  x <- as.character(x)

  x <- x[
    !is.na(x) &
      trimws(x) != ""
  ]

  if (length(x) == 0) {
    return(NA_character_)
  }

  x[1]
}

strict_PSM_evidence <- strict_PSM_evidence_source[
  ,
  .(
    evidence_record_found =
      any(evidence_record_found),

    evidence_sequence_match =
      if (any(evidence_record_found)) {
        all(
          evidence_sequence_match[
            evidence_record_found
          ]
        )
      } else {
        FALSE
      },

    leading_razor_protein =
      first_non_missing_character(
        leading_razor_protein
      ),

    leading_razor_available =
      any(leading_razor_available),

    leading_razor_is_JCAST_AS =
      any(leading_razor_is_JCAST_AS),

    leading_razor_matches_any_strict_JCAST_source =
      any(
        leading_razor_matches_this_JCAST_source
      ),

    evidence_mass_error_ppm =
      first_non_missing_numeric(
        evidence_mass_error_ppm
      ),

    evidence_mass_error_Da =
      first_non_missing_numeric(
        evidence_mass_error_Da
      ),

    evidence_uncalibrated_mass_error_ppm =
      first_non_missing_numeric(
        evidence_uncalibrated_mass_error_ppm
      )
  ),
  by = .(
    PSM_uid,
    peptide,
    Dataset,
    AS_group,
    Comparison,
    Acquisition_type,
    `Raw file`,
    `Scan number`
  )
]

AS_peptide_evidence_summary <-
  strict_PSM_evidence[
    ,
    {

      n_leading_available <-
        sum(
          leading_razor_available
        )

      n_leading_JCAST <-
        sum(
          leading_razor_is_JCAST_AS
        )

      n_leading_exact_match <-
        sum(
          leading_razor_matches_any_strict_JCAST_source
        )

      n_LFQ_PSM <-
        sum(
          Acquisition_type == "LabelFree"
        )

      n_LFQ_mass_available <-
        sum(
          Acquisition_type == "LabelFree" &
            !is.na(
              evidence_mass_error_ppm
            )
        )

      leading_class <-
        if (n_leading_available == 0) {

          "leading_razor_not_available"

        } else if (
          n_leading_JCAST ==
          n_leading_available
        ) {

          "all_available_PSMs_JCAST_leading_razor"

        } else if (
          n_leading_JCAST > 0
        ) {

          "mixed_JCAST_and_non_JCAST_leading_razor"

        } else {

          "all_available_PSMs_non_JCAST_leading_razor"
        }

      mass_status <-
        if (n_LFQ_PSM == 0) {

          "not_reported_TMT_only"

        } else if (
          n_LFQ_mass_available > 0
        ) {

          "LFQ_evidence_mass_error_available"

        } else {

          "LFQ_evidence_mass_error_not_available"
        }

      .(
        n_PSM_evidence =
          uniqueN(PSM_uid),

        n_evidence_records_found =
          sum(
            evidence_record_found
          ),

        n_missing_evidence_records =
          sum(
            !evidence_record_found
          ),

        n_leading_razor_available =
          n_leading_available,

        n_PSM_JCAST_leading_razor =
          n_leading_JCAST,

        n_PSM_leading_razor_matches_strict_JCAST_source =
          n_leading_exact_match,

        any_JCAST_leading_razor =
          n_leading_JCAST > 0,

        any_leading_razor_matches_strict_JCAST_source =
          n_leading_exact_match > 0,

        leading_razor_support_class =
          leading_class,

        n_LabelFree_PSM =
          n_LFQ_PSM,

        n_LFQ_mass_error_available =
          n_LFQ_mass_available,

        median_LFQ_abs_mass_error_ppm =
          if (
            n_LFQ_mass_available > 0
          ) {
            median(
              abs(
                evidence_mass_error_ppm[
                  Acquisition_type ==
                    "LabelFree" &
                    !is.na(
                      evidence_mass_error_ppm
                    )
                ]
              )
            )
          } else {
            NA_real_
          },

        mass_error_reporting_status =
          mass_status
      )
    },
    by = peptide
  ]

AS_peptide_master_with_evidence <- merge(
  AS_peptide_master_base,
  AS_peptide_evidence_summary,
  by = "peptide",
  all.x = TRUE,
  sort = FALSE
)

AS_peptide_master_with_evidence[
  ,
  evidence_PSM_count_matches_master :=
    n_PSM ==
    n_PSM_evidence
]

data.table::fwrite(
  AS_peptide_evidence_summary,
  file.path(OUTDIR, "12b_AS_peptide_evidence_annotation.tsv"),
  sep = "\t", quote = FALSE, na = "NA"
)

data.table::fwrite(
  AS_peptide_master_with_evidence,
  file.path(OUTDIR, "12c_AS_peptide_master_with_evidence.tsv"),
  sep = "\t", quote = FALSE, na = "NA"
)

saveRDS(
  strict_PSM_evidence,
  file.path(OUTDIR, "12d_AS_strict_PSM_evidence_annotation.rds"),
  compress = TRUE
)

saveRDS(
  AS_peptide_master_with_evidence,
  file.path(OUTDIR, "12c_AS_peptide_master_with_evidence.rds"),
  compress = TRUE
)

# ============================================================
# 12. Publication-level classification
# ============================================================

if (!exists("AS_peptide_master_with_evidence")) {
  AS_peptide_master_with_evidence <- readRDS(
    file.path(
      OUTDIR,
      "12c_AS_peptide_master_with_evidence.rds"
    )
  )
}

setDT(AS_peptide_master_with_evidence)

AS_publication_base <- copy(
  AS_peptide_master_with_evidence
)

AS_publication_base[
  ,
  AS_mapping_class :=
    fcase(
      n_AS_events == 1 &
        n_AS_short_IDs == 1,
      "single_AS_source",

      n_AS_events == 1 &
        n_AS_short_IDs > 1,
      "multiple_isoforms_of_one_AS_event",

      n_AS_events > 1,
      "multiple_AS_events",

      default =
        "AS_mapping_unresolved"
    )
]

AS_publication_base[
  ,
  AS_isoform_specificity_class :=
    fcase(
      has_AS_isoform_specific_source_model,
      "has_AS_isoform_specific_source_model",

      n_shared_between_isoform_source_pairs ==
        n_AS_source_peptide_pairs,
      "shared_between_event_isoforms_only",

      n_partner_unavailable_source_pairs ==
        n_AS_source_peptide_pairs,
      "partner_isoform_unavailable_only",

      default =
        "mixed_non_specific_source_context"
    )
]

AS_publication_base[
  ,
  source_uniqueness_class :=
    fcase(
      all_source_pairs_have_non_JCAST_source,
      "shared_with_non_JCAST_search_database_source",

      any_JCAST_only_in_search_database,
      "JCAST_only_in_MaxQuant_Proteins_field",

      default =
        "source_uniqueness_unresolved"
    )
]

AS_publication_base[
  ,
  source_conclusion :=
    fcase(
      all_source_pairs_have_non_JCAST_source &
        has_AS_isoform_specific_source_model,

      paste0(
        "peptide is sequence-specific to at least one modeled ",
        "JCAST AS isoform but is also matched to non-JCAST ",
        "protein sources in MaxQuant"
      ),

      all_source_pairs_have_non_JCAST_source,

      paste0(
        "peptide is mapped to a JCAST AS-derived protein sequence ",
        "but is also matched to non-JCAST protein sources in MaxQuant"
      ),

      any_JCAST_only_in_search_database,

      paste0(
        "at least one JCAST source mapping had no non-JCAST ",
        "co-match in the MaxQuant Proteins field; independent ",
        "proteome-wide uniqueness verification is still required"
      ),

      default =
        "source assignment remains unresolved"
    )
]

AS_publication_base[
  ,
  strict_MS_MHC_support := TRUE
]

AS_publication_base[
  ,
  supports_natural_HLA_presentation := FALSE
]

AS_publication_base[
  ,
  recommended_claim :=
    fifelse(
      has_AS_isoform_specific_source_model,

      paste0(
        "MS-supported and NetMHCpan-predicted HLA-I-binding ",
        "peptide with sequence specificity to at least one ",
        "modeled AS isoform"
      ),

      paste0(
        "MS-supported and NetMHCpan-predicted HLA-I-binding ",
        "peptide mapped to a JCAST AS-derived protein sequence"
      )
    )
]

AS_publication_base[
  ,
  manual_spectrum_review_priority :=
    fcase(
      has_AS_isoform_specific_source_model &
        best_EL_binder_class == "SB",
      "P1_review_AS_isoform_specific_SB",

      has_AS_isoform_specific_source_model &
        best_EL_binder_class == "WB",
      "P2_review_AS_isoform_specific_WB",

      !has_AS_isoform_specific_source_model &
        best_EL_binder_class == "SB",
      "P3_review_AS_associated_SB",

      !has_AS_isoform_specific_source_model &
        best_EL_binder_class == "WB",
      "P4_review_AS_associated_WB",

      default =
        "review_priority_unresolved"
    )
]

# ============================================================
# 13. Publication candidate tables
# ============================================================

if (!exists("AS_publication_base")) {

  if (!exists("AS_peptide_master_with_evidence")) {
    AS_peptide_master_with_evidence <- readRDS(
      file.path(
        OUTDIR,
        "12c_AS_peptide_master_with_evidence.rds"
      )
    )
  }

  stop(
    "AS_publication_base is not available in the current R environment.",
    "Run the publication-classification step before continuing."
  )
}

setDT(AS_publication_base)

required_AS_publication_columns <- c(
  "peptide",
  "AS_short_IDs",
  "gene_ensg",
  "AS_groups",
  "AS_event_keys",
  "event_types",
  "type_orders",
  "tiers",
  "translation_tier_support",
  "AS_mapping_class",
  "AS_isoform_specificity_class",

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

  "HLA_I_alleles",
  "n_HLA_I_alleles",
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

missing_AS_publication_columns <- setdiff(
  required_AS_publication_columns,
  names(AS_publication_base)
)

if (length(missing_AS_publication_columns) > 0) {

  stop(
    "Required columns are missing from the publication table."
  )
}

AS_publication_table <- AS_publication_base[
  ,
  .(
    peptide,

    AS_short_IDs,

    gene_ensg,

    AS_groups,

    AS_event_keys,

    event_types,

    type_orders,

    tiers,

    translation_tier_support,

    AS_mapping_class,

    AS_isoform_specificity_class,

    n_PSM,

    n_raw_files,

    datasets,

    comparisons,

    min_PEP,

    max_Score,

    max_Delta_score,

    max_PIF,

    PIF_availability,

    median_LFQ_abs_mass_error_ppm,

    mass_error_reporting_status,

    HLA_alleles =
      HLA_I_alleles,

    n_HLA_alleles =
      n_HLA_I_alleles,

    best_EL_allele,

    best_EL_rank,

    best_EL_score,

    best_EL_binder_class,

    strict_MS_MHC_support,

    source_uniqueness_class,

    source_conclusion,

    supports_natural_HLA_presentation,

    recommended_claim,

    manual_spectrum_review_priority
  )
]

AS_publication_table[
  ,
  review_priority_order :=
    fcase(
      manual_spectrum_review_priority ==
        "P1_review_AS_isoform_specific_SB",
      1L,

      manual_spectrum_review_priority ==
        "P2_review_AS_isoform_specific_WB",
      2L,

      manual_spectrum_review_priority ==
        "P3_review_AS_associated_SB",
      3L,

      manual_spectrum_review_priority ==
        "P4_review_AS_associated_WB",
      4L,

      default = 99L
    )
]

setorder(
  AS_publication_table,
  review_priority_order,
  best_EL_rank,
  min_PEP,
  -max_Score,
  peptide,
  na.last = TRUE
)

AS_publication_table[
  ,
  review_priority_order := NULL
]

AS_publication_output_file <- file.path(
  OUTDIR,
  "19a_AS_publication_candidate_summary_HLA_fixed.tsv"
)

fwrite(
  AS_publication_table,
  AS_publication_output_file,
  sep = "\t",
  quote = FALSE,
  na = "NA"
)

saveRDS(
  AS_publication_table,
  file.path(
    OUTDIR,
    "19a_AS_publication_candidate_summary_HLA_fixed.rds"
  ),
  compress = TRUE
)

if (!exists("AS_publication_table")) {

  AS_publication_table <- fread(
    file.path(
      OUTDIR,
      "19a_AS_publication_candidate_summary_HLA_fixed.tsv"
    ),
    sep = "\t",
    header = TRUE,
    data.table = TRUE,
    check.names = FALSE,
    na.strings = c(
      "",
      "NA"
    )
  )
}

setDT(AS_publication_table)

AS_isoform_specific_publication <- copy(
  AS_publication_table[
    AS_isoform_specificity_class ==
      "has_AS_isoform_specific_source_model"
  ]
)

AS_isoform_specific_publication[
  ,
  isoform_review_order :=
    fcase(
      manual_spectrum_review_priority ==
        "P1_review_AS_isoform_specific_SB",
      1L,

      manual_spectrum_review_priority ==
        "P2_review_AS_isoform_specific_WB",
      2L,

      default = 99L
    )
]

setorder(
  AS_isoform_specific_publication,
  isoform_review_order,
  best_EL_rank,
  min_PEP,
  -max_Score,
  -n_PSM,
  peptide,
  na.last = TRUE
)

AS_isoform_specific_publication[
  ,
  isoform_review_order := NULL
]

AS_isoform_specific_output <- file.path(
  OUTDIR,
  paste0(
    "19b_AS_isoform_model_specific_",
    "publication_candidate_summary_HLA_fixed.tsv"
  )
)

fwrite(
  AS_isoform_specific_publication,
  AS_isoform_specific_output,
  sep = "\t",
  quote = FALSE,
  na = "NA"
)

saveRDS(
  AS_isoform_specific_publication,
  file.path(
    OUTDIR,
    paste0(
      "19b_AS_isoform_model_specific_",
      "publication_candidate_summary_HLA_fixed.rds"
    )
  ),
  compress = TRUE
)

# ============================================================
# 14. MaxQuant canonical/isoform source classification for model-specific candidates
# ============================================================

if (!exists("AS_publication_table")) {

  publication_file <- file.path(
    OUTDIR,
    "19a_AS_publication_candidate_summary_HLA_fixed.tsv"
  )

  if (!file.exists(publication_file)) {
    stop("Neither AS_publication_table nor the saved 19a publication table is available.")
  }

  AS_publication_table <- fread(
    publication_file,
    sep = "\t",
    header = TRUE,
    check.names = FALSE,
    na.strings = c("", "NA")
  )
}

setDT(AS_publication_table)

AS_315 <- copy(
  AS_publication_table[
    AS_isoform_specificity_class ==
      "has_AS_isoform_specific_source_model"
  ]
)

if (!exists("AS_PSM_high_conf")) {
  stop(
    "AS_PSM_high_conf is not available in the current R environment.",
    "Load the high-confidence PSM object saved in the PSM filtering step."
  )
}

setDT(AS_PSM_high_conf)

peptide_col <- intersect(
  c(
    "peptide",
    "Sequence",
    "sequence"
  ),
  names(AS_PSM_high_conf)
)

protein_col <- intersect(
  c(
    "Proteins",
    "Protein_matches_MaxQuant",
    "protein_matches_MaxQuant",
    "protein_matches"
  ),
  names(AS_PSM_high_conf)
)

if (length(peptide_col) == 0L) {
  stop("No peptide-sequence column was found in AS_PSM_high_conf.")
}

if (length(protein_col) == 0L) {
  stop("The Proteins column was not found in AS_PSM_high_conf.")
}

peptide_col <- peptide_col[1]
protein_col <- protein_col[1]

AS_PSM_high_conf[
  ,
  peptide_for_isoform_audit :=
    toupper(
      trimws(
        as.character(
          get(peptide_col)
        )
      )
    )
]

AS_PSM_high_conf[
  ,
  proteins_for_isoform_audit :=
    as.character(
      get(protein_col)
    )
]

AS_315_PSM <- copy(
  AS_PSM_high_conf[
    peptide_for_isoform_audit %chin%
      AS_315$peptide
  ]
)

AS_315_protein_long <- AS_315_PSM[
  ,
  {

    protein_string <- proteins_for_isoform_audit[1]

    if (
      is.na(protein_string) ||
      trimws(protein_string) == ""
    ) {

      data.table(
        protein_id = NA_character_
      )

    } else {

      protein_ids <- unlist(
        strsplit(
          protein_string,
          split = ";",
          fixed = TRUE
        ),
        use.names = FALSE
      )

      protein_ids <- trimws(protein_ids)

      data.table(
        protein_id = protein_ids
      )
    }
  },
  by = .(
    peptide =
      peptide_for_isoform_audit
  )
]

AS_315_protein_long <- unique(
  AS_315_protein_long[
    !is.na(protein_id) &
      protein_id != ""
  ]
)

classify_uniprot_source <- function(protein_id) {

  if (
    is.na(protein_id) ||
    trimws(protein_id) == ""
  ) {
    return("missing")
  }

  protein_id <- trimws(protein_id)

  fields <- strsplit(
    protein_id,
    split = "|",
    fixed = TRUE
  )[[1]]

  is_jcast <- (
    length(fields) > 3L &&
      any(
        grepl(
          "^ENSG[0-9]+",
          fields
        )
      )
  )

  if (is_jcast) {
    return("JCAST_custom_AS")
  }

  if (
    length(fields) == 3L &&
    fields[1] == "sp"
  ) {

    accession <- fields[2]

    if (
      grepl(
        "-[0-9]+$",
        accession
      )
    ) {
      return("SwissProt_isoform")
    }

    return("SwissProt_canonical")
  }

  if (
    length(fields) == 3L &&
    fields[1] == "tr"
  ) {
    return("TrEMBL")
  }

  return("other_database_source")
}

extract_uniprot_accession <- function(protein_id) {

  if (
    is.na(protein_id) ||
    trimws(protein_id) == ""
  ) {
    return(NA_character_)
  }

  fields <- strsplit(
    trimws(protein_id),
    split = "|",
    fixed = TRUE
  )[[1]]

  if (
    length(fields) >= 2L &&
    fields[1] %in% c("sp", "tr")
  ) {
    return(fields[2])
  }

  NA_character_
}

AS_315_protein_long[
  ,
  protein_source_type :=
    vapply(
      protein_id,
      classify_uniprot_source,
      character(1)
    )
]

AS_315_protein_long[
  ,
  accession :=
    vapply(
      protein_id,
      extract_uniprot_accession,
      character(1)
    )
]

AS_315_protein_long[
  ,
  parent_accession :=
    fifelse(
      protein_source_type ==
        "SwissProt_isoform",

      sub(
        "-[0-9]+$",
        "",
        accession
      ),

      accession
    )
]

collapse_values <- function(x) {

  x <- as.character(x)

  x <- x[
    !is.na(x) &
      x != ""
  ]

  x <- sort(unique(x))

  if (length(x) == 0L) {
    return(NA_character_)
  }

  paste(
    x,
    collapse = ";"
  )
}

AS_315_uniprot_assignment <- AS_315_protein_long[
  ,
  {

    jcast_accessions <- unique(
      accession[
        protein_source_type ==
          "JCAST_custom_AS"
      ]
    )

    sp_canonical_accessions <- unique(
      accession[
        protein_source_type ==
          "SwissProt_canonical"
      ]
    )

    sp_isoform_accessions <- unique(
      accession[
        protein_source_type ==
          "SwissProt_isoform"
      ]
    )

    sp_isoform_parents <- unique(
      parent_accession[
        protein_source_type ==
          "SwissProt_isoform"
      ]
    )

    tr_accessions <- unique(
      accession[
        protein_source_type ==
          "TrEMBL"
      ]
    )

    jcast_accessions <- jcast_accessions[
      !is.na(jcast_accessions)
    ]

    sp_canonical_accessions <-
      sp_canonical_accessions[
        !is.na(sp_canonical_accessions)
      ]

    sp_isoform_accessions <-
      sp_isoform_accessions[
        !is.na(sp_isoform_accessions)
      ]

    sp_isoform_parents <-
      sp_isoform_parents[
        !is.na(sp_isoform_parents)
      ]

    tr_accessions <- tr_accessions[
      !is.na(tr_accessions)
    ]

    same_parent_isoforms <- intersect(
      sp_isoform_parents,
      jcast_accessions
    )

    other_parent_isoforms <- setdiff(
      sp_isoform_parents,
      jcast_accessions
    )

    list(
      n_JCAST_source_hits =
        length(jcast_accessions),

      JCAST_base_accessions =
        collapse_values(jcast_accessions),

      hit_sp_canonical_in_MaxQuant =
        length(sp_canonical_accessions) > 0L,

      n_sp_canonical_hits =
        length(sp_canonical_accessions),

      sp_canonical_accessions =
        collapse_values(
          sp_canonical_accessions
        ),

      hit_sp_isoform_in_MaxQuant =
        length(sp_isoform_accessions) > 0L,

      n_sp_isoform_hits =
        length(sp_isoform_accessions),

      sp_isoform_accessions =
        collapse_values(
          sp_isoform_accessions
        ),

      sp_isoform_parent_accessions =
        collapse_values(
          sp_isoform_parents
        ),

      isoform_same_parent_as_JCAST =
        length(same_parent_isoforms) > 0L,

      same_parent_isoform_accessions =
        collapse_values(
          sp_isoform_accessions[
            sp_isoform_parents %in%
              same_parent_isoforms
          ]
        ),

      isoform_other_parent_than_JCAST =
        length(other_parent_isoforms) > 0L,

      other_parent_isoform_parents =
        collapse_values(
          other_parent_isoforms
        ),

      hit_tr_in_MaxQuant =
        length(tr_accessions) > 0L,

      n_tr_hits =
        length(tr_accessions),

      tr_accessions =
        collapse_values(tr_accessions)
    )
  },
  by = peptide
]
same_parent_isoform_detail <- AS_315_protein_long[
  protein_source_type ==
    "SwissProt_isoform",
  .(
    same_parent_isoform_accessions =
      collapse_values(
        accession[
          parent_accession %chin%
            AS_315_protein_long[
              peptide == .BY$peptide &
                protein_source_type ==
                "JCAST_custom_AS",
              unique(accession)
            ]
        ]
      )
  ),
  by = peptide
]

AS_315_uniprot_assignment[
  same_parent_isoform_detail,
  on = "peptide",
  same_parent_isoform_accessions :=
    i.same_parent_isoform_accessions
]

AS_315_uniprot_assignment[
  ,
  reviewed_protein_assignment_class :=
    fcase(

      isoform_same_parent_as_JCAST &
        !hit_sp_canonical_in_MaxQuant,

      paste0(
        "target_same_parent_isoform_",
        "without_canonical_assignment"
      ),

      isoform_same_parent_as_JCAST &
        hit_sp_canonical_in_MaxQuant,

      paste0(
        "same_parent_isoform_",
        "with_canonical_assignment"
      ),

      hit_sp_isoform_in_MaxQuant &
        !isoform_same_parent_as_JCAST &
        !hit_sp_canonical_in_MaxQuant,

      paste0(
        "other_parent_isoform_",
        "without_canonical_assignment"
      ),

      hit_sp_isoform_in_MaxQuant &
        !isoform_same_parent_as_JCAST &
        hit_sp_canonical_in_MaxQuant,

      paste0(
        "canonical_and_",
        "other_parent_isoform_assignment"
      ),

      !hit_sp_isoform_in_MaxQuant &
        hit_sp_canonical_in_MaxQuant,

      "canonical_only_assignment",

      !hit_sp_isoform_in_MaxQuant &
        !hit_sp_canonical_in_MaxQuant,

      "no_standard_SwissProt_assignment",

      default =
        "unresolved_assignment"
    )
]

AS_315_uniprot_assignment[
  ,
  is_target_isoform_candidate :=
    reviewed_protein_assignment_class ==
    paste0(
      "target_same_parent_isoform_",
      "without_canonical_assignment"
    )
]

AS_315_isoform_classified <- merge(
  AS_315,
  AS_315_uniprot_assignment,
  by = "peptide",
  all.x = TRUE,
  sort = FALSE
)
setorder(
  AS_315_isoform_classified,
  -is_target_isoform_candidate,
  best_EL_rank,
  min_PEP,
  -max_Score,
  peptide,
  na.last = TRUE
)
step15_output <- file.path(
  OUTDIR,
  paste0(
    "20a_315_AS_paired_model_specific_",
    "MaxQuant_canonical_isoform_classification.tsv"
  )
)

fwrite(
  AS_315_isoform_classified,
  step15_output,
  sep = "\t",
  quote = FALSE,
  na = "NA"
)

saveRDS(
  AS_315_isoform_classified,
  file.path(
    OUTDIR,
    paste0(
      "20a_315_AS_paired_model_specific_",
      "MaxQuant_canonical_isoform_classification.rds"
    )
  ),
  compress = TRUE
)

# ============================================================
# 15. Compact model-specific publication table
# ============================================================

publication_compact_cols <- c(
  "peptide",
  "AS_short_IDs",
  "gene_ensg",
  "AS_groups",
  "event_types",
  "type_orders",
  "tiers",
  "AS_mapping_class",
  "AS_isoform_specificity_class",
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
  "reviewed_protein_assignment_class",
  "source_uniqueness_class",
  "source_conclusion",
  "supports_natural_HLA_presentation",
  "recommended_claim",
  "manual_spectrum_review_priority",
  "is_target_isoform_candidate"
)

missing_compact_cols <- setdiff(publication_compact_cols, names(AS_315_isoform_classified))
if (length(missing_compact_cols) > 0) {
  stop("Missing columns for compact publication table: ", paste(missing_compact_cols, collapse = ", "))
}

AS_publication_compact <- AS_315_isoform_classified[, ..publication_compact_cols]

AS_compact_output <- file.path(
  OUTDIR,
  "20b_315_AS_publication_candidate_summary_HLA_fixed.tsv"
)

data.table::fwrite(
  AS_publication_compact,
  AS_compact_output,
  sep = "\t", quote = FALSE, na = "NA"
)

# ============================================================
# 16. Generate UniProt peptide-query batches
# ============================================================
# Submit the generated batch files to the UniProt peptide-search webpage,
# then save the returned Excel files in the same new_315 directory.

query_peptides <- unique(toupper(trimws(AS_publication_compact$peptide)))
query_peptides <- query_peptides[!is.na(query_peptides) & query_peptides != ""]
peptide_batches <- split(query_peptides, ceiling(seq_along(query_peptides) / 100))

uniprot_query_dir <- "./MS_unspecific/AS/result/contain_jcast/IEDB_uniprot_search/new_315"
dir.create(uniprot_query_dir, recursive = TRUE, showWarnings = FALSE)

for (i in seq_along(peptide_batches)) {
  writeLines(
    paste(peptide_batches[[i]], collapse = ","),
    file.path(uniprot_query_dir, sprintf("AS_peptide_batch_%02d.txt", i))
  )
}

# ============================================================
# 17. UniProt peptide-search post-processing
# ============================================================
# The MaxQuant protein assignments above are retained as search-database context.
# Final source attribution is determined from the independent UniProt peptide query.

UNIPROT_RESULT_DIR <- "./MS_unspecific/AS/result/contain_jcast/IEDB_uniprot_search"
SP_TR_MAP_FILE <- "./MS_unspecific/sp_tr_map.txt"

uniprot_batch_names <- c(
  "AS_peptide_batch_01",
  "AS_peptide_batch_02",
  "AS_peptide_batch_03",
  "AS_peptide_batch_04"
)

parse_uniprot_result <- function(input_file, output_dir) {

  x <- readxl::read_excel(input_file) %>%
    dplyr::rename(
      Entry_Name = `Entry Name`,
      Protein_Names = `Protein Names`,
      Gene_Names = `Gene Names`
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

  peptide_class <- x %>%
    dplyr::select(
      Peptide,
      Entry,
      Entry_Name,
      Protein_Names,
      Gene_Names,
      Organism,
      Match_Type
    ) %>%
    dplyr::distinct() %>%
    dplyr::group_by(Peptide) %>%
    dplyr::summarise(
      has_canonical = any(Match_Type == "canonical"),
      has_isoform = any(Match_Type == "isoform"),
      Final_Class = dplyr::case_when(
        has_canonical & !has_isoform ~ "canonical",
        !has_canonical & has_isoform ~ "isoform",
        has_canonical & has_isoform ~ "both",
        TRUE ~ "unknown"
      ),
      Matched_Entries = paste(unique(Entry), collapse = "; "),
      Matched_Genes = paste(unique(Gene_Names), collapse = "; "),
      n_matches = dplyr::n(),
      .groups = "drop"
    )

  sample_name <- tools::file_path_sans_ext(
    basename(input_file)
  )

  readr::write_csv(
    peptide_class,
    file.path(
      output_dir,
      paste0(sample_name, "_output.csv")
    ),
    na = ""
  )

  peptide_class
}

uniprot_result_files <- file.path(
  uniprot_query_dir,
  paste0(uniprot_batch_names, ".xlsx")
)

if (!all(file.exists(uniprot_result_files))) {
  stop(
    "One or more UniProt peptide-search result files are missing: ",
    paste(
      uniprot_result_files[!file.exists(uniprot_result_files)],
      collapse = ", "
    )
  )
}

invisible(
  lapply(
    uniprot_result_files,
    parse_uniprot_result,
    output_dir = UNIPROT_RESULT_DIR
  )
)

AS_uniprot_output <- data.table::rbindlist(
  lapply(
    file.path(
      UNIPROT_RESULT_DIR,
      paste0(uniprot_batch_names, "_output.csv")
    ),
    data.table::fread
  ),
  use.names = TRUE,
  fill = TRUE
)

sp_tr_map <- data.table::fread(
  SP_TR_MAP_FILE
)

AS_uniprot_output[, row_id := .I]

AS_entries_long <- AS_uniprot_output[
  ,
  .(
    accession = trimws(
      unlist(
        strsplit(
          Matched_Entries,
          ";",
          fixed = TRUE
        )
      )
    )
  ),
  by = row_id
]

AS_entries_long[
  ,
  accession := sub(
    "^(sp|tr)\\|",
    "",
    accession
  )
]

AS_entries_long[
  ,
  parent_accession := sub(
    "-[0-9]+$",
    "",
    accession
  )
]

needed_entries <- unique(
  c(
    AS_entries_long$accession,
    AS_entries_long$parent_accession
  )
)

uniprot_map_small <- unique(
  sp_tr_map[
    Entry %chin% needed_entries,
    .(
      Entry,
      database
    )
  ],
  by = "Entry"
)

AS_entries_long[
  uniprot_map_small,
  on = .(accession = Entry),
  database_exact := i.database
]

AS_entries_long[
  uniprot_map_small,
  on = .(parent_accession = Entry),
  database_parent := i.database
]

AS_entries_long[
  ,
  database_final := fifelse(
    !is.na(database_exact),
    database_exact,
    database_parent
  )
]

AS_entries_long[
  ,
  accession_with_db := fifelse(
    !is.na(database_final),
    paste0(database_final, "|", accession),
    accession
  )
]

AS_entries_long[
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

AS_entries_summary <- AS_entries_long[
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
        accession[
          match_type == "not_found"
        ]
      )
      if (length(x) == 0L) {
        NA_character_
      } else {
        paste(x, collapse = "; ")
      }
    }
  ),
  by = row_id
]

all_AS_search_sp_tr <- merge(
  AS_uniprot_output,
  AS_entries_summary,
  by = "row_id",
  all.x = TRUE,
  sort = FALSE
)

setorder(
  all_AS_search_sp_tr,
  row_id
)

all_AS_search_sp_tr[, row_id := NULL]
all_AS_search_sp_tr[, Matched_Entries_original := Matched_Entries]
all_AS_search_sp_tr[, Matched_Entries := Matched_Entries_with_db]
all_AS_search_sp_tr[, Matched_Entries_with_db := NULL]

all_AS_search_sp_tr <- all_AS_search_sp_tr %>%
  dplyr::mutate(
    has_sp_any = stringr::str_detect(
      Matched_Entries,
      "(^|;\\s*)sp\\|"
    ),
    has_tr_any = stringr::str_detect(
      Matched_Entries,
      "(^|;\\s*)tr\\|"
    ),
    sp_tr_status = dplyr::case_when(
      has_tr_any & !has_sp_any ~ "all_tr",
      has_tr_any & has_sp_any ~ "sp_and_tr",
      has_sp_any & !has_tr_any ~ "sp_only",
      TRUE ~ "unknown"
    )
  ) %>%
  dplyr::select(
    Peptide,
    has_sp_any,
    has_tr_any,
    sp_tr_status,
    has_canonical,
    has_isoform,
    Final_Class,
    has_sp_match,
    has_tr_match,
    n_matches,
    all_entries_annotated,
    unmatched_entries,
    Matched_Entries,
    Matched_Entries_original,
    Matched_Genes
  )

write.table(
  all_AS_search_sp_tr,
  "./MS_unspecific/AS/result/contain_jcast/IEDB_uniprot_search/all_AS_search_sp_tr.txt",
  quote = FALSE,
  row.names = FALSE,
  sep = "\t"
)

# ============================================================
# 18. UniProt-based final source attribution
# ============================================================

tidy_uniprot <- function(df) {

  df %>%
    dplyr::select(
      Peptide,
      Matched_Entries,
      Matched_Genes
    ) %>%
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

AS_uniprot_publication <- tidy_uniprot(
  all_AS_search_sp_tr
)

write.table(
  AS_uniprot_publication,
  "./MS_unspecific/AS/result/contain_jcast/IEDB_uniprot_search/new_all_AS_search_sp_tr.txt",
  quote = FALSE,
  row.names = FALSE,
  sep = "\t"
)

write.csv(
  AS_uniprot_publication,
  "./MS_unspecific/AS/result/contain_jcast/IEDB_uniprot_search/new_all_AS_search_sp_tr.csv",
  quote = FALSE,
  row.names = FALSE
)

# Preserve MaxQuant-derived source fields as search-database context only.
# Final source attribution is assigned from the independent UniProt peptide query.
AS_publication_final <- copy(
  AS_publication_compact
)

setnames(
  AS_publication_final,
  old = c(
    "source_uniqueness_class",
    "source_conclusion"
  ),
  new = c(
    "MaxQuant_source_uniqueness_class",
    "MaxQuant_source_conclusion"
  )
)

AS_publication_final <- merge(
  AS_publication_final,
  as.data.table(AS_uniprot_publication),
  by = "peptide",
  all.x = TRUE,
  sort = FALSE
)

# All publication candidates were included in the UniProt query batches.
# Peptides absent from the returned match table are recorded as having no
# detected UniProt sequence match. This does not establish AS-specific origin.
AS_publication_final[
  is.na(UniProt_match_class),
  `:=`(
    UniProt_match_class = "no_UniProt_match",
    has_sp_canonical = FALSE,
    has_sp_isoform = FALSE,
    has_tr_match = FALSE,
    n_sp_canonical_matches = 0L,
    n_sp_isoform_matches = 0L,
    n_tr_matches = 0L,
    n_total_UniProt_matches = 0L,
    sp_canonical_entries = "",
    sp_isoform_entries = "",
    tr_entries = "",
    matched_genes = ""
  )
]

AS_publication_final[
  ,
  source_uniqueness_class :=
    fcase(
      has_sp_canonical,
      "shared_with_SwissProt_canonical",

      !has_sp_canonical &
        has_sp_isoform,
      "shared_with_SwissProt_isoform",

      !has_sp_canonical &
        !has_sp_isoform &
        has_tr_match,
      "shared_with_TrEMBL",

      UniProt_match_class ==
        "no_UniProt_match",
      "no_UniProt_match_detected",

      default =
        "UniProt_source_unresolved"
    )
]

AS_publication_final[
  ,
  source_conclusion :=
    fcase(
      has_sp_canonical,
      paste0(
        "The peptide sequence matches at least one Swiss-Prot canonical protein; ",
        "the sequence evidence does not support unique assignment to the modeled AS source."
      ),

      !has_sp_canonical &
        has_sp_isoform,
      paste0(
        "The peptide sequence matches at least one Swiss-Prot isoform; ",
        "the sequence evidence does not support unique assignment to the modeled AS source."
      ),

      !has_sp_canonical &
        !has_sp_isoform &
        has_tr_match,
      paste0(
        "The peptide sequence matches at least one TrEMBL protein; ",
        "the sequence evidence does not support unique assignment to the modeled AS source."
      ),

      UniProt_match_class ==
        "no_UniProt_match",
      paste0(
        "No Swiss-Prot canonical, Swiss-Prot isoform, or TrEMBL sequence match ",
        "was detected in the UniProt peptide query; AS-source attribution remains provisional."
      ),

      default =
        "UniProt-based source attribution remains unresolved."
    )
]

# Keep the original conservative recommended_claim unchanged.
AS_publication_final[
  ,
  supports_natural_HLA_presentation := FALSE
]

setorder(
  AS_publication_final,
  best_EL_rank,
  min_PEP,
  -max_Score,
  peptide,
  na.last = TRUE
)

AS_final_output <- file.path(
  OUTDIR,
  "20c_315_AS_publication_candidate_summary_UniProt_source_attribution.tsv"
)

data.table::fwrite(
  AS_publication_final,
  AS_final_output,
  sep = "\t",
  quote = FALSE,
  na = "NA"
)

saveRDS(
  AS_publication_final,
  file.path(
    OUTDIR,
    "20c_315_AS_publication_candidate_summary_UniProt_source_attribution.rds"
  ),
  compress = TRUE
)
