# Step 0: Create an independent lncRNA project directory ---------------------------------------------------

lnc_out_dir <- paste0(
  "./MS_unspecific/lncRNA/result/",
  "lncRNA_MS_MHC_integration"
)

dir.create(
  lnc_out_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

setwd(lnc_out_dir)

# Step 1: Build the lncRNA ORF ID dictionary -------------------------------------------------

ori_map_file <- paste0(
  "./circRNA_peptide4/lncRNA_known/netMHCpan/",
  "ori_shortID_map.txt"
)

short_map_file <- paste0(
  "./circRNA_peptide4/lncRNA_known/netMHCpan/",
  "map_2.txt"
)

#
suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(stringr)
  library(tibble)
})

# ============================================================
# 1. Read the two mapping files
# ============================================================

ori_map <- data.table::fread(
  ori_map_file,
  data.table = FALSE,
  check.names = FALSE
)

short_map <- data.table::fread(
  short_map_file,
  data.table = FALSE,
  check.names = FALSE
)

required_ori_columns <- c(
  "ORF_ID",
  "new_ID"
)

required_short_columns <- c(
  "orf_key",
  "new_ID"
)

missing_ori_columns <- setdiff(
  required_ori_columns,
  names(ori_map)
)

missing_short_columns <- setdiff(
  required_short_columns,
  names(short_map)
)

if (length(missing_ori_columns) > 0) {
  stop(
    "ori_shortID_map.txt is missing required columns: ",
    paste(missing_ori_columns, collapse = ", ")
  )
}

if (length(missing_short_columns) > 0) {
  stop(
    "map_2.txt is missing required columns: ",
    paste(missing_short_columns, collapse = ", ")
  )
}


# ============================================================
# 2. Parse information from full FASTA headers
#
# Example:
# lcl|ORF17_ENST00000768434.1:188:66
# unnamed protein product
# ============================================================

parsed_match <- stringr::str_match(
  ori_map$ORF_ID,
  paste0(
    "^lcl\\|",
    "(ORF[0-9]+_",
    "(ENST[0-9]+(?:\\.[0-9]+)?))",
    ":([0-9]+)",
    ":([0-9]+)",
    "(?:\\s|$)"
  )
)

orf_dictionary_full <- ori_map %>%
  dplyr::transmute(
    source_id_full = as.character(ORF_ID),
    
    source_id = parsed_match[, 2],
    
    transcript_id = parsed_match[, 3],
    
    orf_start = suppressWarnings(
      as.integer(parsed_match[, 4])
    ),
    
    orf_end = suppressWarnings(
      as.integer(parsed_match[, 5])
    ),
    
    netmhcpan_id = as.character(new_ID)
  ) %>%
  dplyr::mutate(
    orf_id = stringr::str_extract(
      source_id,
      "^ORF[0-9]+"
    ),
    
    transcript_id_no_version = sub(
      "\\.[0-9]+$",
      "",
      transcript_id
    ),
    
    orf_strand = dplyr::case_when(
      orf_start < orf_end ~ "plus",
      orf_start > orf_end ~ "minus",
      TRUE ~ "unknown"
    ),
    
    orf_coordinate_length_nt = dplyr::if_else(
      !is.na(orf_start) & !is.na(orf_end),
      abs(orf_end - orf_start) + 1L,
      NA_integer_
    ),
    
    orf_partial = stringr::str_detect(
      source_id_full,
      regex(
        "\\bpartial\\b",
        ignore_case = TRUE
      )
    ),
    
    orf_completion_status = dplyr::if_else(
      orf_partial,
      "partial",
      "complete_or_not_marked_partial"
    )
  )

# ============================================================
# 3. Process map_2.txt
# ============================================================

short_map_clean <- short_map %>%
  dplyr::transmute(
    source_id = as.character(orf_key),
    
    netmhcpan_id_from_short_map =
      as.character(new_ID)
  )

# ============================================================
# 4. Merge the two mapping files
# ============================================================

orf_dictionary <- orf_dictionary_full %>%
  dplyr::left_join(
    short_map_clean,
    by = "source_id"
  ) %>%
  dplyr::mutate(
    mapping_consistent = dplyr::case_when(
      is.na(netmhcpan_id) |
        is.na(netmhcpan_id_from_short_map) ~ FALSE,
      
      netmhcpan_id ==
        netmhcpan_id_from_short_map ~ TRUE,
      
      TRUE ~ FALSE
    )
  )

# ============================================================
# 5. Check for ORFs in map_2 that are absent from the full mapping
# ============================================================

short_not_in_full <- short_map_clean %>%
  dplyr::anti_join(
    orf_dictionary_full %>%
      dplyr::select(source_id),
    by = "source_id"
  )

full_not_in_short <- orf_dictionary_full %>%
  dplyr::select(
    source_id,
    netmhcpan_id
  ) %>%
  dplyr::anti_join(
    short_map_clean,
    by = "source_id"
  )

# ============================================================
# 6. Generate QC statistics
# ============================================================

mapping_qc <- tibble::tibble(
  check_item = c(
    "n_full_mapping_rows",
    "n_short_mapping_rows",
    "n_unparsed_full_headers",
    "n_missing_transcript_id",
    "n_duplicate_source_id_full_map",
    "n_duplicate_netmhcpan_id_full_map",
    "n_duplicate_source_id_short_map",
    "n_duplicate_netmhcpan_id_short_map",
    "n_short_not_in_full",
    "n_full_not_in_short",
    "n_inconsistent_new_ID",
    "n_plus_ORF",
    "n_minus_ORF",
    "n_unknown_ORF",
    "n_partial_ORF"
  ),
  
  value = c(
    nrow(orf_dictionary_full),
    
    nrow(short_map_clean),
    
    sum(is.na(orf_dictionary$source_id)),
    
    sum(is.na(orf_dictionary$transcript_id)),
    
    sum(
      duplicated(
        orf_dictionary_full$source_id
      )
    ),
    
    sum(
      duplicated(
        orf_dictionary_full$netmhcpan_id
      )
    ),
    
    sum(
      duplicated(
        short_map_clean$source_id
      )
    ),
    
    sum(
      duplicated(
        short_map_clean$
          netmhcpan_id_from_short_map
      )
    ),
    
    nrow(short_not_in_full),
    
    nrow(full_not_in_short),
    
    sum(
      !orf_dictionary$mapping_consistent,
      na.rm = TRUE
    ),
    
    sum(
      orf_dictionary$orf_strand == "plus",
      na.rm = TRUE
    ),
    
    sum(
      orf_dictionary$orf_strand == "minus",
      na.rm = TRUE
    ),
    
    sum(
      orf_dictionary$orf_strand == "unknown",
      na.rm = TRUE
    ),
    
    sum(
      orf_dictionary$orf_partial,
      na.rm = TRUE
    )
  )
)


# ============================================================
# 7. Perform strict checks for critical mapping errors
# ============================================================

critical_qc_names <- c(
  "n_unparsed_full_headers",
  "n_missing_transcript_id",
  "n_duplicate_source_id_full_map",
  "n_duplicate_netmhcpan_id_full_map",
  "n_duplicate_source_id_short_map",
  "n_duplicate_netmhcpan_id_short_map",
  "n_short_not_in_full",
  "n_full_not_in_short",
  "n_inconsistent_new_ID"
)

critical_failures <- mapping_qc %>%
  dplyr::filter(
    check_item %in% critical_qc_names,
    value > 0
  )

if (nrow(critical_failures) > 0) {
  
  print(
    critical_failures,
    n = Inf
  )
  
  stop(
    "Critical problems were detected in the ORF ID mapping; ",
    "please resolve them before proceeding to the MaxQuant analysis."
  )
}

# ============================================================
# 8. Define the main-analysis and reverse-ORF audit tables
# ============================================================

orf_dictionary_main <- orf_dictionary %>%
  dplyr::filter(
    orf_strand == "plus"
  )

orf_dictionary_reverse <- orf_dictionary %>%
  dplyr::filter(
    orf_strand == "minus"
  )

# ============================================================
# 9. Reorder columns
# ============================================================

dictionary_columns <- c(
  "source_id",
  "source_id_full",
  "orf_id",
  "transcript_id",
  "transcript_id_no_version",
  "orf_start",
  "orf_end",
  "orf_strand",
  "orf_coordinate_length_nt",
  "orf_partial",
  "orf_completion_status",
  "netmhcpan_id",
  "netmhcpan_id_from_short_map",
  "mapping_consistent"
)

orf_dictionary <- orf_dictionary %>%
  dplyr::select(
    dplyr::all_of(
      dictionary_columns
    )
  )

orf_dictionary_main <- orf_dictionary_main %>%
  dplyr::select(
    dplyr::all_of(
      dictionary_columns
    )
  )

orf_dictionary_reverse <- orf_dictionary_reverse %>%
  dplyr::select(
    dplyr::all_of(
      dictionary_columns
    )
  )

# ============================================================
# 10. Save results
# ============================================================

data.table::fwrite(
  orf_dictionary,
  file.path(
    lnc_out_dir,
    "00_lncRNA_ORF_ID_dictionary_all.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  na = "NA"
)

data.table::fwrite(
  orf_dictionary_main,
  file.path(
    lnc_out_dir,
    "00a_lncRNA_ORF_ID_dictionary_main_plus.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  na = "NA"
)

data.table::fwrite(
  orf_dictionary_reverse,
  file.path(
    lnc_out_dir,
    "00b_lncRNA_ORF_ID_dictionary_reverse_audit.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  na = "NA"
)

data.table::fwrite(
  mapping_qc,
  file.path(
    lnc_out_dir,
    "00c_lncRNA_ORF_ID_mapping_QC.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  na = "NA"
)


# Step 2: Inspect lncRNA MaxQuant result directories -----------------------------------------------


suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(stringr)
  library(tibble)
})

mq_result_root <- paste0(
  "./MS_unspecific/lncRNA/result/",
  "contain_lncRNA"
)

# ============================================================
# 1. Locate the three core MaxQuant files
# ============================================================

mq_files <- list.files(
  path = mq_result_root,
  pattern = "^(msms|evidence|peptides)\\.txt$",
  recursive = TRUE,
  full.names = TRUE,
  ignore.case = TRUE
)


if (length(mq_files) == 0) {
  stop("No core MaxQuant files were found. Check mq_result_root.")
}

# ============================================================
# 2. Build the file manifest
# ============================================================

mq_manifest <- tibble::tibble(
  file_path = normalizePath(
    mq_files,
    winslash = "/",
    mustWork = FALSE
  )
) %>%
  dplyr::mutate(
    file_name = basename(file_path),
    
    file_type = stringr::str_remove(
      file_name,
      regex("\\.txt$", ignore_case = TRUE)
    ),
    
    txt_dir = dirname(file_path),
    
    relative_path = stringr::str_remove(
      file_path,
      paste0(
        "^",
        stringr::str_replace_all(
          normalizePath(
            mq_result_root,
            winslash = "/",
            mustWork = FALSE
          ),
          "([.\\^$|()\\[\\]{}*+?\\\\])",
          "\\\\\\1"
        ),
        "/?"
      )
    ),
    
    dataset = stringr::str_match(
      relative_path,
      "^([^/]+)/"
    )[, 2],
    
    comparison_raw = stringr::str_match(
      relative_path,
      "^[^/]+/([^/]+)/combined/txt/"
    )[, 2],
    
    comparison = stringr::str_remove(
      comparison_raw,
      "_result$"
    ),
    
    dataset_type = dplyr::case_when(
      dataset == "PXD037581" ~ "TMT10",
      dataset == "PXD044963" ~ "LabelFree",
      TRUE ~ "unknown"
    ),
    
    file_size_MB = round(
      file.info(file_path)$size / 1024^2,
      3
    )
  ) %>%
  dplyr::arrange(
    dataset,
    comparison,
    file_type
  )


mq_set_summary <- mq_manifest %>%
  dplyr::count(
    dataset,
    dataset_type,
    comparison,
    txt_dir,
    file_type,
    name = "n_files"
  ) %>%
  tidyr::pivot_wider(
    names_from = file_type,
    values_from = n_files,
    values_fill = 0
  )

for (nm in c("msms", "evidence", "peptides")) {
  if (!nm %in% names(mq_set_summary)) {
    mq_set_summary[[nm]] <- 0L
  }
}

mq_set_summary <- mq_set_summary %>%
  dplyr::mutate(
    complete_MaxQuant_set =
      msms == 1L &
      evidence == 1L &
      peptides == 1L
  ) %>%
  dplyr::arrange(
    dataset,
    comparison
  )


# Step 3: Check MaxQuant column names and consistency across result directories -------------------------------------------------

read_header_safely <- function(path) {
  
  tryCatch(
    {
      cols <- names(
        data.table::fread(
          path,
          nrows = 0,
          check.names = FALSE,
          showProgress = FALSE
        )
      )
      
      tibble::tibble(
        file_path = path,
        column_index = seq_along(cols),
        column_name = cols,
        header_read_success = TRUE,
        header_error = NA_character_
      )
    },
    
    error = function(e) {
      tibble::tibble(
        file_path = path,
        column_index = NA_integer_,
        column_name = NA_character_,
        header_read_success = FALSE,
        header_error = conditionMessage(e)
      )
    }
  )
}

mq_column_inventory <- purrr::map_dfr(
  mq_manifest$file_path,
  read_header_safely
) %>%
  dplyr::left_join(
    mq_manifest %>%
      dplyr::select(
        file_path,
        dataset,
        dataset_type,
        comparison,
        txt_dir,
        file_type
      ),
    by = "file_path"
  ) %>%
  dplyr::select(
    dataset,
    dataset_type,
    comparison,
    txt_dir,
    file_type,
    file_path,
    column_index,
    column_name,
    header_read_success,
    header_error
  ) %>%
  dplyr::arrange(
    dataset,
    comparison,
    file_type,
    column_index
  )

# Check for header-reading failures
header_failures <- mq_column_inventory %>%
  dplyr::filter(
    !header_read_success
  )


# 3.4 Check whether the six result sets have consistent column structures

column_signatures <- mq_column_inventory %>%
  dplyr::filter(
    header_read_success
  ) %>%
  dplyr::group_by(
    dataset,
    dataset_type,
    comparison,
    txt_dir,
    file_type,
    file_path
  ) %>%
  dplyr::summarise(
    n_columns = dplyr::n(),
    column_signature = paste(
      column_name,
      collapse = " || "
    ),
    .groups = "drop"
  )

# Count distinct column structures for each file type
column_structure_qc <- column_signatures %>%
  dplyr::group_by(
    file_type
  ) %>%
  dplyr::summarise(
    n_files = dplyr::n(),
    n_distinct_column_structures =
      dplyr::n_distinct(column_signature),
    min_n_columns = min(n_columns),
    max_n_columns = max(n_columns),
    .groups = "drop"
  )


#
column_structure_by_dataset_type <- column_signatures %>%
  dplyr::group_by(
    dataset_type,
    file_type
  ) %>%
  dplyr::summarise(
    n_files = dplyr::n(),
    n_distinct_column_structures =
      dplyr::n_distinct(column_signature),
    min_n_columns = min(n_columns),
    max_n_columns = max(n_columns),
    .groups = "drop"
  )


# Check required columns
required_columns <- tibble::tribble(
  ~file_type,  ~required_column,
  "msms",      "Sequence",
  "msms",      "Proteins",
  "msms",      "PEP",
  "msms",      "Score",
  "msms",      "Delta score",
  "msms",      "Evidence ID",
  "msms",      "Peptide ID",
  "msms",      "Raw file",
  "msms",      "Charge",
  "evidence",  "id",
  "evidence",  "Leading razor protein",
  "peptides",  "id",
  "peptides",  "Leading razor protein"
)

available_columns <- mq_column_inventory %>%
  dplyr::filter(
    header_read_success
  ) %>%
  dplyr::distinct(
    dataset,
    comparison,
    file_type,
    column_name
  )

required_column_qc <- mq_set_summary %>%
  dplyr::filter(
    complete_MaxQuant_set
  ) %>%
  dplyr::select(
    dataset,
    dataset_type,
    comparison
  ) %>%
  tidyr::crossing(
    required_columns
  ) %>%
  dplyr::left_join(
    available_columns %>%
      dplyr::mutate(
        column_present = TRUE
      ),
    by = c(
      "dataset",
      "comparison",
      "file_type",
      "required_column" = "column_name"
    )
  ) %>%
  dplyr::mutate(
    column_present = tidyr::replace_na(
      column_present,
      FALSE
    )
  )

missing_required_columns <- required_column_qc %>%
  dplyr::filter(
    !column_present
  )

# Save QC results
data.table::fwrite(
  mq_column_inventory,
  file.path(
    lnc_out_dir,
    "01c_MaxQuant_column_inventory.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  na = "NA"
)

data.table::fwrite(
  column_structure_qc,
  file.path(
    lnc_out_dir,
    "01d_MaxQuant_column_structure_QC.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  na = "NA"
)

data.table::fwrite(
  column_structure_by_dataset_type,
  file.path(
    lnc_out_dir,
    "01e_MaxQuant_column_structure_by_dataset_type.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  na = "NA"
)

data.table::fwrite(
  missing_required_columns,
  file.path(
    lnc_out_dir,
    "01f_MaxQuant_missing_required_columns.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  na = "NA"
)


# Step 4A: Read and filter high-confidence PSMs -------------------------------------------------------
# 4A.1 Define the processing workflow


# ============================================================
# PSM filtering parameters
# ============================================================

PSM_MIN_LENGTH <- 7L
PSM_MAX_PEP <- 0.01
PSM_MIN_SCORE <- 40
PSM_MIN_DELTA_SCORE <- 15

# ============================================================
# Function for processing msms.txt from one MaxQuant result directory
# ============================================================

process_one_msms <- function(
    dataset,
    dataset_type,
    comparison,
    txt_dir
) {
  
  msms_file <- file.path(
    txt_dir,
    "msms.txt"
  )
  
  if (!file.exists(msms_file)) {
    stop(
      "msms.txt was not found: ",
      msms_file
    )
  }
  
  required_columns <- c(
    "id",
    "Raw file",
    "Scan number",
    "Scan index",
    "Sequence",
    "Length",
    "Modified sequence",
    "Proteins",
    "Charge",
    "Mass error [ppm]",
    "PEP",
    "Score",
    "Delta score",
    "PIF",
    "Reverse",
    "Contaminant",
    "Peptide ID",
    "Evidence ID"
  )
  
  # Read the header first to ensure columns are selected by name rather than position
  available_columns <- names(
    data.table::fread(
      msms_file,
      nrows = 0,
      check.names = FALSE,
      showProgress = FALSE
    )
  )
  
  missing_columns <- setdiff(
    required_columns,
    available_columns
  )
  
  if (length(missing_columns) > 0) {
    stop(
      dataset,
      " / ",
      comparison,
      " msms.txt is missing required columns: ",
      paste(
        missing_columns,
        collapse = ", "
      )
    )
  }
  
  message(
    "Reading: ",
    dataset,
    " / ",
    comparison
  )
  
  msms_raw <- data.table::fread(
    msms_file,
    select = required_columns,
    data.table = FALSE,
    check.names = FALSE,
    showProgress = TRUE
  )
  
  # ==========================================================
  # Standardize column names and data types
  # ==========================================================
  
  msms <- msms_raw %>%
    dplyr::transmute(
      dataset = dataset,
      dataset_type = dataset_type,
      comparison = comparison,
      
      msms_id = suppressWarnings(
        as.integer(id)
      ),
      
      raw_file = as.character(
        `Raw file`
      ),
      
      scan_number = suppressWarnings(
        as.integer(`Scan number`)
      ),
      
      scan_index = suppressWarnings(
        as.integer(`Scan index`)
      ),
      
      peptide = toupper(
        as.character(Sequence)
      ),
      
      peptide_length_reported =
        suppressWarnings(
          as.integer(Length)
        ),
      
      modified_sequence =
        as.character(
          `Modified sequence`
        ),
      
      proteins =
        as.character(Proteins),
      
      charge = suppressWarnings(
        as.integer(Charge)
      ),
      
      mass_error_ppm =
        suppressWarnings(
          as.numeric(
            `Mass error [ppm]`
          )
        ),
      
      PEP = suppressWarnings(
        as.numeric(PEP)
      ),
      
      Score = suppressWarnings(
        as.numeric(Score)
      ),
      
      Delta_score =
        suppressWarnings(
          as.numeric(
            `Delta score`
          )
        ),
      
      PIF = suppressWarnings(
        as.numeric(PIF)
      ),
      
      Reverse = as.character(
        Reverse
      ),
      
      Contaminant = as.character(
        Contaminant
      ),
      
      peptide_id = suppressWarnings(
        as.integer(`Peptide ID`)
      ),
      
      evidence_id = suppressWarnings(
        as.integer(`Evidence ID`)
      )
    ) %>%
    dplyr::mutate(
      peptide_length = dplyr::coalesce(
        peptide_length_reported,
        nchar(peptide)
      ),
      
      Reverse = dplyr::coalesce(
        Reverse,
        ""
      ),
      
      Contaminant = dplyr::coalesce(
        Contaminant,
        ""
      ),
      
      proteins = dplyr::coalesce(
        proteins,
        ""
      ),
      
      # The id field is unique within each msms.txt file;
      # adding dataset and comparison creates a project-wide unique PSM key
      psm_key = paste(
        dataset,
        comparison,
        msms_id,
        sep = "|"
      ),
      
      contains_lncRNA_candidate =
        stringr::str_detect(
          proteins,
          stringr::regex(
            paste0(
              "lcl\\|",
              "ORF[0-9]+_",
              "ENST[0-9]+",
              "(?:\\.[0-9]+)?",
              ":[0-9]+:[0-9]+"
            )
          )
        ),
      
      contains_SwissProt =
        stringr::str_detect(
          proteins,
          "(^|;)sp\\|"
        ),
      
      contains_TrEMBL =
        stringr::str_detect(
          proteins,
          "(^|;)tr\\|"
        )
    )
  
  # ==========================================================
  # Apply sequential filters for waterfall QC
  # ==========================================================
  
  pass_0 <- rep(
    TRUE,
    nrow(msms)
  )
  
  pass_1 <- pass_0 &
    !is.na(msms$peptide_length) &
    msms$peptide_length >=
    PSM_MIN_LENGTH
  
  pass_2 <- pass_1 &
    !is.na(msms$PEP) &
    msms$PEP < PSM_MAX_PEP
  
  pass_3 <- pass_2 &
    !is.na(msms$Score) &
    msms$Score > PSM_MIN_SCORE
  
  # Strict policy:
  # PSMs with missing Delta score do not pass
  pass_4 <- pass_3 &
    !is.na(msms$Delta_score) &
    msms$Delta_score >=
    PSM_MIN_DELTA_SCORE
  
  pass_5 <- pass_4 &
    msms$Reverse != "+"
  
  pass_6 <- pass_5 &
    msms$Contaminant != "+"
  
  pass_7 <- pass_6 &
    msms$contains_lncRNA_candidate
  
  qc <- tibble::tibble(
    dataset = dataset,
    dataset_type = dataset_type,
    comparison = comparison,
    
    filter_step = c(
      "00_raw_msms",
      "01_length_ge_7",
      "02_PEP_lt_0.01",
      "03_Score_gt_40",
      "04_Delta_score_ge_15_strict",
      "05_remove_reverse",
      "06_remove_contaminant",
      "07_contains_lncRNA_ORF"
    ),
    
    n_PSM = c(
      sum(pass_0),
      sum(pass_1),
      sum(pass_2),
      sum(pass_3),
      sum(pass_4),
      sum(pass_5),
      sum(pass_6),
      sum(pass_7)
    )
  ) %>%
    dplyr::mutate(
      n_removed_at_step =
        dplyr::lag(n_PSM) - n_PSM,
      
      n_removed_at_step =
        dplyr::if_else(
          is.na(n_removed_at_step),
          0L,
          as.integer(
            n_removed_at_step
          )
        ),
      
      fraction_of_raw =
        n_PSM / dplyr::first(n_PSM)
    )
  
  # All PSMs passing the quality filters
  high_confidence_all <- msms[
    pass_6,
    ,
    drop = FALSE
  ]
  
  # PSMs that also map to at least one lncRNA ORF
  high_confidence_candidate <- msms[
    pass_7,
    ,
    drop = FALSE
  ]
  
  list(
    high_confidence_all =
      high_confidence_all,
    
    high_confidence_candidate =
      high_confidence_candidate,
    
    qc = qc
  )
}


# 4A.2 Process the six result sets -----------------------------------------------------------
# ============================================================
# Build the formal MaxQuant analysis manifest
# ============================================================

if (!exists("mq_set_summary")) {
  stop(
    "mq_set_summary is not available in the current environment; ",
    "rerun the MaxQuant file scan and directory summary code first."
  )
}

mq_analysis_manifest <- mq_set_summary %>%
  dplyr::transmute(
    dataset = as.character(dataset),
    dataset_type = as.character(dataset_type),
    comparison = as.character(comparison),
    txt_dir = as.character(txt_dir),
    include_in_analysis = as.logical(
      complete_MaxQuant_set
    )
  ) %>%
  dplyr::arrange(
    dataset,
    comparison
  )

# Save the analysis manifest
data.table::fwrite(
  mq_analysis_manifest,
  file.path(
    lnc_out_dir,
    "01_MaxQuant_analysis_manifest.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  na = "NA"
)
#
mq_sets_to_process <- mq_analysis_manifest %>%
  dplyr::filter(
    include_in_analysis
  ) %>%
  dplyr::select(
    dataset,
    dataset_type,
    comparison,
    txt_dir
  ) %>%
  dplyr::arrange(
    dataset,
    comparison
  )


if (nrow(mq_sets_to_process) != 6) {
  stop(
    "Expected six MaxQuant result sets, but found: ",
    nrow(mq_sets_to_process)
  )
}

msms_results <- purrr::pmap(
  mq_sets_to_process,
  process_one_msms
)

names(msms_results) <- paste(
  mq_sets_to_process$dataset,
  mq_sets_to_process$comparison,
  sep = "__"
)

#


# Step 4A continued: Combine PSM results from the six result sets --------------------------------------------------------
high_confidence_psm_all <- purrr::map_dfr(
  msms_results,
  "high_confidence_all"
)

candidate_psm_preliminary <- purrr::map_dfr(
  msms_results,
  "high_confidence_candidate"
)


# 1. Check uniqueness of PSM keys

psm_key_qc <- tibble::tibble(
  object = c(
    "high_confidence_psm_all",
    "candidate_psm_preliminary"
  ),
  
  n_rows = c(
    nrow(high_confidence_psm_all),
    nrow(candidate_psm_preliminary)
  ),
  
  n_unique_psm_key = c(
    dplyr::n_distinct(
      high_confidence_psm_all$psm_key
    ),
    
    dplyr::n_distinct(
      candidate_psm_preliminary$psm_key
    )
  )
) %>%
  dplyr::mutate(
    psm_key_unique =
      n_rows == n_unique_psm_key
  )


if (any(!psm_key_qc$psm_key_unique)) {
  stop(
    "Duplicate psm_key values were detected; Proteins cannot be expanded safely."
  )
}


#


# 2. Summarize by dataset and comparison


# Step 4B: Expand all protein mappings for each PSM
# 4B.1 Define the protein-accession parser
parse_mq_protein_accession <- function(x) {
  
  x <- trimws(as.character(x))
  
  is_lncRNA <- stringr::str_detect(
    x,
    paste0(
      "^lcl\\|",
      "ORF[0-9]+_",
      "ENST[0-9]+",
      "(?:\\.[0-9]+)?",
      ":[0-9]+:[0-9]+$"
    )
  )
  
  source_id <- dplyr::if_else(
    is_lncRNA,
    x %>%
      stringr::str_remove("^lcl\\|") %>%
      stringr::str_remove(
        ":[0-9]+:[0-9]+$"
      ),
    NA_character_
  )
  
  protein_class <- dplyr::case_when(
    is_lncRNA ~ "lncRNA_ORF",
    stringr::str_detect(x, "^sp\\|") ~ "SwissProt",
    stringr::str_detect(x, "^tr\\|") ~ "TrEMBL",
    stringr::str_detect(
      x,
      "^REV__|^REV_|^reverse_"
    ) ~ "reverse_database",
    stringr::str_detect(
      x,
      "^CON__|^CON_"
    ) ~ "contaminant_database",
    TRUE ~ "other"
  )
  
  tibble::tibble(
    protein_accession = x,
    protein_class = protein_class,
    source_id = source_id
  )
}

# 4B.2 Expand Proteins by semicolon
candidate_protein_mapping_all <- candidate_psm_preliminary %>%
  dplyr::select(
    dataset,
    dataset_type,
    comparison,
    psm_key,
    msms_id,
    raw_file,
    scan_number,
    scan_index,
    peptide,
    peptide_length,
    modified_sequence,
    proteins,
    charge,
    mass_error_ppm,
    PEP,
    Score,
    Delta_score,
    PIF,
    peptide_id,
    evidence_id
  ) %>%
  tidyr::separate_rows(
    proteins,
    sep = ";"
  ) %>%
  dplyr::mutate(
    proteins = trimws(proteins)
  ) %>%
  dplyr::filter(
    !is.na(proteins),
    proteins != ""
  ) %>%
  dplyr::rename(
    protein_accession_raw = proteins
  )

# Parse protein accessions
parsed_proteins <- purrr::map_dfr(
  candidate_protein_mapping_all$
    protein_accession_raw,
  parse_mq_protein_accession
)

candidate_protein_mapping_all <-
  dplyr::bind_cols(
    candidate_protein_mapping_all,
    parsed_proteins %>%
      dplyr::select(
        protein_class,
        source_id
      )
  )


# 4B.3 Retain lncRNA ORF mappings and join the ORF dictionary ---------------------------------------------
candidate_lncRNA_mapping_all <-
  candidate_protein_mapping_all %>%
  dplyr::filter(
    protein_class == "lncRNA_ORF"
  ) %>%
  dplyr::left_join(
    orf_dictionary,
    by = "source_id"
  )


# 4B.4 Check ORFs that could not be mapped to the dictionary ---------------------------------------------------------
lncRNA_mapping_join_qc <-
  candidate_lncRNA_mapping_all %>%
  dplyr::summarise(
    n_mapping_rows = dplyr::n(),
    
    n_unique_PSM =
      dplyr::n_distinct(psm_key),
    
    n_unique_source_id =
      dplyr::n_distinct(source_id),
    
    n_missing_dictionary_match =
      sum(
        is.na(transcript_id)
      ),
    
    n_inconsistent_mapping =
      sum(
        !mapping_consistent,
        na.rm = TRUE
      ),
    
    n_plus_mapping =
      sum(
        orf_strand == "plus",
        na.rm = TRUE
      ),
    
    n_minus_mapping =
      sum(
        orf_strand == "minus",
        na.rm = TRUE
      )
  )


# 4B.6 Build the formal plus-strand mapping table for the main analysis -----------------------------------------------
candidate_lncRNA_mapping_main <-
  candidate_lncRNA_mapping_all %>%
  dplyr::filter(
    orf_strand == "plus"
  )

# 4B.7 Identify PSMs with minus-strand mappings only ----------------------------------------------------

psm_strand_status <-
  candidate_lncRNA_mapping_all %>%
  dplyr::group_by(
    dataset,
    dataset_type,
    comparison,
    psm_key,
    peptide
  ) %>%
  dplyr::summarise(
    has_plus_ORF =
      any(
        orf_strand == "plus",
        na.rm = TRUE
      ),
    
    has_minus_ORF =
      any(
        orf_strand == "minus",
        na.rm = TRUE
      ),
    
    n_plus_ORF =
      dplyr::n_distinct(
        source_id[
          orf_strand == "plus"
        ]
      ),
    
    n_minus_ORF =
      dplyr::n_distinct(
        source_id[
          orf_strand == "minus"
        ]
      ),
    
    .groups = "drop"
  ) %>%
  dplyr::mutate(
    strand_mapping_class =
      dplyr::case_when(
        has_plus_ORF &
          has_minus_ORF ~
          "plus_and_minus",
        
        has_plus_ORF ~
          "plus_only",
        
        has_minus_ORF ~
          "minus_only",
        
        TRUE ~
          "unknown"
      )
  )

# Number of PSMs retained in the main analysis

n_main_plus_psm <-
  candidate_lncRNA_mapping_main %>%
  dplyr::summarise(
    n = dplyr::n_distinct(psm_key)
  ) %>%
  dplyr::pull(n)


# 4B.8 Save results ---------------------------------------------------------------
data.table::fwrite(
  candidate_protein_mapping_all,
  file.path(
    lnc_out_dir,
    "03_all_protein_mappings_from_candidate_PSM.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  na = "NA"
)

data.table::fwrite(
  candidate_lncRNA_mapping_all,
  file.path(
    lnc_out_dir,
    "03a_lncRNA_ORF_mapping_all_strands.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  na = "NA"
)

data.table::fwrite(
  candidate_lncRNA_mapping_main,
  file.path(
    lnc_out_dir,
    "03b_lncRNA_ORF_mapping_main_plus.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  na = "NA"
)

data.table::fwrite(
  psm_strand_status,
  file.path(
    lnc_out_dir,
    "03c_PSM_ORF_strand_mapping_status.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  na = "NA"
)

data.table::fwrite(
  lncRNA_mapping_join_qc,
  file.path(
    lnc_out_dir,
    "03d_lncRNA_ORF_dictionary_join_QC.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  na = "NA"
)


# Step 5: Add evidence, peptide, and protein-mapping context ---------------------------------------
# 5.1 Check uniqueness of PSM + source_id
mapping_key_duplicates <- candidate_lncRNA_mapping_main %>%
  dplyr::count(
    dataset,
    comparison,
    psm_key,
    source_id,
    name = "n"
  ) %>%
  dplyr::filter(
    n > 1
  )


if (nrow(mapping_key_duplicates) > 0) {
  stop(
    "Duplicate ",
    "PSM + source_id keys were detected in candidate_lncRNA_mapping_main."
  )
}

# 5.2 Summarize the complete protein-mapping context for each PSM
psm_protein_context <- candidate_protein_mapping_all %>%
  dplyr::group_by(
    dataset,
    dataset_type,
    comparison,
    psm_key
  ) %>%
  dplyr::summarise(
    n_total_protein_mappings =
      dplyr::n_distinct(
        protein_accession_raw
      ),
    
    n_lncRNA_ORF_mappings_all =
      dplyr::n_distinct(
        protein_accession_raw[
          protein_class == "lncRNA_ORF"
        ]
      ),
    
    n_SwissProt_mappings =
      dplyr::n_distinct(
        protein_accession_raw[
          protein_class == "SwissProt"
        ]
      ),
    
    n_TrEMBL_mappings =
      dplyr::n_distinct(
        protein_accession_raw[
          protein_class == "TrEMBL"
        ]
      ),
    
    n_other_protein_mappings =
      dplyr::n_distinct(
        protein_accession_raw[
          protein_class == "other"
        ]
      ),
    
    has_SwissProt_mapping =
      any(
        protein_class == "SwissProt"
      ),
    
    has_TrEMBL_mapping =
      any(
        protein_class == "TrEMBL"
      ),
    
    SwissProt_accessions =
      paste(
        sort(
          unique(
            protein_accession_raw[
              protein_class == "SwissProt"
            ]
          )
        ),
        collapse = ";"
      ),
    
    TrEMBL_accessions =
      paste(
        sort(
          unique(
            protein_accession_raw[
              protein_class == "TrEMBL"
            ]
          )
        ),
        collapse = ";"
      ),
    
    all_protein_accessions =
      paste(
        sort(
          unique(protein_accession_raw)
        ),
        collapse = ";"
      ),
    
    .groups = "drop"
  ) %>%
  dplyr::mutate(
    SwissProt_accessions =
      dplyr::na_if(
        SwissProt_accessions,
        ""
      ),
    
    TrEMBL_accessions =
      dplyr::na_if(
        TrEMBL_accessions,
        ""
      )
  )
# 5.3 Summarize lncRNA ORF strand context for each PSM

psm_lncRNA_context <- candidate_lncRNA_mapping_all %>%
  dplyr::group_by(
    dataset,
    dataset_type,
    comparison,
    psm_key
  ) %>%
  dplyr::summarise(
    n_lncRNA_source_ids_all =
      dplyr::n_distinct(source_id),
    
    n_plus_source_ids =
      dplyr::n_distinct(
        source_id[
          orf_strand == "plus"
        ]
      ),
    
    n_minus_source_ids =
      dplyr::n_distinct(
        source_id[
          orf_strand == "minus"
        ]
      ),
    
    n_transcripts_all =
      dplyr::n_distinct(transcript_id),
    
    n_plus_transcripts =
      dplyr::n_distinct(
        transcript_id[
          orf_strand == "plus"
        ]
      ),
    
    plus_source_ids =
      paste(
        sort(
          unique(
            source_id[
              orf_strand == "plus"
            ]
          )
        ),
        collapse = ";"
      ),
    
    minus_source_ids =
      paste(
        sort(
          unique(
            source_id[
              orf_strand == "minus"
            ]
          )
        ),
        collapse = ";"
      ),
    
    .groups = "drop"
  ) %>%
  dplyr::mutate(
    plus_source_ids =
      dplyr::na_if(
        plus_source_ids,
        ""
      ),
    
    minus_source_ids =
      dplyr::na_if(
        minus_source_ids,
        ""
      ),
    
    strand_mapping_class =
      dplyr::case_when(
        n_plus_source_ids > 0 &
          n_minus_source_ids > 0 ~
          "plus_and_minus",
        
        n_plus_source_ids > 0 ~
          "plus_only",
        
        n_minus_source_ids > 0 ~
          "minus_only",
        
        TRUE ~
          "unknown"
      )
  )

# 5.4 Read evidence.txt and peptides.txt from the six result sets
read_one_mq_context <- function(
    dataset,
    dataset_type,
    comparison,
    txt_dir
) {
  
  evidence_file <- file.path(
    txt_dir,
    "evidence.txt"
  )
  
  peptides_file <- file.path(
    txt_dir,
    "peptides.txt"
  )
  
  evidence_required <- c(
    "id",
    "Sequence",
    "Proteins",
    "Leading razor protein"
  )
  
  peptides_required <- c(
    "id",
    "Sequence",
    "Proteins",
    "Leading razor protein",
    "Unique (Proteins)"
  )
  
  evidence_available <- names(
    data.table::fread(
      evidence_file,
      nrows = 0,
      check.names = FALSE,
      showProgress = FALSE
    )
  )
  
  peptides_available <- names(
    data.table::fread(
      peptides_file,
      nrows = 0,
      check.names = FALSE,
      showProgress = FALSE
    )
  )
  
  evidence_missing <- setdiff(
    evidence_required,
    evidence_available
  )
  
  peptides_missing <- setdiff(
    peptides_required,
    peptides_available
  )
  
  if (length(evidence_missing) > 0) {
    stop(
      dataset,
      " / ",
      comparison,
      " evidence.txt is missing required columns: ",
      paste(
        evidence_missing,
        collapse = ", "
      )
    )
  }
  
  if (length(peptides_missing) > 0) {
    stop(
      dataset,
      " / ",
      comparison,
      " peptides.txt is missing required columns: ",
      paste(
        peptides_missing,
        collapse = ", "
      )
    )
  }
  
  message(
    "Reading evidence/peptides: ",
    dataset,
    " / ",
    comparison
  )
  
  evidence_raw <- data.table::fread(
    evidence_file,
    select = evidence_required,
    data.table = FALSE,
    check.names = FALSE,
    showProgress = TRUE
  )
  
  peptides_raw <- data.table::fread(
    peptides_file,
    select = peptides_required,
    data.table = FALSE,
    check.names = FALSE,
    showProgress = TRUE
  )
  
  evidence_context <- evidence_raw %>%
    dplyr::transmute(
      dataset = dataset,
      dataset_type = dataset_type,
      comparison = comparison,
      
      evidence_id =
        suppressWarnings(
          as.integer(id)
        ),
      
      evidence_sequence =
        toupper(
          as.character(Sequence)
        ),
      
      evidence_proteins =
        as.character(Proteins),
      
      evidence_leading_razor_protein =
        as.character(
          `Leading razor protein`
        )
    )
  
  peptide_context <- peptides_raw %>%
    dplyr::transmute(
      dataset = dataset,
      dataset_type = dataset_type,
      comparison = comparison,
      
      peptide_id =
        suppressWarnings(
          as.integer(id)
        ),
      
      peptide_table_sequence =
        toupper(
          as.character(Sequence)
        ),
      
      peptide_table_proteins =
        as.character(Proteins),
      
      peptide_leading_razor_protein =
        as.character(
          `Leading razor protein`
        ),
      
      unique_proteins_status =
        as.character(
          `Unique (Proteins)`
        )
    )
  
  list(
    evidence = evidence_context,
    peptides = peptide_context
  )
}


# Read all result sets --------------------------------------------------------------------
mq_context_results <- purrr::pmap(
  mq_sets_to_process,
  read_one_mq_context
)

names(mq_context_results) <- paste(
  mq_sets_to_process$dataset,
  mq_sets_to_process$comparison,
  sep = "__"
)

evidence_context_all <- purrr::map_dfr(
  mq_context_results,
  "evidence"
)

peptide_context_all <- purrr::map_dfr(
  mq_context_results,
  "peptides"
)


# 5.5 Check uniqueness of evidence and peptide IDs ----------------------------------------
evidence_id_qc <- evidence_context_all %>%
  dplyr::count(
    dataset,
    comparison,
    evidence_id,
    name = "n"
  ) %>%
  dplyr::filter(
    n > 1
  )

peptide_id_qc <- peptide_context_all %>%
  dplyr::count(
    dataset,
    comparison,
    peptide_id,
    name = "n"
  ) %>%
  dplyr::filter(
    n > 1
  )


if (
  nrow(evidence_id_qc) > 0 ||
  nrow(peptide_id_qc) > 0
) {
  stop(
    "Duplicate evidence_id or peptide_id values were detected."
  )
}


# 5.6 Build the formal PSM x plus-strand ORF detail table ---------------------------------------------
main_plus_psm_source_detail <-
  candidate_lncRNA_mapping_main %>%
  dplyr::left_join(
    evidence_context_all,
    by = c(
      "dataset",
      "dataset_type",
      "comparison",
      "evidence_id"
    )
  ) %>%
  dplyr::left_join(
    peptide_context_all,
    by = c(
      "dataset",
      "dataset_type",
      "comparison",
      "peptide_id"
    )
  ) %>%
  dplyr::left_join(
    psm_protein_context,
    by = c(
      "dataset",
      "dataset_type",
      "comparison",
      "psm_key"
    )
  ) %>%
  dplyr::left_join(
    psm_lncRNA_context,
    by = c(
      "dataset",
      "dataset_type",
      "comparison",
      "psm_key"
    )
  )


# 5.7 Standardize the leading razor protein --------------------------------------------
main_plus_psm_source_detail <-
  main_plus_psm_source_detail %>%
  dplyr::mutate(
    evidence_leading_razor_protein =
      dplyr::na_if(
        trimws(
          evidence_leading_razor_protein
        ),
        ""
      ),
    
    peptide_leading_razor_protein =
      dplyr::na_if(
        trimws(
          peptide_leading_razor_protein
        ),
        ""
      ),
    
    leading_razor_protein =
      dplyr::coalesce(
        evidence_leading_razor_protein,
        peptide_leading_razor_protein
      ),
    
    leading_razor_class =
      dplyr::case_when(
        stringr::str_detect(
          leading_razor_protein,
          "^lcl\\|ORF[0-9]+_ENST"
        ) ~
          "lncRNA_ORF",
        
        stringr::str_detect(
          leading_razor_protein,
          "^sp\\|"
        ) ~
          "SwissProt",
        
        stringr::str_detect(
          leading_razor_protein,
          "^tr\\|"
        ) ~
          "TrEMBL",
        
        is.na(leading_razor_protein) ~
          "missing",
        
        TRUE ~
          "other"
      ),
    
    leading_razor_source_id =
      dplyr::if_else(
        leading_razor_class ==
          "lncRNA_ORF",
        
        leading_razor_protein %>%
          stringr::str_remove(
            "^lcl\\|"
          ) %>%
          stringr::str_remove(
            ":[0-9]+:[0-9]+$"
          ),
        
        NA_character_
      ),
    
    current_source_is_leading_razor =
      !is.na(
        leading_razor_source_id
      ) &
      source_id ==
      leading_razor_source_id,
    
    evidence_sequence_match =
      !is.na(evidence_sequence) &
      peptide ==
      evidence_sequence,
    
    peptide_table_sequence_match =
      !is.na(
        peptide_table_sequence
      ) &
      peptide ==
      peptide_table_sequence
  )


# 5.8 Generate join QC --------------------------------------------------------------

main_detail_qc <-
  main_plus_psm_source_detail %>%
  dplyr::summarise(
    n_detail_rows =
      dplyr::n(),
    
    n_unique_PSM =
      dplyr::n_distinct(psm_key),
    
    n_unique_peptides =
      dplyr::n_distinct(peptide),
    
    n_unique_source_id =
      dplyr::n_distinct(source_id),
    
    n_unique_transcripts =
      dplyr::n_distinct(transcript_id),
    
    n_missing_evidence_join =
      sum(
        is.na(evidence_sequence)
      ),
    
    n_missing_peptide_join =
      sum(
        is.na(peptide_table_sequence)
      ),
    
    n_evidence_sequence_mismatch =
      sum(
        !evidence_sequence_match,
        na.rm = TRUE
      ),
    
    n_peptide_sequence_mismatch =
      sum(
        !peptide_table_sequence_match,
        na.rm = TRUE
      ),
    
    n_missing_leading_razor =
      sum(
        is.na(leading_razor_protein)
      ),
    
    n_current_source_is_leading =
      sum(
        current_source_is_leading_razor,
        na.rm = TRUE
      ),
    
    n_has_SwissProt_mapping =
      sum(
        has_SwissProt_mapping,
        na.rm = TRUE
      ),
    
    n_has_TrEMBL_mapping =
      sum(
        has_TrEMBL_mapping,
        na.rm = TRUE
      )
  )


# 5.9 Audit PSMs with both plus- and minus-strand mappings ---------------------------------------------
mixed_strand_psm_audit <-
  main_plus_psm_source_detail %>%
  dplyr::filter(
    strand_mapping_class ==
      "plus_and_minus"
  ) %>%
  dplyr::select(
    dataset,
    dataset_type,
    comparison,
    psm_key,
    raw_file,
    scan_number,
    peptide,
    source_id,
    transcript_id,
    orf_start,
    orf_end,
    plus_source_ids,
    minus_source_ids,
    PEP,
    Score,
    Delta_score,
    PIF,
    leading_razor_protein,
    leading_razor_class,
    has_SwissProt_mapping,
    has_TrEMBL_mapping
  ) %>%
  dplyr::arrange(
    dataset,
    comparison,
    psm_key,
    source_id
  )


# 5.10 Save results ---------------------------------------------------------------
data.table::fwrite(
  main_plus_psm_source_detail,
  file.path(
    lnc_out_dir,
    "04_main_plus_PSM_source_detail.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  na = "NA"
)

data.table::fwrite(
  psm_protein_context,
  file.path(
    lnc_out_dir,
    "04a_PSM_all_protein_mapping_context.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  na = "NA"
)

data.table::fwrite(
  psm_lncRNA_context,
  file.path(
    lnc_out_dir,
    "04b_PSM_lncRNA_ORF_mapping_context.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  na = "NA"
)

data.table::fwrite(
  main_detail_qc,
  file.path(
    lnc_out_dir,
    "04c_main_plus_PSM_source_detail_QC.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  na = "NA"
)

data.table::fwrite(
  mixed_strand_psm_audit,
  file.path(
    lnc_out_dir,
    "04d_plus_and_minus_PSM_audit.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  na = "NA"
)


# Step 6: Summarize by dataset + comparison + source_id + peptide ----------------------

collapse_unique <- function(x, sep = ";") {
  
  x <- as.character(x)
  
  x <- x[
    !is.na(x) &
      trimws(x) != ""
  ]
  
  x <- sort(unique(x))
  
  if (length(x) == 0) {
    return(NA_character_)
  }
  
  paste(x, collapse = sep)
}


# 6.2 Summarize PSM-level mapping context ------------------------------------------------------

main_plus_psm_level <- main_plus_psm_source_detail %>%
  dplyr::select(
    dataset,
    dataset_type,
    comparison,
    psm_key,
    peptide,
    n_plus_source_ids,
    n_minus_source_ids,
    strand_mapping_class,
    has_SwissProt_mapping,
    has_TrEMBL_mapping,
    leading_razor_protein,
    leading_razor_class
  ) %>%
  dplyr::distinct(
    dataset,
    comparison,
    psm_key,
    .keep_all = TRUE
  )
# QC check
psm_level_context_qc <- main_plus_psm_level %>%
  dplyr::summarise(
    n_PSM = dplyr::n(),
    
    n_unique_peptides =
      dplyr::n_distinct(peptide),
    
    n_PSM_one_plus_ORF =
      sum(n_plus_source_ids == 1),
    
    n_PSM_two_plus_ORF =
      sum(n_plus_source_ids == 2),
    
    n_PSM_more_than_two_plus_ORF =
      sum(n_plus_source_ids > 2),
    
    n_PSM_with_SwissProt =
      sum(
        has_SwissProt_mapping,
        na.rm = TRUE
      ),
    
    n_PSM_without_SwissProt =
      sum(
        !has_SwissProt_mapping,
        na.rm = TRUE
      ),
    
    n_PSM_with_TrEMBL =
      sum(
        has_TrEMBL_mapping,
        na.rm = TRUE
      ),
    
    n_leading_SwissProt =
      sum(
        leading_razor_class == "SwissProt",
        na.rm = TRUE
      ),
    
    n_leading_lncRNA_ORF =
      sum(
        leading_razor_class == "lncRNA_ORF",
        na.rm = TRUE
      )
  )


# 6.3 Select a representative PSM for each group ---------------------------------------------------------

representative_psm_by_source_peptide <-
  main_plus_psm_source_detail %>%
  dplyr::mutate(
    abs_mass_error_ppm =
      abs(mass_error_ppm),
    
    PIF_for_sort =
      dplyr::coalesce(
        PIF,
        -Inf
      ),
    
    abs_mass_error_for_sort =
      dplyr::coalesce(
        abs_mass_error_ppm,
        Inf
      )
  ) %>%
  dplyr::arrange(
    dataset,
    comparison,
    source_id,
    peptide,
    
    PEP,
    dplyr::desc(Score),
    dplyr::desc(Delta_score),
    dplyr::desc(PIF_for_sort),
    abs_mass_error_for_sort,
    psm_key
  ) %>%
  dplyr::group_by(
    dataset,
    dataset_type,
    comparison,
    source_id,
    peptide
  ) %>%
  dplyr::slice_head(
    n = 1
  ) %>%
  dplyr::ungroup() %>%
  dplyr::transmute(
    dataset,
    dataset_type,
    comparison,
    source_id,
    peptide,
    
    representative_psm_key =
      psm_key,
    
    representative_msms_id =
      msms_id,
    
    representative_raw_file =
      raw_file,
    
    representative_scan_number =
      scan_number,
    
    representative_scan_index =
      scan_index,
    
    representative_modified_sequence =
      modified_sequence,
    
    representative_charge =
      charge,
    
    representative_PEP =
      PEP,
    
    representative_Score =
      Score,
    
    representative_Delta_score =
      Delta_score,
    
    representative_PIF =
      PIF,
    
    representative_mass_error_ppm =
      mass_error_ppm,
    
    representative_abs_mass_error_ppm =
      abs_mass_error_ppm,
    
    representative_evidence_id =
      evidence_id,
    
    representative_peptide_id =
      peptide_id,
    
    representative_leading_razor_protein =
      leading_razor_protein,
    
    representative_leading_razor_class =
      leading_razor_class,
    
    representative_source_is_leading_razor =
      current_source_is_leading_razor,
    
    representative_has_SwissProt_mapping =
      has_SwissProt_mapping,
    
    representative_has_TrEMBL_mapping =
      has_TrEMBL_mapping
  )


# 6.4 Generate the formal MS summary table -----------------------------------------------------------
ms_source_peptide_summary <-
  main_plus_psm_source_detail %>%
  dplyr::group_by(
    dataset,
    dataset_type,
    comparison,
    source_id,
    orf_id,
    transcript_id,
    transcript_id_no_version,
    orf_start,
    orf_end,
    orf_strand,
    orf_partial,
    orf_completion_status,
    netmhcpan_id,
    peptide
  ) %>%
  dplyr::summarise(
    peptide_length =
      dplyr::first(peptide_length),
    
    n_PSM =
      dplyr::n_distinct(psm_key),
    
    n_raw_files =
      dplyr::n_distinct(raw_file),
    
    n_charge_states =
      dplyr::n_distinct(
        charge[
          !is.na(charge)
        ]
      ),
    
    n_modified_sequences =
      dplyr::n_distinct(
        modified_sequence[
          !is.na(modified_sequence)
        ]
      ),
    
    n_evidence_ids =
      dplyr::n_distinct(
        evidence_id[
          !is.na(evidence_id)
        ]
      ),
    
    min_PEP =
      min(
        PEP,
        na.rm = TRUE
      ),
    
    max_Score =
      max(
        Score,
        na.rm = TRUE
      ),
    
    max_Delta_score =
      max(
        Delta_score,
        na.rm = TRUE
      ),
    
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
    
    n_PSM_source_is_leading_razor =
      dplyr::n_distinct(
        psm_key[
          current_source_is_leading_razor
        ]
      ),
    
    source_is_leading_razor_in_any_PSM =
      any(
        current_source_is_leading_razor,
        na.rm = TRUE
      ),
    
    source_is_leading_razor_in_all_PSM =
      all(
        current_source_is_leading_razor
      ),
    
    n_PSM_with_SwissProt_mapping =
      dplyr::n_distinct(
        psm_key[
          has_SwissProt_mapping
        ]
      ),
    
    n_PSM_with_TrEMBL_mapping =
      dplyr::n_distinct(
        psm_key[
          has_TrEMBL_mapping
        ]
      ),
    
    any_SwissProt_mapping =
      any(
        has_SwissProt_mapping,
        na.rm = TRUE
      ),
    
    all_PSM_have_SwissProt_mapping =
      all(
        has_SwissProt_mapping
      ),
    
    any_TrEMBL_mapping =
      any(
        has_TrEMBL_mapping,
        na.rm = TRUE
      ),
    
    SwissProt_accessions =
      collapse_unique(
        SwissProt_accessions
      ),
    
    TrEMBL_accessions =
      collapse_unique(
        TrEMBL_accessions
      ),
    
    leading_razor_proteins =
      collapse_unique(
        leading_razor_protein
      ),
    
    leading_razor_classes =
      collapse_unique(
        leading_razor_class
      ),
    
    unique_proteins_statuses =
      collapse_unique(
        unique_proteins_status
      ),
    
    strand_mapping_classes =
      collapse_unique(
        strand_mapping_class
      ),
    
    n_PSM_also_mapped_to_minus_ORF =
      dplyr::n_distinct(
        psm_key[
          n_minus_source_ids > 0
        ]
      ),
    
    .groups = "drop"
  ) %>%
  dplyr::left_join(
    representative_psm_by_source_peptide,
    by = c(
      "dataset",
      "dataset_type",
      "comparison",
      "source_id",
      "peptide"
    )
  ) %>%
  dplyr::arrange(
    dataset,
    comparison,
    source_id,
    peptide
  )


# 6.5 Add preliminary source-context classification ----------------------------------------------------------

ms_source_peptide_summary <-
  ms_source_peptide_summary %>%
  dplyr::mutate(
    MS_database_mapping_class =
      dplyr::case_when(
        any_SwissProt_mapping &
          any_TrEMBL_mapping ~
          "lncRNA_ORF_plus_SwissProt_plus_TrEMBL",
        
        any_SwissProt_mapping ~
          "lncRNA_ORF_plus_SwissProt",
        
        any_TrEMBL_mapping ~
          "lncRNA_ORF_plus_TrEMBL",
        
        TRUE ~
          "lncRNA_ORF_only_in_MaxQuant_Proteins"
      ),
    
    leading_razor_support_class =
      dplyr::case_when(
        source_is_leading_razor_in_all_PSM ~
          "source_is_leading_in_all_PSM",
        
        source_is_leading_razor_in_any_PSM ~
          "source_is_leading_in_some_PSM",
        
        TRUE ~
          "source_is_never_leading"
      )
  )


# 6.6 QC for the MS summary table ---------------------------------------------------------------

ms_summary_duplicate_keys <-
  ms_source_peptide_summary %>%
  dplyr::count(
    dataset,
    comparison,
    source_id,
    peptide,
    name = "n"
  ) %>%
  dplyr::filter(
    n > 1
  )

ms_summary_qc <- tibble::tibble(
  check_item = c(
    "n_summary_rows",
    "n_unique_summary_keys",
    "n_duplicate_summary_keys",
    "n_unique_peptides",
    "n_unique_source_ids",
    "n_unique_transcripts",
    "n_missing_representative_psm",
    "n_invalid_n_PSM",
    "n_minus_ORF_in_main_summary",
    "n_missing_transcript_id",
    "n_missing_netmhcpan_id"
  ),
  
  value = c(
    nrow(ms_source_peptide_summary),
    
    ms_source_peptide_summary %>%
      dplyr::distinct(
        dataset,
        comparison,
        source_id,
        peptide
      ) %>%
      nrow(),
    
    nrow(ms_summary_duplicate_keys),
    
    dplyr::n_distinct(
      ms_source_peptide_summary$peptide
    ),
    
    dplyr::n_distinct(
      ms_source_peptide_summary$source_id
    ),
    
    dplyr::n_distinct(
      ms_source_peptide_summary$transcript_id
    ),
    
    sum(
      is.na(
        ms_source_peptide_summary$
          representative_psm_key
      )
    ),
    
    sum(
      is.na(ms_source_peptide_summary$n_PSM) |
        ms_source_peptide_summary$n_PSM < 1
    ),
    
    sum(
      ms_source_peptide_summary$orf_strand !=
        "plus"
    ),
    
    sum(
      is.na(
        ms_source_peptide_summary$transcript_id
      )
    ),
    
    sum(
      is.na(
        ms_source_peptide_summary$netmhcpan_id
      )
    )
  )
)


# 6.7 Summarize MS results by dataset and comparison ----------------------------------------------------
ms_summary_by_dataset <-
  ms_source_peptide_summary %>%
  dplyr::group_by(
    dataset,
    dataset_type,
    comparison
  ) %>%
  dplyr::summarise(
    n_source_peptide_rows =
      dplyr::n(),
    
    n_unique_PSM =
      dplyr::n_distinct(
        representative_psm_key
      ),
    
    n_unique_peptides =
      dplyr::n_distinct(peptide),
    
    n_unique_source_ids =
      dplyr::n_distinct(source_id),
    
    n_unique_transcripts =
      dplyr::n_distinct(transcript_id),
    
    n_rows_with_SwissProt =
      sum(
        any_SwissProt_mapping
      ),
    
    n_rows_lncRNA_only_in_MQ =
      sum(
        !any_SwissProt_mapping &
          !any_TrEMBL_mapping
      ),
    
    .groups = "drop"
  )


# 6.9 Save results ----------------------------------------------------------------

data.table::fwrite(
  ms_source_peptide_summary,
  file.path(
    lnc_out_dir,
    "05_MS_source_peptide_summary.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  na = "NA"
)

data.table::fwrite(
  representative_psm_by_source_peptide,
  file.path(
    lnc_out_dir,
    "05a_representative_PSM_by_source_peptide.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  na = "NA"
)

data.table::fwrite(
  ms_summary_qc,
  file.path(
    lnc_out_dir,
    "05b_MS_source_peptide_summary_QC.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  na = "NA"
)

data.table::fwrite(
  ms_summary_by_dataset,
  file.path(
    lnc_out_dir,
    "05c_MS_summary_by_dataset_comparison.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  na = "NA"
)

data.table::fwrite(
  psm_level_context_qc,
  file.path(
    lnc_out_dir,
    "05d_PSM_level_mapping_context_QC.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  na = "NA"
)

# 


# Step 7: Process lncRNA NetMHCpan results ------------------------------------------------


# 7.1 Set the NetMHCpan result path -----------------------------------------------------
netmhcpan_file <- paste0(
  "./circRNA_peptide4/",
  "lncRNA_known/netMHCpan/HLA_1/combined_f.txt"
)

if (!file.exists(netmhcpan_file)) {
  stop(
    "NetMHCpan result file was not found: ",
    netmhcpan_file
  )
}


# 7.2 Read the NetMHCpan file and validate columns -----------------------------------------------------------
netmhcpan_raw <- data.table::fread(
  netmhcpan_file,
  data.table = FALSE,
  check.names = FALSE,
  showProgress = TRUE
)


# Check required columns
required_netmhc_columns <- c(
  "Pos",
  "Peptide",
  "ID",
  "type",
  "_core",
  "_icore",
  "_EL-score",
  "_EL_Rank",
  "_BA-score",
  "_BA_Rank"
)

missing_netmhc_columns <- setdiff(
  required_netmhc_columns,
  names(netmhcpan_raw)
)

if (length(missing_netmhc_columns) > 0) {
  stop(
    "NetMHCpan file is missing required columns: ",
    paste(
      missing_netmhc_columns,
      collapse = ", "
    )
  )
}


# 7.3 Standardize NetMHCpan columns ------------------------------------------------------
netmhcpan_all <- netmhcpan_raw %>%
  dplyr::transmute(
    netmhcpan_id =
      as.character(.data[["ID"]]),
    
    peptide =
      toupper(
        as.character(
          .data[["Peptide"]]
        )
      ),
    
    peptide_position =
      suppressWarnings(
        as.integer(
          .data[["Pos"]]
        )
      ),
    
    HLA_allele =
      as.character(
        .data[["type"]]
      ),
    
    binding_core =
      as.character(
        .data[["_core"]]
      ),
    
    binding_icore =
      as.character(
        .data[["_icore"]]
      ),
    
    EL_score =
      suppressWarnings(
        as.numeric(
          .data[["_EL-score"]]
        )
      ),
    
    EL_rank =
      suppressWarnings(
        as.numeric(
          .data[["_EL_Rank"]]
        )
      ),
    
    BA_score =
      suppressWarnings(
        as.numeric(
          .data[["_BA-score"]]
        )
      ),
    
    BA_rank =
      suppressWarnings(
        as.numeric(
          .data[["_BA_Rank"]]
        )
      )
  ) %>%
  dplyr::mutate(
    peptide_length = nchar(peptide),
    
    EL_binder_class =
      dplyr::case_when(
        is.na(EL_rank) ~
          "missing",
        
        EL_rank < 0.5 ~
          "SB",
        
        EL_rank < 2 ~
          "WB",
        
        TRUE ~
          "NB"
      ),
    
    is_EL_binder =
      EL_binder_class %in%
      c("SB", "WB")
  )

# 7.4 Join the ORF dictionary using the NetMHCpan ID ------------------------------------------------------
netmhcpan_all <- netmhcpan_all %>%
  dplyr::left_join(
    orf_dictionary %>%
      dplyr::select(
        source_id,
        orf_id,
        transcript_id,
        transcript_id_no_version,
        orf_start,
        orf_end,
        orf_strand,
        orf_partial,
        orf_completion_status,
        netmhcpan_id,
        mapping_consistent
      ),
    by = "netmhcpan_id"
  )


# 7.5 Check ID mapping and ORF strand -----------------------------------------------------------
netmhcpan_mapping_qc <- netmhcpan_all %>%
  dplyr::summarise(
    n_raw_rows =
      dplyr::n(),
    
    n_unique_netmhcpan_id =
      dplyr::n_distinct(
        netmhcpan_id
      ),
    
    n_unique_peptides =
      dplyr::n_distinct(
        peptide
      ),
    
    n_unique_HLA_alleles =
      dplyr::n_distinct(
        HLA_allele
      ),
    
    n_missing_ORF_mapping =
      sum(
        is.na(source_id)
      ),
    
    n_inconsistent_mapping =
      sum(
        !mapping_consistent,
        na.rm = TRUE
      ),
    
    n_plus_rows =
      sum(
        orf_strand == "plus",
        na.rm = TRUE
      ),
    
    n_minus_rows =
      sum(
        orf_strand == "minus",
        na.rm = TRUE
      ),
    
    n_SB_rows =
      sum(
        EL_binder_class == "SB",
        na.rm = TRUE
      ),
    
    n_WB_rows =
      sum(
        EL_binder_class == "WB",
        na.rm = TRUE
      ),
    
    n_NB_rows =
      sum(
        EL_binder_class == "NB",
        na.rm = TRUE
      ),
    
    n_missing_EL_rank =
      sum(
        is.na(EL_rank)
      )
  )


# 7.6 Generate binder audit tables for all strands ------------------------------------------------
netmhcpan_binders_all_strands <-
  netmhcpan_all %>%
  dplyr::filter(
    is_EL_binder
  )

# Only plus-strand ORFs are retained in the formal main analysis

netmhcpan_binders_plus <-
  netmhcpan_binders_all_strands %>%
  dplyr::filter(
    orf_strand == "plus"
  )

netmhcpan_binders_minus_audit <-
  netmhcpan_binders_all_strands %>%
  dplyr::filter(
    orf_strand == "minus"
  )


# Step 8: Deduplicate by ORF, peptide, and HLA allele --------------------------------------------------


# 8.1 Helper function ----------------------------------------------------------------
collapse_unique_numeric <- function(x) {
  
  x <- suppressWarnings(
    as.numeric(x)
  )
  
  x <- sort(
    unique(
      x[
        !is.na(x)
      ]
    )
  )
  
  if (length(x) == 0) {
    return(NA_character_)
  }
  
  paste(
    x,
    collapse = ";"
  )
}


# 8.2 Count predicted positions for each allele ---------------------------------------------------

netmhcpan_position_summary <-
  netmhcpan_binders_plus %>%
  dplyr::group_by(
    source_id,
    orf_id,
    transcript_id,
    netmhcpan_id,
    peptide,
    HLA_allele
  ) %>%
  dplyr::summarise(
    n_predicted_positions =
      dplyr::n_distinct(
        peptide_position
      ),
    
    predicted_positions =
      collapse_unique_numeric(
        peptide_position
      ),
    
    .groups = "drop"
  )


# 8.3 Select one row for each ORF, peptide, and allele --------------------------------------------
netmhcpan_allele_dedup <-
  netmhcpan_binders_plus %>%
  dplyr::mutate(
    EL_rank_sort =
      dplyr::coalesce(
        EL_rank,
        Inf
      ),
    
    EL_score_sort =
      dplyr::coalesce(
        EL_score,
        -Inf
      ),
    
    BA_rank_sort =
      dplyr::coalesce(
        BA_rank,
        Inf
      ),
    
    BA_score_sort =
      dplyr::coalesce(
        BA_score,
        -Inf
      )
  ) %>%
  dplyr::arrange(
    source_id,
    peptide,
    HLA_allele,
    EL_rank_sort,
    dplyr::desc(
      EL_score_sort
    ),
    BA_rank_sort,
    dplyr::desc(
      BA_score_sort
    ),
    peptide_position
  ) %>%
  dplyr::group_by(
    source_id,
    orf_id,
    transcript_id,
    netmhcpan_id,
    peptide,
    HLA_allele
  ) %>%
  dplyr::slice_head(
    n = 1
  ) %>%
  dplyr::ungroup() %>%
  dplyr::select(
    -EL_rank_sort,
    -EL_score_sort,
    -BA_rank_sort,
    -BA_score_sort
  ) %>%
  dplyr::left_join(
    netmhcpan_position_summary,
    by = c(
      "source_id",
      "orf_id",
      "transcript_id",
      "netmhcpan_id",
      "peptide",
      "HLA_allele"
    )
  )


# 8.4 Validate deduplication keys ---------------------------------------------------------------

netmhcpan_allele_duplicate_qc <-
  netmhcpan_allele_dedup %>%
  dplyr::count(
    source_id,
    peptide,
    HLA_allele,
    name = "n"
  ) %>%
  dplyr::filter(
    n > 1
  )


if (nrow(netmhcpan_allele_duplicate_qc) > 0) {
  stop(
    "Duplicate source_id + peptide + HLA allele keys remain after NetMHCpan deduplication."
  )
}


# Step 9: Summarize by source_id + peptide ---------------------------------------------

# 9.1 Select the best EL allele -------------------------------------------------------

best_EL_allele_table <-
  netmhcpan_allele_dedup %>%
  dplyr::mutate(
    EL_rank_sort =
      dplyr::coalesce(
        EL_rank,
        Inf
      ),
    
    EL_score_sort =
      dplyr::coalesce(
        EL_score,
        -Inf
      )
  ) %>%
  dplyr::arrange(
    source_id,
    peptide,
    EL_rank_sort,
    dplyr::desc(
      EL_score_sort
    ),
    HLA_allele
  ) %>%
  dplyr::group_by(
    source_id,
    peptide
  ) %>%
  dplyr::slice_head(
    n = 1
  ) %>%
  dplyr::ungroup() %>%
  dplyr::transmute(
    source_id,
    peptide,
    
    best_EL_allele =
      HLA_allele,
    
    best_EL_rank =
      EL_rank,
    
    best_EL_score =
      EL_score,
    
    best_EL_binder_class =
      EL_binder_class,
    
    best_EL_position =
      peptide_position,
    
    best_EL_binding_core =
      binding_core,
    
    best_EL_binding_icore =
      binding_icore
  )


# 9.2 Select the best BA allele -------------------------------------------------------
best_BA_allele_table <-
  netmhcpan_allele_dedup %>%
  dplyr::mutate(
    BA_rank_sort =
      dplyr::coalesce(
        BA_rank,
        Inf
      ),
    
    BA_score_sort =
      dplyr::coalesce(
        BA_score,
        -Inf
      )
  ) %>%
  dplyr::arrange(
    source_id,
    peptide,
    BA_rank_sort,
    dplyr::desc(
      BA_score_sort
    ),
    HLA_allele
  ) %>%
  dplyr::group_by(
    source_id,
    peptide
  ) %>%
  dplyr::slice_head(
    n = 1
  ) %>%
  dplyr::ungroup() %>%
  dplyr::transmute(
    source_id,
    peptide,
    
    best_BA_allele =
      HLA_allele,
    
    best_BA_rank =
      BA_rank,
    
    best_BA_score =
      BA_score,
    
    best_BA_position =
      peptide_position
  )


# 9.3 Generate the source-peptide-level MHC summary ----------------------------------------------
netmhcpan_source_peptide_summary <-
  netmhcpan_allele_dedup %>%
  dplyr::group_by(
    source_id,
    orf_id,
    transcript_id,
    transcript_id_no_version,
    orf_start,
    orf_end,
    orf_strand,
    orf_partial,
    orf_completion_status,
    netmhcpan_id,
    peptide
  ) %>%
  dplyr::summarise(
    peptide_length =
      dplyr::first(
        peptide_length
      ),
    
    n_HLA_alleles =
      dplyr::n_distinct(
        HLA_allele
      ),
    
    HLA_alleles =
      collapse_unique(
        HLA_allele
      ),
    
    n_SB_alleles =
      sum(
        EL_binder_class == "SB"
      ),
    
    n_WB_alleles =
      sum(
        EL_binder_class == "WB"
      ),
    
    predicted_EL_binder_classes =
      collapse_unique(
        EL_binder_class
      ),
    
    max_EL_score_across_alleles =
      if (
        all(is.na(EL_score))
      ) {
        NA_real_
      } else {
        max(
          EL_score,
          na.rm = TRUE
        )
      },
    
    min_EL_rank_across_alleles =
      if (
        all(is.na(EL_rank))
      ) {
        NA_real_
      } else {
        min(
          EL_rank,
          na.rm = TRUE
        )
      },
    
    max_BA_score_across_alleles =
      if (
        all(is.na(BA_score))
      ) {
        NA_real_
      } else {
        max(
          BA_score,
          na.rm = TRUE
        )
      },
    
    min_BA_rank_across_alleles =
      if (
        all(is.na(BA_rank))
      ) {
        NA_real_
      } else {
        min(
          BA_rank,
          na.rm = TRUE
        )
      },
    
    predicted_positions =
      collapse_unique(
        predicted_positions
      ),
    
    .groups = "drop"
  ) %>%
  dplyr::left_join(
    best_EL_allele_table,
    by = c(
      "source_id",
      "peptide"
    )
  ) %>%
  dplyr::left_join(
    best_BA_allele_table,
    by = c(
      "source_id",
      "peptide"
    )
  ) %>%
  dplyr::arrange(
    source_id,
    peptide
  )

# Step 10: Strictly integrate MS and MHC results ----------------------------------------------------------

netmhcpan_annotation_for_join <-
  netmhcpan_source_peptide_summary %>%
  dplyr::select(
    source_id,
    peptide,
    n_HLA_alleles,
    HLA_alleles,
    n_SB_alleles,
    n_WB_alleles,
    predicted_EL_binder_classes,
    best_EL_allele,
    best_EL_rank,
    best_EL_score,
    best_EL_binder_class,
    best_EL_position,
    best_EL_binding_core,
    best_EL_binding_icore,
    max_EL_score_across_alleles,
    min_EL_rank_across_alleles,
    best_BA_allele,
    best_BA_rank,
    best_BA_score,
    best_BA_position,
    max_BA_score_across_alleles,
    min_BA_rank_across_alleles,
    predicted_positions
  )

# Strict join
ms_mhc_strict_overlap <-
  ms_source_peptide_summary %>%
  dplyr::inner_join(
    netmhcpan_annotation_for_join,
    by = c(
      "source_id",
      "peptide"
    )
  ) %>%
  dplyr::arrange(
    dataset,
    comparison,
    source_id,
    peptide
  )


# 10.1 Retain MS rows without MHC binder support -----------------------------------------------
ms_without_mhc_binder <-
  ms_source_peptide_summary %>%
  dplyr::anti_join(
    netmhcpan_annotation_for_join,
    by = c(
      "source_id",
      "peptide"
    )
  )


# 10.2 QC for strict MS-MHC integration -------------------------------------------------------------
strict_integration_qc <- tibble::tibble(
  check_item = c(
    "n_MS_source_peptide_rows",
    "n_MS_unique_peptides",
    "n_MHC_plus_binder_source_peptide_rows",
    "n_MHC_plus_binder_unique_peptides",
    "n_strict_overlap_rows",
    "n_strict_overlap_unique_peptides",
    "n_strict_overlap_unique_source_ids",
    "n_MS_rows_without_MHC_binder",
    "n_duplicate_strict_overlap_keys"
  ),
  
  value = c(
    nrow(
      ms_source_peptide_summary
    ),
    
    dplyr::n_distinct(
      ms_source_peptide_summary$
        peptide
    ),
    
    nrow(
      netmhcpan_source_peptide_summary
    ),
    
    dplyr::n_distinct(
      netmhcpan_source_peptide_summary$
        peptide
    ),
    
    nrow(
      ms_mhc_strict_overlap
    ),
    
    dplyr::n_distinct(
      ms_mhc_strict_overlap$
        peptide
    ),
    
    dplyr::n_distinct(
      ms_mhc_strict_overlap$
        source_id
    ),
    
    nrow(
      ms_without_mhc_binder
    ),
    
    ms_mhc_strict_overlap %>%
      dplyr::count(
        dataset,
        comparison,
        source_id,
        peptide
      ) %>%
      dplyr::filter(
        n > 1
      ) %>%
      nrow()
  )
)


# Save NetMHCpan and strict-integration results ------------------------------------------------------

data.table::fwrite(
  netmhcpan_all,
  file.path(
    lnc_out_dir,
    "06_NetMHCpan_all_rows_mapped.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  na = "NA"
)

data.table::fwrite(
  netmhcpan_mapping_qc,
  file.path(
    lnc_out_dir,
    "06a_NetMHCpan_mapping_QC.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  na = "NA"
)

data.table::fwrite(
  netmhcpan_binders_plus,
  file.path(
    lnc_out_dir,
    "06b_NetMHCpan_plus_binder_detail.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  na = "NA"
)

data.table::fwrite(
  netmhcpan_binders_minus_audit,
  file.path(
    lnc_out_dir,
    "06c_NetMHCpan_minus_binder_audit.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  na = "NA"
)

data.table::fwrite(
  netmhcpan_allele_dedup,
  file.path(
    lnc_out_dir,
    "06d_NetMHCpan_plus_binder_by_allele.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  na = "NA"
)

data.table::fwrite(
  netmhcpan_source_peptide_summary,
  file.path(
    lnc_out_dir,
    "06e_NetMHCpan_source_peptide_summary.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  na = "NA"
)

data.table::fwrite(
  ms_mhc_strict_overlap,
  file.path(
    lnc_out_dir,
    "07_MS_NetMHCpan_strict_source_peptide_overlap.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  na = "NA"
)

data.table::fwrite(
  ms_without_mhc_binder,
  file.path(
    lnc_out_dir,
    "07a_MS_source_peptide_without_MHC_binder.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  na = "NA"
)

data.table::fwrite(
  strict_integration_qc,
  file.path(
    lnc_out_dir,
    "07b_MS_NetMHCpan_strict_integration_QC.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  na = "NA"
)


# Step 7 correction: Identify and remove malformed NetMHCpan rows -------------------------------------------------

# 1. Inspect malformed rows

netmhcpan_invalid_rows <- netmhcpan_all %>%
  dplyr::filter(
    is.na(source_id) |
      is.na(EL_rank) |
      !stringr::str_detect(
        netmhcpan_id,
        "^LncID_[0-9]+$"
      )
  )


# Check whether any valid lncRNA IDs truly failed to map

valid_lncid_but_unmapped <- netmhcpan_all %>%
  dplyr::filter(
    stringr::str_detect(
      netmhcpan_id,
      "^LncID_[0-9]+$"
    ),
    is.na(source_id)
  ) %>%
  dplyr::distinct(
    netmhcpan_id
  ) %>%
  dplyr::arrange(
    netmhcpan_id
  )

if (nrow(valid_lncid_but_unmapped) > 0) {
  stop(
    "Validly formatted LncID entries failed to map; ",
    "the NetMHCpan file and ID mapping table may come from different versions."
  )
}

# Remove duplicated header rows
netmhcpan_clean <- netmhcpan_all %>%
  dplyr::filter(
    netmhcpan_id != "ID",
    peptide != "PEPTIDE",
    HLA_allele != "type",
    stringr::str_detect(
      netmhcpan_id,
      "^LncID_[0-9]+$"
    ),
    !is.na(source_id),
    !is.na(EL_rank),
    mapping_consistent
  )

netmhcpan_cleaning_qc <- tibble::tibble(
  check_item = c(
    "n_rows_before_cleaning",
    "n_repeated_header_rows",
    "n_rows_after_cleaning",
    "n_missing_source_id_after_cleaning",
    "n_missing_EL_rank_after_cleaning",
    "n_invalid_netmhcpan_id_after_cleaning"
  ),
  value = c(
    nrow(netmhcpan_all),
    nrow(netmhcpan_invalid_rows),
    nrow(netmhcpan_clean),
    sum(is.na(netmhcpan_clean$source_id)),
    sum(is.na(netmhcpan_clean$EL_rank)),
    sum(
      !stringr::str_detect(
        netmhcpan_clean$netmhcpan_id,
        "^LncID_[0-9]+$"
      )
    )
  )
)


# Save malformed-row records and cleanup QC
data.table::fwrite(
  netmhcpan_invalid_rows,
  file.path(
    lnc_out_dir,
    "06a1_NetMHCpan_repeated_header_rows.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  na = "NA"
)

data.table::fwrite(
  netmhcpan_cleaning_qc,
  file.path(
    lnc_out_dir,
    "06a2_NetMHCpan_cleaning_QC.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  na = "NA"
)

# Rebuild the formal binder objects ----------------------------------------------------------

netmhcpan_binders_all_strands <- netmhcpan_clean %>%
  dplyr::filter(
    is_EL_binder
  )

netmhcpan_binders_plus <- netmhcpan_binders_all_strands %>%
  dplyr::filter(
    orf_strand == "plus"
  )

netmhcpan_binders_minus_audit <- netmhcpan_binders_all_strands %>%
  dplyr::filter(
    orf_strand == "minus"
  )


# Step 11: Peptide-only overlap and transcript-consistency diagnostics -------------------------------------


# 11.1 Ensure helper functions are available -----------------------------------------------------------

if (!exists("collapse_unique")) {
  
  collapse_unique <- function(x, sep = ";") {
    
    x <- as.character(x)
    
    x <- x[
      !is.na(x) &
        trimws(x) != ""
    ]
    
    x <- sort(unique(x))
    
    if (length(x) == 0) {
      return(NA_character_)
    }
    
    paste(
      x,
      collapse = sep
    )
  }
}

# 11.2 Prepare the MS-side data ------------------------------------------------------------

ms_overlap_input <-
  ms_source_peptide_summary %>%
  dplyr::transmute(
    dataset,
    dataset_type,
    comparison,
    peptide,
    
    ms_source_id =
      source_id,
    
    ms_orf_id =
      orf_id,
    
    ms_transcript_id =
      transcript_id,
    
    ms_netmhcpan_id =
      netmhcpan_id,
    
    n_PSM,
    min_PEP,
    max_Score,
    max_Delta_score,
    max_PIF,
    
    MS_database_mapping_class,
    leading_razor_support_class
  ) %>%
  dplyr::distinct()
# 11.3 Prepare the NetMHCpan-side data -----------------------------------------------------

mhc_overlap_input <-
  netmhcpan_source_peptide_summary %>%
  dplyr::transmute(
    peptide,
    
    mhc_source_id =
      source_id,
    
    mhc_orf_id =
      orf_id,
    
    mhc_transcript_id =
      transcript_id,
    
    mhc_netmhcpan_id =
      netmhcpan_id,
    
    n_HLA_alleles,
    HLA_alleles,
    n_SB_alleles,
    n_WB_alleles,
    
    best_EL_allele,
    best_EL_rank,
    best_EL_score,
    best_EL_binder_class,
    
    best_BA_allele,
    best_BA_rank,
    best_BA_score
  ) %>%
  dplyr::distinct()


# 11.4 Perform a many-to-many join by peptide only ---------------------------------------------------
peptide_only_overlap_mapping <-
  ms_overlap_input %>%
  dplyr::inner_join(
    mhc_overlap_input,
    by = "peptide",
    relationship = "many-to-many"
  ) %>%
  dplyr::mutate(
    same_source_id =
      ms_source_id ==
      mhc_source_id,
    
    same_transcript =
      ms_transcript_id ==
      mhc_transcript_id,
    
    same_transcript_different_ORF =
      same_transcript &
      !same_source_id,
    
    different_transcript =
      !same_transcript,
    
    source_relationship =
      dplyr::case_when(
        same_source_id ~
          "same_source_id",
        
        same_transcript_different_ORF ~
          "same_transcript_different_ORF",
        
        TRUE ~
          "different_transcript"
      )
  ) %>%
  dplyr::arrange(
    dataset,
    comparison,
    peptide,
    ms_transcript_id,
    ms_source_id,
    mhc_transcript_id,
    mhc_source_id
  )


# 11.5 Summarize matches for each MS source-peptide pair -----------------------------------------
ms_overlap_annotation <-
  peptide_only_overlap_mapping %>%
  dplyr::group_by(
    dataset,
    dataset_type,
    comparison,
    peptide,
    ms_source_id,
    ms_orf_id,
    ms_transcript_id,
    ms_netmhcpan_id
  ) %>%
  dplyr::summarise(
    n_all_MHC_source_matches =
      dplyr::n_distinct(
        mhc_source_id
      ),
    
    n_all_MHC_transcript_matches =
      dplyr::n_distinct(
        mhc_transcript_id
      ),
    
    n_exact_source_matches =
      dplyr::n_distinct(
        mhc_source_id[
          same_source_id
        ]
      ),
    
    n_same_transcript_different_ORF_matches =
      dplyr::n_distinct(
        mhc_source_id[
          same_transcript_different_ORF
        ]
      ),
    
    n_different_transcript_matches =
      dplyr::n_distinct(
        mhc_source_id[
          different_transcript
        ]
      ),
    
    any_same_source_id =
      any(
        same_source_id,
        na.rm = TRUE
      ),
    
    any_same_transcript_different_ORF =
      any(
        same_transcript_different_ORF,
        na.rm = TRUE
      ),
    
    any_same_transcript =
      any(
        same_transcript,
        na.rm = TRUE
      ),
    
    exact_MHC_source_ids =
      collapse_unique(
        mhc_source_id[
          same_source_id
        ]
      ),
    
    same_transcript_different_ORF_MHC_source_ids =
      collapse_unique(
        mhc_source_id[
          same_transcript_different_ORF
        ]
      ),
    
    different_transcript_MHC_source_ids =
      collapse_unique(
        mhc_source_id[
          different_transcript
        ]
      ),
    
    all_peptide_match_MHC_source_ids =
      collapse_unique(
        mhc_source_id
      ),
    
    all_peptide_match_MHC_transcripts =
      collapse_unique(
        mhc_transcript_id
      ),
    
    .groups = "drop"
  )


# 11.6 Restore MS rows with no peptide-level match ----------------------------------------------

ms_source_peptide_overlap_status <-
  ms_overlap_input %>%
  dplyr::left_join(
    ms_overlap_annotation,
    by = c(
      "dataset",
      "dataset_type",
      "comparison",
      "peptide",
      "ms_source_id",
      "ms_orf_id",
      "ms_transcript_id",
      "ms_netmhcpan_id"
    )
  ) %>%
  dplyr::mutate(
    n_all_MHC_source_matches =
      tidyr::replace_na(
        n_all_MHC_source_matches,
        0L
      ),
    
    n_all_MHC_transcript_matches =
      tidyr::replace_na(
        n_all_MHC_transcript_matches,
        0L
      ),
    
    n_exact_source_matches =
      tidyr::replace_na(
        n_exact_source_matches,
        0L
      ),
    
    n_same_transcript_different_ORF_matches =
      tidyr::replace_na(
        n_same_transcript_different_ORF_matches,
        0L
      ),
    
    n_different_transcript_matches =
      tidyr::replace_na(
        n_different_transcript_matches,
        0L
      ),
    
    any_same_source_id =
      tidyr::replace_na(
        any_same_source_id,
        FALSE
      ),
    
    any_same_transcript_different_ORF =
      tidyr::replace_na(
        any_same_transcript_different_ORF,
        FALSE
      ),
    
    any_same_transcript =
      tidyr::replace_na(
        any_same_transcript,
        FALSE
      ),
    
    MS_MHC_overlap_class =
      dplyr::case_when(
        any_same_source_id ~
          "exact_source_id_and_peptide",
        
        any_same_transcript_different_ORF ~
          "same_transcript_different_ORF",
        
        n_all_MHC_source_matches > 0 ~
          "same_peptide_different_transcript_only",
        
        TRUE ~
          "no_MHC_binder_for_peptide"
      )
  ) %>%
  dplyr::arrange(
    dataset,
    comparison,
    peptide,
    ms_source_id
  )


# 11.7 Check diagnostic classifications -------------------------------------------------------------
overlap_class_qc <-
  ms_source_peptide_overlap_status %>%
  dplyr::group_by(
    MS_MHC_overlap_class
  ) %>%
  dplyr::summarise(
    n_MS_source_peptide_rows =
      dplyr::n(),
    
    n_unique_peptides =
      dplyr::n_distinct(
        peptide
      ),
    
    n_unique_MS_source_ids =
      dplyr::n_distinct(
        ms_source_id
      ),
    
    n_unique_MS_transcripts =
      dplyr::n_distinct(
        ms_transcript_id
      ),
    
    .groups = "drop"
  )


# exact_source_id_and_peptide
# Exact source_id and peptide match: strict support
# 
# same_transcript_different_ORF
# Peptide and transcript match, but the ORF differs
# 
# same_peptide_different_transcript_only
# Peptide sequence matches, but the transcript differs
# 
# no_MHC_binder_for_peptide
# The MS peptide has no plus-strand binder prediction


# 11.8 Check consistency with the strict-join results -----------------------------------------------------

strict_overlap_keys <-
  ms_mhc_strict_overlap %>%
  dplyr::distinct(
    dataset,
    dataset_type,
    comparison,
    peptide,
    source_id
  ) %>%
  dplyr::rename(
    ms_source_id =
      source_id
  )

diagnostic_exact_keys <-
  ms_source_peptide_overlap_status %>%
  dplyr::filter(
    MS_MHC_overlap_class ==
      "exact_source_id_and_peptide"
  ) %>%
  dplyr::distinct(
    dataset,
    dataset_type,
    comparison,
    peptide,
    ms_source_id
  )

# Compare both result sets
strict_not_in_diagnostic <-
  strict_overlap_keys %>%
  dplyr::anti_join(
    diagnostic_exact_keys,
    by = c(
      "dataset",
      "dataset_type",
      "comparison",
      "peptide",
      "ms_source_id"
    )
  )

diagnostic_not_in_strict <-
  diagnostic_exact_keys %>%
  dplyr::anti_join(
    strict_overlap_keys,
    by = c(
      "dataset",
      "dataset_type",
      "comparison",
      "peptide",
      "ms_source_id"
    )
  )

# Generate QC:


overlap_consistency_qc <-
  tibble::tibble(
    check_item = c(
      "n_MS_input_rows",
      "n_strict_overlap_keys",
      "n_diagnostic_exact_keys",
      "n_strict_not_in_diagnostic",
      "n_diagnostic_not_in_strict"
    ),
    
    value = c(
      nrow(ms_overlap_input),
      nrow(strict_overlap_keys),
      nrow(diagnostic_exact_keys),
      nrow(strict_not_in_diagnostic),
      nrow(diagnostic_not_in_strict)
    )
  )


# Add hard checks

if (
  nrow(strict_not_in_diagnostic) > 0 ||
  nrow(diagnostic_not_in_strict) > 0
) {
  stop(
    "Peptide-only diagnostics are inconsistent with the strict source_id join results."
  )
}


# 11.9 Count peptides supported at different mapping levels

overlap_peptide_counts <-
  tibble::tibble(
    support_level = c(
      "any_peptide_sequence_overlap",
      "same_transcript_or_exact_source",
      "exact_source_id_and_peptide",
      "same_transcript_different_ORF_only",
      "different_transcript_only",
      "no_MHC_binder"
    ),
    
    n_unique_peptides = c(
      ms_source_peptide_overlap_status %>%
        dplyr::filter(
          n_all_MHC_source_matches > 0
        ) %>%
        dplyr::summarise(
          n = dplyr::n_distinct(peptide)
        ) %>%
        dplyr::pull(n),
      
      ms_source_peptide_overlap_status %>%
        dplyr::filter(
          any_same_transcript
        ) %>%
        dplyr::summarise(
          n = dplyr::n_distinct(peptide)
        ) %>%
        dplyr::pull(n),
      
      ms_source_peptide_overlap_status %>%
        dplyr::filter(
          any_same_source_id
        ) %>%
        dplyr::summarise(
          n = dplyr::n_distinct(peptide)
        ) %>%
        dplyr::pull(n),
      
      ms_source_peptide_overlap_status %>%
        dplyr::filter(
          MS_MHC_overlap_class ==
            "same_transcript_different_ORF"
        ) %>%
        dplyr::summarise(
          n = dplyr::n_distinct(peptide)
        ) %>%
        dplyr::pull(n),
      
      ms_source_peptide_overlap_status %>%
        dplyr::filter(
          MS_MHC_overlap_class ==
            "same_peptide_different_transcript_only"
        ) %>%
        dplyr::summarise(
          n = dplyr::n_distinct(peptide)
        ) %>%
        dplyr::pull(n),
      
      ms_source_peptide_overlap_status %>%
        dplyr::filter(
          MS_MHC_overlap_class ==
            "no_MHC_binder_for_peptide"
        ) %>%
        dplyr::summarise(
          n = dplyr::n_distinct(peptide)
        ) %>%
        dplyr::pull(n)
    )
  )


# 11.10 Save results --------------------------------------------------------------

data.table::fwrite(
  peptide_only_overlap_mapping,
  file.path(
    lnc_out_dir,
    "08_peptide_only_MS_MHC_overlap_mapping.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  na = "NA"
)

data.table::fwrite(
  ms_source_peptide_overlap_status,
  file.path(
    lnc_out_dir,
    "08a_MS_source_peptide_overlap_classification.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  na = "NA"
)

data.table::fwrite(
  overlap_class_qc,
  file.path(
    lnc_out_dir,
    "08b_MS_MHC_overlap_class_QC.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  na = "NA"
)

data.table::fwrite(
  overlap_consistency_qc,
  file.path(
    lnc_out_dir,
    "08c_strict_vs_peptide_overlap_consistency_QC.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  na = "NA"
)

data.table::fwrite(
  overlap_peptide_counts,
  file.path(
    lnc_out_dir,
    "08d_MS_MHC_overlap_peptide_counts.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  na = "NA"
)


# Step 12: Build strict candidates and a one-row-per-peptide master table -----------------------------------------------

# 12.1 Build strict candidate ORF-peptide keys -------------------------------------------------

exact_candidate_pairs <-
  ms_mhc_strict_overlap %>%
  dplyr::select(
    source_id,
    orf_id,
    transcript_id,
    transcript_id_no_version,
    orf_start,
    orf_end,
    orf_strand,
    orf_partial,
    orf_completion_status,
    netmhcpan_id,
    peptide
  ) %>%
  dplyr::distinct() %>%
  dplyr::arrange(
    peptide,
    source_id
  )


# Check how many ORFs map to each peptide:
candidate_source_multiplicity <-
  exact_candidate_pairs %>%
  dplyr::group_by(
    peptide
  ) %>%
  dplyr::summarise(
    n_exact_source_ids =
      dplyr::n_distinct(source_id),
    
    n_exact_transcripts =
      dplyr::n_distinct(transcript_id),
    
    exact_source_ids =
      collapse_unique(source_id),
    
    exact_transcript_ids =
      collapse_unique(transcript_id),
    
    exact_netmhcpan_ids =
      collapse_unique(netmhcpan_id),
    
    .groups = "drop"
  )


# 12.2 Extract all candidate PSM mapping details

exact_candidate_psm_mapping_detail <-
  main_plus_psm_source_detail %>%
  dplyr::semi_join(
    exact_candidate_pairs %>%
      dplyr::select(
        source_id,
        peptide
      ),
    by = c(
      "source_id",
      "peptide"
    )
  ) %>%
  dplyr::arrange(
    peptide,
    source_id,
    dataset,
    comparison,
    psm_key
  )


# 12.3 Keep one row per PSM for peptide-level PSM counting -----------------------------------------

exact_candidate_psm_unique <-
  exact_candidate_psm_mapping_detail %>%
  dplyr::arrange(
    peptide,
    psm_key,
    source_id
  ) %>%
  dplyr::group_by(
    peptide,
    psm_key
  ) %>%
  dplyr::slice_head(
    n = 1
  ) %>%
  dplyr::ungroup()


# Step 13: Build a one-row-per-source_id + peptide candidate table ------------------------------------

# 13.1 Select a representative PSM at the ORF-peptide level ---------------------------------------------
candidate_representative_psm_source_peptide <-
  exact_candidate_psm_mapping_detail %>%
  dplyr::mutate(
    abs_mass_error_ppm =
      abs(mass_error_ppm),
    
    PIF_sort =
      dplyr::coalesce(
        PIF,
        -Inf
      ),
    
    abs_mass_error_sort =
      dplyr::coalesce(
        abs_mass_error_ppm,
        Inf
      )
  ) %>%
  dplyr::arrange(
    source_id,
    peptide,
    PEP,
    dplyr::desc(Score),
    dplyr::desc(Delta_score),
    dplyr::desc(PIF_sort),
    abs_mass_error_sort,
    dataset,
    comparison,
    psm_key
  ) %>%
  dplyr::group_by(
    source_id,
    peptide
  ) %>%
  dplyr::slice_head(
    n = 1
  ) %>%
  dplyr::ungroup() %>%
  dplyr::transmute(
    source_id,
    peptide,
    
    representative_dataset =
      dataset,
    
    representative_dataset_type =
      dataset_type,
    
    representative_comparison =
      comparison,
    
    representative_psm_key =
      psm_key,
    
    representative_raw_file =
      raw_file,
    
    representative_scan_number =
      scan_number,
    
    representative_modified_sequence =
      modified_sequence,
    
    representative_charge =
      charge,
    
    representative_PEP =
      PEP,
    
    representative_Score =
      Score,
    
    representative_Delta_score =
      Delta_score,
    
    representative_PIF =
      PIF,
    
    representative_mass_error_ppm =
      mass_error_ppm,
    
    representative_abs_mass_error_ppm =
      abs_mass_error_ppm,
    
    representative_leading_razor_protein =
      leading_razor_protein,
    
    representative_leading_razor_class =
      leading_razor_class,
    
    representative_source_is_leading_razor =
      current_source_is_leading_razor
  )


# 13.2 Summarize MS evidence at the ORF-peptide level -----------------------------------------------
candidate_source_peptide_MS <-
  exact_candidate_psm_mapping_detail %>%
  dplyr::group_by(
    source_id,
    orf_id,
    transcript_id,
    transcript_id_no_version,
    orf_start,
    orf_end,
    orf_strand,
    orf_partial,
    orf_completion_status,
    netmhcpan_id,
    peptide
  ) %>%
  dplyr::summarise(
    peptide_length =
      dplyr::first(peptide_length),
    
    n_PSM =
      dplyr::n_distinct(psm_key),
    
    n_PSM_TMT10 =
      dplyr::n_distinct(
        psm_key[
          dataset_type == "TMT10"
        ]
      ),
    
    n_PSM_LabelFree =
      dplyr::n_distinct(
        psm_key[
          dataset_type == "LabelFree"
        ]
      ),
    
    n_datasets =
      dplyr::n_distinct(dataset),
    
    datasets =
      collapse_unique(dataset),
    
    n_dataset_comparison_groups =
      dplyr::n_distinct(
        paste(
          dataset,
          comparison,
          sep = "|"
        )
      ),
    
    dataset_comparison_groups =
      collapse_unique(
        paste(
          dataset,
          comparison,
          sep = "|"
        )
      ),
    
    n_raw_files =
      dplyr::n_distinct(raw_file),
    
    raw_files =
      collapse_unique(raw_file),
    
    n_charge_states =
      dplyr::n_distinct(
        charge[
          !is.na(charge)
        ]
      ),
    
    charge_states =
      collapse_unique(
        charge
      ),
    
    min_PEP =
      min(
        PEP,
        na.rm = TRUE
      ),
    
    max_Score =
      max(
        Score,
        na.rm = TRUE
      ),
    
    max_Delta_score =
      max(
        Delta_score,
        na.rm = TRUE
      ),
    
    max_PIF =
      if (all(is.na(PIF))) {
        NA_real_
      } else {
        max(
          PIF,
          na.rm = TRUE
        )
      },
    
    n_PSM_source_is_leading =
      dplyr::n_distinct(
        psm_key[
          current_source_is_leading_razor
        ]
      ),
    
    source_is_leading_in_any_PSM =
      any(
        current_source_is_leading_razor,
        na.rm = TRUE
      ),
    
    source_is_leading_in_all_PSM =
      all(
        current_source_is_leading_razor
      ),
    
    n_PSM_with_SwissProt_mapping =
      dplyr::n_distinct(
        psm_key[
          has_SwissProt_mapping
        ]
      ),
    
    n_PSM_with_TrEMBL_mapping =
      dplyr::n_distinct(
        psm_key[
          has_TrEMBL_mapping
        ]
      ),
    
    any_SwissProt_mapping =
      any(
        has_SwissProt_mapping,
        na.rm = TRUE
      ),
    
    all_PSM_have_SwissProt_mapping =
      all(
        has_SwissProt_mapping
      ),
    
    any_TrEMBL_mapping =
      any(
        has_TrEMBL_mapping,
        na.rm = TRUE
      ),
    
    SwissProt_accessions =
      collapse_unique(
        SwissProt_accessions
      ),
    
    TrEMBL_accessions =
      collapse_unique(
        TrEMBL_accessions
      ),
    
    leading_razor_proteins =
      collapse_unique(
        leading_razor_protein
      ),
    
    leading_razor_classes =
      collapse_unique(
        leading_razor_class
      ),
    
    .groups = "drop"
  ) %>%
  dplyr::left_join(
    candidate_representative_psm_source_peptide,
    by = c(
      "source_id",
      "peptide"
    )
  )


# 13.3 Add NetMHCpan annotations ------------------------------------------------------

candidate_source_peptide_core <-
  candidate_source_peptide_MS %>%
  dplyr::left_join(
    netmhcpan_annotation_for_join,
    by = c(
      "source_id",
      "peptide"
    )
  ) %>%
  dplyr::mutate(
    MS_database_mapping_class =
      dplyr::case_when(
        any_SwissProt_mapping &
          any_TrEMBL_mapping ~
          "lncRNA_ORF_plus_SwissProt_plus_TrEMBL",
        
        any_SwissProt_mapping ~
          "lncRNA_ORF_plus_SwissProt",
        
        any_TrEMBL_mapping ~
          "lncRNA_ORF_plus_TrEMBL",
        
        TRUE ~
          "lncRNA_ORF_only_in_MaxQuant_Proteins"
      ),
    
    leading_razor_support_class =
      dplyr::case_when(
        source_is_leading_in_all_PSM ~
          "source_is_leading_in_all_PSM",
        
        source_is_leading_in_any_PSM ~
          "source_is_leading_in_some_PSM",
        
        TRUE ~
          "source_is_never_leading"
      )
  ) %>%
  dplyr::arrange(
    peptide,
    source_id
  )

# Step 14: Build a one-row-per-peptide master table ----------------------------------------------------

# 14.1 Summarize ORF and transcript information for candidate peptides ----------------------------------------------

candidate_peptide_source_summary <-
  exact_candidate_pairs %>%
  dplyr::group_by(
    peptide
  ) %>%
  dplyr::summarise(
    peptide_length =
      nchar(
        dplyr::first(peptide)
      ),
    
    n_exact_source_ids =
      dplyr::n_distinct(source_id),
    
    exact_source_ids =
      collapse_unique(source_id),
    
    n_exact_transcripts =
      dplyr::n_distinct(transcript_id),
    
    exact_transcript_ids =
      collapse_unique(transcript_id),
    
    exact_ORF_ids =
      collapse_unique(orf_id),
    
    exact_netmhcpan_ids =
      collapse_unique(netmhcpan_id),
    
    ORF_partial_statuses =
      collapse_unique(
        as.character(orf_partial)
      ),
    
    ORF_completion_statuses =
      collapse_unique(
        orf_completion_status
      ),
    
    .groups = "drop"
  )


# 14.2 Summarize PSM-level MS evidence for candidate peptides ------------------------------------------------------
candidate_peptide_MS_summary <-
  exact_candidate_psm_unique %>%
  dplyr::group_by(
    peptide
  ) %>%
  dplyr::summarise(
    n_PSM_total =
      dplyr::n_distinct(psm_key),
    
    n_PSM_TMT10 =
      dplyr::n_distinct(
        psm_key[
          dataset_type == "TMT10"
        ]
      ),
    
    n_PSM_LabelFree =
      dplyr::n_distinct(
        psm_key[
          dataset_type == "LabelFree"
        ]
      ),
    
    n_datasets =
      dplyr::n_distinct(dataset),
    
    datasets =
      collapse_unique(dataset),
    
    n_dataset_comparison_groups =
      dplyr::n_distinct(
        paste(
          dataset,
          comparison,
          sep = "|"
        )
      ),
    
    dataset_comparison_groups =
      collapse_unique(
        paste(
          dataset,
          comparison,
          sep = "|"
        )
      ),
    
    n_raw_files =
      dplyr::n_distinct(raw_file),
    
    raw_files =
      collapse_unique(raw_file),
    
    n_charge_states =
      dplyr::n_distinct(
        charge[
          !is.na(charge)
        ]
      ),
    
    charge_states =
      collapse_unique(charge),
    
    min_PEP =
      min(
        PEP,
        na.rm = TRUE
      ),
    
    max_Score =
      max(
        Score,
        na.rm = TRUE
      ),
    
    max_Delta_score =
      max(
        Delta_score,
        na.rm = TRUE
      ),
    
    max_PIF =
      if (all(is.na(PIF))) {
        NA_real_
      } else {
        max(
          PIF,
          na.rm = TRUE
        )
      },
    
    n_PSM_with_SwissProt_mapping =
      dplyr::n_distinct(
        psm_key[
          has_SwissProt_mapping
        ]
      ),
    
    n_PSM_with_TrEMBL_mapping =
      dplyr::n_distinct(
        psm_key[
          has_TrEMBL_mapping
        ]
      ),
    
    any_SwissProt_mapping =
      any(
        has_SwissProt_mapping,
        na.rm = TRUE
      ),
    
    all_PSM_have_SwissProt_mapping =
      all(
        has_SwissProt_mapping
      ),
    
    any_TrEMBL_mapping =
      any(
        has_TrEMBL_mapping,
        na.rm = TRUE
      ),
    
    leading_razor_proteins =
      collapse_unique(
        leading_razor_protein
      ),
    
    leading_razor_classes =
      collapse_unique(
        leading_razor_class
      ),
    
    .groups = "drop"
  )

# Whether a candidate source is the leading razor protein must be calculated from the full mapping detail:

candidate_peptide_source_leading_summary <-
  exact_candidate_psm_mapping_detail %>%
  dplyr::group_by(
    peptide
  ) %>%
  dplyr::summarise(
    n_PSM_candidate_source_is_leading =
      dplyr::n_distinct(
        psm_key[
          current_source_is_leading_razor
        ]
      ),
    
    any_candidate_source_is_leading =
      any(
        current_source_is_leading_razor,
        na.rm = TRUE
      ),
    
    .groups = "drop"
  )


# 14.3 Select a representative PSM at the peptide level ------------------------------------------------------


candidate_representative_psm_peptide <-
  exact_candidate_psm_unique %>%
  dplyr::mutate(
    abs_mass_error_ppm =
      abs(mass_error_ppm),
    
    PIF_sort =
      dplyr::coalesce(
        PIF,
        -Inf
      ),
    
    abs_mass_error_sort =
      dplyr::coalesce(
        abs_mass_error_ppm,
        Inf
      )
  ) %>%
  dplyr::arrange(
    peptide,
    PEP,
    dplyr::desc(Score),
    dplyr::desc(Delta_score),
    dplyr::desc(PIF_sort),
    abs_mass_error_sort,
    dataset,
    comparison,
    psm_key
  ) %>%
  dplyr::group_by(
    peptide
  ) %>%
  dplyr::slice_head(
    n = 1
  ) %>%
  dplyr::ungroup() %>%
  dplyr::transmute(
    peptide,
    
    representative_dataset =
      dataset,
    
    representative_dataset_type =
      dataset_type,
    
    representative_comparison =
      comparison,
    
    representative_psm_key =
      psm_key,
    
    representative_raw_file =
      raw_file,
    
    representative_scan_number =
      scan_number,
    
    representative_modified_sequence =
      modified_sequence,
    
    representative_charge =
      charge,
    
    representative_PEP =
      PEP,
    
    representative_Score =
      Score,
    
    representative_Delta_score =
      Delta_score,
    
    representative_PIF =
      PIF,
    
    representative_mass_error_ppm =
      mass_error_ppm,
    
    representative_abs_mass_error_ppm =
      abs_mass_error_ppm,
    
    representative_leading_razor_protein =
      leading_razor_protein,
    
    representative_leading_razor_class =
      leading_razor_class
  )


# 14.4 Build the peptide-level MHC summary -------------------------------------------------------

candidate_MHC_allele_detail <-
  netmhcpan_allele_dedup %>%
  dplyr::semi_join(
    exact_candidate_pairs %>%
      dplyr::select(
        source_id,
        peptide
      ),
    by = c(
      "source_id",
      "peptide"
    )
  )

# Select the best peptide-level EL result ------------------------------------------------------
candidate_best_EL_peptide <-
  candidate_MHC_allele_detail %>%
  dplyr::mutate(
    EL_rank_sort =
      dplyr::coalesce(
        EL_rank,
        Inf
      ),
    
    EL_score_sort =
      dplyr::coalesce(
        EL_score,
        -Inf
      )
  ) %>%
  dplyr::arrange(
    peptide,
    EL_rank_sort,
    dplyr::desc(EL_score_sort),
    HLA_allele,
    source_id
  ) %>%
  dplyr::group_by(
    peptide
  ) %>%
  dplyr::slice_head(
    n = 1
  ) %>%
  dplyr::ungroup() %>%
  dplyr::transmute(
    peptide,
    
    best_EL_source_id =
      source_id,
    
    best_EL_allele =
      HLA_allele,
    
    best_EL_rank =
      EL_rank,
    
    best_EL_score =
      EL_score,
    
    best_EL_binder_class =
      EL_binder_class,
    
    best_EL_position =
      peptide_position,
    
    best_EL_binding_core =
      binding_core,
    
    best_EL_binding_icore =
      binding_icore
  )


# Select the best BA result ----------------------------------------------------------------

candidate_best_BA_peptide <-
  candidate_MHC_allele_detail %>%
  dplyr::mutate(
    BA_rank_sort =
      dplyr::coalesce(
        BA_rank,
        Inf
      ),
    
    BA_score_sort =
      dplyr::coalesce(
        BA_score,
        -Inf
      )
  ) %>%
  dplyr::arrange(
    peptide,
    BA_rank_sort,
    dplyr::desc(BA_score_sort),
    HLA_allele,
    source_id
  ) %>%
  dplyr::group_by(
    peptide
  ) %>%
  dplyr::slice_head(
    n = 1
  ) %>%
  dplyr::ungroup() %>%
  dplyr::transmute(
    peptide,
    
    best_BA_source_id =
      source_id,
    
    best_BA_allele =
      HLA_allele,
    
    best_BA_rank =
      BA_rank,
    
    best_BA_score =
      BA_score,
    
    best_BA_position =
      peptide_position
  )

# Summarize all HLA alleles: ----------------------------------------------------------------

candidate_peptide_MHC_summary <-
  candidate_MHC_allele_detail %>%
  dplyr::group_by(
    peptide
  ) %>%
  dplyr::summarise(
    n_MHC_source_ids =
      dplyr::n_distinct(source_id),
    
    MHC_source_ids =
      collapse_unique(source_id),
    
    n_HLA_alleles =
      dplyr::n_distinct(HLA_allele),
    
    HLA_alleles =
      collapse_unique(HLA_allele),
    
    n_source_HLA_pairs =
      dplyr::n(),
    
    n_SB_source_HLA_pairs =
      sum(
        EL_binder_class == "SB"
      ),
    
    n_WB_source_HLA_pairs =
      sum(
        EL_binder_class == "WB"
      ),
    
    min_EL_rank_across_all_sources_alleles =
      min(
        EL_rank,
        na.rm = TRUE
      ),
    
    max_EL_score_across_all_sources_alleles =
      max(
        EL_score,
        na.rm = TRUE
      ),
    
    min_BA_rank_across_all_sources_alleles =
      if (all(is.na(BA_rank))) {
        NA_real_
      } else {
        min(
          BA_rank,
          na.rm = TRUE
        )
      },
    
    max_BA_score_across_all_sources_alleles =
      if (all(is.na(BA_score))) {
        NA_real_
      } else {
        max(
          BA_score,
          na.rm = TRUE
        )
      },
    
    .groups = "drop"
  ) %>%
  dplyr::left_join(
    candidate_best_EL_peptide,
    by = "peptide"
  ) %>%
  dplyr::left_join(
    candidate_best_BA_peptide,
    by = "peptide"
  )


# 14.5 Generate the formal one-row-per-peptide master table --------------------------------------------------

candidate_peptide_master <-
  candidate_peptide_source_summary %>%
  dplyr::left_join(
    candidate_peptide_MS_summary,
    by = "peptide"
  ) %>%
  dplyr::left_join(
    candidate_peptide_source_leading_summary,
    by = "peptide"
  ) %>%
  dplyr::left_join(
    candidate_representative_psm_peptide,
    by = "peptide"
  ) %>%
  dplyr::left_join(
    candidate_peptide_MHC_summary,
    by = "peptide"
  ) %>%
  dplyr::mutate(
    strict_MS_MHC_support =
      TRUE,
    
    source_mapping_uniqueness =
      dplyr::case_when(
        n_exact_source_ids == 1 ~
          "single_exact_ORF",
        
        n_exact_transcripts == 1 ~
          "multiple_ORFs_same_transcript",
        
        TRUE ~
          "multiple_transcripts"
      ),
    
    MS_database_mapping_class =
      dplyr::case_when(
        any_SwissProt_mapping &
          any_TrEMBL_mapping ~
          "lncRNA_ORF_plus_SwissProt_plus_TrEMBL",
        
        any_SwissProt_mapping ~
          "lncRNA_ORF_plus_SwissProt",
        
        any_TrEMBL_mapping ~
          "lncRNA_ORF_plus_TrEMBL",
        
        TRUE ~
          "lncRNA_ORF_only_in_MaxQuant_Proteins"
      ),
    
    source_interpretation =
      dplyr::case_when(
        !any_SwissProt_mapping &
          !any_TrEMBL_mapping ~
          "no_canonical_protein_mapping_in_MaxQuant",
        
        any_candidate_source_is_leading ~
          "candidate_ORF_leading_but_shared_protein_mapping",
        
        TRUE ~
          "shared_with_canonical_protein_and_candidate_not_leading"
      )
  ) %>%
  dplyr::arrange(
    best_EL_rank,
    peptide
  )


# 14.6 Hard QC for the one-row-per-peptide table ----------------------------------------------------

candidate_peptide_duplicate_qc <-
  candidate_peptide_master %>%
  dplyr::count(
    peptide,
    name = "n"
  ) %>%
  dplyr::filter(
    n > 1
  )

#
candidate_peptide_master_qc <-
  tibble::tibble(
    check_item = c(
      "n_candidate_peptide_rows",
      "n_unique_candidate_peptides",
      "n_duplicate_peptide_keys",
      "n_missing_MS_evidence",
      "n_missing_MHC_evidence",
      "n_missing_source_id",
      "n_missing_transcript_id",
      "n_missing_representative_PSM",
      "n_non_strict_support",
      "n_exact_source_peptide_pairs"
    ),
    
    value = c(
      nrow(
        candidate_peptide_master
      ),
      
      dplyr::n_distinct(
        candidate_peptide_master$peptide
      ),
      
      nrow(
        candidate_peptide_duplicate_qc
      ),
      
      sum(
        is.na(
          candidate_peptide_master$n_PSM_total
        ) |
          candidate_peptide_master$n_PSM_total < 1
      ),
      
      sum(
        is.na(
          candidate_peptide_master$best_EL_allele
        )
      ),
      
      sum(
        is.na(
          candidate_peptide_master$exact_source_ids
        )
      ),
      
      sum(
        is.na(
          candidate_peptide_master$exact_transcript_ids
        )
      ),
      
      sum(
        is.na(
          candidate_peptide_master$
            representative_psm_key
        )
      ),
      
      sum(
        !candidate_peptide_master$
          strict_MS_MHC_support
      ),
      
      nrow(
        exact_candidate_pairs
      )
    )
  )


# 14.8 Save core candidate results -----------------------------------------------------------

data.table::fwrite(
  exact_candidate_pairs,
  file.path(
    lnc_out_dir,
    "09_exact_candidate_source_peptide_pairs.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  na = "NA"
)

data.table::fwrite(
  exact_candidate_psm_mapping_detail,
  file.path(
    lnc_out_dir,
    "09a_exact_candidate_PSM_source_detail.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  na = "NA"
)

data.table::fwrite(
  candidate_source_peptide_core,
  file.path(
    lnc_out_dir,
    "09b_exact_candidate_source_peptide_core.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  na = "NA"
)

data.table::fwrite(
  candidate_peptide_master,
  file.path(
    lnc_out_dir,
    "09c_exact_candidate_one_row_per_peptide_core.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  na = "NA"
)

data.table::fwrite(
  candidate_peptide_master_qc,
  file.path(
    lnc_out_dir,
    "09d_exact_candidate_one_row_per_peptide_QC.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  na = "NA"
)

data.table::fwrite(
  candidate_source_multiplicity,
  file.path(
    lnc_out_dir,
    "09e_exact_candidate_source_multiplicity.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  na = "NA"
)


# Correction: retain the nine-peptide strict MHC subset

strict_MHC_candidate_peptide_subset <-
  candidate_peptide_master

rm(candidate_peptide_master)


# Save as an explicit intermediate result
data.table::fwrite(
  strict_MHC_candidate_peptide_subset,
  file.path(
    lnc_out_dir,
    "09c_strict_MHC_candidate_peptide_subset.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  na = "NA"
)

# Build the complete 143-peptide master analysis table --------------------------------------------------


# 1. Count each PSM only once -----------------------------------------------------------

all_MS_PSM_unique <- main_plus_psm_source_detail %>%
  dplyr::arrange(
    peptide,
    psm_key,
    source_id
  ) %>%
  dplyr::group_by(
    peptide,
    psm_key
  ) %>%
  dplyr::slice_head(n = 1) %>%
  dplyr::ungroup()

# 2. Summarize ORFs and transcripts for each peptide -----------------------------------------
all_MS_peptide_source_summary <-
  main_plus_psm_source_detail %>%
  dplyr::group_by(
    peptide
  ) %>%
  dplyr::summarise(
    peptide_length =
      dplyr::first(peptide_length),
    
    n_MS_source_ids =
      dplyr::n_distinct(source_id),
    
    MS_source_ids =
      collapse_unique(source_id),
    
    n_MS_transcripts =
      dplyr::n_distinct(transcript_id),
    
    MS_transcript_ids =
      collapse_unique(transcript_id),
    
    MS_ORF_ids =
      collapse_unique(orf_id),
    
    MS_netmhcpan_ids =
      collapse_unique(netmhcpan_id),
    
    ORF_strands =
      collapse_unique(orf_strand),
    
    ORF_partial_statuses =
      collapse_unique(
        as.character(orf_partial)
      ),
    
    ORF_completion_statuses =
      collapse_unique(
        orf_completion_status
      ),
    
    n_partial_ORFs =
      dplyr::n_distinct(
        source_id[
          !is.na(orf_partial) &
            orf_partial
        ]
      ),
    
    .groups = "drop"
  )


# 3. Summarize MS evidence for all 143 peptides -------------------------------------------------

all_MS_peptide_evidence_summary <-
  all_MS_PSM_unique %>%
  dplyr::group_by(
    peptide
  ) %>%
  dplyr::summarise(
    n_PSM_total =
      dplyr::n_distinct(psm_key),
    
    n_PSM_TMT10 =
      dplyr::n_distinct(
        psm_key[
          dataset_type == "TMT10"
        ]
      ),
    
    n_PSM_LabelFree =
      dplyr::n_distinct(
        psm_key[
          dataset_type == "LabelFree"
        ]
      ),
    
    n_datasets =
      dplyr::n_distinct(dataset),
    
    datasets =
      collapse_unique(dataset),
    
    n_dataset_comparison_groups =
      dplyr::n_distinct(
        paste(
          dataset,
          comparison,
          sep = "|"
        )
      ),
    
    dataset_comparison_groups =
      collapse_unique(
        paste(
          dataset,
          comparison,
          sep = "|"
        )
      ),
    
    n_raw_files =
      dplyr::n_distinct(raw_file),
    
    raw_files =
      collapse_unique(raw_file),
    
    n_charge_states =
      dplyr::n_distinct(
        charge[
          !is.na(charge)
        ]
      ),
    
    charge_states =
      collapse_unique(charge),
    
    min_PEP =
      min(
        PEP,
        na.rm = TRUE
      ),
    
    max_Score =
      max(
        Score,
        na.rm = TRUE
      ),
    
    max_Delta_score =
      max(
        Delta_score,
        na.rm = TRUE
      ),
    
    max_PIF =
      if (all(is.na(PIF))) {
        NA_real_
      } else {
        max(PIF, na.rm = TRUE)
      },
    
    n_PSM_with_SwissProt_mapping =
      dplyr::n_distinct(
        psm_key[
          has_SwissProt_mapping
        ]
      ),
    
    n_PSM_with_TrEMBL_mapping =
      dplyr::n_distinct(
        psm_key[
          has_TrEMBL_mapping
        ]
      ),
    
    any_SwissProt_mapping =
      any(
        has_SwissProt_mapping,
        na.rm = TRUE
      ),
    
    all_PSM_have_SwissProt_mapping =
      all(
        has_SwissProt_mapping
      ),
    
    any_TrEMBL_mapping =
      any(
        has_TrEMBL_mapping,
        na.rm = TRUE
      ),
    
    any_lncRNA_ORF_is_leading =
      any(
        leading_razor_class ==
          "lncRNA_ORF",
        na.rm = TRUE
      ),
    
    n_PSM_lncRNA_ORF_is_leading =
      dplyr::n_distinct(
        psm_key[
          leading_razor_class ==
            "lncRNA_ORF"
        ]
      ),
    
    leading_razor_proteins =
      collapse_unique(
        leading_razor_protein
      ),
    
    leading_razor_classes =
      collapse_unique(
        leading_razor_class
      ),
    
    .groups = "drop"
  )


# 4. Select a representative PSM for each of the 143 peptides -----------------------------------------------
all_MS_representative_PSM <-
  all_MS_PSM_unique %>%
  dplyr::mutate(
    abs_mass_error_ppm =
      abs(mass_error_ppm),
    
    PIF_sort =
      dplyr::coalesce(
        PIF,
        -Inf
      ),
    
    abs_mass_error_sort =
      dplyr::coalesce(
        abs_mass_error_ppm,
        Inf
      )
  ) %>%
  dplyr::arrange(
    peptide,
    PEP,
    dplyr::desc(Score),
    dplyr::desc(Delta_score),
    dplyr::desc(PIF_sort),
    abs_mass_error_sort,
    dataset,
    comparison,
    psm_key
  ) %>%
  dplyr::group_by(
    peptide
  ) %>%
  dplyr::slice_head(n = 1) %>%
  dplyr::ungroup() %>%
  dplyr::transmute(
    peptide,
    
    representative_dataset =
      dataset,
    
    representative_dataset_type =
      dataset_type,
    
    representative_comparison =
      comparison,
    
    representative_psm_key =
      psm_key,
    
    representative_raw_file =
      raw_file,
    
    representative_scan_number =
      scan_number,
    
    representative_scan_index =
      scan_index,
    
    representative_modified_sequence =
      modified_sequence,
    
    representative_charge =
      charge,
    
    representative_PEP =
      PEP,
    
    representative_Score =
      Score,
    
    representative_Delta_score =
      Delta_score,
    
    representative_PIF =
      PIF,
    
    representative_mass_error_ppm =
      mass_error_ppm,
    
    representative_abs_mass_error_ppm =
      abs_mass_error_ppm,
    
    representative_leading_razor_protein =
      leading_razor_protein,
    
    representative_leading_razor_class =
      leading_razor_class
  )

# 5. Prepare annotations for the nine peptides with strict MHC support ------------------------------------------------
strict_MHC_annotation <-
  candidate_peptide_MHC_summary %>%
  dplyr::mutate(
    strict_MS_MHC_support = TRUE
  )


# 6. Build the complete base master table ----------------------------------------------------------
lncRNA_peptide_master_base <-
  all_MS_peptide_source_summary %>%
  dplyr::left_join(
    all_MS_peptide_evidence_summary,
    by = "peptide"
  ) %>%
  dplyr::left_join(
    all_MS_representative_PSM,
    by = "peptide"
  ) %>%
  dplyr::left_join(
    strict_MHC_annotation,
    by = "peptide"
  ) %>%
  dplyr::mutate(
    strict_MS_MHC_support =
      tidyr::replace_na(
        strict_MS_MHC_support,
        FALSE
      ),
    
    MS_MHC_support_class =
      dplyr::if_else(
        strict_MS_MHC_support,
        "exact_source_id_and_peptide_binder",
        "no_exact_plus_strand_MHC_binder"
      ),
    
    source_mapping_class =
      dplyr::case_when(
        n_MS_source_ids == 1 ~
          "single_plus_ORF",
        
        n_MS_transcripts == 1 ~
          "multiple_plus_ORFs_same_transcript",
        
        TRUE ~
          "multiple_transcripts"
      ),
    
    MS_database_mapping_class =
      dplyr::case_when(
        any_SwissProt_mapping &
          any_TrEMBL_mapping ~
          "lncRNA_ORF_plus_SwissProt_plus_TrEMBL",
        
        any_SwissProt_mapping ~
          "lncRNA_ORF_plus_SwissProt",
        
        any_TrEMBL_mapping ~
          "lncRNA_ORF_plus_TrEMBL",
        
        TRUE ~
          "lncRNA_ORF_only_in_MaxQuant_Proteins"
      ),
    
    preliminary_source_interpretation =
      dplyr::case_when(
        !any_SwissProt_mapping &
          !any_TrEMBL_mapping ~
          "no_canonical_mapping_in_MaxQuant_Proteins",
        
        any_lncRNA_ORF_is_leading ~
          "lncRNA_ORF_leading_but_shared_with_canonical_protein",
        
        TRUE ~
          "shared_with_canonical_protein_and_lncRNA_ORF_not_leading"
      )
  ) %>%
  dplyr::arrange(
    dplyr::desc(strict_MS_MHC_support),
    best_EL_rank,
    min_PEP,
    dplyr::desc(max_Score),
    peptide
  )

# 7. QC for the corrected master table -------------------------------------------------------------
lncRNA_master_base_duplicate_qc <-
  lncRNA_peptide_master_base %>%
  dplyr::count(
    peptide,
    name = "n"
  ) %>%
  dplyr::filter(
    n > 1
  )

expected_strict_MHC_peptides <-
  dplyr::n_distinct(
    ms_mhc_strict_overlap$peptide
  )

lncRNA_master_base_qc <-
  tibble::tibble(
    check_item = c(
      "n_master_rows",
      "n_unique_peptides",
      "n_duplicate_peptide_keys",
      "n_strict_MHC_supported_peptides",
      "n_without_strict_MHC_support",
      "expected_strict_MHC_peptides",
      "n_missing_MS_evidence",
      "n_missing_source_mapping",
      "n_missing_representative_PSM",
      "n_non_plus_ORF_status"
    ),
    
    value = c(
      nrow(
        lncRNA_peptide_master_base
      ),
      
      dplyr::n_distinct(
        lncRNA_peptide_master_base$peptide
      ),
      
      nrow(
        lncRNA_master_base_duplicate_qc
      ),
      
      sum(
        lncRNA_peptide_master_base$
          strict_MS_MHC_support
      ),
      
      sum(
        !lncRNA_peptide_master_base$
          strict_MS_MHC_support
      ),
      
      expected_strict_MHC_peptides,
      
      sum(
        is.na(
          lncRNA_peptide_master_base$
            n_PSM_total
        ) |
          lncRNA_peptide_master_base$
          n_PSM_total < 1
      ),
      
      sum(
        is.na(
          lncRNA_peptide_master_base$
            MS_source_ids
        )
      ),
      
      sum(
        is.na(
          lncRNA_peptide_master_base$
            representative_psm_key
        )
      ),
      
      sum(
        lncRNA_peptide_master_base$
          ORF_strands != "plus"
      )
    )
  )


# Add hard checks： -----------------------------------------------------------------
stopifnot(
  nrow(lncRNA_peptide_master_base) ==
    dplyr::n_distinct(
      lncRNA_peptide_master_base$peptide
    ),
  
  nrow(lncRNA_master_base_duplicate_qc) == 0,
  
  sum(
    lncRNA_peptide_master_base$
      strict_MS_MHC_support
  ) ==
    expected_strict_MHC_peptides
)

# Save results ----------------------------------------------------------------------

data.table::fwrite(
  lncRNA_peptide_master_base,
  file.path(
    lnc_out_dir,
    "10_all_high_confidence_lncRNA_peptide_master_base.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  na = "NA"
)

data.table::fwrite(
  lncRNA_master_base_qc,
  file.path(
    lnc_out_dir,
    "10a_all_high_confidence_lncRNA_peptide_master_base_QC.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  na = "NA"
)


# Prepare publication-ready results -----------------------------------------------------------------

lnc_master <- lncRNA_peptide_master_base


#
comparison_summary <-
  all_MS_PSM_unique %>%
  dplyr::group_by(peptide) %>%
  dplyr::summarise(
    comparisons =
      collapse_unique(comparison),
    .groups = "drop"
  )


# 3. Complete mass-error fields -------------------------------------------------------------
if (
  !"median_LFQ_abs_mass_error_ppm" %in%
  names(lnc_master)
) {
  lnc_master$
    median_LFQ_abs_mass_error_ppm <-
    NA_real_
}

if (
  !"mass_error_reporting_status" %in%
  names(lnc_master)
) {
  lnc_master <-
    lnc_master %>%
    dplyr::mutate(
      mass_error_reporting_status =
        dplyr::case_when(
          n_PSM_LabelFree == 0 ~
            "not_reported_TMT_only",
          
          n_PSM_LabelFree > 0 &
            is.na(
              median_LFQ_abs_mass_error_ppm
            ) ~
            "LabelFree_mass_error_not_yet_calculated",
          
          TRUE ~
            "reported_from_LabelFree_PSM"
        )
    )
}

# 4. Generate lncRNA source classifications ---------------------------------------------------------

lnc_master <-
  lnc_master %>%
  dplyr::mutate(
    source_uniqueness_class =
      dplyr::case_when(
        stringr::str_detect(
          MS_database_mapping_class,
          "SwissProt"
        ) &
          stringr::str_detect(
            MS_database_mapping_class,
            "TrEMBL"
          ) ~
          "shared_with_SwissProt_and_TrEMBL",
        
        stringr::str_detect(
          MS_database_mapping_class,
          "SwissProt"
        ) ~
          "shared_with_SwissProt",
        
        stringr::str_detect(
          MS_database_mapping_class,
          "TrEMBL"
        ) ~
          "shared_with_TrEMBL",
        
        MS_database_mapping_class ==
          "lncRNA_ORF_only_in_MaxQuant_Proteins" ~
          paste0(
            "no_SwissProt_or_TrEMBL_mapping_",
            "in_MaxQuant_database"
          ),
        
        TRUE ~
          "source_unresolved"
      )
  )


# 5. Generate source conclusions and manual-review priorities -------------------------------------------------------
lnc_master <-
  lnc_master %>%
  dplyr::mutate(
    PIF_availability =
      dplyr::case_when(
        is.na(max_PIF) ~
          "not_available",
        
        TRUE ~
          "available"
      ),
    
    has_TMT_PSM =
      n_PSM_TMT10 > 0,
    
    has_LFQ_PSM =
      n_PSM_LabelFree > 0,
    
    supports_natural_HLA_presentation =
      FALSE,
    
    source_conclusion =
      dplyr::case_when(
        source_uniqueness_class %in%
          c(
            "shared_with_SwissProt",
            "shared_with_TrEMBL",
            "shared_with_SwissProt_and_TrEMBL"
          ) ~
          paste0(
            "peptide sequence shared with a known ",
            "protein database entry; shotgun MS ",
            "cannot establish lncRNA ORF-specific origin"
          ),
        
        source_uniqueness_class ==
          paste0(
            "no_SwissProt_or_TrEMBL_mapping_",
            "in_MaxQuant_database"
          ) ~
          paste0(
            "no Swiss-Prot or TrEMBL mapping was ",
            "observed in the searched MaxQuant ",
            "database; lncRNA ORF-specific origin ",
            "is not yet proven"
          ),
        
        TRUE ~
          "peptide source remains unresolved"
      ),
    
    recommended_claim =
      dplyr::case_when(
        strict_MS_MHC_support ~
          paste0(
            "MS-supported and NetMHCpan-predicted ",
            "HLA-binding peptide"
          ),
        
        TRUE ~
          paste0(
            "MS-supported peptide without strict ",
            "source-consistent HLA binder prediction"
          )
      ),
    
    manual_spectrum_review_priority =
      dplyr::case_when(
        strict_MS_MHC_support &
          source_uniqueness_class ==
          paste0(
            "no_SwissProt_or_TrEMBL_mapping_",
            "in_MaxQuant_database"
          ) &
          best_EL_binder_class == "SB" ~
          "P1_review_noncanonical_SB",
        
        strict_MS_MHC_support &
          source_uniqueness_class ==
          paste0(
            "no_SwissProt_or_TrEMBL_mapping_",
            "in_MaxQuant_database"
          ) ~
          "P2_review_noncanonical_WB",
        
        strict_MS_MHC_support &
          best_EL_binder_class == "SB" ~
          "P3_review_shared_source_SB",
        
        strict_MS_MHC_support &
          best_EL_binder_class == "WB" ~
          "P4_review_shared_source_WB",
        
        TRUE ~
          "P5_review_MS_only"
      )
  )

# 6. Generate the final publication-ready lncRNA summary table -----------------------------------------------------
lncRNA_publication_table <-
  lnc_master %>%
  dplyr::left_join(
    comparison_summary,
    by = "peptide"
  ) %>%
  dplyr::filter(
    strict_MS_MHC_support
  ) %>%
  dplyr::transmute(
    peptide,
    
    # lncRNA and ORF source information
    lncRNA_transcript_ids =
      MS_transcript_ids,
    
    lncRNA_ORF_ids =
      MS_ORF_ids,
    
    lncRNA_mapping_class =
      source_mapping_class,
    
    ORF_completion_statuses,
    
    # MS evidence: column names aligned with the circRNA workflow
    n_PSM =
      n_PSM_total,
    
    n_raw_files,
    datasets,
    comparisons,
    
    min_PEP,
    max_Score,
    max_Delta_score,
    max_PIF,
    PIF_availability,
    
    # Mass error: column names aligned with the circRNA workflow
    median_LFQ_abs_mass_error_ppm,
    mass_error_reporting_status,
    
    # HLA prediction: column names aligned with the circRNA workflow
    HLA_alleles,
    n_HLA_alleles,
    best_EL_allele,
    best_EL_rank,
    best_EL_score,
    best_EL_binder_class,
    
    # Source attribution and interpretation
    strict_MS_MHC_support,
    source_uniqueness_class,
    source_conclusion,
    
    supports_natural_HLA_presentation,
    recommended_claim,
    manual_spectrum_review_priority
  ) %>%
  dplyr::arrange(
    factor(
      manual_spectrum_review_priority,
      levels = c(
        "P1_review_noncanonical_SB",
        "P2_review_noncanonical_WB",
        "P3_review_shared_source_SB",
        "P4_review_shared_source_WB",
        "P5_review_MS_only"
      )
    ),
    best_EL_rank,
    min_PEP,
    dplyr::desc(max_Score),
    peptide
  )

# 7. Final QC -----------------------------------------------------------------

lncRNA_publication_qc <-
  tibble::tibble(
    check_item = c(
      "n_rows",
      "n_unique_peptides",
      "n_duplicate_peptides",
      "n_missing_transcript_ids",
      "n_missing_ORF_ids",
      "n_missing_best_EL_allele",
      "n_missing_best_EL_rank",
      "n_non_strict_MHC_rows"
    ),
    
    value = c(
      nrow(
        lncRNA_publication_table
      ),
      
      dplyr::n_distinct(
        lncRNA_publication_table$peptide
      ),
      
      nrow(
        lncRNA_publication_table %>%
          dplyr::count(
            peptide,
            name = "n"
          ) %>%
          dplyr::filter(n > 1)
      ),
      
      sum(
        is.na(
          lncRNA_publication_table$
            lncRNA_transcript_ids
        )
      ),
      
      sum(
        is.na(
          lncRNA_publication_table$
            lncRNA_ORF_ids
        )
      ),
      
      sum(
        is.na(
          lncRNA_publication_table$
            best_EL_allele
        )
      ),
      
      sum(
        is.na(
          lncRNA_publication_table$
            best_EL_rank
        )
      ),
      
      sum(
        !lncRNA_publication_table$
          strict_MS_MHC_support
      )
    )
  )


# Add hard checks ------------------------------------------------------------------
stopifnot(
  nrow(lncRNA_publication_table) ==
    dplyr::n_distinct(
      lncRNA_publication_table$peptide
    ),
  
  all(
    lncRNA_publication_table$
      strict_MS_MHC_support
  ),
  
  all(
    !is.na(
      lncRNA_publication_table$
        best_EL_allele
    )
  ),
  
  all(
    !is.na(
      lncRNA_publication_table$
        best_EL_rank
    )
  )
)

# 8. Save the publication-ready summary table --------------------------------------------------------------

data.table::fwrite(
  lncRNA_publication_table,
  file.path(
    lnc_out_dir,
    paste0(
      "19a_lncRNA_publication_",
      "candidate_summary_HLA_fixed.tsv"
    )
  ),
  sep = "\t",
  quote = FALSE,
  na = "NA"
)


# Save resultsQC： -------------------------------------------------------------------

data.table::fwrite(
  lncRNA_publication_qc,
  file.path(
    lnc_out_dir,
    paste0(
      "19a_lncRNA_publication_",
      "candidate_summary_HLA_fixed_QC.tsv"
    )
  ),
  sep = "\t",
  quote = FALSE,
  na = "NA"
)


# Prepare peptide lists for UniProt queries and downstream analyses -----------------------------------------------------
X19a_lncRNA_publication_candidate_summary_HLA_fixed <- read_delim("19a_lncRNA_publication_candidate_summary_HLA_fixed.tsv",
                                                                  delim = "\t", escape_double = FALSE,
                                                                  trim_ws = TRUE)


peptide <- paste(X19a_lncRNA_publication_candidate_summary_HLA_fixed$peptide, collapse = ",")


# Process downloaded web-query result files ---------------------------------------------------------------

# Process UniProt query results ----------------------------------------------------------


uniprot_result <- function(input_file, output_dir) {
  
  # Read the current Excel input file
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
    " ; unique peptides: ", dplyr::n_distinct(peptide_detail$Peptide)
  )
  
  return(peptide_class)
}
# Run
input_dir <- "./MS_unspecific/lncRNA/result/IEDB_uniprot_search"


output_dir <- "./MS_unspecific/lncRNA/result/IEDB_uniprot_search"

sample_names <- "lncRNA_uniprot"

full_file_paths <- file.path(
  input_dir,
  paste0(sample_names, ".xlsx")
)

results <- lapply(
  full_file_paths,
  uniprot_result,
  output_dir = output_dir
)

names(results) <- sample_names

# Add Swiss-Prot, TrEMBL, or isoform annotations ---------------------------------------------------------

# Annotate Swiss-Prot versus TrEMBL ---------------------------------------------------------------

# Match Swiss-Prot and TrEMBL IDs --------------------------------------------------------------

library(dplyr)
library(tidyr)
library(stringr)

library(data.table)

# Convert directly to data.table to minimize unnecessary object copying

lncRNA_uniprot_result <- read_csv("./MS_unspecific/lncRNA/result/IEDB_uniprot_search/lncRNA_uniprot_output.csv")
setDT(lncRNA_uniprot_result)

sp_tr_map <- read.delim("./MS_unspecific/sp_tr_map.txt")

setDT(sp_tr_map)

# Add row identifiers to the original table
lncRNA_uniprot_result[, row_id := .I]

# 1. Split Matched_Entries
lncRNA_entries_long <- lncRNA_uniprot_result[
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

lncRNA_entries_long[
  ,
  accession := sub(
    "^(sp|tr)\\|",
    "",
    accession
  )
]

lncRNA_entries_long[
  ,
  parent_accession := sub(
    "-[0-9]+$",
    "",
    accession
  )
]

# 2. Extract accessions that require lookup
needed_entries <- unique(
  c(
    lncRNA_entries_long$accession,
    lncRNA_entries_long$parent_accession
  )
)


# 3. Filter the large mapping table to relevant entries only

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

# 4. Perform exact accession matching
lncRNA_entries_long[
  uniprot_map_small,
  on = .(accession = Entry),
  database_exact := i.database
]

# 5. Use canonical accessions to annotate isoforms

lncRNA_entries_long[
  uniprot_map_small,
  on = .(parent_accession = Entry),
  database_parent := i.database
]

lncRNA_entries_long[
  ,
  database_final := fifelse(
    !is.na(database_exact),
    database_exact,
    database_parent
  )
]
# Add sp| or tr| prefixes:
lncRNA_entries_long[
  ,
  accession_with_db := fifelse(
    !is.na(database_final),
    paste0(database_final, "|", accession),
    accession
  )
]

# 6. Annotate match type
lncRNA_entries_long[
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
# 7. Collapse matches back to semicolon-separated format

lncRNA_entries_summary <- lncRNA_entries_long[
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

# 8. Merge annotations back into the original table
all_lncRNA_search_with_db <- merge(
  lncRNA_uniprot_result,
  lncRNA_entries_summary,
  by = "row_id",
  all.x = TRUE,
  sort = FALSE
)

setorder(
  all_lncRNA_search_with_db,
  row_id
)

all_lncRNA_search_with_db[
  ,
  row_id := NULL
]


# Replace the original column
all_lncRNA_search_with_db[
  ,
  Matched_Entries_original := Matched_Entries
]

all_lncRNA_search_with_db[
  ,
  Matched_Entries := Matched_Entries_with_db
]

all_lncRNA_search_with_db[
  ,
  Matched_Entries_with_db := NULL
]


write.table(all_lncRNA_search_with_db,"./MS_unspecific/lncRNA/result/IEDB_uniprot_search/all_lncRNA_search_sp_tr.txt",quote = FALSE,row.names = FALSE,sep = "\t")


# Process UniProt query results ----------------------------------------------------------

# Prepare the publication-ready table ---------------------------------------------------------------

tidy_uniprot <- function(df) {
  
  df %>%
    dplyr::select(
      Peptide,
      Matched_Entries,
      Matched_Genes
    ) %>%
    
    # One UniProt entry per row
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

#

# Run ----------------------------------------------------------------------
all_lncRNA_search_sp_tr <- read.delim2("./MS_unspecific/lncRNA/result/IEDB_uniprot_search/all_lncRNA_search_sp_tr.txt")

lncRNA_uniprot_publication  <- tidy_uniprot(
  all_lncRNA_search_sp_tr
)


write.table(lncRNA_uniprot_publication ,"./MS_unspecific/lncRNA/result/IEDB_uniprot_search/new_all_lncRNA_search_sp_tr.txt",quote = FALSE,row.names = FALSE,sep = "\t")

write.csv(lncRNA_uniprot_publication,"./MS_unspecific/lncRNA/result/IEDB_uniprot_search/new_all_lncRNA_search_sp_tr.csv",quote = FALSE,row.names = FALSE,sep = ",")
