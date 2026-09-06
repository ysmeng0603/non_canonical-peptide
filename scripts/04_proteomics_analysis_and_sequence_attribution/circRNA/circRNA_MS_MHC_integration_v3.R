#!/usr/bin/env Rscript

# circRNA_MS_MHC_integration.R
#
# Purpose
#   1. Read MaxQuant msms/evidence/peptides files for PXD037581 (TMT10)
#      and PXD044963 (label-free).
#   2. Apply PSM-level quality filters.
#   3. Retain every lcl|ORF... candidate mapping while preserving shared
#      Swiss-Prot/TrEMBL matches and the MaxQuant leading razor protein.
#   4. Normalize and deduplicate NetMHCpan predictions.
#   5. Join MS and MHC evidence by exact source_id + peptide.
#   6. Export a compact candidate table for subsequent BSJ annotation.
#
# This script does NOT perform BSJ, UniProt, or IEDB annotation.

options(stringsAsFactors = FALSE)

required_packages <- c("data.table", "dplyr", "tidyr", "stringr", "purrr", "tibble")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages) > 0L) {
  stop(
    "Missing R packages: ", paste(missing_packages, collapse = ", "),
    "\nInstall them before running this script."
  )
}

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(purrr)
  library(tibble)
})

# ============================================================
# 0. Configuration: edit this section when paths differ
# ============================================================
CFG <- list(
  mq_result_root = "./MS_unspecific/circRNA/result",
  search_group = "contain_circRNA",

  dataset_design = tibble::tribble(
    ~dataset,     ~dataset_type,
    "PXD037581", "TMT10",
    "PXD044963", "LabelFree"
  ),

  comparisons = c("HC_eRA", "HC_RA", "RA_eRA"),

  # Change this path if combined_f.txt is stored elsewhere.
  netmhcpan_file = "./circRNA_peptide3/part4_netMHCpan/MHC_1/combined_f.txt",

  out_dir = "./MS_unspecific/circRNA/result/circRNA_MS_MHC_integration",

  # MaxQuant PSM filters.
  min_peptide_length = 7L,
  max_peptide_length = Inf,
  pep_threshold = 0.01,
  score_threshold = 40,
  delta_score_threshold = 15,

  # "strict": Delta score must be present and >= threshold.
  # "allow_na": missing Delta score is retained; observed values must pass.
  # "ignore": Delta score is not used for filtering.
  delta_score_policy = "strict",

  # The circRNA input FASTA is expected to be in transcript orientation.
  # ORFfinder also predicts reverse-complement ORFs; the main analysis keeps plus strand.
  main_orf_strands = c("plus"),
  keep_unknown_strand_in_main = FALSE,

  # Writing every high-confidence non-candidate PSM can create a very large file.
  write_all_high_confidence_psm = FALSE,

  # NetMHCpan thresholds based on EL rank.
  mhc_strong_rank = 0.5,
  mhc_weak_rank = 2.0,
  keep_only_mhc_binders = TRUE
)

# ============================================================
# 1. General helpers
# ============================================================
dir.create(CFG$out_dir, recursive = TRUE, showWarnings = FALSE)

write_tsv <- function(x, filename) {
  path <- file.path(CFG$out_dir, filename)
  data.table::fwrite(
    data.table::as.data.table(x),
    file = path,
    sep = "\t",
    quote = FALSE,
    na = "NA"
  )
  message("Wrote: ", path)
  invisible(path)
}

normalize_peptide <- function(x) {
  x <- as.character(x)
  x <- stringr::str_replace_all(x, "\\s+", "")
  stringr::str_to_upper(x)
}

normalize_source_id <- function(x) {
  x <- stringr::str_trim(as.character(x))
  x <- stringr::str_remove(x, "^lcl\\|")
  x <- stringr::str_remove(x, ":[0-9]+:[0-9]+.*$")
  x <- stringr::str_replace(
    x,
    "^ORF([0-9]+)_([0-9]+)$",
    "ORF\\1_hsa_circ_\\2"
  )
  x <- stringr::str_replace(
    x,
    "^ORF([0-9]+)_circ_([0-9]+)$",
    "ORF\\1_hsa_circ_\\2"
  )
  x <- stringr::str_replace(
    x,
    "^ORF([0-9]+)_hsa_circ_hsa_circ_([0-9]+)$",
    "ORF\\1_hsa_circ_\\2"
  )
  x
}

collapse_unique <- function(x, sep = ";") {
  x <- unique(as.character(x[!is.na(x) & as.character(x) != ""]))
  if (length(x) == 0L) NA_character_ else paste(x, collapse = sep)
}

collapse_sorted_numeric <- function(x, sep = ";") {
  x <- suppressWarnings(as.numeric(x))
  x <- sort(unique(x[is.finite(x)]))
  if (length(x) == 0L) NA_character_ else paste(x, collapse = sep)
}

safe_min <- function(x) {
  x <- suppressWarnings(as.numeric(x))
  if (length(x) == 0L || all(is.na(x))) NA_real_ else min(x, na.rm = TRUE)
}

safe_max <- function(x) {
  x <- suppressWarnings(as.numeric(x))
  if (length(x) == 0L || all(is.na(x))) NA_real_ else max(x, na.rm = TRUE)
}

get_col <- function(df, name, default = NA) {
  if (name %in% names(df)) df[[name]] else rep(default, nrow(df))
}

require_columns <- function(df, columns, label) {
  missing <- setdiff(columns, names(df))
  if (length(missing) > 0L) {
    stop(label, " is missing required columns: ", paste(missing, collapse = ", "))
  }
}

flag_is_plus <- function(x) {
  !is.na(x) & trimws(as.character(x)) == "+"
}

# Read only selected MaxQuant columns. Names are matched after make.names().
read_selected_mq <- function(path, exact_clean_names, regex_clean_names = character()) {
  if (!file.exists(path)) stop("File does not exist: ", path)

  raw_names <- names(data.table::fread(path, nrows = 0L, check.names = FALSE))
  clean_names <- make.names(raw_names, unique = TRUE)

  keep <- clean_names %in% exact_clean_names
  if (length(regex_clean_names) > 0L) {
    for (pattern in regex_clean_names) {
      keep <- keep | grepl(pattern, clean_names, perl = TRUE)
    }
  }

  if (!any(keep)) stop("No requested columns were found in: ", path)

  selected_raw_names <- raw_names[keep]
  data.table::fread(
    path,
    select = selected_raw_names,
    data.table = FALSE,
    check.names = TRUE,
    na.strings = c("", "NA", "NaN")
  )
}

classify_protein_string <- function(proteins) {
  proteins <- ifelse(is.na(proteins), "", proteins)
  tibble(
    hit_any_candidate = str_detect(proteins, "(^|;)lcl\\|"),
    hit_any_sp = str_detect(proteins, "(^|;)sp\\|"),
    hit_any_tr = str_detect(proteins, "(^|;)tr\\|"),
    hit_any_reverse = str_detect(proteins, "(^|;)REV__"),
    hit_any_contaminant = str_detect(proteins, "(^|;)(CON__|[^;]*CON_)")
  )
}

parse_candidate_accession <- function(candidate_accession) {
  source_id_raw <- str_match(candidate_accession, "^lcl\\|([^:[:space:]]+)")[, 2]
  source_id <- normalize_source_id(source_id_raw)
  coords <- str_match(candidate_accession, "^lcl\\|[^:[:space:]]+:([0-9]+):([0-9]+)")
  start0 <- suppressWarnings(as.integer(coords[, 2]))
  end0 <- suppressWarnings(as.integer(coords[, 3]))

  strand <- case_when(
    !is.na(start0) & !is.na(end0) & start0 < end0 ~ "plus",
    !is.na(start0) & !is.na(end0) & start0 > end0 ~ "minus",
    TRUE ~ "unknown"
  )

  tibble(
    source_id = source_id,
    orf_id = str_extract(source_id, "^ORF[0-9]+"),
    circRNA_id = str_extract(source_id, "hsa_circ_[0-9]+"),
    orf_start_0based = start0,
    orf_end_0based = end0,
    orf_strand = strand
  )
}

is_main_strand <- function(strand) {
  strand %in% CFG$main_orf_strands |
    (CFG$keep_unknown_strand_in_main & strand == "unknown")
}

# ============================================================
# 2. Read one MaxQuant result and create PSM-level candidate mappings
# ============================================================
read_maxquant_one <- function(dataset, dataset_type, comparison) {
  txt_dir <- file.path(
    CFG$mq_result_root,
    CFG$search_group,
    dataset,
    paste0(comparison, "_result"),
    "combined",
    "txt"
  )

  msms_file <- file.path(txt_dir, "msms.txt")
  evidence_file <- file.path(txt_dir, "evidence.txt")
  peptides_file <- file.path(txt_dir, "peptides.txt")

  message("Reading MaxQuant: ", dataset, " / ", comparison)

  msms_exact <- c(
    "Raw.file", "Scan.number", "Scan.index", "Sequence", "Length",
    "Modifications", "Modified.sequence", "Proteins", "Charge",
    "Fragmentation", "Mass.analyzer", "m.z", "Mass",
    "Mass.error..ppm.", "Mass.error..Da.", "Retention.time",
    "PEP", "Score", "Delta.score", "PIF", "Fraction.of.total.spectrum",
    "Base.peak.fraction", "Precursor.Intensity", "Matches",
    "Number.of.matches", "Intensity.coverage", "Peak.coverage",
    "Reverse", "Contaminant", "id", "Peptide.ID", "Evidence.ID"
  )

  msms <- read_selected_mq(
    msms_file,
    exact_clean_names = msms_exact,
    regex_clean_names = c("^Reporter\\.intensity\\.corrected\\.[0-9]+$")
  )

  require_columns(
    msms,
    c("Raw.file", "Scan.number", "Sequence", "Proteins", "PEP", "Score", "Delta.score"),
    paste0("msms.txt [", dataset, "/", comparison, "]")
  )

  evidence_exact <- c(
    "id", "Leading.proteins", "Leading.razor.protein", "Intensity",
    "Experiment", "Raw.file", "Fraction", "PIF", "Mass.error..ppm.",
    "Potential.contaminant", "Reverse"
  )
  evidence <- read_selected_mq(evidence_file, evidence_exact)

  peptides_exact <- c(
    "id", "Leading.razor.protein", "Unique..Groups.", "Unique..Proteins.",
    "MS.MS.Count", "PEP", "Score", "Proteins"
  )
  peptides <- read_selected_mq(peptides_file, peptides_exact)

  # Standard PSM table.
  psm <- tibble(
    dataset = dataset,
    dataset_type = dataset_type,
    comparison = comparison,
    search_group = CFG$search_group,
    msms_id = suppressWarnings(as.integer(get_col(msms, "id"))),
    evidence_id = suppressWarnings(as.integer(get_col(msms, "Evidence.ID"))),
    peptide_id = suppressWarnings(as.integer(get_col(msms, "Peptide.ID"))),
    raw_file = as.character(get_col(msms, "Raw.file")),
    scan_number = suppressWarnings(as.integer(get_col(msms, "Scan.number"))),
    scan_index = suppressWarnings(as.integer(get_col(msms, "Scan.index"))),
    peptide = normalize_peptide(get_col(msms, "Sequence")),
    modified_sequence = as.character(get_col(msms, "Modified.sequence")),
    modifications = as.character(get_col(msms, "Modifications")),
    peptide_length = suppressWarnings(as.integer(get_col(msms, "Length"))),
    charge = suppressWarnings(as.integer(get_col(msms, "Charge"))),
    fragmentation = as.character(get_col(msms, "Fragmentation")),
    mass_analyzer = as.character(get_col(msms, "Mass.analyzer")),
    precursor_mz = suppressWarnings(as.numeric(get_col(msms, "m.z"))),
    peptide_mass = suppressWarnings(as.numeric(get_col(msms, "Mass"))),
    mass_error_ppm = suppressWarnings(as.numeric(get_col(msms, "Mass.error..ppm."))),
    mass_error_da = suppressWarnings(as.numeric(get_col(msms, "Mass.error..Da."))),
    retention_time = suppressWarnings(as.numeric(get_col(msms, "Retention.time"))),
    PEP = suppressWarnings(as.numeric(get_col(msms, "PEP"))),
    Score = suppressWarnings(as.numeric(get_col(msms, "Score"))),
    Delta_score = suppressWarnings(as.numeric(get_col(msms, "Delta.score"))),
    PIF = suppressWarnings(as.numeric(get_col(msms, "PIF"))),
    fraction_total_spectrum = suppressWarnings(as.numeric(get_col(msms, "Fraction.of.total.spectrum"))),
    base_peak_fraction = suppressWarnings(as.numeric(get_col(msms, "Base.peak.fraction"))),
    precursor_intensity = suppressWarnings(as.numeric(get_col(msms, "Precursor.Intensity"))),
    matched_ions = as.character(get_col(msms, "Matches")),
    number_of_matches = suppressWarnings(as.integer(get_col(msms, "Number.of.matches"))),
    intensity_coverage = suppressWarnings(as.numeric(get_col(msms, "Intensity.coverage"))),
    peak_coverage = suppressWarnings(as.numeric(get_col(msms, "Peak.coverage"))),
    proteins_original = as.character(get_col(msms, "Proteins")),
    reverse_flag = as.character(get_col(msms, "Reverse")),
    contaminant_flag = as.character(get_col(msms, "Contaminant"))
  )

  reporter_cols <- grep(
    "^Reporter\\.intensity\\.corrected\\.[0-9]+$",
    names(msms),
    value = TRUE
  )
  if (length(reporter_cols) > 0L) {
    psm <- dplyr::bind_cols(psm, msms[, reporter_cols, drop = FALSE])
  }

  # Join evidence metadata.
  evidence_std <- tibble(
    evidence_id = suppressWarnings(as.integer(get_col(evidence, "id"))),
    evidence_leading_proteins = as.character(get_col(evidence, "Leading.proteins")),
    evidence_leading_razor = as.character(get_col(evidence, "Leading.razor.protein")),
    evidence_intensity = suppressWarnings(as.numeric(get_col(evidence, "Intensity"))),
    experiment = as.character(get_col(evidence, "Experiment")),
    evidence_raw_file = as.character(get_col(evidence, "Raw.file")),
    fraction = as.character(get_col(evidence, "Fraction")),
    evidence_PIF = suppressWarnings(as.numeric(get_col(evidence, "PIF")))
  ) %>%
    dplyr::distinct(evidence_id, .keep_all = TRUE)

  peptide_std <- tibble(
    peptide_id = suppressWarnings(as.integer(get_col(peptides, "id"))),
    peptide_leading_razor = as.character(get_col(peptides, "Leading.razor.protein")),
    unique_groups = as.character(get_col(peptides, "Unique..Groups.")),
    unique_proteins = as.character(get_col(peptides, "Unique..Proteins.")),
    peptide_table_msms_count = suppressWarnings(as.integer(get_col(peptides, "MS.MS.Count")))
  ) %>%
    dplyr::distinct(peptide_id, .keep_all = TRUE)

  psm <- psm %>%
    dplyr::left_join(evidence_std, by = "evidence_id") %>%
    dplyr::left_join(peptide_std, by = "peptide_id") %>%
    dplyr::mutate(
      leading_razor_protein = coalesce(evidence_leading_razor, peptide_leading_razor),
      psm_key = paste(dataset, comparison, raw_file, scan_number, msms_id, sep = "|")
    )

  protein_flags <- classify_protein_string(psm$proteins_original)
  psm <- dplyr::bind_cols(psm, protein_flags)

  delta_pass <- switch(
    CFG$delta_score_policy,
    strict = !is.na(psm$Delta_score) & psm$Delta_score >= CFG$delta_score_threshold,
    allow_na = is.na(psm$Delta_score) | psm$Delta_score >= CFG$delta_score_threshold,
    ignore = rep(TRUE, nrow(psm)),
    stop("Unknown CFG$delta_score_policy: ", CFG$delta_score_policy)
  )

  psm <- psm %>%
    dplyr::mutate(
      pass_not_reverse = !flag_is_plus(reverse_flag) & !hit_any_reverse,
      pass_not_contaminant = !flag_is_plus(contaminant_flag) & !hit_any_contaminant,
      pass_sequence = !is.na(peptide) & peptide != "",
      pass_length = !is.na(peptide_length) &
        peptide_length >= CFG$min_peptide_length &
        peptide_length <= CFG$max_peptide_length,
      pass_pep = !is.na(PEP) & PEP < CFG$pep_threshold,
      pass_score = !is.na(Score) & Score > CFG$score_threshold,
      pass_delta = delta_pass,
      pass_all = pass_not_reverse & pass_not_contaminant & pass_sequence &
        pass_length & pass_pep & pass_score & pass_delta
    )

  audit <- tibble(
    dataset = dataset,
    dataset_type = dataset_type,
    comparison = comparison,
    n_msms_rows = nrow(psm),
    n_not_reverse = sum(psm$pass_not_reverse, na.rm = TRUE),
    n_not_contaminant = sum(psm$pass_not_contaminant, na.rm = TRUE),
    n_length_pass = sum(psm$pass_length, na.rm = TRUE),
    n_pep_pass = sum(psm$pass_pep, na.rm = TRUE),
    n_score_pass = sum(psm$pass_score, na.rm = TRUE),
    n_delta_pass = sum(psm$pass_delta, na.rm = TRUE),
    n_high_confidence_psm = sum(psm$pass_all, na.rm = TRUE)
  )

  high_conf <- psm %>% dplyr::filter(pass_all)

  # One row per PSM x protein mapping.
  protein_long <- high_conf %>%
    dplyr::select(-starts_with("pass_")) %>%
    dplyr::mutate(
      all_protein_hits = proteins_original,
      protein_hit = proteins_original
    ) %>%
    tidyr::separate_rows(protein_hit, sep = ";") %>%
    dplyr::mutate(
      protein_hit = str_trim(protein_hit),
      protein_class = case_when(
        str_starts(protein_hit, fixed("lcl|")) ~ "candidate",
        str_starts(protein_hit, fixed("sp|")) ~ "SwissProt",
        str_starts(protein_hit, fixed("tr|")) ~ "TrEMBL",
        str_starts(protein_hit, fixed("REV__")) ~ "reverse",
        str_detect(protein_hit, "CON_") ~ "contaminant",
        TRUE ~ "other"
      )
    )

  candidate_all <- protein_long %>%
    dplyr::filter(protein_class == "candidate") %>%
    dplyr::rename(candidate_accession = protein_hit)

  parsed <- parse_candidate_accession(candidate_all$candidate_accession)
  candidate_all <- dplyr::bind_cols(candidate_all, parsed) %>%
    dplyr::filter(!is.na(source_id), source_id != "") %>%
    dplyr::mutate(
      candidate_is_leading = if_else(
        is.na(leading_razor_protein),
        FALSE,
        str_detect(leading_razor_protein, fixed(paste0("lcl|", source_id)))
      ),
      leading_is_candidate = !is.na(leading_razor_protein) & str_starts(leading_razor_protein, fixed("lcl|")),
      leading_is_sp = !is.na(leading_razor_protein) & str_starts(leading_razor_protein, fixed("sp|")),
      leading_is_tr = !is.na(leading_razor_protein) & str_starts(leading_razor_protein, fixed("tr|")),
      mq_source_class = case_when(
        hit_any_sp & hit_any_tr ~ "candidate_shared_with_sp_and_tr",
        hit_any_sp ~ "candidate_shared_with_sp",
        hit_any_tr ~ "candidate_shared_with_tr",
        TRUE ~ "candidate_only_in_search_database"
      ),
      main_strand_pass = is_main_strand(orf_strand)
    )

  candidate_main <- candidate_all %>% dplyr::filter(main_strand_pass)

  audit <- audit %>%
    dplyr::mutate(
      n_high_confidence_psm_with_candidate = n_distinct(candidate_all$psm_key),
      n_candidate_psm_source_rows_all_strands = nrow(candidate_all),
      n_candidate_psm_source_rows_main_strand = nrow(candidate_main),
      n_candidate_source_peptide_all_strands = n_distinct(paste(candidate_all$source_id, candidate_all$peptide)),
      n_candidate_source_peptide_main_strand = n_distinct(paste(candidate_main$source_id, candidate_main$peptide))
    )

  list(
    audit = audit,
    high_conf = if (CFG$write_all_high_confidence_psm) high_conf else NULL,
    candidate_all = candidate_all,
    candidate_main = candidate_main
  )
}

# ============================================================
# 3. Run all MaxQuant dataset/comparison combinations
# ============================================================
run_grid <- tidyr::crossing(
  dataset = CFG$dataset_design$dataset,
  comparison = CFG$comparisons
) %>%
  dplyr::left_join(CFG$dataset_design, by = "dataset") %>%
  dplyr::select(dataset, dataset_type, comparison)

mq_results <- purrr::pmap(
  run_grid,
  function(dataset, dataset_type, comparison) {
    read_maxquant_one(dataset, dataset_type, comparison)
  }
)

mq_audit <- dplyr::bind_rows(map(mq_results, "audit"))
mq_high_conf <- if (CFG$write_all_high_confidence_psm) {
  dplyr::bind_rows(map(mq_results, "high_conf"))
} else {
  tibble()
}
mq_candidate_all <- dplyr::bind_rows(map(mq_results, "candidate_all"))
mq_candidate_main <- dplyr::bind_rows(map(mq_results, "candidate_main"))

write_tsv(mq_audit, "00_MaxQuant_filter_audit.tsv")
if (CFG$write_all_high_confidence_psm) {
  write_tsv(mq_high_conf, "01_high_confidence_PSM_all_proteins.tsv")
}
write_tsv(mq_candidate_all, "02_candidate_PSM_source_all_strands.tsv")
write_tsv(mq_candidate_main, "03_candidate_PSM_source_main_strand.tsv")

# ============================================================
# 4. Summarize MS evidence by dataset + comparison + source_id + peptide
# ============================================================
ms_source_summary <- mq_candidate_main %>%
  dplyr::group_by(dataset, dataset_type, comparison, source_id, circRNA_id, orf_id, peptide) %>%
  dplyr::summarise(
    orf_strand = collapse_unique(orf_strand),
    orf_start_0based = safe_min(orf_start_0based),
    orf_end_0based = safe_max(orf_end_0based),
    candidate_accessions = collapse_unique(candidate_accession),
    n_candidate_accessions = n_distinct(candidate_accession),
    n_PSM = n_distinct(psm_key),
    n_evidence = n_distinct(evidence_id[!is.na(evidence_id)]),
    n_raw_files = n_distinct(raw_file[!is.na(raw_file)]),
    raw_files = collapse_unique(raw_file),
    scan_numbers = collapse_unique(paste0(raw_file, ":", scan_number)),
    charges = collapse_unique(charge),
    best_PEP = safe_min(PEP),
    best_Score = safe_max(Score),
    best_Delta_score = safe_max(Delta_score),
    best_abs_mass_error_ppm = safe_min(abs(mass_error_ppm)),
    best_PIF = safe_max(PIF),
    max_number_of_matches = safe_max(number_of_matches),
    max_intensity_coverage = safe_max(intensity_coverage),
    max_peak_coverage = safe_max(peak_coverage),
    hit_any_sp = any(hit_any_sp, na.rm = TRUE),
    hit_any_tr = any(hit_any_tr, na.rm = TRUE),
    candidate_is_leading_any = any(candidate_is_leading, na.rm = TRUE),
    leading_is_candidate_any = any(leading_is_candidate, na.rm = TRUE),
    leading_is_sp_any = any(leading_is_sp, na.rm = TRUE),
    leading_is_tr_any = any(leading_is_tr, na.rm = TRUE),
    leading_razor_proteins = collapse_unique(leading_razor_protein),
    mq_source_class = case_when(
      any(hit_any_sp, na.rm = TRUE) & any(hit_any_tr, na.rm = TRUE) ~ "candidate_shared_with_sp_and_tr",
      any(hit_any_sp, na.rm = TRUE) ~ "candidate_shared_with_sp",
      any(hit_any_tr, na.rm = TRUE) ~ "candidate_shared_with_tr",
      TRUE ~ "candidate_only_in_search_database"
    ),
    .groups = "drop"
  )

representative_psm <- mq_candidate_main %>%
  dplyr::arrange(PEP, desc(Score), desc(Delta_score), desc(PIF), abs(mass_error_ppm)) %>%
  dplyr::group_by(dataset, comparison, source_id, peptide) %>%
  slice_head(n = 1L) %>%
  dplyr::ungroup() %>%
  transmute(
    dataset, comparison, source_id, peptide,
    representative_raw_file = raw_file,
    representative_scan_number = scan_number,
    representative_msms_id = msms_id,
    representative_charge = charge,
    representative_modified_sequence = modified_sequence,
    representative_PEP = PEP,
    representative_Score = Score,
    representative_Delta_score = Delta_score,
    representative_mass_error_ppm = mass_error_ppm,
    representative_PIF = PIF,
    representative_number_of_matches = number_of_matches,
    representative_intensity_coverage = intensity_coverage,
    representative_peak_coverage = peak_coverage,
    representative_matched_ions = matched_ions
  )

ms_source_summary <- ms_source_summary %>%
  dplyr::left_join(
    representative_psm,
    by = c("dataset", "comparison", "source_id", "peptide")
  )

write_tsv(ms_source_summary, "04_MS_source_peptide_summary.tsv")
write_tsv(representative_psm, "04b_representative_PSM.tsv")

# TMT reporter table. Values are exported, not normalized or statistically compared here.
reporter_cols_main <- grep(
  "^Reporter\\.intensity\\.corrected\\.[0-9]+$",
  names(mq_candidate_main),
  value = TRUE
)

if (length(reporter_cols_main) > 0L) {
  tmt_reporter_long <- mq_candidate_main %>%
    dplyr::filter(dataset_type == "TMT10") %>%
    dplyr::select(
      dataset, comparison, source_id, circRNA_id, orf_id, peptide,
      psm_key, raw_file, scan_number, evidence_id,
      all_of(reporter_cols_main)
    ) %>%
    pivot_longer(
      cols = all_of(reporter_cols_main),
      names_to = "reporter_column",
      values_to = "reporter_intensity_corrected"
    ) %>%
    dplyr::mutate(
      reporter_channel = suppressWarnings(as.integer(str_extract(reporter_column, "[0-9]+$")))
    )

  write_tsv(tmt_reporter_long, "04c_TMT_reporter_intensity_long.tsv")
}

# Label-free evidence intensity by raw file. Sample/group mapping can be added later.
lfq_evidence <- mq_candidate_main %>%
  dplyr::filter(dataset_type == "LabelFree") %>%
  dplyr::distinct(
    dataset, comparison, source_id, circRNA_id, orf_id, peptide,
    evidence_id, raw_file, evidence_intensity, .keep_all = FALSE
  ) %>%
  dplyr::group_by(dataset, comparison, source_id, circRNA_id, orf_id, peptide, raw_file) %>%
  dplyr::summarise(
    n_evidence_features = n_distinct(evidence_id[!is.na(evidence_id)]),
    summed_evidence_intensity = if (all(is.na(evidence_intensity))) NA_real_ else sum(evidence_intensity, na.rm = TRUE),
    max_evidence_intensity = safe_max(evidence_intensity),
    .groups = "drop"
  )

write_tsv(lfq_evidence, "04d_LabelFree_evidence_intensity_by_raw.tsv")

# ============================================================
# 5. Read and normalize NetMHCpan
# ============================================================
read_netmhcpan <- function(path) {
  if (!file.exists(path)) {
    stop(
      "NetMHCpan file does not exist: ", path,
      "\nEdit CFG$netmhcpan_file at the top of the script."
    )
  }

  x <- data.table::fread(
    path,
    data.table = FALSE,
    check.names = TRUE,
    na.strings = c("", "NA", "NaN")
  )

  require_columns(x, c("Peptide", "ID", "X_EL_Rank"), "NetMHCpan table")

  hla_col <- if ("type" %in% names(x)) {
    "type"
  } else if ("HLA" %in% names(x)) {
    "HLA"
  } else if ("Allele" %in% names(x)) {
    "Allele"
  } else {
    stop("NetMHCpan table is missing the HLA allele column (expected type/HLA/Allele).")
  }

  pos_col <- if ("Pos" %in% names(x)) "Pos" else if ("Position" %in% names(x)) "Position" else NULL

  raw_source_id <- str_trim(as.character(x$ID))
  normalized_source_id <- normalize_source_id(raw_source_id)

  out <- tibble(
    source_id_raw = raw_source_id,
    source_id = normalized_source_id,
    circRNA_id = str_extract(normalized_source_id, "hsa_circ_[0-9]+"),
    orf_id = str_extract(normalized_source_id, "^ORF[0-9]+"),
    peptide = normalize_peptide(x$Peptide),
    HLA_allele = str_trim(as.character(x[[hla_col]])),
    position = if (is.null(pos_col)) NA_integer_ else suppressWarnings(as.integer(x[[pos_col]])),
    core = as.character(get_col(x, "X_core")),
    icore = as.character(get_col(x, "X_icore")),
    EL_score = suppressWarnings(as.numeric(get_col(x, "X_EL.score"))),
    EL_rank = suppressWarnings(as.numeric(get_col(x, "X_EL_Rank"))),
    BA_score = suppressWarnings(as.numeric(get_col(x, "X_BA.score"))),
    BA_rank = suppressWarnings(as.numeric(get_col(x, "X_BA_Rank")))
  ) %>%
    dplyr::filter(
      !is.na(source_id), source_id != "",
      !is.na(peptide), peptide != "",
      !is.na(HLA_allele), HLA_allele != "",
      !is.na(EL_rank)
    ) %>%
    dplyr::mutate(
      binder_class = case_when(
        EL_rank < CFG$mhc_strong_rank ~ "SB",
        EL_rank < CFG$mhc_weak_rank ~ "WB",
        TRUE ~ "NB"
      )
    )

  if (CFG$keep_only_mhc_binders) {
    out <- out %>% dplyr::filter(binder_class %in% c("SB", "WB"))
  }

  out
}

mhc_detail <- read_netmhcpan(CFG$netmhcpan_file)
write_tsv(mhc_detail, "05_NetMHCpan_binder_detail.tsv")

mhc_by_allele <- mhc_detail %>%
  dplyr::group_by(source_id, circRNA_id, orf_id, peptide, HLA_allele) %>%
  dplyr::summarise(
    n_predicted_positions = n_distinct(position[!is.na(position)]),
    all_positions = collapse_sorted_numeric(position),
    best_EL_rank = safe_min(EL_rank),
    best_EL_score = safe_max(EL_score),
    best_BA_rank = safe_min(BA_rank),
    best_BA_score = safe_max(BA_score),
    cores = collapse_unique(core),
    icores = collapse_unique(icore),
    binder_class = case_when(
      best_EL_rank < CFG$mhc_strong_rank ~ "SB",
      best_EL_rank < CFG$mhc_weak_rank ~ "WB",
      TRUE ~ "NB"
    ),
    .groups = "drop"
  )

mhc_source_peptide <- mhc_by_allele %>%
  dplyr::group_by(source_id, circRNA_id, orf_id, peptide) %>%
  dplyr::summarise(
    n_HLA_alleles = n_distinct(HLA_allele),
    HLA_alleles = collapse_unique(HLA_allele),
    n_total_predicted_positions = sum(n_predicted_positions, na.rm = TRUE),
    best_EL_rank = safe_min(best_EL_rank),
    best_EL_score = safe_max(best_EL_score),
    best_BA_rank = safe_min(best_BA_rank),
    best_BA_score = safe_max(best_BA_score),
    best_binder_class = case_when(
      best_EL_rank < CFG$mhc_strong_rank ~ "SB",
      best_EL_rank < CFG$mhc_weak_rank ~ "WB",
      TRUE ~ "NB"
    ),
    .groups = "drop"
  )

write_tsv(mhc_by_allele, "05b_NetMHCpan_deduplicated_by_allele.tsv")
write_tsv(mhc_source_peptide, "05c_NetMHCpan_source_peptide_summary.tsv")

# ============================================================
# 6. Exact MS + MHC integration by source_id + peptide
# ============================================================
ms_mhc_by_allele <- ms_source_summary %>%
  dplyr::inner_join(
    mhc_by_allele,
    by = c("source_id", "peptide"),
    suffix = c("_MS", "_MHC")
  ) %>%
  dplyr::mutate(
    source_annotation_consistent =
      (is.na(circRNA_id_MHC) | is.na(circRNA_id_MS) | circRNA_id_MHC == circRNA_id_MS) &
      (is.na(orf_id_MHC) | is.na(orf_id_MS) | orf_id_MHC == orf_id_MS),
    circRNA_id = coalesce(circRNA_id_MS, circRNA_id_MHC),
    orf_id = coalesce(orf_id_MS, orf_id_MHC)
  )

ms_mhc_candidates <- ms_source_summary %>%
  dplyr::inner_join(
    mhc_source_peptide,
    by = c("source_id", "peptide"),
    suffix = c("_MS", "_MHC")
  ) %>%
  dplyr::mutate(
    source_annotation_consistent =
      (is.na(circRNA_id_MHC) | is.na(circRNA_id_MS) | circRNA_id_MHC == circRNA_id_MS) &
      (is.na(orf_id_MHC) | is.na(orf_id_MS) | orf_id_MHC == orf_id_MS),
    circRNA_id = coalesce(circRNA_id_MS, circRNA_id_MHC),
    orf_id = coalesce(orf_id_MS, orf_id_MHC),
    observed_in_TMT = dataset_type == "TMT10",
    observed_in_LabelFree = dataset_type == "LabelFree",
    source_ambiguity = case_when(
      hit_any_sp & leading_is_sp_any ~ "high: shared_with_sp_and_sp_is_leading",
      hit_any_sp ~ "high: shared_with_sp",
      hit_any_tr ~ "moderate: shared_with_tr",
      !candidate_is_leading_any & leading_is_candidate_any ~ "moderate: another_candidate_is_leading",
      candidate_is_leading_any ~ "lower: candidate_is_leading",
      TRUE ~ "undetermined"
    )
  )

write_tsv(ms_mhc_by_allele, "06_MS_MHC_overlap_by_HLA_allele.tsv")
write_tsv(ms_mhc_candidates, "07_MS_MHC_candidate_summary.tsv")

# Collapse across datasets/comparisons to create the next-step BSJ input.
bsj_input <- ms_mhc_candidates %>%
  dplyr::group_by(source_id, circRNA_id, orf_id, peptide) %>%
  dplyr::summarise(
    datasets = collapse_unique(dataset),
    dataset_types = collapse_unique(dataset_type),
    comparisons = collapse_unique(comparison),
    n_dataset_comparison_hits = n_distinct(paste(dataset, comparison)),
    total_PSM = sum(n_PSM, na.rm = TRUE),
    total_raw_files = sum(n_raw_files, na.rm = TRUE),
    best_MS_PEP = safe_min(best_PEP),
    best_MS_Score = safe_max(best_Score),
    best_MS_Delta_score = safe_max(best_Delta_score),
    hit_any_sp = any(hit_any_sp, na.rm = TRUE),
    hit_any_tr = any(hit_any_tr, na.rm = TRUE),
    candidate_is_leading_any = any(candidate_is_leading_any, na.rm = TRUE),
    leading_is_sp_any = any(leading_is_sp_any, na.rm = TRUE),
    candidate_accessions = collapse_unique(candidate_accessions),
    orf_strand = collapse_unique(orf_strand),
    orf_start_0based = safe_min(orf_start_0based),
    orf_end_0based = safe_max(orf_end_0based),
    n_HLA_alleles = safe_max(n_HLA_alleles),
    HLA_alleles = collapse_unique(HLA_alleles),
    best_EL_rank = safe_min(best_EL_rank),
    best_EL_score = safe_max(best_EL_score),
    best_BA_rank = safe_min(best_BA_rank),
    best_BA_score = safe_max(best_BA_score),
    best_binder_class = case_when(
      best_EL_rank < CFG$mhc_strong_rank ~ "SB",
      best_EL_rank < CFG$mhc_weak_rank ~ "WB",
      TRUE ~ "NB"
    ),
    source_ambiguity = case_when(
      any(hit_any_sp, na.rm = TRUE) & any(leading_is_sp_any, na.rm = TRUE) ~ "high: shared_with_sp_and_sp_is_leading",
      any(hit_any_sp, na.rm = TRUE) ~ "high: shared_with_sp",
      any(hit_any_tr, na.rm = TRUE) ~ "moderate: shared_with_tr",
      any(candidate_is_leading_any, na.rm = TRUE) ~ "lower: candidate_is_leading",
      TRUE ~ "undetermined"
    ),
    .groups = "drop"
  ) %>%
  dplyr::arrange(best_EL_rank, best_MS_PEP, desc(best_MS_Score))

write_tsv(bsj_input, "08_BSJ_input_unique_source_peptide.tsv")

# ============================================================
# 7. Diagnostics
# ============================================================
ms_key <- ms_source_summary %>% dplyr::distinct(source_id, circRNA_id, orf_id, peptide)
mhc_key <- mhc_source_peptide %>% dplyr::distinct(source_id, circRNA_id, orf_id, peptide)
exact_key_overlap <- dplyr::inner_join(ms_key, mhc_key, by = c("source_id", "peptide"))
circ_peptide_overlap <- dplyr::inner_join(
  ms_key %>% dplyr::distinct(circRNA_id, peptide),
  mhc_key %>% dplyr::distinct(circRNA_id, peptide),
  by = c("circRNA_id", "peptide")
)

ms_peptides <- unique(ms_key$peptide)
mhc_peptides <- unique(mhc_key$peptide)
ms_sources <- unique(ms_key$source_id)
mhc_sources <- unique(mhc_key$source_id)
ms_circs <- unique(stats::na.omit(ms_key$circRNA_id))
mhc_circs <- unique(stats::na.omit(mhc_key$circRNA_id))

id_diagnostics <- tibble(
  metric = c(
    "MS unique source_id + peptide",
    "MHC unique source_id + peptide",
    "Exact source_id + peptide overlap",
    "MS unique peptides",
    "MHC unique peptides",
    "Peptide-only overlap",
    "MS unique source IDs",
    "MHC unique source IDs",
    "Source-ID-only overlap",
    "MS unique circRNA IDs",
    "MHC unique circRNA IDs",
    "circRNA-ID-only overlap",
    "circRNA ID + peptide overlap",
    "Final unique BSJ-input candidates"
  ),
  value = c(
    nrow(ms_key),
    nrow(mhc_key),
    nrow(exact_key_overlap),
    length(ms_peptides),
    length(mhc_peptides),
    length(intersect(ms_peptides, mhc_peptides)),
    length(ms_sources),
    length(mhc_sources),
    length(intersect(ms_sources, mhc_sources)),
    length(ms_circs),
    length(mhc_circs),
    length(intersect(ms_circs, mhc_circs)),
    nrow(circ_peptide_overlap),
    nrow(bsj_input)
  )
)

write_tsv(id_diagnostics, "09_ID_and_overlap_diagnostics.tsv")
#
# Export peptide-only overlap and assess circRNA consistency
overlap_185 <- ms_source_summary %>%
  dplyr::select(
    peptide,
    MS_source_id = source_id,
    MS_circRNA_id = circRNA_id
  ) %>%
  distinct() %>%
  inner_join(
    mhc_source_peptide %>%
      dplyr::select(
        peptide,
        MHC_source_id = source_id,
        MHC_circRNA_id = circRNA_id,
        HLA_alleles,
        best_EL_rank,
        best_binder_class
      ) %>%
      distinct(),
    by = "peptide",
    relationship = "many-to-many"
  ) %>%
  mutate(
    same_circRNA = MS_circRNA_id == MHC_circRNA_id,
    same_source_id = MS_source_id == MHC_source_id
  )

write_tsv(
  overlap_185,
  "09d_peptide_overlap_mapping.tsv"
)


#

write_tsv(
  mhc_detail %>% dplyr::distinct(source_id_raw, source_id, circRNA_id, orf_id),
  "09a_NetMHCpan_source_ID_normalization.tsv"
)

write_tsv(
  dplyr::inner_join(
    ms_key %>% dplyr::distinct(circRNA_id),
    mhc_key %>% dplyr::distinct(circRNA_id),
    by = "circRNA_id"
  ),
  "09a2_shared_circRNA_IDs.tsv"
)

write_tsv(
  dplyr::anti_join(ms_key, mhc_key, by = c("source_id", "peptide")),
  "09b_MS_source_peptide_without_MHC_match.tsv"
)

write_tsv(
  dplyr::anti_join(mhc_key, ms_key, by = c("source_id", "peptide")),
  "09c_MHC_source_peptide_without_MS_match.tsv"
)

comparison_summary <- ms_source_summary %>%
  dplyr::group_by(dataset, dataset_type, comparison) %>%
  dplyr::summarise(
    n_MS_source_peptide = n(),
    n_MS_unique_peptides = n_distinct(peptide),
    n_MS_unique_sources = n_distinct(source_id),
    .groups = "drop"
  ) %>%
  dplyr::left_join(
    ms_mhc_candidates %>%
      dplyr::group_by(dataset, dataset_type, comparison) %>%
      dplyr::summarise(
        n_MS_MHC_source_peptide = n(),
        n_MS_MHC_unique_peptides = n_distinct(peptide),
        n_MS_MHC_unique_sources = n_distinct(source_id),
        .groups = "drop"
      ),
    by = c("dataset", "dataset_type", "comparison")
  ) %>%
  dplyr::mutate(across(starts_with("n_MS_MHC"), ~tidyr::replace_na(.x, 0L)))

write_tsv(comparison_summary, "10_dataset_comparison_summary.tsv")

# Save configuration and session information for reproducibility.
config_lines <- c(
  paste0("run_time=", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
  paste0("mq_result_root=", CFG$mq_result_root),
  paste0("search_group=", CFG$search_group),
  paste0("netmhcpan_file=", CFG$netmhcpan_file),
  paste0("pep_threshold=", CFG$pep_threshold),
  paste0("score_threshold=", CFG$score_threshold),
  paste0("delta_score_threshold=", CFG$delta_score_threshold),
  paste0("delta_score_policy=", CFG$delta_score_policy),
  paste0("main_orf_strands=", paste(CFG$main_orf_strands, collapse = ",")),
  paste0("mhc_strong_rank=", CFG$mhc_strong_rank),
  paste0("mhc_weak_rank=", CFG$mhc_weak_rank),
  capture.output(sessionInfo())
)
writeLines(config_lines, file.path(CFG$out_dir, "RUN_INFO.txt"))

message("\nPipeline completed successfully.")
message("Next-step BSJ input: ", file.path(CFG$out_dir, "08_BSJ_input_unique_source_peptide.tsv"))


# ============================================================
# Generate a one-row-per-peptide integrated table for the 185 overlapping peptides
# Run directly in R / RStudio
# ============================================================

# 1. Set working directory
setwd(
  "./MS_unspecific/circRNA/result/circRNA_MS_MHC_integration"
)

# 2. Load packages
suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
})

# ============================================================
# Helper functions
# ============================================================

safe_min <- function(x) {
  x <- suppressWarnings(as.numeric(x))
  
  if (length(x) == 0 || all(is.na(x))) {
    return(NA_real_)
  }
  
  min(x, na.rm = TRUE)
}

safe_max <- function(x) {
  x <- suppressWarnings(as.numeric(x))
  
  if (length(x) == 0 || all(is.na(x))) {
    return(NA_real_)
  }
  
  max(x, na.rm = TRUE)
}

collapse_unique <- function(x) {
  x <- as.character(x)
  x <- x[!is.na(x) & x != ""]
  x <- sort(unique(x))
  
  if (length(x) == 0) {
    return(NA_character_)
  }
  
  paste(x, collapse = ";")
}

any_true <- function(x) {
  x <- toupper(as.character(x))
  any(x == "TRUE", na.rm = TRUE)
}

# ============================================================
# Check input files
# ============================================================

required_files <- c(
  "09d_peptide_overlap_mapping.tsv",
  "09g_185_peptide_circRNA_summary.tsv",
  "03_candidate_PSM_source_main_strand.tsv",
  "05b_NetMHCpan_deduplicated_by_allele.tsv"
)

missing_files <- required_files[!file.exists(required_files)]

if (length(missing_files) > 0) {
  stop(
    "The following input files are missing:\n",
    paste(missing_files, collapse = "\n")
  )
}

# ============================================================
# 3. Extract valid peptide-circRNA pairs
# ============================================================

overlap_mapping <- fread(
  "09d_peptide_overlap_mapping.tsv",
  data.table = FALSE
)

# Coerce to character to avoid logical/character join conflicts
overlap_mapping <- overlap_mapping %>%
  mutate(
    peptide = as.character(peptide),
    MS_circRNA_id = as.character(MS_circRNA_id),
    same_circRNA = as.character(same_circRNA)
  )

valid_pairs <- overlap_mapping %>%
  filter(toupper(same_circRNA) == "TRUE") %>%
  transmute(
    peptide = peptide,
    circRNA_id = MS_circRNA_id
  ) %>%
  filter(
    !is.na(peptide),
    !is.na(circRNA_id),
    peptide != "",
    circRNA_id != ""
  ) %>%
  distinct()


# ============================================================
# 4. Read the one-row-per-peptide circRNA mapping table
# ============================================================

circ_summary <- fread(
  "09g_185_peptide_circRNA_summary.tsv",
  data.table = FALSE
) %>%
  mutate(
    peptide = as.character(peptide),
    circRNA_ids = as.character(circRNA_ids),
    circRNA_mapping_class =
      as.character(circRNA_mapping_class)
  )


# ============================================================
# 5. Summarize MaxQuant PSM evidence
# ============================================================

ms_raw <- fread(
  "03_candidate_PSM_source_main_strand.tsv",
  data.table = FALSE
) %>%
  mutate(
    peptide = as.character(peptide),
    circRNA_id = as.character(circRNA_id),
    source_id = as.character(source_id)
  )


ms_detail <- ms_raw %>%
  inner_join(
    valid_pairs,
    by = c("peptide", "circRNA_id")
  )


ms_summary <- ms_detail %>%
  group_by(peptide) %>%
  summarise(
    datasets = collapse_unique(dataset),
    comparisons = collapse_unique(comparison),
    
    MS_source_ids = collapse_unique(source_id),
    
    n_MS_source_ids = n_distinct(
      source_id[
        !is.na(source_id) &
          source_id != ""
      ]
    ),
    
    n_PSM = n_distinct(psm_key),
    
    n_raw_files = n_distinct(
      paste(
        dataset,
        comparison,
        raw_file,
        sep = "|"
      )[
        !is.na(raw_file) &
          raw_file != ""
      ]
    ),
    
    n_charge_states = n_distinct(
      charge[!is.na(charge)]
    ),
    
    best_PEP = safe_min(PEP),
    best_Score = safe_max(Score),
    best_Delta_score = safe_max(Delta_score),
    
    best_abs_mass_error_ppm =
      safe_min(abs(mass_error_ppm)),
    
    best_PIF = safe_max(PIF),
    
    hit_any_sp = any_true(hit_any_sp),
    hit_any_tr = any_true(hit_any_tr),
    
    candidate_is_leading_any =
      any_true(candidate_is_leading),
    
    leading_is_sp_any =
      any_true(leading_is_sp),
    
    leading_razor_proteins =
      collapse_unique(leading_razor_protein),
    
    .groups = "drop"
  )

# ============================================================
# 6. Summarize NetMHCpan evidence
# ============================================================

mhc_raw <- fread(
  "05b_NetMHCpan_deduplicated_by_allele.tsv",
  data.table = FALSE
) %>%
  mutate(
    peptide = as.character(peptide),
    circRNA_id = as.character(circRNA_id),
    source_id = as.character(source_id),
    HLA_allele = as.character(HLA_allele)
  )


mhc_detail <- mhc_raw %>%
  inner_join(
    valid_pairs,
    by = c("peptide", "circRNA_id")
  )


mhc_summary <- mhc_detail %>%
  group_by(peptide) %>%
  summarise(
    MHC_source_ids = collapse_unique(source_id),
    
    n_MHC_source_ids = n_distinct(
      source_id[
        !is.na(source_id) &
          source_id != ""
      ]
    ),
    
    HLA_alleles = collapse_unique(HLA_allele),
    
    n_HLA_alleles = n_distinct(
      HLA_allele[
        !is.na(HLA_allele) &
          HLA_allele != ""
      ]
    ),
    
    best_EL_rank = safe_min(best_EL_rank),
    best_EL_score = safe_max(best_EL_score),
    best_BA_rank = safe_min(best_BA_rank),
    
    .groups = "drop"
  ) %>%
  mutate(
    best_binder_class = case_when(
      is.na(best_EL_rank) ~ NA_character_,
      best_EL_rank < 0.5 ~ "SB",
      best_EL_rank < 2 ~ "WB",
      TRUE ~ "NB"
    )
  )

# ============================================================
# 7. Merge into a one-row-per-peptide integrated table
# ============================================================

integrated <- circ_summary %>%
  left_join(
    ms_summary,
    by = "peptide"
  ) %>%
  left_join(
    mhc_summary,
    by = "peptide"
  ) %>%
  mutate(
    MQ_source_class = case_when(
      hit_any_sp %in% TRUE &
        leading_is_sp_any %in% TRUE ~
        "shared_with_SwissProt_and_SP_is_leading",
      
      hit_any_sp %in% TRUE ~
        "shared_with_SwissProt",
      
      hit_any_tr %in% TRUE ~
        "shared_with_TrEMBL",
      
      hit_any_sp %in% FALSE &
        hit_any_tr %in% FALSE ~
        "candidate_only_in_current_search_database",
      
      TRUE ~
        "undetermined"
    ),
    
    current_candidate_class = case_when(
      circRNA_mapping_class == "single_circRNA" &
        MQ_source_class ==
        "candidate_only_in_current_search_database" ~
        "higher_priority_for_follow_up",
      
      circRNA_mapping_class == "multiple_circRNAs" ~
        "multiple_circRNA_mapping",
      
      grepl(
        "SwissProt",
        MQ_source_class,
        fixed = TRUE
      ) ~
        "shared_with_known_protein",
      
      MQ_source_class == "shared_with_TrEMBL" ~
        "shared_with_TrEMBL",
      
      TRUE ~
        "requires_review"
    )
  ) %>%
  arrange(
    factor(
      current_candidate_class,
      levels = c(
        "higher_priority_for_follow_up",
        "requires_review",
        "multiple_circRNA_mapping",
        "shared_with_TrEMBL",
        "shared_with_known_protein"
      )
    ),
    best_EL_rank,
    best_PEP,
    desc(best_Score)
  )

# ============================================================
# 8. Export results
# ============================================================

output_file <- "10_185_integrated_peptide_summary.tsv"

fwrite(
  integrated,
  file = output_file,
  sep = "\t",
  quote = FALSE,
  na = "NA"
)


# Key check
if (nrow(integrated) != 185) {
  warning(
    "Output does not contain 185 rows; check for duplicated or missing peptides."
  )
}


# ============================================================
# Recover NetMHCpan circRNA_id values and update the 185-peptide integrated table
# Run directly in R / RStudio
# ============================================================

setwd(
  "./MS_unspecific/circRNA/result/circRNA_MS_MHC_integration"
)

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
})

# ------------------------------------------------------------
# Helper functions
# ------------------------------------------------------------

safe_min <- function(x) {
  x <- suppressWarnings(as.numeric(x))
  
  if (length(x) == 0 || all(is.na(x))) {
    return(NA_real_)
  }
  
  min(x, na.rm = TRUE)
}

safe_max <- function(x) {
  x <- suppressWarnings(as.numeric(x))
  
  if (length(x) == 0 || all(is.na(x))) {
    return(NA_real_)
  }
  
  max(x, na.rm = TRUE)
}

collapse_unique <- function(x) {
  x <- as.character(x)
  x <- x[!is.na(x) & x != ""]
  x <- sort(unique(x))
  
  if (length(x) == 0) {
    return(NA_character_)
  }
  
  paste(x, collapse = ";")
}

# ============================================================
# 1. Read valid peptide-circRNA pairs
# ============================================================

overlap_mapping <- fread(
  "09d_peptide_overlap_mapping.tsv",
  data.table = FALSE
) %>%
  mutate(
    peptide = as.character(peptide),
    MS_circRNA_id = as.character(MS_circRNA_id),
    same_circRNA = as.character(same_circRNA)
  )

valid_pairs <- overlap_mapping %>%
  filter(toupper(same_circRNA) == "TRUE") %>%
  transmute(
    peptide = peptide,
    circRNA_id = MS_circRNA_id
  ) %>%
  filter(
    !is.na(peptide),
    !is.na(circRNA_id),
    peptide != "",
    circRNA_id != ""
  ) %>%
  distinct()


# ============================================================
# 2. Read only the required columns from the NetMHCpan table
# ============================================================

mhc_raw <- fread(
  "05b_NetMHCpan_deduplicated_by_allele.tsv",
  select = c(
    "source_id",
    "peptide",
    "HLA_allele",
    "best_EL_rank",
    "best_EL_score",
    "best_BA_rank",
    "best_BA_score",
    "binder_class"
  ),
  data.table = FALSE
)

mhc_raw <- mhc_raw %>%
  mutate(
    source_id = as.character(source_id),
    peptide = as.character(peptide),
    HLA_allele = as.character(HLA_allele)
  )

# ============================================================
# 3. Recover circRNA_id from source_id
# ============================================================

has_short_id <- grepl(
  "^ORF[0-9]+_[0-9]+$",
  mhc_raw$source_id
)

has_full_id <- grepl(
  "^ORF[0-9]+_hsa_circ_[0-9]+",
  mhc_raw$source_id
)

mhc_raw$circRNA_id <- NA_character_
mhc_raw$source_id_full <- NA_character_

# Current format, e.g. ORF101_0126885
mhc_raw$circRNA_id[has_short_id] <- paste0(
  "hsa_circ_",
  sub(
    "^ORF[0-9]+_",
    "",
    mhc_raw$source_id[has_short_id]
  )
)

mhc_raw$source_id_full[has_short_id] <- sub(
  "^(ORF[0-9]+)_([0-9]+)$",
  "\\1_hsa_circ_\\2",
  mhc_raw$source_id[has_short_id]
)

# Also support source_id values already in full format
mhc_raw$circRNA_id[has_full_id] <- sub(
  "^ORF[0-9]+_(hsa_circ_[0-9]+).*$",
  "\\1",
  mhc_raw$source_id[has_full_id]
)

mhc_raw$source_id_full[has_full_id] <-
  mhc_raw$source_id[has_full_id]


# ============================================================
# 4. Retain only MHC results consistent with the sources of the 185 peptides
# ============================================================

mhc_detail <- mhc_raw %>%
  inner_join(
    valid_pairs,
    by = c(
      "peptide",
      "circRNA_id"
    )
  )


# ============================================================
# 5. Summarize MHC results to one row per peptide
# ============================================================

mhc_summary_fixed <- mhc_detail %>%
  group_by(peptide) %>%
  summarise(
    MHC_source_ids = collapse_unique(
      source_id_full
    ),
    
    n_MHC_source_ids = n_distinct(
      source_id_full[
        !is.na(source_id_full) &
          source_id_full != ""
      ]
    ),
    
    HLA_alleles = collapse_unique(
      HLA_allele
    ),
    
    n_HLA_alleles = n_distinct(
      HLA_allele[
        !is.na(HLA_allele) &
          HLA_allele != ""
      ]
    ),
    
    best_EL_rank = safe_min(
      best_EL_rank
    ),
    
    best_EL_score = safe_max(
      best_EL_score
    ),
    
    best_BA_rank = safe_min(
      best_BA_rank
    ),
    
    best_BA_score = safe_max(
      best_BA_score
    ),
    
    predicted_binder_classes =
      collapse_unique(binder_class),
    
    .groups = "drop"
  ) %>%
  mutate(
    best_binder_class = case_when(
      is.na(best_EL_rank) ~
        NA_character_,
      
      best_EL_rank < 0.5 ~
        "SB",
      
      best_EL_rank < 2 ~
        "WB",
      
      TRUE ~
        "NB"
    )
  )


# ============================================================
# 6. Update the previously generated integrated table
# ============================================================

integrated_old <- fread(
  "10_185_integrated_peptide_summary.tsv",
  data.table = FALSE
)

# Remove the old MHC columns that contain only NA values
old_mhc_columns <- c(
  "MHC_source_ids",
  "n_MHC_source_ids",
  "HLA_alleles",
  "n_HLA_alleles",
  "best_EL_rank",
  "best_EL_score",
  "best_BA_rank",
  "best_BA_score",
  "predicted_binder_classes",
  "best_binder_class"
)

integrated_fixed <- integrated_old %>%
  dplyr::select(
    -any_of(old_mhc_columns)
  ) %>%
  left_join(
    mhc_summary_fixed,
    by = "peptide"
  ) %>%
  arrange(
    circRNA_mapping_class,
    best_EL_rank,
    best_PEP,
    desc(best_Score)
  )

# ============================================================
# 7. Save corrected results
# ============================================================

fwrite(
  integrated_fixed,
  "10_185_integrated_peptide_summary_fixed.tsv",
  sep = "\t",
  quote = FALSE,
  na = "NA"
)


# Prepare files required for BSJ-crossing analysis --------------------------------------------------------------------

# ============================================================
# Merge and deduplicate five single-circle circRNA nucleotide FASTA files
# ============================================================

nucleotide_files <- c(
  "./circRNA_peptide3/part1_coding_prob/IRES_finder/merge_fasta/RA_eRA_down_merge.fasta",
  "./circRNA_peptide3/part1_coding_prob/IRES_finder/merge_fasta/H_RA_down_merge.fasta",
  "./circRNA_peptide3/part1_coding_prob/IRES_finder/merge_fasta/H_eRA_down_merge.fasta",
  "./circRNA_peptide3/part1_coding_prob/IRES_finder/merge_fasta/H_eRA_up_merge.fasta",
  "./circRNA_peptide3/part1_coding_prob/IRES_finder/merge_fasta/H_RA_up_merge.fasta"
)

output_file <- paste0(
  "./MS_unspecific/circRNA/result/",
  "circRNA_MS_MHC_integration/",
  "11_merged_single_circle_circRNA.fasta"
)

# Read FASTA and return IDs and sequences
read_fasta_base <- function(file) {
  
  lines <- readLines(file, warn = FALSE)
  header_index <- which(grepl("^>", lines))
  
  if (length(header_index) == 0) {
    stop("No FASTA header found in file: ", file)
  }
  
  end_index <- c(
    header_index[-1] - 1,
    length(lines)
  )
  
  ids <- sub(
    "^>(\\S+).*$",
    "\\1",
    lines[header_index]
  )
  
  sequences <- vapply(
    seq_along(header_index),
    function(i) {
      sequence_lines <- lines[
        (header_index[i] + 1):end_index[i]
      ]
      
      sequence_lines <- sequence_lines[
        !grepl("^>", sequence_lines)
      ]
      
      toupper(
        gsub(
          "\\s+",
          "",
          paste(sequence_lines, collapse = "")
        )
      )
    },
    character(1)
  )
  
  data.frame(
    circRNA_id = ids,
    sequence = sequences,
    source_file = basename(file),
    stringsAsFactors = FALSE
  )
}

# Read all files
circ_list <- lapply(
  nucleotide_files,
  read_fasta_base
)

circ_all <- do.call(
  rbind,
  circ_list
)


# Check whether the same ID is associated with different sequences
conflict_check <- aggregate(
  sequence ~ circRNA_id,
  data = circ_all,
  FUN = function(x) length(unique(x))
)

conflicting_ids <- conflict_check[
  conflict_check$sequence > 1,
  "circRNA_id"
]


if (length(conflicting_ids) > 0) {
  
  conflict_detail <- circ_all[
    circ_all$circRNA_id %in% conflicting_ids,
  ]
  
  write.table(
    conflict_detail,
    file = paste0(
      "./MS_unspecific/circRNA/result/",
      "circRNA_MS_MHC_integration/",
      "11a_conflicting_circRNA_sequences.tsv"
    ),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )
  
  stop(
    "The same circRNA ID is associated with different sequences; ",
    "11a_conflicting_circRNA_sequences.tsv has been written; ",
    "the merged FASTA was not generated."
  )
}

# Keep only one record when the same ID has an identical sequence
circ_unique <- circ_all[
  !duplicated(circ_all$circRNA_id),
  c("circRNA_id", "sequence")
]

# Check for empty sequences
if (any(nchar(circ_unique$sequence) == 0)) {
  stop("Empty circRNA sequences are present.")
}

# Write the merged FASTA with 80 nt per line
connection <- file(
  output_file,
  open = "w"
)

for (i in seq_len(nrow(circ_unique))) {
  
  writeLines(
    paste0(">", circ_unique$circRNA_id[i]),
    connection
  )
  
  sequence_i <- circ_unique$sequence[i]
  
  starts <- seq(
    1,
    nchar(sequence_i),
    by = 80
  )
  
  sequence_lines <- substring(
    sequence_i,
    starts,
    pmin(
      starts + 79,
      nchar(sequence_i)
    )
  )
  
  writeLines(
    sequence_lines,
    connection
  )
}

close(connection)


# ============================================================
# Merge three MaxQuant circRNA ORF protein FASTA files
# ============================================================

orf_files <- c(
  "./circRNA_peptide4/circRNA/MaxQuant/need_file/only_circRNA/merge_circ_RA_eRA.fasta",
  "./circRNA_peptide4/circRNA/MaxQuant/need_file/only_circRNA/merge_circ_HC_eRA.fasta",
  "./circRNA_peptide4/circRNA/MaxQuant/need_file/only_circRNA/merge_circ_HC_RA.fasta"
)

output_dir <- paste0(
  "./MS_unspecific/circRNA/result/",
  "circRNA_MS_MHC_integration"
)

output_fasta <- file.path(
  output_dir,
  "12_merged_MaxQuant_circRNA_ORF.fasta"
)

# ------------------------------------------------------------
# FASTA reader
# ------------------------------------------------------------

read_protein_fasta <- function(file) {
  
  lines <- readLines(file, warn = FALSE)
  
  header_index <- which(grepl("^>", lines))
  
  if (length(header_index) == 0) {
    stop("No FASTA header found in file: ", file)
  }
  
  end_index <- c(
    header_index[-1] - 1,
    length(lines)
  )
  
  headers <- sub(
    "^>",
    "",
    lines[header_index]
  )
  
  # Retain only the first field of each FASTA header
  # Example:
  # lcl|ORF4_hsa_circ_0000798:2490:4208
  protein_ids <- sub(
    "\\s+.*$",
    "",
    headers
  )
  
  sequences <- vapply(
    seq_along(header_index),
    function(i) {
      
      start_line <- header_index[i] + 1
      end_line <- end_index[i]
      
      if (start_line > end_line) {
        return("")
      }
      
      sequence_lines <- lines[start_line:end_line]
      
      sequence_lines <- sequence_lines[
        !grepl("^>", sequence_lines)
      ]
      
      toupper(
        gsub(
          "\\s+",
          "",
          paste(sequence_lines, collapse = "")
        )
      )
    },
    character(1)
  )
  
  data.frame(
    protein_id = protein_ids,
    original_header = headers,
    protein_sequence = sequences,
    source_file = basename(file),
    stringsAsFactors = FALSE
  )
}

# ------------------------------------------------------------
# Read the three files
# ------------------------------------------------------------

orf_list <- lapply(
  orf_files,
  read_protein_fasta
)

orf_all <- do.call(
  rbind,
  orf_list
)


# ------------------------------------------------------------
# Check whether the same protein ID is associated with different protein sequences
# ------------------------------------------------------------

conflict_check <- aggregate(
  protein_sequence ~ protein_id,
  data = orf_all,
  FUN = function(x) length(unique(x))
)

conflicting_ids <- conflict_check$protein_id[
  conflict_check$protein_sequence > 1
]


if (length(conflicting_ids) > 0) {
  
  conflict_detail <- orf_all[
    orf_all$protein_id %in% conflicting_ids,
  ]
  
  write.table(
    conflict_detail,
    file = file.path(
      output_dir,
      "12a_conflicting_ORF_sequences.tsv"
    ),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )
  
  stop(
    "The same protein ID is associated with different protein sequences; ",
    "12a_conflicting_ORF_sequences.tsv has been generated."
  )
}

# ------------------------------------------------------------
# Remove duplicated protein IDs
# ------------------------------------------------------------

orf_unique <- orf_all[
  !duplicated(orf_all$protein_id),
]

if (any(nchar(orf_unique$protein_sequence) == 0)) {
  stop("Empty ORF protein sequences are present.")
}

# ------------------------------------------------------------
# Write the merged FASTA
# ------------------------------------------------------------

connection <- file(
  output_fasta,
  open = "w"
)

for (i in seq_len(nrow(orf_unique))) {
  
  writeLines(
    paste0(">", orf_unique$original_header[i]),
    connection
  )
  
  sequence_i <- orf_unique$protein_sequence[i]
  
  starts <- seq(
    1,
    nchar(sequence_i),
    by = 80
  )
  
  sequence_lines <- substring(
    sequence_i,
    starts,
    pmin(
      starts + 79,
      nchar(sequence_i)
    )
  )
  
  writeLines(
    sequence_lines,
    connection
  )
}

close(connection)

# Also save a TSV manifest for inspection
write.table(
  orf_unique,
  file = file.path(
    output_dir,
    "12_merged_MaxQuant_circRNA_ORF_manifest.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)


# ============================================================
# Check ORF coordinates against single-circle circRNA lengths
# ============================================================

setwd(
  "./MS_unspecific/circRNA/result/circRNA_MS_MHC_integration"
)

library(data.table)
library(dplyr)

# ------------------------------------------------------------
# Read the single-circle circRNA FASTA
# ------------------------------------------------------------

read_fasta_simple <- function(file) {
  
  lines <- readLines(file, warn = FALSE)
  header_index <- which(grepl("^>", lines))
  end_index <- c(header_index[-1] - 1, length(lines))
  
  ids <- sub(
    "^>(\\S+).*$",
    "\\1",
    lines[header_index]
  )
  
  sequences <- vapply(
    seq_along(header_index),
    function(i) {
      
      start_line <- header_index[i] + 1
      end_line <- end_index[i]
      
      if (start_line > end_line) {
        return("")
      }
      
      paste(
        lines[start_line:end_line],
        collapse = ""
      ) |>
        gsub("\\s+", "", x = _) |>
        toupper()
    },
    character(1)
  )
  
  data.frame(
    circRNA_id = ids,
    circRNA_sequence = sequences,
    circRNA_length = nchar(sequences),
    stringsAsFactors = FALSE
  )
}

circ_info <- read_fasta_simple(
  "11_merged_single_circle_circRNA.fasta"
)


# ------------------------------------------------------------
# Read the ORF manifest
# ------------------------------------------------------------

orf_manifest <- fread(
  "12_merged_MaxQuant_circRNA_ORF_manifest.tsv",
  data.table = FALSE
)

# ------------------------------------------------------------
# Parse ORF, circRNA, and coordinate information from the header
# ------------------------------------------------------------

orf_parsed <- orf_manifest %>%
  mutate(
    protein_id = as.character(protein_id),
    
    clean_id = sub(
      "^lcl\\|",
      "",
      protein_id
    ),
    
    orf_id = sub(
      "_hsa_circ_.*$",
      "",
      clean_id
    ),
    
    circRNA_id = sub(
      "^ORF[0-9]+_(hsa_circ_[0-9]+):.*$",
      "\\1",
      clean_id
    ),
    
    orf_start = as.integer(
      sub(
        "^.*:([0-9]+):([0-9]+)$",
        "\\1",
        clean_id
      )
    ),
    
    orf_end = as.integer(
      sub(
        "^.*:([0-9]+):([0-9]+)$",
        "\\2",
        clean_id
      )
    ),
    
    protein_length_aa = nchar(
      protein_sequence
    ),
    
    coordinate_length_nt =
      orf_end - orf_start + 1,
    
    expected_length_nt =
      protein_length_aa * 3,
    
    coordinate_difference_nt =
      coordinate_length_nt -
      expected_length_nt
  ) %>%
  left_join(
    circ_info %>%
      dplyr::select(
        circRNA_id,
        circRNA_length
      ),
    by = "circRNA_id"
  )

# ------------------------------------------------------------
# Determine whether each ORF crosses one or more BSJs
# ------------------------------------------------------------

orf_parsed <- orf_parsed %>%
  mutate(
    start_circle_index = floor(
      (orf_start - 1) / circRNA_length
    ),
    
    end_circle_index = floor(
      (orf_end - 1) / circRNA_length
    ),
    
    orf_crosses_BSJ =
      end_circle_index >
      start_circle_index,
    
    n_BSJ_crossed =
      end_circle_index -
      start_circle_index
  )

# ------------------------------------------------------------
# Save diagnostic results
# ------------------------------------------------------------

fwrite(
  orf_parsed,
  "13_ORF_coordinate_and_BSJ_diagnostic.tsv",
  sep = "\t",
  quote = FALSE,
  na = "NA"
)


# ============================================================
# Peptide-level BSJ analysis for the 185 candidate peptides
# ============================================================

setwd(
  "./MS_unspecific/circRNA/result/circRNA_MS_MHC_integration"
)

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
})

# ------------------------------------------------------------
# 1. Normalize ORF ID format
# ------------------------------------------------------------

normalize_source_id <- function(x) {
  
  x <- as.character(x)
  
  # Remove lcl|
  x <- sub("^lcl\\|", "", x)
  
  # Remove header annotations
  x <- sub("\\s+.*$", "", x)
  
  # Remove ORF coordinates
  x <- sub(
    ":[0-9]+:[0-9]+$",
    "",
    x
  )
  
  # Support the ORF101_0126885 format
  short_format <- grepl(
    "^ORF[0-9]+_[0-9]+$",
    x
  )
  
  x[short_format] <- sub(
    "^(ORF[0-9]+)_([0-9]+)$",
    "\\1_hsa_circ_\\2",
    x[short_format]
  )
  
  x
}

# ------------------------------------------------------------
# 2. Read ORF coordinates and sequences
# ------------------------------------------------------------

orf_info <- fread(
  "13_ORF_coordinate_and_BSJ_diagnostic.tsv",
  data.table = FALSE
) %>%
  mutate(
    protein_id = as.character(protein_id),
    source_id = normalize_source_id(protein_id),
    circRNA_id = as.character(circRNA_id),
    protein_sequence = as.character(protein_sequence),
    orf_start = as.integer(orf_start),
    orf_end = as.integer(orf_end),
    circRNA_length = as.integer(circRNA_length)
  )


# ------------------------------------------------------------
# 3. Read mappings between the 185 peptides and MaxQuant ORFs
# ------------------------------------------------------------

mapping_raw <- fread(
  "09d_peptide_overlap_mapping.tsv",
  data.table = FALSE
)

candidate_mapping <- mapping_raw %>%
  mutate(
    peptide = as.character(peptide),
    MS_source_id = as.character(MS_source_id),
    MS_circRNA_id = as.character(MS_circRNA_id),
    same_circRNA = as.character(same_circRNA)
  ) %>%
  filter(
    toupper(same_circRNA) == "TRUE"
  ) %>%
  transmute(
    peptide,
    source_id = normalize_source_id(
      MS_source_id
    ),
    circRNA_id = MS_circRNA_id
  ) %>%
  filter(
    !is.na(peptide),
    peptide != "",
    !is.na(source_id),
    source_id != ""
  ) %>%
  distinct()


# ------------------------------------------------------------
# 4. Join ORF coordinates, sequences, and circRNA lengths
# ------------------------------------------------------------

candidate_orf <- candidate_mapping %>%
  left_join(
    orf_info %>%
      dplyr::select(
        source_id,
        protein_id,
        circRNA_id,
        protein_sequence,
        orf_start,
        orf_end,
        protein_length_aa,
        circRNA_length
      ),
    by = c(
      "source_id",
      "circRNA_id"
    )
  )

# Records with missing ORF information or single-circle length
unresolved_mapping <- candidate_orf %>%
  filter(
    is.na(protein_sequence) |
      is.na(circRNA_length) |
      is.na(orf_start)
  )

fwrite(
  unresolved_mapping,
  "14a_unresolved_candidate_ORF_mapping.tsv",
  sep = "\t",
  quote = FALSE,
  na = "NA"
)


candidate_orf_resolved <- candidate_orf %>%
  filter(
    !is.na(protein_sequence),
    !is.na(circRNA_length),
    !is.na(orf_start)
  )

# ------------------------------------------------------------
# 5. BSJ analysis function for one peptide-ORF mapping
# ------------------------------------------------------------

analyse_one_mapping <- function(i) {
  
  peptide_i <- candidate_orf_resolved$peptide[i]
  protein_i <- candidate_orf_resolved$protein_sequence[i]
  
  source_i <- candidate_orf_resolved$source_id[i]
  protein_id_i <- candidate_orf_resolved$protein_id[i]
  circRNA_i <- candidate_orf_resolved$circRNA_id[i]
  
  orf_start_i <- candidate_orf_resolved$orf_start[i]
  orf_end_i <- candidate_orf_resolved$orf_end[i]
  circ_length_i <- candidate_orf_resolved$circRNA_length[i]
  
  peptide_length <- nchar(peptide_i)
  
  # All peptide occurrence positions in the protein sequence
  positions <- gregexpr(
    peptide_i,
    protein_i,
    fixed = TRUE
  )[[1]]
  
  # Peptide not found in the corresponding ORF protein
  if (
    length(positions) == 1 &&
    positions[1] == -1
  ) {
    
    return(
      data.frame(
        peptide = peptide_i,
        source_id = source_i,
        protein_id = protein_id_i,
        circRNA_id = circRNA_i,
        occurrence_found = FALSE,
        peptide_occurrence = NA_integer_,
        aa_start = NA_integer_,
        aa_end = NA_integer_,
        peptide_length_aa = peptide_length,
        peptide_nt_start = NA_integer_,
        peptide_nt_end = NA_integer_,
        BSJ_boundary_nt = NA_integer_,
        peptide_crosses_BSJ = NA,
        aa_left_of_BSJ = NA_integer_,
        aa_right_of_BSJ = NA_integer_,
        junction_codon_count = NA_integer_,
        min_complete_aa_each_side = NA_integer_,
        BSJ_overlap_class =
          "peptide_not_found_in_ORF",
        stringsAsFactors = FALSE
      )
    )
  }
  
  occurrence_results <- lapply(
    seq_along(positions),
    function(j) {
      
      aa_start_i <- positions[j]
      aa_end_i <-
        aa_start_i +
        peptide_length -
        1
      
      # Nucleotide coordinates of the peptide on the repeated circRNA sequence
      peptide_nt_start_i <-
        orf_start_i +
        (aa_start_i - 1) * 3
      
      peptide_nt_end_i <-
        peptide_nt_start_i +
        peptide_length * 3 -
        1
      
      # Biological BSJs are located at:
      # circRNA_length、2*circRNA_length……
      first_boundary <- ceiling(
        peptide_nt_start_i /
          circ_length_i
      ) * circ_length_i
      
      last_boundary <- floor(
        (peptide_nt_end_i - 1) /
          circ_length_i
      ) * circ_length_i
      
      if (first_boundary > last_boundary) {
        
        boundaries <- integer(0)
        
      } else {
        
        boundaries <- seq(
          first_boundary,
          last_boundary,
          by = circ_length_i
        )
      }
      
      # Peptide does not cross a BSJ
      if (length(boundaries) == 0) {
        
        return(
          data.frame(
            peptide = peptide_i,
            source_id = source_i,
            protein_id = protein_id_i,
            circRNA_id = circRNA_i,
            occurrence_found = TRUE,
            peptide_occurrence = j,
            aa_start = aa_start_i,
            aa_end = aa_end_i,
            peptide_length_aa = peptide_length,
            peptide_nt_start =
              peptide_nt_start_i,
            peptide_nt_end =
              peptide_nt_end_i,
            BSJ_boundary_nt = NA_integer_,
            peptide_crosses_BSJ = FALSE,
            aa_left_of_BSJ = NA_integer_,
            aa_right_of_BSJ = NA_integer_,
            junction_codon_count = 0L,
            min_complete_aa_each_side =
              NA_integer_,
            BSJ_overlap_class = "non_BSJ",
            stringsAsFactors = FALSE
          )
        )
      }
      
      # A short peptide will usually encounter only one BSJ,
      # but the code supports multiple BSJs
      boundary_results <- lapply(
        boundaries,
        function(boundary_i) {
          
          codon_starts <-
            peptide_nt_start_i +
            (seq_len(peptide_length) - 1) * 3
          
          codon_ends <-
            codon_starts + 2
          
          # Amino acids located completely to the left of the BSJ
          aa_left <- sum(
            codon_ends <= boundary_i
          )
          
          # Amino acids located completely to the right of the BSJ
          aa_right <- sum(
            codon_starts > boundary_i
          )
          
          # Codons that themselves cross the BSJ
          junction_codons <- sum(
            codon_starts <= boundary_i &
              codon_ends > boundary_i
          )
          
          min_each_side <- min(
            aa_left,
            aa_right
          )
          
          overlap_class <- dplyr::case_when(
            min_each_side == 0 ~
              "junction_codon_only",
            
            min_each_side == 1 ~
              "marginal_BSJ_overlap",
            
            min_each_side >= 2 ~
              "robust_BSJ_overlap",
            
            TRUE ~
              "undetermined"
          )
          
          data.frame(
            peptide = peptide_i,
            source_id = source_i,
            protein_id = protein_id_i,
            circRNA_id = circRNA_i,
            occurrence_found = TRUE,
            peptide_occurrence = j,
            aa_start = aa_start_i,
            aa_end = aa_end_i,
            peptide_length_aa = peptide_length,
            peptide_nt_start =
              peptide_nt_start_i,
            peptide_nt_end =
              peptide_nt_end_i,
            BSJ_boundary_nt =
              boundary_i,
            peptide_crosses_BSJ = TRUE,
            aa_left_of_BSJ = aa_left,
            aa_right_of_BSJ = aa_right,
            junction_codon_count =
              junction_codons,
            min_complete_aa_each_side =
              min_each_side,
            BSJ_overlap_class =
              overlap_class,
            stringsAsFactors = FALSE
          )
        }
      )
      
      do.call(
        rbind,
        boundary_results
      )
    }
  )
  
  do.call(
    rbind,
    occurrence_results
  )
}

# ------------------------------------------------------------
# 6. Analyze all candidates
# ------------------------------------------------------------

bsj_detail_list <- lapply(
  seq_len(nrow(candidate_orf_resolved)),
  analyse_one_mapping
)

bsj_detail <- bind_rows(
  bsj_detail_list
)

fwrite(
  bsj_detail,
  "14_peptide_BSJ_mapping_detail.tsv",
  sep = "\t",
  quote = FALSE,
  na = "NA"
)


# ------------------------------------------------------------
# 7. Summarize to one row per peptide
# ------------------------------------------------------------

bsj_found <- bsj_detail %>%
  filter(
    occurrence_found
  ) %>%
  mutate(
    occurrence_key = paste(
      source_id,
      aa_start,
      aa_end,
      sep = "|"
    )
  )

bsj_summary_observed <- bsj_found %>%
  group_by(peptide) %>%
  summarise(
    n_ORF_mappings =
      n_distinct(source_id),
    
    n_peptide_occurrences =
      n_distinct(occurrence_key),
    
    n_BSJ_crossing_occurrences =
      n_distinct(
        occurrence_key[
          peptide_crosses_BSJ %in% TRUE
        ]
      ),
    
    n_non_BSJ_occurrences =
      n_distinct(
        occurrence_key[
          peptide_crosses_BSJ %in% FALSE
        ]
      ),
    
    any_junction_codon_only =
      any(
        BSJ_overlap_class ==
          "junction_codon_only"
      ),
    
    any_marginal_BSJ =
      any(
        BSJ_overlap_class ==
          "marginal_BSJ_overlap"
      ),
    
    any_robust_BSJ =
      any(
        BSJ_overlap_class ==
          "robust_BSJ_overlap"
      ),
    
    max_complete_aa_each_side =
      if (
        all(
          is.na(
            min_complete_aa_each_side
          )
        )
      ) {
        NA_integer_
      } else {
        max(
          min_complete_aa_each_side,
          na.rm = TRUE
        )
      },
    
    .groups = "drop"
  ) %>%
  mutate(
    best_observed_BSJ_class =
      case_when(
        any_robust_BSJ ~
          "robust_BSJ_overlap",
        
        any_marginal_BSJ ~
          "marginal_BSJ_overlap",
        
        any_junction_codon_only ~
          "junction_codon_only",
        
        TRUE ~
          "non_BSJ"
      ),
    
    BSJ_position_ambiguity =
      case_when(
        n_BSJ_crossing_occurrences > 0 &
          n_non_BSJ_occurrences > 0 ~
          "mixed_crossing_and_non_crossing_positions",
        
        n_peptide_occurrences > 1 ~
          "multiple_mapped_positions",
        
        TRUE ~
          "single_mapped_position"
      )
  )

# ------------------------------------------------------------
# 8. Ensure that the final set still contains 185 peptides
# ------------------------------------------------------------

integrated <- fread(
  "10_185_integrated_peptide_summary_fixed.tsv",
  data.table = FALSE
) %>%
  mutate(
    peptide = as.character(peptide)
  )

unresolved_summary <- candidate_orf %>%
  group_by(peptide) %>%
  summarise(
    n_total_candidate_ORF_mappings = n(),
    
    n_unresolved_ORF_mappings = sum(
      is.na(protein_sequence) |
        is.na(circRNA_length) |
        is.na(orf_start)
    ),
    
    .groups = "drop"
  )

bsj_summary <- integrated %>%
  dplyr::select(
    peptide,
    hit_any_sp,
    hit_any_tr
  ) %>%
  left_join(
    unresolved_summary,
    by = "peptide"
  ) %>%
  left_join(
    bsj_summary_observed,
    by = "peptide"
  ) %>%
  mutate(
    # All 185 peptides are currently shared with Swiss-Prot,
    # so all values here are expected to be FALSE
    BSJ_sequence_unique =
      !(hit_any_sp %in% TRUE |
          hit_any_tr %in% TRUE),
    
    BSJ_interpretation =
      case_when(
        best_observed_BSJ_class ==
          "robust_BSJ_overlap" &
          hit_any_sp %in% TRUE ~
          paste0(
            "robust positional BSJ overlap, ",
            "but peptide is shared with Swiss-Prot"
          ),
        
        best_observed_BSJ_class ==
          "marginal_BSJ_overlap" ~
          paste0(
            "only one complete amino acid ",
            "on one side of BSJ"
          ),
        
        best_observed_BSJ_class ==
          "junction_codon_only" ~
          paste0(
            "BSJ is crossed only within ",
            "a junction codon"
          ),
        
        best_observed_BSJ_class ==
          "non_BSJ" ~
          "peptide does not cross BSJ",
        
        TRUE ~
          "BSJ status unresolved"
      )
  )

fwrite(
  bsj_summary,
  "14b_185_peptide_BSJ_summary.tsv",
  sep = "\t",
  quote = FALSE,
  na = "NA"
)


# ============================================================
# Merge the revised peptide-level BSJ annotations into the 185-peptide integrated master table
# Run directly in R / RStudio
# ============================================================

setwd(
  "./MS_unspecific/circRNA/result/circRNA_MS_MHC_integration"
)

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
})

master_file <- "10_185_integrated_peptide_summary_fixed.tsv"
bsj_file <- "14b_185_peptide_BSJ_summary_revised.tsv"

required_files <- c(master_file, bsj_file)
missing_files <- required_files[!file.exists(required_files)]

if (length(missing_files) > 0) {
  stop(
    "The following files are missing:\n",
    paste(missing_files, collapse = "\n"),
    "\nPlace the revised 14b file in the current working directory first."
  )
}

master <- fread(master_file, data.table = FALSE) %>%
  mutate(peptide = as.character(peptide))

bsj <- fread(bsj_file, data.table = FALSE) %>%
  mutate(peptide = as.character(peptide))

# ------------------------------------------------------------
# Input checks
# ------------------------------------------------------------

if (nrow(master) != 185 || n_distinct(master$peptide) != 185) {
  stop("The integrated master table does not contain 185 unique peptides; check the input file first.")
}

if (nrow(bsj) != 185 || n_distinct(bsj$peptide) != 185) {
  stop("The revised BSJ table does not contain 185 unique peptides; check the input file first.")
}

missing_in_bsj <- setdiff(master$peptide, bsj$peptide)
extra_in_bsj <- setdiff(bsj$peptide, master$peptide)

if (length(missing_in_bsj) > 0 || length(extra_in_bsj) > 0) {
  stop(
    "The peptide sets in the two tables are inconsistent.\n",
    "Present in the master table but missing from the BSJ table: ", length(missing_in_bsj), "\n",
    "Present in the BSJ table but missing from the master table: ", length(extra_in_bsj)
  )
}

# ------------------------------------------------------------
# Remove existing BSJ columns from the master table to avoid duplicate fields
# ------------------------------------------------------------

old_bsj_columns <- intersect(
  names(master),
  setdiff(names(bsj), "peptide")
)

master_clean <- master %>%
  dplyr::select(-any_of(old_bsj_columns))

# ------------------------------------------------------------
# Merge
# ------------------------------------------------------------

integrated_bsj <- master_clean %>%
  left_join(
    bsj,
    by = "peptide"
  ) %>%
  mutate(
    # Under the current evidence, none of the 185 peptides supports a circRNA-specific origin
    evidence_statement = case_when(
      best_observed_BSJ_class == "marginal_BSJ_overlap" ~
        paste0(
          "shotgun-MS-supported and HLA-predicted; ",
          "weak positional BSJ overlap; sequence shared with Swiss-Prot"
        ),
      
      best_observed_BSJ_class == "junction_codon_only" ~
        paste0(
          "shotgun-MS-supported and HLA-predicted; ",
          "BSJ lies within one codon; sequence shared with Swiss-Prot"
        ),
      
      TRUE ~
        paste0(
          "shotgun-MS-supported and HLA-predicted; ",
          "non-BSJ peptide shared with Swiss-Prot"
        )
    ),
    
    # This is an operational priority for manual spectrum review, not a measure of biological confidence
    manual_spectrum_review_priority = case_when(
      best_observed_BSJ_class == "marginal_BSJ_overlap" ~
        "P1_review_BSJ_marginal",
      
      best_observed_BSJ_class == "junction_codon_only" ~
        "P2_review_BSJ_codon_only",
      
      best_binder_class == "SB" ~
        "P3_review_non_BSJ_SB",
      
      TRUE ~
        "P4_review_non_BSJ_WB"
    ),
    
    recommended_claim =
      "source-ambiguous MS-supported predicted HLA-binding peptide",
    
    supports_natural_HLA_presentation = FALSE,
    supports_circRNA_specific_origin = FALSE
  ) %>%
  arrange(
    factor(
      manual_spectrum_review_priority,
      levels = c(
        "P1_review_BSJ_marginal",
        "P2_review_BSJ_codon_only",
        "P3_review_non_BSJ_SB",
        "P4_review_non_BSJ_WB"
      )
    ),
    best_EL_rank,
    best_PEP,
    desc(best_Score)
  )

# ------------------------------------------------------------
# Export the integrated master table
# ------------------------------------------------------------

fwrite(
  integrated_bsj,
  "15_185_integrated_with_revised_BSJ.tsv",
  sep = "\t",
  quote = FALSE,
  na = "NA"
)

# Export only the six peptides with positional BSJ overlap
bsj_positional_6 <- integrated_bsj %>%
  filter(peptide_crosses_BSJ_any %in% TRUE)

fwrite(
  bsj_positional_6,
  "15a_6_BSJ_positional_candidates.tsv",
  sep = "\t",
  quote = FALSE,
  na = "NA"
)

# ------------------------------------------------------------
# Summary statistics
# Avoid count() to prevent function conflicts with other packages
# ------------------------------------------------------------

summary_table <- dplyr::bind_rows(
  
  # BSJ classification
  integrated_bsj %>%
    dplyr::transmute(
      peptide = peptide,
      category = best_observed_BSJ_class
    ) %>%
    dplyr::group_by(category) %>%
    dplyr::summarise(
      n_unique_peptides = dplyr::n_distinct(peptide),
      .groups = "drop"
    ) %>%
    dplyr::mutate(
      summary_type = "BSJ_class"
    ),
  
  # Manual spectrum review priority
  integrated_bsj %>%
    dplyr::transmute(
      peptide = peptide,
      category = manual_spectrum_review_priority
    ) %>%
    dplyr::group_by(category) %>%
    dplyr::summarise(
      n_unique_peptides = dplyr::n_distinct(peptide),
      .groups = "drop"
    ) %>%
    dplyr::mutate(
      summary_type = "manual_review_priority"
    ),
  
  # HLA binder classification
  integrated_bsj %>%
    dplyr::transmute(
      peptide = peptide,
      category = best_binder_class
    ) %>%
    dplyr::group_by(category) %>%
    dplyr::summarise(
      n_unique_peptides = dplyr::n_distinct(peptide),
      .groups = "drop"
    ) %>%
    dplyr::mutate(
      summary_type = "HLA_binder_class"
    )
  
) %>%
  dplyr::select(
    summary_type,
    category,
    n_unique_peptides
  ) %>%
  dplyr::arrange(
    summary_type,
    category
  )

data.table::fwrite(
  summary_table,
  "15b_integrated_candidate_summary_counts.tsv",
  sep = "\t",
  quote = FALSE,
  na = "NA"
)
# ------------------------------------------------------------
# Final check
# ------------------------------------------------------------


###
# Check that each summary contains all 185 unique peptides
summary_check <- summary_table %>%
  dplyr::group_by(summary_type) %>%
  dplyr::summarise(
    total_unique_peptides = sum(n_unique_peptides),
    .groups = "drop"
  )


if (any(summary_check$total_unique_peptides != 185)) {
  warning("At least one summary category does not total 185 peptides; check missing values or classification fields.")
}

###

# ============================================================
# Generate the manual MS/MS spectrum review queue
# Automatically identify PEP / Score columns in the representative PSM file
# Run directly in R / RStudio
# ============================================================

setwd(
  "./MS_unspecific/circRNA/result/circRNA_MS_MHC_integration"
)

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
})

integrated_file <- "15_185_integrated_with_revised_BSJ.tsv"

# ============================================================
# 1. Check the integrated master table
# ============================================================

if (!file.exists(integrated_file)) {
  stop(
    "Integrated master table not found: ",
    integrated_file
  )
}

# ============================================================
# 2. Automatically locate the representative PSM file
# ============================================================

representative_psm_candidates <- c(
  "04b_representative_PSM.tsv",
  "04b_representative_PSMs.tsv",
  "04b_MS_representative_PSM.tsv"
)

existing_psm_files <- representative_psm_candidates[
  file.exists(representative_psm_candidates)
]

if (length(existing_psm_files) == 0) {
  stop(
    "No representative PSM file was found. Check whether the following files exist:\n",
    paste(
      representative_psm_candidates,
      collapse = "\n"
    )
  )
}

representative_psm_file <- existing_psm_files[1]


# ============================================================
# 3. Read the integrated master table and representative PSM table
# ============================================================

integrated_bsj <- data.table::fread(
  integrated_file,
  data.table = FALSE
) %>%
  dplyr::mutate(
    peptide = as.character(peptide)
  )

representative_psm <- data.table::fread(
  representative_psm_file,
  data.table = FALSE
) %>%
  dplyr::mutate(
    peptide = as.character(peptide)
  )


if (
  nrow(integrated_bsj) != 185 ||
  dplyr::n_distinct(integrated_bsj$peptide) != 185
) {
  stop(
    "The integrated master table does not contain 185 unique peptides; check the input file first."
  )
}

if (!"peptide" %in% names(representative_psm)) {
  stop(
    "The representative PSM file does not contain a peptide column."
  )
}

# ============================================================
# 4. Automatically identify quality fields
# ============================================================

find_first_column <- function(data, candidates) {
  
  matched <- candidates[
    candidates %in% names(data)
  ]
  
  if (length(matched) == 0) {
    return(NA_character_)
  }
  
  matched[1]
}

pep_column <- find_first_column(
  representative_psm,
  c(
    "PEP",
    "pep",
    "best_PEP",
    "Best_PEP",
    "Posterior.Error.Probability"
  )
)

score_column <- find_first_column(
  representative_psm,
  c(
    "Score",
    "score",
    "best_Score",
    "Best_Score",
    "Andromeda.Score"
  )
)

delta_score_column <- find_first_column(
  representative_psm,
  c(
    "Delta_score",
    "Delta.Score",
    "delta_score",
    "best_Delta_score",
    "Best_Delta_score"
  )
)

mass_error_column <- find_first_column(
  representative_psm,
  c(
    "mass_error_ppm",
    "Mass.error..ppm.",
    "Mass.error.ppm",
    "best_abs_mass_error_ppm"
  )
)

pif_column <- find_first_column(
  representative_psm,
  c(
    "PIF",
    "pif",
    "best_PIF"
  )
)


# ============================================================
# 5. Create standardized representative PSM quality fields
# ============================================================

representative_psm_tmp <- representative_psm %>%
  dplyr::mutate(
    representative_PEP = if (!is.na(pep_column)) {
      suppressWarnings(
        as.numeric(.data[[pep_column]])
      )
    } else {
      rep(NA_real_, dplyr::n())
    },
    
    representative_Score = if (!is.na(score_column)) {
      suppressWarnings(
        as.numeric(.data[[score_column]])
      )
    } else {
      rep(NA_real_, dplyr::n())
    },
    
    representative_Delta_score =
      if (!is.na(delta_score_column)) {
        suppressWarnings(
          as.numeric(.data[[delta_score_column]])
        )
      } else {
        rep(NA_real_, dplyr::n())
      },
    
    representative_mass_error_ppm =
      if (!is.na(mass_error_column)) {
        suppressWarnings(
          as.numeric(.data[[mass_error_column]])
        )
      } else {
        rep(NA_real_, dplyr::n())
      },
    
    representative_PIF =
      if (!is.na(pif_column)) {
        suppressWarnings(
          as.numeric(.data[[pif_column]])
        )
      } else {
        rep(NA_real_, dplyr::n())
      },
    
    # Used for sorting; missing PEP values are placed last
    PEP_sort = dplyr::if_else(
      is.na(representative_PEP),
      Inf,
      representative_PEP
    ),
    
    # Used for sorting; missing Score values are placed last
    Score_sort = dplyr::if_else(
      is.na(representative_Score),
      -Inf,
      representative_Score
    )
  )

# ============================================================
# 6. Select the best representative PSM for each peptide
# Prioritize the lowest PEP, then the highest Score
# ============================================================

representative_psm_best <- representative_psm_tmp %>%
  dplyr::arrange(
    peptide,
    PEP_sort,
    dplyr::desc(Score_sort)
  ) %>%
  dplyr::group_by(peptide) %>%
  dplyr::slice_head(n = 1) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(
    representative_PSM_found = TRUE
  ) %>%
  dplyr::select(
    -dplyr::any_of(
      c(
        "PEP_sort",
        "Score_sort"
      )
    )
  )


# ============================================================
# 7. Retain only representative PSM fields needed for manual review
# ============================================================

preferred_psm_columns <- c(
  "peptide",
  
  "representative_PSM_found",
  "representative_PEP",
  "representative_Score",
  "representative_Delta_score",
  "representative_mass_error_ppm",
  "representative_PIF",
  
  "dataset",
  "dataset_type",
  "comparison",
  
  "raw_file",
  "Raw.file",
  
  "scan_number",
  "scan",
  "Scan.number",
  "MS.MS.scan.number",
  
  "charge",
  "Charge",
  
  "intensity_coverage",
  "peak_coverage",
  
  "leading_razor_protein",
  "Leading.razor.protein",
  
  "proteins",
  "Proteins",
  
  "psm_key",
  "id"
)

available_psm_columns <- intersect(
  preferred_psm_columns,
  names(representative_psm_best)
)

representative_psm_best <- representative_psm_best %>%
  dplyr::select(
    dplyr::all_of(
      available_psm_columns
    )
  )


# ============================================================
# 8. Add the representative_ prefix to non-standard fields
# Avoid name conflicts with fields in the integrated master table
# ============================================================

columns_to_prefix <- setdiff(
  names(representative_psm_best),
  c(
    "peptide",
    "representative_PSM_found",
    "representative_PEP",
    "representative_Score",
    "representative_Delta_score",
    "representative_mass_error_ppm",
    "representative_PIF"
  )
)

representative_psm_best <- representative_psm_best %>%
  dplyr::rename_with(
    .fn = function(x) {
      paste0(
        "representative_",
        x
      )
    },
    .cols = dplyr::all_of(
      columns_to_prefix
    )
  )

# ============================================================
# 9. Merge candidate information with representative PSMs
# ============================================================

review_queue <- integrated_bsj %>%
  dplyr::left_join(
    representative_psm_best,
    by = "peptide"
  ) %>%
  dplyr::mutate(
    representative_PSM_found =
      representative_PSM_found %in% TRUE,
    
    manual_review_status =
      "not_reviewed",
    
    peptide_sequence_confirmed =
      NA,
    
    major_b_y_ions_supported =
      NA,
    
    unexplained_major_peaks =
      NA,
    
    spectrum_quality =
      NA_character_,
    
    spectrum_exported =
      NA,
    
    annotated_spectrum_file =
      NA_character_,
    
    manual_reviewer =
      NA_character_,
    
    manual_review_date =
      NA_character_,
    
    reviewer_notes =
      NA_character_
  ) %>%
  dplyr::arrange(
    factor(
      manual_spectrum_review_priority,
      levels = c(
        "P1_review_BSJ_marginal",
        "P2_review_BSJ_codon_only",
        "P3_review_non_BSJ_SB",
        "P4_review_non_BSJ_WB"
      )
    ),
    best_PEP,
    dplyr::desc(best_Score)
  )

# ============================================================
# 10. Export the complete manual review queue
# ============================================================

data.table::fwrite(
  review_queue,
  "16_manual_spectrum_review_queue.tsv",
  sep = "\t",
  quote = FALSE,
  na = "NA"
)

# ============================================================
# 11. First batch: six positional BSJ candidates
# ============================================================

review_batch_1 <- review_queue %>%
  dplyr::filter(
    manual_spectrum_review_priority %in%
      c(
        "P1_review_BSJ_marginal",
        "P2_review_BSJ_codon_only"
      )
  )

data.table::fwrite(
  review_batch_1,
  "16a_manual_review_batch1_6_BSJ_peptides.tsv",
  sep = "\t",
  quote = FALSE,
  na = "NA"
)

# ============================================================
# 12. Export candidates lacking a representative PSM
# ============================================================

missing_representative_psm <- review_queue %>%
  dplyr::filter(
    !representative_PSM_found
  )

data.table::fwrite(
  missing_representative_psm,
  "16b_candidates_missing_representative_PSM.tsv",
  sep = "\t",
  quote = FALSE,
  na = "NA"
)


# ============================================================
# Correct mass-error fields in the integrated master table
# Use only reliable mass-error values from PXD044963 LabelFree data
# ============================================================

setwd(
  "./MS_unspecific/circRNA/result/circRNA_MS_MHC_integration"
)

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
})

safe_min <- function(x) {
  x <- suppressWarnings(as.numeric(x))
  x <- x[!is.na(x)]
  
  if (length(x) == 0) {
    return(NA_real_)
  }
  
  min(x)
}

safe_median <- function(x) {
  x <- suppressWarnings(as.numeric(x))
  x <- x[!is.na(x)]
  
  if (length(x) == 0) {
    return(NA_real_)
  }
  
  median(x)
}

# ------------------------------------------------------------
# 1. Read PSM-level data
# ------------------------------------------------------------

ms_raw <- fread(
  "03_candidate_PSM_source_main_strand.tsv",
  data.table = FALSE
) %>%
  mutate(
    peptide = as.character(peptide),
    dataset = as.character(dataset),
    dataset_type = as.character(dataset_type),
    mass_error_ppm =
      suppressWarnings(as.numeric(mass_error_ppm))
  )

# ------------------------------------------------------------
# 2. Remove duplicated rows caused by one PSM mapping to multiple ORFs
# ------------------------------------------------------------

ms_psm_unique <- ms_raw %>%
  arrange(
    peptide,
    psm_key
  ) %>%
  distinct(
    peptide,
    psm_key,
    .keep_all = TRUE
  )


# ------------------------------------------------------------
# 3. Use only plausible mass-error values from LabelFree data
# ------------------------------------------------------------

lfq_mass_error <- ms_psm_unique %>%
  filter(
    dataset == "PXD044963",
    dataset_type == "LabelFree"
  ) %>%
  group_by(peptide) %>%
  summarise(
    n_LFQ_PSM_total = n(),
    
    n_LFQ_PSM_with_mass_error = sum(
      !is.na(mass_error_ppm) &
        abs(mass_error_ppm) <= 50
    ),
    
    n_LFQ_PSM_missing_mass_error = sum(
      is.na(mass_error_ppm)
    ),
    
    n_LFQ_PSM_mass_error_outlier = sum(
      !is.na(mass_error_ppm) &
        abs(mass_error_ppm) > 50
    ),
    
    best_LFQ_abs_mass_error_ppm = safe_min(
      abs(
        mass_error_ppm[
          !is.na(mass_error_ppm) &
            abs(mass_error_ppm) <= 50
        ]
      )
    ),
    
    median_LFQ_abs_mass_error_ppm = safe_median(
      abs(
        mass_error_ppm[
          !is.na(mass_error_ppm) &
            abs(mass_error_ppm) <= 50
        ]
      )
    ),
    
    .groups = "drop"
  )

# ------------------------------------------------------------
# 4. Determine mass-error reporting status for each peptide
# ------------------------------------------------------------

dataset_support <- ms_psm_unique %>%
  group_by(peptide) %>%
  summarise(
    has_TMT_PSM = any(
      dataset == "PXD037581"
    ),
    
    has_LFQ_PSM = any(
      dataset == "PXD044963"
    ),
    
    .groups = "drop"
  )

mass_error_annotation <- dataset_support %>%
  left_join(
    lfq_mass_error,
    by = "peptide"
  ) %>%
  mutate(
    mass_error_reporting_status = case_when(
      has_LFQ_PSM &
        n_LFQ_PSM_with_mass_error > 0 ~
        "reported_from_LabelFree_PSM",
      
      has_LFQ_PSM &
        (
          is.na(n_LFQ_PSM_with_mass_error) |
            n_LFQ_PSM_with_mass_error == 0
        ) ~
        "LabelFree_mass_error_missing",
      
      !has_LFQ_PSM &
        has_TMT_PSM ~
        "not_reported_TMT_only",
      
      TRUE ~
        "mass_error_unavailable"
    ),
    
    mass_error_note = case_when(
      mass_error_reporting_status ==
        "reported_from_LabelFree_PSM" ~
        paste0(
          "Mass error calculated only from ",
          "PXD044963 label-free PSMs"
        ),
      
      mass_error_reporting_status ==
        "LabelFree_mass_error_missing" ~
        paste0(
          "Label-free PSM available, but ",
          "mass error was missing"
        ),
      
      mass_error_reporting_status ==
        "not_reported_TMT_only" ~
        paste0(
          "Not reported because available ",
          "PSMs were from TMT10 data"
        ),
      
      TRUE ~
        "Mass error unavailable"
    )
  )

# ------------------------------------------------------------
# 5. Update the integrated master table
# ------------------------------------------------------------

integrated_old <- fread(
  "15_185_integrated_with_revised_BSJ.tsv",
  data.table = FALSE
) %>%
  mutate(
    peptide = as.character(peptide)
  )

integrated_mass_fixed <- integrated_old %>%
  dplyr::select(
    -any_of(
      c(
        "best_abs_mass_error_ppm",
        "best_LFQ_abs_mass_error_ppm",
        "median_LFQ_abs_mass_error_ppm",
        "n_LFQ_PSM_total",
        "n_LFQ_PSM_with_mass_error",
        "n_LFQ_PSM_missing_mass_error",
        "n_LFQ_PSM_mass_error_outlier",
        "has_TMT_PSM",
        "has_LFQ_PSM",
        "mass_error_reporting_status",
        "mass_error_note"
      )
    )
  ) %>%
  left_join(
    mass_error_annotation,
    by = "peptide"
  )

# ------------------------------------------------------------
# 6. Save the corrected master table
# ------------------------------------------------------------

fwrite(
  integrated_mass_fixed,
  "17_185_integrated_with_corrected_mass_error.tsv",
  sep = "\t",
  quote = FALSE,
  na = "NA"
)

# ------------------------------------------------------------

# 18 ----------------------------------------------------------------------


# Final table cleanup -------------------------------------------------------------------

# ============================================================
# Final cleanup of the integrated master table:
# 1. Correct candidate classification
# 2. Replace NA LFQ counts with 0
# 3. Clarify the EL binder field
# 4. Clarify MS summary field names
# 5. Add the best EL/BA allele
# 6. Add PIF availability
# 7. Export the full audit table and compact publication table
# ============================================================

setwd(
  "./MS_unspecific/circRNA/result/circRNA_MS_MHC_integration"
)

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
})

master_file <- "17_185_integrated_with_corrected_mass_error.tsv"
mhc_file <- "05b_NetMHCpan_deduplicated_by_allele.tsv"
mapping_file <- "09d_peptide_overlap_mapping.tsv"

required_files <- c(
  master_file,
  mhc_file,
  mapping_file
)

missing_files <- required_files[
  !file.exists(required_files)
]

if (length(missing_files) > 0) {
  stop(
    "The following input files are missing:\n",
    paste(missing_files, collapse = "\n")
  )
}

# ============================================================
# 1. Read the current official master table
# ============================================================

master <- data.table::fread(
  master_file,
  data.table = FALSE
) %>%
  dplyr::mutate(
    peptide = as.character(peptide)
  )

if (
  nrow(master) != 185 ||
  dplyr::n_distinct(master$peptide) != 185
) {
  stop(
    "The input master table does not contain 185 unique peptides."
  )
}


# ============================================================
# 2. Correct LFQ count columns
#
# For TMT-only peptides:
# the LFQ PSM count should be 0 rather than NA.
#
# However, LFQ mass-error measurements remain NA.
# ============================================================

lfq_count_columns <- intersect(
  c(
    "n_LFQ_PSM_total",
    "n_LFQ_PSM_with_mass_error",
    "n_LFQ_PSM_missing_mass_error",
    "n_LFQ_PSM_mass_error_outlier"
  ),
  names(master)
)

master_clean <- master %>%
  dplyr::mutate(
    dplyr::across(
      dplyr::all_of(lfq_count_columns),
      ~ dplyr::coalesce(
        suppressWarnings(as.integer(.x)),
        0L
      )
    )
  )

# ============================================================
# 3. Correct candidate classification
#
# Do not let multiple_circRNAs mask sharing with Swiss-Prot.
# Treat source uniqueness and circRNA mapping as two independent dimensions.
# ============================================================

master_clean <- master_clean %>%
  dplyr::mutate(
    candidate_summary_class = paste(
      source_uniqueness_class,
      circRNA_mapping_class,
      sep = ";"
    )
  ) %>%
  dplyr::select(
    -dplyr::any_of(
      "current_candidate_class"
    )
  )

# ============================================================
# 4. Rename fields for clearer interpretation
# ============================================================

rename_map <- c(
  # MS summary values are aggregated across different PSMs, so avoid the prefix best
  "min_PEP" = "best_PEP",
  "max_Score" = "best_Score",
  "max_Delta_score" = "best_Delta_score",
  "max_PIF" = "best_PIF",
  
  # HLA binder class is defined by EL rank
  "best_EL_binder_class" = "best_binder_class",
  "predicted_EL_binder_classes" =
    "predicted_binder_classes"
)

for (new_name in names(rename_map)) {
  
  old_name <- rename_map[[new_name]]
  
  if (
    old_name %in% names(master_clean) &&
    !new_name %in% names(master_clean)
  ) {
    names(master_clean)[
      names(master_clean) == old_name
    ] <- new_name
  }
}

# ============================================================
# 5. Add the PIF availability field
# ============================================================

if ("max_PIF" %in% names(master_clean)) {
  
  master_clean <- master_clean %>%
    dplyr::mutate(
      PIF_availability = dplyr::case_when(
        is.na(max_PIF) ~ "not_available",
        TRUE ~ "available"
      )
    )
}

# ============================================================
# 6. Read valid peptide-circRNA pairs
# ============================================================

valid_pairs <- data.table::fread(
  mapping_file,
  data.table = FALSE
) %>%
  dplyr::mutate(
    peptide = as.character(peptide),
    MS_circRNA_id =
      as.character(MS_circRNA_id),
    same_circRNA =
      as.character(same_circRNA)
  ) %>%
  dplyr::filter(
    toupper(same_circRNA) == "TRUE"
  ) %>%
  dplyr::transmute(
    peptide,
    circRNA_id = MS_circRNA_id
  ) %>%
  dplyr::filter(
    !is.na(peptide),
    peptide != "",
    !is.na(circRNA_id),
    circRNA_id != ""
  ) %>%
  dplyr::distinct()


# ============================================================
# 7. Extract the best allele from the NetMHCpan table
#
# Read only required fields to avoid loading irrelevant large columns.
# ============================================================

mhc <- data.table::fread(
  mhc_file,
  select = c(
    "source_id",
    "peptide",
    "HLA_allele",
    "best_EL_rank",
    "best_EL_score",
    "best_BA_rank",
    "best_BA_score"
  ),
  data.table = FALSE
) %>%
  dplyr::mutate(
    source_id = as.character(source_id),
    peptide = as.character(peptide),
    HLA_allele = as.character(HLA_allele),
    
    best_EL_rank =
      suppressWarnings(as.numeric(best_EL_rank)),
    
    best_EL_score =
      suppressWarnings(as.numeric(best_EL_score)),
    
    best_BA_rank =
      suppressWarnings(as.numeric(best_BA_rank)),
    
    best_BA_score =
      suppressWarnings(as.numeric(best_BA_score))
  ) %>%
  # Filter by peptide first to substantially reduce downstream data volume
  dplyr::filter(
    peptide %in% master_clean$peptide
  )


# ============================================================
# 8. Recover circRNA_id from short source_id values
#
# Example:
# ORF101_0126885 -> hsa_circ_0126885
# ============================================================

mhc <- mhc %>%
  dplyr::mutate(
    circRNA_id = dplyr::case_when(
      
      grepl(
        "^ORF[0-9]+_[0-9]+$",
        source_id
      ) ~ paste0(
        "hsa_circ_",
        sub(
          "^ORF[0-9]+_",
          "",
          source_id
        )
      ),
      
      grepl(
        "^ORF[0-9]+_hsa_circ_[0-9]+",
        source_id
      ) ~ sub(
        "^ORF[0-9]+_(hsa_circ_[0-9]+).*$",
        "\\1",
        source_id
      ),
      
      TRUE ~ NA_character_
    )
  )

# Retain only source-consistent MHC records
mhc_valid <- mhc %>%
  dplyr::inner_join(
    valid_pairs,
    by = c(
      "peptide",
      "circRNA_id"
    )
  )


if (
  dplyr::n_distinct(mhc_valid$peptide) != 185
) {
  warning(
    "Valid MHC records do not cover all 185 peptides."
  )
}

# ============================================================
# 9. Select the best EL allele for each peptide
#
# Sorting rules:
# 1. Lowest EL rank
# 2. Highest EL score
# 3. Alphabetical allele order to ensure reproducible tie-breaking
# ============================================================

best_el <- mhc_valid %>%
  dplyr::filter(
    !is.na(best_EL_rank)
  ) %>%
  dplyr::arrange(
    peptide,
    best_EL_rank,
    dplyr::desc(best_EL_score),
    HLA_allele
  ) %>%
  dplyr::group_by(peptide) %>%
  dplyr::slice_head(n = 1) %>%
  dplyr::ungroup() %>%
  dplyr::transmute(
    peptide,
    best_EL_allele = HLA_allele,
    selected_best_EL_rank = best_EL_rank,
    selected_best_EL_score = best_EL_score
  )

# ============================================================
# 10. Select the best BA allele for each peptide
# ============================================================

best_ba <- mhc_valid %>%
  dplyr::filter(
    !is.na(best_BA_rank)
  ) %>%
  dplyr::arrange(
    peptide,
    best_BA_rank,
    dplyr::desc(best_BA_score),
    HLA_allele
  ) %>%
  dplyr::group_by(peptide) %>%
  dplyr::slice_head(n = 1) %>%
  dplyr::ungroup() %>%
  dplyr::transmute(
    peptide,
    best_BA_allele = HLA_allele,
    selected_best_BA_rank = best_BA_rank,
    selected_best_BA_score = best_BA_score
  )

# ============================================================
# 11. Merge the best allele assignments back into the master table
# ============================================================

master_final <- master_clean %>%
  dplyr::select(
    -dplyr::any_of(
      c(
        "best_EL_allele",
        "best_BA_allele",
        "selected_best_EL_rank",
        "selected_best_EL_score",
        "selected_best_BA_rank",
        "selected_best_BA_score"
      )
    )
  ) %>%
  dplyr::left_join(
    best_el,
    by = "peptide"
  ) %>%
  dplyr::left_join(
    best_ba,
    by = "peptide"
  )

# ============================================================
# 12. Check consistency between original summary optima and the selected allele rows
# ============================================================

master_final <- master_final %>%
  dplyr::mutate(
    EL_rank_consistency = dplyr::case_when(
      is.na(best_EL_rank) &
        is.na(selected_best_EL_rank) ~ TRUE,
      
      abs(
        best_EL_rank -
          selected_best_EL_rank
      ) < 1e-10 ~ TRUE,
      
      TRUE ~ FALSE
    ),
    
    EL_score_consistency = dplyr::case_when(
      is.na(best_EL_score) &
        is.na(selected_best_EL_score) ~ TRUE,
      
      abs(
        best_EL_score -
          selected_best_EL_score
      ) < 1e-10 ~ TRUE,
      
      TRUE ~ FALSE
    ),
    
    BA_rank_consistency = dplyr::case_when(
      is.na(best_BA_rank) &
        is.na(selected_best_BA_rank) ~ TRUE,
      
      abs(
        best_BA_rank -
          selected_best_BA_rank
      ) < 1e-10 ~ TRUE,
      
      TRUE ~ FALSE
    )
  )


# ============================================================
# 13. Add table-note-level fields
# ============================================================

master_final <- master_final %>%
  dplyr::mutate(
    HLA_binder_definition =
      paste0(
        "NetMHCpan EL rank: ",
        "SB < 0.5%; WB < 2%"
      ),
    
    MS_summary_note =
      paste0(
        "min_PEP, max_Score, max_Delta_score ",
        "and max_PIF may originate from different PSMs"
      ),
    
    source_conclusion =
      paste0(
        "peptide sequence shared with Swiss-Prot; ",
        "shotgun MS cannot establish circRNA-specific origin"
      )
  )

# ============================================================
# 14. Reorder key columns
# ============================================================

priority_columns <- c(
  "peptide",
  
  "circRNA_ids",
  "n_circRNA",
  "circRNA_mapping_class",
  "candidate_summary_class",
  
  "MS_source_ids",
  "n_PSM",
  "n_raw_files",
  "datasets",
  "comparisons",
  
  "min_PEP",
  "max_Score",
  "max_Delta_score",
  "max_PIF",
  "PIF_availability",
  
  "has_TMT_PSM",
  "has_LFQ_PSM",
  "n_LFQ_PSM_total",
  "n_LFQ_PSM_with_mass_error",
  "median_LFQ_abs_mass_error_ppm",
  "mass_error_reporting_status",
  
  "HLA_alleles",
  "n_HLA_alleles",
  "best_EL_allele",
  "best_EL_rank",
  "best_EL_score",
  "best_EL_binder_class",
  
  "best_BA_allele",
  "best_BA_rank",
  "best_BA_score",
  
  "best_observed_BSJ_class",
  "peptide_crosses_BSJ_any",
  "max_complete_aa_each_side",
  "occurrence_pattern",
  
  "source_uniqueness_class",
  "BSJ_sequence_unique",
  "supports_circRNA_specific_origin",
  "supports_natural_HLA_presentation",
  
  "manual_spectrum_review_priority",
  "recommended_claim",
  "source_conclusion"
)

priority_columns <- intersect(
  priority_columns,
  names(master_final)
)

remaining_columns <- setdiff(
  names(master_final),
  priority_columns
)

master_final <- master_final %>%
  dplyr::select(
    dplyr::all_of(priority_columns),
    dplyr::all_of(remaining_columns)
  )

# ============================================================
# 15. Export the full audit master table
# ============================================================

audit_output <- paste0(
  "18_185_integrated_final_audit_table.tsv"
)

data.table::fwrite(
  master_final,
  audit_output,
  sep = "\t",
  quote = FALSE,
  na = "NA"
)

# ============================================================
# 16. Generate the compact publication table
# ============================================================

publication_columns <- c(
  "peptide",
  
  "circRNA_ids",
  "circRNA_mapping_class",
  
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
  
  "best_EL_allele",
  "best_EL_rank",
  "best_EL_score",
  "best_EL_binder_class",
  "n_HLA_alleles",
  
  "best_observed_BSJ_class",
  "max_complete_aa_each_side",
  "occurrence_pattern",
  
  "source_uniqueness_class",
  "source_conclusion",
  
  "manual_spectrum_review_priority"
)

publication_columns <- intersect(
  publication_columns,
  names(master_final)
)

publication_table <- master_final %>%
  dplyr::select(
    dplyr::all_of(publication_columns)
  )

data.table::fwrite(
  publication_table,
  "18a_185_publication_candidate_summary.tsv",
  sep = "\t",
  quote = FALSE,
  na = "NA"
)

# ============================================================
# 17. Export the field dictionary
# ============================================================

field_dictionary <- data.frame(
  field = c(
    "min_PEP",
    "max_Score",
    "max_Delta_score",
    "max_PIF",
    "median_LFQ_abs_mass_error_ppm",
    "best_EL_allele",
    "best_EL_rank",
    "best_EL_binder_class",
    "best_observed_BSJ_class",
    "max_complete_aa_each_side",
    "source_uniqueness_class",
    "supports_circRNA_specific_origin"
  ),
  
  description = c(
    "Minimum PEP observed across supporting PSMs",
    "Maximum Andromeda score across supporting PSMs",
    "Maximum delta score across supporting PSMs",
    "Maximum PIF across supporting PSMs",
    "Median absolute precursor mass error using label-free PSMs only",
    "HLA allele associated with the lowest NetMHCpan EL rank",
    "Lowest NetMHCpan EL rank across valid source-consistent predictions",
    "Binder class based on EL rank: SB <0.5%, WB <2%",
    "Best observed peptide-level BSJ positional class",
    "Minimum number of complete amino acids on either side of BSJ; NA for non-BSJ peptides",
    "Whether peptide sequence is shared with a known protein database entry",
    "Whether current evidence supports a circRNA-specific peptide origin"
  ),
  
  stringsAsFactors = FALSE
)

data.table::fwrite(
  field_dictionary,
  "18b_final_table_field_dictionary.tsv",
  sep = "\t",
  quote = FALSE,
  na = "NA"
)


# Final reorganization --------------------------------------------------------------------
library(data.table)
library(dplyr)

# Re-read the original table 18 to avoid using modified objects from the current environment
audit18 <- data.table::fread(
  "18_185_integrated_final_audit_table.tsv",
  data.table = FALSE
)

# Confirm that all required fields are present
required_HLA_columns <- c(
  "best_EL_score",
  "best_BA_score",
  "selected_best_EL_score",
  "selected_best_BA_score"
)

missing_HLA_columns <- setdiff(
  required_HLA_columns,
  names(audit18)
)

if (length(missing_HLA_columns) > 0) {
  stop(
    "Missing fields: ",
    paste(missing_HLA_columns, collapse = ", ")
  )
}

# Call dplyr::rename explicitly to avoid function conflicts
audit18_fixed <- audit18 %>%
  dplyr::rename(
    max_EL_score_across_alleles = best_EL_score,
    max_BA_score_across_alleles = best_BA_score
  ) %>%
  dplyr::mutate(
    # Score from the same row as the best rank and best allele
    best_EL_score = selected_best_EL_score,
    best_BA_score = selected_best_BA_score,
    
    EL_score_consistency = dplyr::case_when(
      is.na(best_EL_score) &
        is.na(selected_best_EL_score) ~ TRUE,
      
      !is.na(best_EL_score) &
        !is.na(selected_best_EL_score) &
        abs(
          best_EL_score -
            selected_best_EL_score
        ) < 1e-10 ~ TRUE,
      
      TRUE ~ FALSE
    ),
    
    BA_score_consistency = dplyr::case_when(
      is.na(best_BA_score) &
        is.na(selected_best_BA_score) ~ TRUE,
      
      !is.na(best_BA_score) &
        !is.na(selected_best_BA_score) &
        abs(
          best_BA_score -
            selected_best_BA_score
        ) < 1e-10 ~ TRUE,
      
      TRUE ~ FALSE
    )
  )

audit18_fixed <- audit18_fixed %>%
  dplyr::relocate(
    best_EL_allele,
    best_EL_rank,
    best_EL_score,
    best_EL_binder_class,
    max_EL_score_across_alleles,
    .after = n_HLA_alleles
  ) %>%
  dplyr::relocate(
    best_BA_allele,
    best_BA_rank,
    best_BA_score,
    max_BA_score_across_alleles,
    .after = max_EL_score_across_alleles
  )
#
data.table::fwrite(
  audit18_fixed,
  "19_185_integrated_final_audit_table_HLA_fixed.tsv",
  sep = "\t",
  quote = FALSE,
  na = "NA"
)

# ============================================================
# Generate the 19a compact publication table
# Extract core fields from the corrected audit18_fixed table
# ============================================================

publication_columns_19a <- c(
  "peptide",
  
  # circRNA mapping
  "circRNA_ids",
  "circRNA_mapping_class",
  
  # MS evidence
  "n_PSM",
  "n_raw_files",
  "datasets",
  "comparisons",
  "min_PEP",
  "max_Score",
  "max_Delta_score",
  "max_PIF",
  "PIF_availability",
  
  # Mass error
  "median_LFQ_abs_mass_error_ppm",
  "mass_error_reporting_status",
  
  # HLA prediction
  "best_EL_allele",
  "best_EL_rank",
  "best_EL_score",
  "best_EL_binder_class",
  "n_HLA_alleles",
  
  # BSJ information
  "best_observed_BSJ_class",
  "max_complete_aa_each_side",
  "occurrence_pattern",
  
  # Source attribution conclusion
  "source_uniqueness_class",
  "source_conclusion",
  
  # Manual spectrum review
  "manual_spectrum_review_priority"
)

# Check whether any planned fields are missing
missing_publication_columns <- setdiff(
  publication_columns_19a,
  names(audit18_fixed)
)

if (length(missing_publication_columns) > 0) {
  warning(
    "The following fields were not found in audit18_fixed and will be skipped automatically:\n",
    paste(
      missing_publication_columns,
      collapse = "\n"
    )
  )
}

# Retain only fields that are actually present
publication19 <- audit18_fixed %>%
  dplyr::select(
    dplyr::all_of(
      intersect(
        publication_columns_19a,
        names(audit18_fixed)
      )
    )
  )

# Save the compact publication table
data.table::fwrite(
  publication19,
  "19a_185_publication_candidate_summary_HLA_fixed.tsv",
  sep = "\t",
  quote = FALSE,
  na = "NA"
)


# Prepare peptides for UniProt query -----------------------------------------------------------
# Query only peptides not overlapping IEDB entries
peptide <- paste(publication19$peptide, collapse = ",")

writeLines(peptide, "./MS_unspecific/circRNA/result/IEDB_uniprot_search_new/need_file/circRNA_final.txt")

# Process UniProt query results ----------------------------------------------------------


uniprot_result <- function(input_file, output_dir) {
  
  # Read the current input Excel file
  circ_uniprot <- readxl::read_excel(input_file)
  
  table_copy <- circ_uniprot %>%
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
    "; peptide count: ", dplyr::n_distinct(peptide_detail$Peptide)
  )
  
  return(peptide_class)
}
# Run
input_dir <- paste0(
  "./MS_unspecific/circRNA/result/",
  "IEDB_uniprot_search_new/need_file"
)

output_dir <- paste0(
  "./MS_unspecific/circRNA/result/",
  "IEDB_uniprot_search_new"
)

sample_names <- c(
  "circRNA_final_1",
  "circRNA_final_2",
  "circRNA_final_3"
)

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

##
circRNA_final_1_output <- read.csv("./MS_unspecific/circRNA/result/IEDB_uniprot_search_new/circRNA_final_1_output.csv")
circRNA_final_2_output <- read.csv("./MS_unspecific/circRNA/result/IEDB_uniprot_search_new/circRNA_final_2_output.csv")
circRNA_final_3_output <- read.csv("./MS_unspecific/circRNA/result/IEDB_uniprot_search_new/circRNA_final_3_output.csv")

all_uniprot_search <- rbind(circRNA_final_1_output,circRNA_final_2_output,circRNA_final_3_output)


write.table(all_uniprot_search,"./MS_unspecific/circRNA/result/IEDB_uniprot_search_new/all_uniprot_search.txt",quote = FALSE,row.names = FALSE,sep = "\t")
###
# Match Swiss-Prot and TrEMBL IDs --------------------------------------------------------------

library(dplyr)
library(tidyr)
library(stringr)

library(data.table)

# Convert directly to data.table to avoid unnecessary object copies
setDT(all_uniprot_search)

sp_tr_map <- read.delim("./MS_unspecific/sp_tr_map.txt")

setDT(sp_tr_map)

# Add row IDs to the original table
all_uniprot_search[, row_id := .I]

# 1. Split Matched_Entries
circRNA_entries_long <- all_uniprot_search[
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

# 2. Extract accessions that actually need to be queried
needed_entries <- unique(
  c(
    circRNA_entries_long$accession,
    circRNA_entries_long$parent_accession
  )
)


# 3. Filter only relevant entries from the large mapping table

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

# 4. Match accessions exactly
circRNA_entries_long[
  uniprot_map_small,
  on = .(accession = Entry),
  database_exact := i.database
]

# 5. Supplement isoform entries using the parent accession

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
# Add sp| or tr| prefixes:
circRNA_entries_long[
  ,
  accession_with_db := fifelse(
    !is.na(database_final),
    paste0(database_final, "|", accession),
    accession
  )
]

# 6. Classify match types
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
# 7. Collapse entries back into semicolon-separated strings

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

# 8. Merge annotations back into the original table
all_circrNA_search_with_db <- merge(
  all_uniprot_search,
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


# Replace the original column
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


write.table(all_circrNA_search_with_db,"./MS_unspecific/circRNA/result/IEDB_uniprot_search_new/all_circrNA_search_sp_tr.txt",quote = FALSE,row.names = FALSE,sep = "\t")

# Process UniProt query output ----------------------------------------------------------

# Format for the publication table ---------------------------------------------------------------

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
all_circRNA_search_sp_tr <- read.delim2("./MS_unspecific/circRNA/result/IEDB_uniprot_search/all_circrNA_search_sp_tr.txt")

circRNA_uniprot_publication  <- tidy_uniprot(
  all_circRNA_search_sp_tr
)


write.table(circRNA_uniprot_publication ,"./MS_unspecific/circRNA/result/IEDB_uniprot_search/new_all_circRNA_search_sp_tr.txt",quote = FALSE,row.names = FALSE,sep = "\t")

write.csv(circRNA_uniprot_publication,"./MS_unspecific/circRNA/result/IEDB_uniprot_search/new_all_circRNA_search_sp_tr.csv",quote = FALSE,row.names = FALSE,sep = ",")
