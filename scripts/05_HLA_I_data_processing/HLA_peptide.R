# ============================================================
# HLA-I prediction landscape across non-canonical peptide sources
# Sources: AS, circRNA, lncRNA, and pseudogene
# ============================================================

# HLA-I prediction landscape across four non-canonical sources

library(data.table)
library(ggplot2)
library(dplyr)
library(stringr)
library(patchwork)

outdir <- "./circpeptide4/new_HLA_part/result"

dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

paths <- list(

  AS = "./circpeptide4/AS_rMATS_jcast/netMHCpan/HLA_1/combined_f_sigAS.txt",

  circRNA = "./circpeptide4/circRNA/MHC/combined_f.txt",

  lncRNA = "./circpeptide4/lncRNA_known/combined_f.txt",

  Pseudogene = "./circpeptide4/pseudogenes/result/step8_netMHCpan/result/HLA_1/combined_f_new.txt"

)

# 1. Standardize NetMHCpan output column names ----------------

# Different source files may use different column names.

# Standardize them before downstream integration.

# AS/lncRNA/pseudogene may use X_EL_Rank.

# circRNA may use _EL_Rank, _EL-score, and related variants.

# Target standardized fields:

# Pos, Peptide, ID, type, X_core, X_icore, X_EL_score, X_EL_Rank, X_BA_score, X_BA_Rank

standardize_netmhc_cols <- function(dt, source_name) {

  # Normalize periods, hyphens, and spaces in column names.

  old_names <- names(dt)

  new_names <- old_names

  new_names <- gsub("\\.", "_", new_names)

  new_names <- gsub("-", "_", new_names)

  new_names <- gsub(" ", "_", new_names)

  new_names <- gsub("__+", "_", new_names)

  setnames(dt, old = old_names, new = new_names)

  # Handle circRNA columns that begin with an underscore.

  if ("_core" %in% names(dt)) {

    setnames(dt, "_core", "X_core")

  }

  if ("_icore" %in% names(dt)) {

    setnames(dt, "_icore", "X_icore")

  }

  if ("_EL_score" %in% names(dt)) {

    setnames(dt, "_EL_score", "X_EL_score")

  }

  if ("_EL_Rank" %in% names(dt)) {

    setnames(dt, "_EL_Rank", "X_EL_Rank")

  }

  if ("_BA_score" %in% names(dt)) {

    setnames(dt, "_BA_score", "X_BA_score")

  }

  if ("_BA_Rank" %in% names(dt)) {

    setnames(dt, "_BA_Rank", "X_BA_Rank")

  }

  # Handle lowercase names and other common aliases.

  if ("peptide" %in% names(dt)) {

    setnames(dt, "peptide", "Peptide")

  }

  if ("id" %in% names(dt)) {

    setnames(dt, "id", "ID")

  }

  if ("allele" %in% names(dt)) {

    setnames(dt, "allele", "type")

  }

  if ("HLA" %in% names(dt)) {

    setnames(dt, "HLA", "type")

  }

  # If X_EL_Rank is still absent, detect an EL + Rank column automatically.

  if (!"X_EL_Rank" %in% names(dt)) {

    rank_candidates <- names(dt)[

      grepl("EL", names(dt), ignore.case = TRUE) &

        grepl("Rank", names(dt), ignore.case = TRUE)

    ]

    if (length(rank_candidates) > 0) {

      setnames(dt, rank_candidates[1], "X_EL_Rank")

    } else {

      message("Columns in ", source_name, ":")

      print(names(dt))

      stop("No EL Rank-like column found in: ", source_name)

    }

  }

  # If X_EL_score is absent, detect an EL + score column automatically.

  if (!"X_EL_score" %in% names(dt)) {

    elscore_candidates <- names(dt)[

      grepl("EL", names(dt), ignore.case = TRUE) &

        grepl("score", names(dt), ignore.case = TRUE)

    ]

    if (length(elscore_candidates) > 0) {

      setnames(dt, elscore_candidates[1], "X_EL_score")

    }

  }

  # If X_BA_Rank is absent, detect a BA + Rank column automatically.

  if (!"X_BA_Rank" %in% names(dt)) {

    barank_candidates <- names(dt)[

      grepl("BA", names(dt), ignore.case = TRUE) &

        grepl("Rank", names(dt), ignore.case = TRUE)

    ]

    if (length(barank_candidates) > 0) {

      setnames(dt, barank_candidates[1], "X_BA_Rank")

    }

  }

  # If X_BA_score is absent, detect a BA + score column automatically.

  if (!"X_BA_score" %in% names(dt)) {

    bascore_candidates <- names(dt)[

      grepl("BA", names(dt), ignore.case = TRUE) &

        grepl("score", names(dt), ignore.case = TRUE)

    ]

    if (length(bascore_candidates) > 0) {

      setnames(dt, bascore_candidates[1], "X_BA_score")

    }

  }

  return(dt)

}

# 2. Read and standardize each NetMHCpan file -----------------

read_netmhc_file <- function(path, source_name) {

  message("Reading: ", source_name)

  dt <- fread(path)

  # Standardize column names.

  dt <- standardize_netmhc_cols(dt, source_name)

  # Check required columns.

  required_cols <- c("Peptide", "ID", "type", "X_EL_Rank")

  missing_cols <- setdiff(required_cols, names(dt))

  if (length(missing_cols) > 0) {

    message("Columns in ", source_name, ":")

    print(names(dt))

    stop("Missing columns in ", source_name, ": ", paste(missing_cols, collapse = ", "))

  }

  # Coerce prediction fields to numeric.

  dt[, X_EL_Rank := suppressWarnings(as.numeric(X_EL_Rank))]

  if ("X_EL_score" %in% names(dt)) {

    dt[, X_EL_score := suppressWarnings(as.numeric(X_EL_score))]

  }

  if ("X_BA_score" %in% names(dt)) {

    dt[, X_BA_score := suppressWarnings(as.numeric(X_BA_score))]

  }

  if ("X_BA_Rank" %in% names(dt)) {

    dt[, X_BA_Rank := suppressWarnings(as.numeric(X_BA_Rank))]

  }

  # Remove rows with non-numeric EL rank values, such as repeated headers.

  n_bad <- sum(is.na(dt$X_EL_Rank))

  if (n_bad > 0) {

    message(source_name, ": removing ", n_bad, " rows with non-numeric X_EL_Rank")

    dt <- dt[!is.na(X_EL_Rank)]

  }

  # Correct the pseudogene ID prefix.

  if (source_name == "Pseudogene") {

    dt[, ID := gsub("^LncID_", "pseudo_", ID)]

  }

  # Add source label.

  dt[, source := source_name]

  # Extract HLA locus.

  dt[, HLA_locus := dplyr::case_when(

    stringr::str_detect(type, "^HLA-A") ~ "HLA-A",

    stringr::str_detect(type, "^HLA-B") ~ "HLA-B",

    stringr::str_detect(type, "^HLA-C") ~ "HLA-C",

    TRUE ~ "Other"

  )]

  # Define binder class.

  dt[, binder_class := dplyr::case_when(

    X_EL_Rank < 0.5 ~ "Strong binder",

    X_EL_Rank >= 0.5 & X_EL_Rank < 2 ~ "Weak binder",

    TRUE ~ "Non-binder"

  )]

  return(dt)

}

# 3. Audit EL-rank ranges in the four combined_f files --------

rank_check <- lapply(names(paths), function(src) {

  dt_tmp <- read_netmhc_file(paths[[src]], src)

  data.frame(

    source = src,

    min_EL_rank = min(dt_tmp$X_EL_Rank, na.rm = TRUE),

    max_EL_rank = max(dt_tmp$X_EL_Rank, na.rm = TRUE),

    n_rows = nrow(dt_tmp),

    n_unique_peptides = data.table::uniqueN(dt_tmp$Peptide),

    n_NA_EL_rank = sum(is.na(dt_tmp$X_EL_Rank))

  )

})

rank_check <- dplyr::bind_rows(rank_check)

fwrite(

  rank_check,

  file.path(outdir, "rank_check_summary.csv")

)

# 4. Summarize HLA-I predictions for each source --------------

summarise_netmhc_source <- function(path, source_name, density_n = 50000) {

  message("Summarising: ", source_name)

  dt <- read_netmhc_file(path, source_name)

  # Retain predicted HLA-A/B/C binders with EL rank < 2.

  binders <- dt[

    X_EL_Rank < 2 & HLA_locus %in% c("HLA-A", "HLA-B", "HLA-C")

  ]

  # 1) Best EL rank at the unique-peptide level.

  best_rank <- binders[

    ,

    .(

      best_EL_Rank = min(X_EL_Rank, na.rm = TRUE),

      best_HLA = type[which.min(X_EL_Rank)],

      n_binding_alleles = uniqueN(type),

      n_binding_loci = uniqueN(HLA_locus)

    ),

    by = .(source, Peptide)

  ]

  best_rank[

    ,

    best_binder_class := dplyr::case_when(

      best_EL_Rank < 0.5 ~ "Strong binder",

      best_EL_Rank >= 0.5 & best_EL_Rank < 2 ~ "Weak binder",

      TRUE ~ "Non-binder"

    )

  ]

  # 2) Fig. 4A / Fig. S4B: unique peptide counts by best binder class.

  unique_binder_summary <- best_rank[

    ,

    .N,

    by = .(source, best_binder_class)

  ]

  # 3) Fig. 4B: peptide-HLA predicted binding events.

  event_summary <- binders[

    ,

    .N,

    by = .(source, HLA_locus, binder_class)

  ]

  # 4) Sample data for the optional Fig. S4A EL-rank distribution.

  set.seed(123)

  density_sample <- binders[

    ,

    .SD[sample(.N, min(.N, density_n))],

    by = .(source, HLA_locus)

  ][

    ,

    .(source, HLA_locus, X_EL_Rank)

  ]

  # 5) Fig. S4C: predicted binding events contributed by each HLA allele.

  allele_event_summary <- binders[

    ,

    .N,

    by = .(source, HLA_locus, type)

  ]

  # 6) Fig. S4D: defer extraction of representative peptide-HLA data.

  # Extract only top peptides later to avoid storing an unnecessarily large matrix.

  # 7) Save the best-rank summary for downstream figures.

  fwrite(

    best_rank,

    file.path(outdir, paste0(source_name, "_best_EL_rank_by_peptide.csv"))

  )

  return(list(

    best_rank = best_rank,

    unique_binder_summary = unique_binder_summary,

    event_summary = event_summary,

    density_sample = density_sample,

    allele_event_summary = allele_event_summary

  ))

}

summary_list <- lapply(names(paths), function(src) {

  summarise_netmhc_source(paths[[src]], src)

})

names(summary_list) <- names(paths)

best_rank_all <- rbindlist(lapply(summary_list, `[[`, "best_rank"))

unique_binder_summary <- rbindlist(lapply(summary_list, `[[`, "unique_binder_summary"))

event_summary <- rbindlist(lapply(summary_list, `[[`, "event_summary"))

density_sample_all <- rbindlist(lapply(summary_list, `[[`, "density_sample"))

allele_event_summary <- rbindlist(lapply(summary_list, `[[`, "allele_event_summary"))

# Save combined summary tables.

fwrite(best_rank_all, file.path(outdir, "all_sources_best_EL_rank_by_peptide.csv"))

fwrite(unique_binder_summary, file.path(outdir, "Fig4A_unique_binder_summary.csv"))

fwrite(event_summary, file.path(outdir, "Fig4B_event_summary.csv"))

fwrite(density_sample_all, file.path(outdir, "FigS4A_density_sample.csv"))

fwrite(allele_event_summary, file.path(outdir, "FigS4C_allele_event_summary.csv"))

# 5. Fig. 4A: strong/weak composition of unique predicted binders ----

comma_format2 <- function(x) format(x, big.mark = ",", scientific = FALSE, trim = TRUE)

fig4A_df <- unique_binder_summary %>%

  filter(best_binder_class %in% c("Strong binder", "Weak binder")) %>%

  mutate(

    source = factor(source, levels = c("AS", "circRNA", "lncRNA", "Pseudogene")),

    best_binder_class = factor(

      best_binder_class,

      levels = c("Strong binder", "Weak binder")

    )

  )

p4A <- ggplot(fig4A_df, aes(x = source, y = N, fill = best_binder_class)) +

  geom_col(width = 0.7) +

  scale_fill_manual(

    values = c(

      "Strong binder" = "#ea8f98",

      "Weak binder" = "#9692af"

    )

  ) +

  scale_y_continuous(labels = comma_format2) +

  labs(

    x = NULL,

    y = "Number of unique predicted binders",

    fill = NULL,

    title = "Unique candidate peptides with predicted HLA-I compatibility"

  ) +

  theme_classic(base_size = 13) +

  theme(

    legend.position = "top",

    plot.title = element_text(face = "bold")

  )

ggsave(

  file.path(outdir, "Fig4A_unique_predicted_binders.pdf"),

  p4A, width = 8, height = 6

)

ggsave(

  file.path(outdir, "Fig4A_unique_predicted_binders.png"),

  p4A, width = 8, height = 6, dpi = 300

)

# 6. Fig. 4B: source x HLA locus x binder-class events -------------

fig4B_df <- event_summary %>%

  filter(

    HLA_locus %in% c("HLA-A", "HLA-B", "HLA-C"),

    binder_class %in% c("Strong binder", "Weak binder")

  ) %>%

  mutate(

    source = factor(source, levels = c("AS", "circRNA", "lncRNA", "Pseudogene")),

    HLA_locus = factor(HLA_locus, levels = c("HLA-A", "HLA-B", "HLA-C")),

    binder_class = factor(binder_class, levels = c("Weak binder", "Strong binder"))

  )

p4B <- ggplot(fig4B_df, aes(x = HLA_locus, y = N, fill = binder_class)) +

  geom_col(width = 0.7) +

  facet_wrap(~ source, nrow = 1, scales = "free_y") +

  scale_fill_manual(

    values = c(

      "Weak binder" = "#9692af",

      "Strong binder" = "#ea8f98"

    )

  ) +

  scale_y_continuous(labels = comma_format2) +

  labs(

    x = NULL,

    y = "Predicted peptide-HLA binding events",

    fill = NULL,

    title = "Predicted HLA-I binding events across HLA loci"

  ) +

  theme_classic(base_size = 13) +

  theme(

    legend.position = "top",

    strip.background = element_blank(),

    strip.text = element_text(face = "bold"),

    plot.title = element_text(face = "bold")

  )

ggsave(

  file.path(outdir, "Fig4B_HLA_locus_binding_events.pdf"),

  p4B, width = 12, height = 5

)

ggsave(

  file.path(outdir, "Fig4B_HLA_locus_binding_events.png"),

  p4B, width = 12, height = 5, dpi = 300

)

# 7. Fig. 4C: best EL-rank violin/boxplot --------------------------

set.seed(123)

fig4C_dt <- as.data.table(best_rank_all)

fig4C_dt <- fig4C_dt[best_EL_Rank < 2]

fig4C_dt <- fig4C_dt[

  ,

  .SD[sample(.N, min(.N, 50000))],

  by = source

]

fig4C_dt[

  ,

  source := factor(source, levels = c("AS", "circRNA", "lncRNA", "Pseudogene"))

]

fig4C_df <- as.data.frame(fig4C_dt)

p4C <- ggplot(fig4C_df, aes(x = source, y = best_EL_Rank, fill = source)) +

  geom_violin(

    trim = TRUE,

    scale = "width",

    alpha = 0.85,

    color = "grey35",

    linewidth = 0.4

  ) +

  geom_boxplot(

    width = 0.13,

    outlier.shape = NA,

    fill = "white",

    color = "grey20",

    linewidth = 0.45

  ) +

  geom_hline(

    yintercept = 0.5,

    linetype = "dashed",

    linewidth = 0.45,

    color = "grey20"

  ) +

  geom_hline(

    yintercept = 2.0,

    linetype = "dashed",

    linewidth = 0.45,

    color = "grey20"

  ) +

  scale_fill_manual(

    values = c(

      "AS" = "#9e9e9e",

      "circRNA" = "#6baed6",

      "lncRNA" = "#9ecae1",

      "Pseudogene" = "#c7c7c7"

    )

  ) +

  scale_y_continuous(

    limits = c(0, 2.15),

    breaks = c(0, 0.5, 1.0, 1.5, 2.0)

  ) +

  labs(

    x = NULL,

    y = "Best EL rank per unique peptide",

    fill = NULL,

    title = "Best predicted HLA-I EL rank of unique candidate peptides"

  ) +

  theme_classic(base_size = 13) +

  theme(

    legend.position = "none",

    plot.title = element_text(face = "bold", hjust = 0),

    axis.text.x = element_text(size = 12),

    axis.text.y = element_text(size = 11),

    axis.title.y = element_text(size = 13),

    plot.margin = margin(8, 12, 8, 8)

  )

ggsave(

  file.path(outdir, "Fig4C_best_EL_rank_violin_boxplot_clean.pdf"),

  p4C, width = 8, height = 6

)

ggsave(

  file.path(outdir, "Fig4C_best_EL_rank_violin_boxplot_clean.png"),

  p4C, width = 8, height = 6, dpi = 300

)

# 8. Fig. 4D: HLA allele breadth stacked bar -----------------------

fig4D_df <- best_rank_all %>%

  dplyr::filter(best_EL_Rank < 2) %>%

  dplyr::mutate(

    source = factor(source, levels = c("AS", "circRNA", "lncRNA", "Pseudogene")),

    breadth_group = cut(

      n_binding_alleles,

      breaks = c(0, 1, 3, 6, Inf),

      labels = c("1 allele", "2-3 alleles", "4-6 alleles", ">6 alleles"),

      include.lowest = TRUE

    )

  ) %>%

  dplyr::count(source, breadth_group) %>%

  dplyr::group_by(source) %>%

  dplyr::mutate(percent = n / sum(n) * 100) %>%

  dplyr::ungroup()

p4D <- ggplot(fig4D_df, aes(x = source, y = percent, fill = breadth_group)) +

  geom_col(width = 0.7) +

  scale_fill_manual(

    values = c(

      "1 allele" = "#d9d9d9",

      "2-3 alleles" = "#bdbdbd",

      "4-6 alleles" = "#9e9ac8",

      ">6 alleles" = "#5e62a9"

    )

  ) +

  labs(

    x = NULL,

    y = "Unique predicted binders (%)",

    fill = "Binding HLA alleles",

    title = "HLA allele breadth of predicted binders"

  ) +

  theme_classic(base_size = 13) +

  theme(

    legend.position = "right",

    plot.title = element_text(face = "bold", hjust = 0),

    axis.text.x = element_text(size = 12),

    axis.text.y = element_text(size = 11),

    axis.title.y = element_text(size = 13),

    plot.margin = margin(8, 12, 8, 8)

  )

ggsave(

  file.path(outdir, "Fig4D_HLA_allele_breadth.pdf"),

  p4D, width = 8, height = 6

)

ggsave(

  file.path(outdir, "Fig4D_HLA_allele_breadth.png"),

  p4D, width = 8, height = 6, dpi = 300

)

# 9. Assemble main-text Fig. 4 -------------------------------

fig4 <- (p4A | p4B) / (p4C | p4D) +

  plot_annotation(

    title = "Predicted HLA-I compatibility landscape of candidate non-canonical peptides"

  ) &

  theme(

    plot.title = element_text(face = "bold", size = 18, hjust = 0.5)

  )

ggsave(

  file.path(outdir, "Fig4_HLA_I_compatibility_landscape.pdf"),

  fig4, width = 16, height = 11

)

ggsave(

  file.path(outdir, "Fig4_HLA_I_compatibility_landscape.png"),

  fig4, width = 16, height = 11, dpi = 300

)

# 10. Fig. S4B: strong/weak binder proportion ----------------------

figS4B_df <- best_rank_all %>%

  dplyr::filter(best_EL_Rank < 2) %>%

  dplyr::count(source, best_binder_class) %>%

  dplyr::group_by(source) %>%

  dplyr::mutate(percent = n / sum(n) * 100) %>%

  dplyr::ungroup() %>%

  dplyr::mutate(

    source = factor(source, levels = c("AS", "circRNA", "lncRNA", "Pseudogene")),

    best_binder_class = factor(

      best_binder_class,

      levels = c("Weak binder", "Strong binder")

    )

  )

pS4B <- ggplot(figS4B_df, aes(x = source, y = percent, fill = best_binder_class)) +

  geom_col(width = 0.7) +

  scale_fill_manual(

    values = c(

      "Weak binder" = "#5e62a9",

      "Strong binder" = "#c45161"

    )

  ) +

  labs(

    x = NULL,

    y = "Unique predicted binders (%)",

    fill = NULL,

    title = "Proportion of strong and weak predicted HLA-I binders"

  ) +

  theme_classic(base_size = 13) +

  theme(

    legend.position = "top",

    plot.title = element_text(face = "bold", hjust = 0),

    axis.text.x = element_text(size = 12),

    axis.text.y = element_text(size = 11),

    axis.title.y = element_text(size = 13)

  )

ggsave(

  file.path(outdir, "FigS4B_strong_weak_binder_proportion.pdf"),

  pS4B, width = 8, height = 6

)

ggsave(

  file.path(outdir, "FigS4B_strong_weak_binder_proportion.png"),

  pS4B, width = 8, height = 6, dpi = 300

)

# Save the corresponding summary table.

fwrite(

  figS4B_df,

  file.path(outdir, "FigS4B_strong_weak_binder_proportion_summary.csv")

)

# 11. Fig. S4C: top HLA allele heatmap -----------------------------

# Select the top 20 HLA alleles by total predicted binding events.

top_hla <- allele_event_summary %>%

  dplyr::group_by(type) %>%

  dplyr::summarise(total_events = sum(N), .groups = "drop") %>%

  dplyr::arrange(dplyr::desc(total_events)) %>%

  dplyr::slice_head(n = 20) %>%

  dplyr::pull(type)

figS4C_df <- allele_event_summary %>%

  dplyr::filter(type %in% top_hla) %>%

  dplyr::mutate(

    source = factor(source, levels = c("AS", "circRNA", "lncRNA", "Pseudogene")),

    type = factor(type, levels = rev(top_hla))

  )

pS4C <- ggplot(figS4C_df, aes(x = source, y = type, fill = N)) +

  geom_tile(

    color = "white",

    linewidth = 0.4,

    width = 0.65,

    height = 0.9

  ) +

  scale_fill_gradient(

    low = "#f2f2f2",

    high = "#5e62a9",

    labels = comma_format2

  ) +

  labs(

    x = NULL,

    y = "HLA allele",

    fill = "Events",

    title = "Top HLA-I alleles contributing predicted binding events"

  ) +

  theme_classic(base_size = 13) +

  theme(

    plot.title = element_text(face = "bold", hjust = 0),

    axis.text.x = element_text(size = 12),

    axis.text.y = element_text(size = 10),

    axis.title.y = element_text(size = 13),

    legend.title = element_text(size = 11),

    legend.text = element_text(size = 10)

  )

ggsave(

  file.path(outdir, "FigS4C_top_HLA_allele_heatmap.pdf"),

  pS4C, width = 7.5, height = 8

)

ggsave(

  file.path(outdir, "FigS4C_top_HLA_allele_heatmap.png"),

  pS4C, width = 7.5, height = 8, dpi = 300

)

# 12. Select representative peptides -------------------------

top_peptides <- best_rank_all %>%

  dplyr::filter(best_EL_Rank < 0.5) %>%   # retain strong binders only

  dplyr::group_by(source) %>%

  dplyr::arrange(best_EL_Rank, .by_group = TRUE) %>%

  dplyr::slice_head(n = 6) %>%

  dplyr::ungroup() %>%

  dplyr::select(source, Peptide, best_EL_Rank, best_HLA, n_binding_alleles)

fwrite(

  top_peptides,

  file.path(outdir, "FigS4D_top_representative_peptides.csv")

)

# 13. Extract representative peptide-HLA data ----------------

extract_top_peptide_hla <- function(path, source_name, top_peptides) {

  message("Extracting peptide-HLA data: ", source_name)

  dt <- read_netmhc_file(path, source_name)

  pep_vec <- top_peptides %>%

    dplyr::filter(source == source_name) %>%

    dplyr::pull(Peptide)

  dt_sub <- dt[

    Peptide %in% pep_vec &

      X_EL_Rank < 2 &

      HLA_locus %in% c("HLA-A", "HLA-B", "HLA-C"),

    .(source, Peptide, type, HLA_locus, X_EL_Rank)

  ]

  return(dt_sub)

}

top_peptide_hla <- rbindlist(lapply(names(paths), function(src) {

  extract_top_peptide_hla(paths[[src]], src, top_peptides)

}))

fwrite(

  top_peptide_hla,

  file.path(outdir, "FigS4D_top_peptide_HLA_long.csv")

)

# 14. Select top HLA alleles ---------------------------------

top_hla_for_heatmap <- top_peptide_hla %>%

  dplyr::count(type, sort = TRUE) %>%

  dplyr::slice_head(n = 15) %>%

  dplyr::pull(type)

# 15. Prepare heatmap data -----------------------------------

figS4D_df <- top_peptide_hla %>%

  dplyr::filter(type %in% top_hla_for_heatmap) %>%

  dplyr::group_by(source, Peptide, type) %>%

  dplyr::summarise(best_EL_Rank = min(X_EL_Rank), .groups = "drop") %>%

  dplyr::left_join(

    top_peptides %>% dplyr::select(source, Peptide, best_EL_Rank_overall = best_EL_Rank),

    by = c("source", "Peptide")

  ) %>%

  dplyr::mutate(

    source = factor(source, levels = c("AS", "circRNA", "lncRNA", "Pseudogene")),

    peptide_label = paste(source, Peptide, sep = " | ")

  )

# Order peptides by source and then by overall best EL rank within each source.

peptide_order <- top_peptides %>%

  dplyr::mutate(

    source = factor(source, levels = c("AS", "circRNA", "lncRNA", "Pseudogene")),

    peptide_label = paste(source, Peptide, sep = " | ")

  ) %>%

  dplyr::arrange(source, best_EL_Rank) %>%

  dplyr::pull(peptide_label)

figS4D_df <- figS4D_df %>%

  dplyr::mutate(

    peptide_label = factor(peptide_label, levels = rev(unique(peptide_order))),

    type = factor(type, levels = top_hla_for_heatmap)

  )

# 16. Representative peptide-HLA heatmap ---------------------------------------

pS4D <- ggplot(figS4D_df, aes(x = type, y = peptide_label, fill = best_EL_Rank)) +

  geom_tile(color = "white", linewidth = 0.35, width = 0.85, height = 0.85) +

  scale_fill_gradient(

    low = "#c45161",     # stronger binding, darker color

    high = "#f3eaea",    # weaker binding, lighter color

    limits = c(0, 2),

    name = "Best EL rank"

  ) +

  labs(

    x = "HLA allele",

    y = NULL,

    title = "Representative peptide-HLA predicted binding heatmap"

  ) +

  theme_classic(base_size = 12) +

  theme(

    plot.title = element_text(face = "bold", hjust = 0, size = 15),

    axis.text.x = element_text(angle = 45, hjust = 1, size = 10),

    axis.text.y = element_text(size = 9),

    legend.title = element_text(size = 10),

    legend.text = element_text(size = 9),

    plot.margin = margin(8, 12, 8, 8)

  )

ggsave(

  file.path(outdir, "FigS4D_representative_peptide_HLA_heatmap.pdf"),

  pS4D, width = 10, height = 8

)

ggsave(

  file.path(outdir, "FigS4D_representative_peptide_HLA_heatmap.png"),

  pS4D, width = 10, height = 8, dpi = 300

)
