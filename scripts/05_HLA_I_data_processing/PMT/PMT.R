library(readr)
library(Biostrings)
library(dplyr)
library(tidyr)
library(stringr)
library(networkD3)
library(htmlwidgets)
library(RColorBrewer)

# ============================================================
# 1. Prepare circRNA ORF sequences for NetChop
# ============================================================

X19_185_integrated_final_audit_table_HLA_fixed <- read_delim(
  "19_185_integrated_final_audit_table_HLA_fixed.tsv",
  delim = "\t",
  escape_double = FALSE,
  trim_ws = TRUE
)

fa1 <- readAAStringSet(
  "./circRNA_peptide4/circRNA/MaxQuant/need_file/only_circRNA/merge_circ_HC_eRA.fasta"
)

fa2 <- readAAStringSet(
  "./circRNA_peptide4/circRNA/MaxQuant/need_file/only_circRNA/merge_circ_HC_RA.fasta"
)

fa3 <- readAAStringSet(
  "./circRNA_peptide4/circRNA/MaxQuant/need_file/only_circRNA/merge_circ_RA_eRA.fasta"
)

all_fa <- c(fa1, fa2, fa3)
all_fa_dedup <- all_fa[!duplicated(names(all_fa))]

names(all_fa_dedup) <- sub(" .*", "", names(all_fa_dedup))
names(all_fa_dedup) <- str_extract(
  names(all_fa_dedup),
  "(?<=lcl\\|).+?(?=:\\d+)"
)

need_185_ID <- X19_185_integrated_final_audit_table_HLA_fixed %>%
  dplyr::select(MS_source_ids) %>%
  tidyr::separate_rows(MS_source_ids, sep = ";") %>%
  unique()

netchop_need_fasta <- all_fa_dedup[
  names(all_fa_dedup) %in% need_185_ID$MS_source_ids
]

netchop_need_fasta <- netchop_need_fasta[
  !duplicated(names(netchop_need_fasta))
]

output_dir <- "./MS_unspecific/circRNA/result/PMT/netchop"

dir.create(
  output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

stopifnot(length(netchop_need_fasta) > 0)

original_ids <- names(netchop_need_fasta)
all_sequences <- as.character(netchop_need_fasta)

# Deduplicate by amino acid sequence
keep <- !duplicated(all_sequences)
netchop_need_fasta_uniq <- netchop_need_fasta[keep]
unique_sequences <- all_sequences[keep]

# Generate short IDs for NetChop
new_ids <- sprintf(
  "nCp_%03d",
  seq_along(netchop_need_fasta_uniq)
)

names(netchop_need_fasta_uniq) <- new_ids

id_map <- data.frame(
  original_id = original_ids,
  sequence = all_sequences,
  stringsAsFactors = FALSE
) %>%
  dplyr::mutate(
    short_id = new_ids[
      match(sequence, unique_sequences)
    ]
  ) %>%
  dplyr::select(
    short_id,
    original_id
  ) %>%
  dplyr::arrange(
    match(short_id, new_ids),
    original_id
  )

write.csv(
  id_map,
  file.path(
    output_dir,
    "NetChop_ID_mapping.csv"
  ),
  row.names = FALSE,
  quote = FALSE
)

# Split NetChop input into batches of at most 100 sequences
chunk_size <- 100L
n <- length(netchop_need_fasta_uniq)

index_list <- split(
  seq_len(n),
  ceiling(seq_len(n) / chunk_size)
)

for (i in seq_along(index_list)) {
  fasta_part <- netchop_need_fasta_uniq[
    index_list[[i]]
  ]

  output_file <- file.path(
    output_dir,
    sprintf(
      "netChop_input_part%d.fasta",
      i
    )
  )

  Biostrings::writeXStringSet(
    fasta_part,
    filepath = output_file,
    format = "fasta"
  )
}

# NetChop result processing is performed in:
# ./MS_unspecific/circRNA/result/PMT/netchop/netChop.ipynb


# ============================================================
# 2. Prepare disease-associated HLA and TCR information
# ============================================================

SampleGroup2 <- read.delim("~/circRNA/a_info/SampleGroup2.txt")

need_sample <- SampleGroup2 %>%
  dplyr::filter(Group != "healthy")

sample_circ_HLA_1 <- read.delim2("~/TCR/need/sample_circ_HLA_1.txt")

df_split <- sample_circ_HLA_1 %>%
  tidyr::separate_rows(HLA_1, sep = ",") %>%
  unique()

# Convert HLA allele names to the format required by TABR-BERT
df_split <- df_split %>%
  dplyr::mutate(
    allele = HLA_1 %>%
      stringr::str_replace("HLA-A", "HLA-A*") %>%
      stringr::str_replace("HLA-B", "HLA-B*") %>%
      stringr::str_replace("HLA-C", "HLA-C*")
  ) %>%
  dplyr::select(name, allele) %>%
  unique() %>%
  dplyr::filter(name %in% need_sample$Sample)

all_SRR_TCR <- read.delim("~/TCR/all_SRR_TCR.txt")

all_SRR_TCR <- all_SRR_TCR %>%
  dplyr::select(SRR, CDR3_amino_acids) %>%
  unique()

colnames(SampleGroup2)[1] <- colnames(all_SRR_TCR)[1]

all_SRR_TCR <- left_join(
  all_SRR_TCR,
  SampleGroup2,
  by = "SRR"
)

tcr_hc <- all_SRR_TCR %>%
  dplyr::filter(Group == "healthy")

tcr_eRA_RA <- all_SRR_TCR %>%
  dplyr::filter(Group != "healthy")

IS_in_hc <- tcr_eRA_RA$CDR3_amino_acids %in%
  tcr_hc$CDR3_amino_acids

tcr_eRA_RA$IS_in_hc <- IS_in_hc

tcr_eRA_RA <- tcr_eRA_RA %>%
  dplyr::filter(IS_in_hc == "FALSE")


# ============================================================
# 3. Prepare TABR-BERT input
# ============================================================

HLA_TAP <- read.delim(
  "./MS_unspecific/circRNA/PMT/HLA_TAP.txt"
)

HLA_TAP <- HLA_TAP %>%
  dplyr::filter(TAP == "TAP") %>%
  dplyr::select(peptide, HLA_type) %>%
  dplyr::mutate(
    allele = HLA_type %>%
      stringr::str_replace("HLA-A", "HLA-A*") %>%
      stringr::str_replace("HLA-B", "HLA-B*") %>%
      stringr::str_replace("HLA-C", "HLA-C*")
  ) %>%
  dplyr::select(peptide, allele) %>%
  unique()

for_tabr_bert <- left_join(
  HLA_TAP,
  df_split,
  by = "allele",
  relationship = "many-to-many"
)

colnames(for_tabr_bert)[3] <- "SRR"

for_tabr_bert <- left_join(
  for_tabr_bert,
  tcr_eRA_RA,
  by = "SRR",
  relationship = "many-to-many"
)

for_tabr_bert <- na.omit(for_tabr_bert)
colnames(for_tabr_bert)[4] <- "cdr3"

# Remove CDR3 sequences containing underscores from TRUST4 output
for_tabr_bert <- for_tabr_bert %>%
  dplyr::filter(!grepl("_", cdr3))

write.csv(
  for_tabr_bert,
  "./MS_unspecific/circRNA/PMT/TABR_BERT/need_file/for_sofwer_srr.csv"
)

final <- for_tabr_bert %>%
  dplyr::select(peptide, allele, cdr3)

write.table(
  final,
  "./MS_unspecific/circRNA/PMT/TABR_BERT/need_file/input.txt",
  quote = FALSE,
  row.names = FALSE,
  sep = "\t"
)

# Remove CDR3 sequences containing question marks for CSV input
final <- final %>%
  dplyr::filter(!grepl("?", cdr3, fixed = TRUE))

write.csv(
  x = final,
  file = "./MS_unspecific/circRNA/PMT/TABR_BERT/need_file/input.csv"
)


# ============================================================
# 4. Generate peptide-HLA Sankey diagram
# ============================================================

# Retain the 30 most frequent peptides
top_pep <- HLA_TAP %>%
  dplyr::count(peptide, name = "total") %>%
  dplyr::slice_max(
    total,
    n = 30,
    with_ties = FALSE
  ) %>%
  dplyr::pull(peptide)

df <- HLA_TAP %>%
  dplyr::filter(peptide %in% top_pep) %>%
  dplyr::select(peptide, allele) %>%
  dplyr::filter(
    !is.na(peptide),
    !is.na(allele)
  )

links_raw <- df %>%
  dplyr::count(
    peptide,
    allele,
    name = "value"
  ) %>%
  dplyr::rename(
    source_name = peptide,
    target_name = allele
  )

pep_levels <- unique(links_raw$source_name)
allele_levels <- unique(links_raw$target_name)

nodes <- data.frame(
  name = c(
    pep_levels,
    allele_levels
  ),
  stringsAsFactors = FALSE
)

nodes$group <- nodes$name

links <- links_raw %>%
  dplyr::mutate(
    source = match(source_name, nodes$name) - 1,
    target = match(target_name, nodes$name) - 1
  ) %>%
  dplyr::select(
    source,
    target,
    value
  )

soft_pal <- c(
  brewer.pal(12, "Set3"),
  brewer.pal(9, "Pastel1"),
  brewer.pal(8, "Pastel2")
)

mycol2 <- rep(
  soft_pal,
  length.out = nrow(nodes)
)

colourScale <- sprintf(
  'd3.scaleOrdinal().domain(["%s"]).range(["%s"])',
  paste(nodes$group, collapse = '\",\"'),
  paste(mycol2, collapse = '\",\"')
)

p <- sankeyNetwork(
  Links = links,
  Nodes = nodes,
  Source = "source",
  Target = "target",
  Value = "value",
  NodeID = "name",
  NodeGroup = "group",
  fontSize = 15,
  fontFamily = "Times New Roman",
  nodeWidth = 36,
  nodePadding = 8,
  iterations = 0,
  width = 720,
  height = 800,
  colourScale = colourScale
)

print(p)

saveWidget(
  p,
  file = "./MS_unspecific/circRNA/PMT/sankey_diagram.html",
  selfcontained = TRUE
)


# ============================================================
# 5. Filter TABR-BERT prediction results
# ============================================================

output <- read.csv(
  "./MS_unspecific/circRNA/PMT/TABR_BERT/result/output.csv"
)

output <- output %>%
  dplyr::filter(rank >= 0.901) %>%
  dplyr::select(
    peptide,
    allele,
    cdr3
  ) %>%
  unique()
