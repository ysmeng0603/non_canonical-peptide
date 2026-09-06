# 1. Process the count matrix -------------------------------------------------
pseudogene.featureCounts.strict <- read.delim("./pseudogenes/result/step2_feature_count/pseudogene.featureCounts.strict.txt", comment.char="#")
summary(colSums(pseudogene.featureCounts.strict[,7:ncol(pseudogene.featureCounts.strict)]))

new_colnames <- colnames(pseudogene.featureCounts.strict)
new_colnames

new_colnames[7:length(new_colnames)] <- gsub(
  ".*(SRR[0-9]+).*",
  "\\1",
  new_colnames[7:length(new_colnames)]
)
new_colnames

colnames(pseudogene.featureCounts.strict) <- new_colnames
colnames(pseudogene.featureCounts.strict)

# Differential expression analysis -------------------------------------------
count_matrix <- pseudogene.featureCounts.strict[, 7:ncol(pseudogene.featureCounts.strict)]  
rownames(count_matrix) <- pseudogene.featureCounts.strict$Geneid
# Prepare group information
SampleGroup2 <- read.delim("~/A_info/SampleGroup2.txt") # 
SampleGroup2 <- SampleGroup2 %>% 
  dplyr::filter(Sample %in% colnames(count_matrix)) # 
SampleGroup2$Group <- as.factor(SampleGroup2$Group)
all(colnames(count_matrix) == SampleGroup2$Sample) # T

###
dds <- DESeqDataSetFromMatrix(
  countData = count_matrix,
  colData   = SampleGroup2,
  design    = ~ Group
)
dds
#dds <- dds[rowSums(counts(dds)) >= 10, ]
keep <- rowSums(counts(dds) >= 10) >= 0.15*nrow(SampleGroup2) #
dds <- dds[keep, ]

#
dds$Group <- relevel(dds$Group, ref = "healthy")
dds <- DESeq(
  dds,
  minReplicatesForReplace = Inf  #
)
#
res <- results(
  dds,
  contrast = c("Group", "eRA", "healthy"),
  independentFiltering = FALSE   # 
)

res <- res[order(res$padj), ]
summary(res)
res_df <- as.data.frame(res)
res_df$pseudogene_id <- rownames(res_df) # 

res_df <- na.omit(res_df) # 

res_df$significance <- ifelse(res_df$padj < 0.05 & abs(res_df$log2FoldChange) > 1, 
                              ifelse(res_df$log2FoldChange > 0, "Up", "Down"), 
                              "Stable")
write.table(res_df,"./pseudogenes/result/step3_deseq2/pseudogene_diff_H_eRA.txt",quote = FALSE,sep = "\t")


# Hc-RA -------------------------------------------------------------------
res <- results(
  dds,
  contrast = c("Group", "RA", "healthy"),
  independentFiltering = FALSE   # 
)

res <- res[order(res$padj), ]
summary(res)
res_df <- as.data.frame(res)
res_df$pseudogene_id <- rownames(res_df) # 
res_df <- na.omit(res_df) # 
# 添加一列用于火山图的显著性标记
res_df$significance <- ifelse(res_df$padj < 0.05 & abs(res_df$log2FoldChange) > 1, 
                              ifelse(res_df$log2FoldChange > 0, "Up", "Down"), 
                              "Stable")
write.table(res_df,"./pseudogenes/result/step3_deseq2/pseudogene_diff_H_RA.txt",quote = FALSE,sep = "\t")
# RA_eRA ------------------------------------------------------------------
#
res <- results(
  dds,
  contrast = c("Group", "eRA", "RA"),
  independentFiltering = FALSE   #
)

res <- res[order(res$padj), ]
summary(res)
res_df <- as.data.frame(res)
res_df$pseudogene_id <- rownames(res_df) # 


res_df <- na.omit(res_df) # 
res_df$significance <- ifelse(res_df$padj < 0.05 & abs(res_df$log2FoldChange) > 1, 
                              ifelse(res_df$log2FoldChange > 0, "Up", "Down"), 
                              "Stable")
write.table(res_df,"./pseudogenes/result/step3_deseq2/pseudogene_diff_RA_eRA.txt",quote = FALSE,sep = "\t")

# Import differential expression results --------------------------------------
setwd("./pseudogenes/result/step3_deseq2")
#
pseudogene_diff_H_eRA <- read.delim("./pseudogenes/result/step3_deseq2/pseudogene_diff_H_eRA.txt")
pseudogene_diff_H_RA <- read.delim("./pseudogenes/result/step3_deseq2/pseudogene_diff_H_RA.txt")
pseudogene_diff_RA_eRA <- read.delim("./pseudogenes/result/step3_deseq2/pseudogene_diff_RA_eRA.txt")
#

head(pseudogene_diff_H_eRA)
table(pseudogene_diff_H_eRA$significance)

table(pseudogene_diff_H_RA$significance)

#
pdf("./pseudogenes/result/step3_deseq2/plot/pseudogenes_H_eRA.pdf")
ggplot(pseudogene_diff_H_eRA, aes(x = log2FoldChange, y = -log10(padj), color = significance)) +
  geom_point(alpha = 0.6, size = 1.5) +
  scale_color_manual(values = c("Stable" = "gray", 
                                "Up" = "red", 
                                "Down" = "blue")) +  
  labs(title = "HC-vs-eRA volcano plot",
       x = "log2FoldChange",
       y = "-log10(padj)") +
  geom_vline(xintercept = c(-1, 1), lty = 4, col = "black", lwd = 0.5) +         
  geom_hline(yintercept = -log10(0.05), lty = 4, col = "black", lwd = 0.5) +     
  theme_minimal() +
  theme(
    text = element_text(size = 14),
    panel.grid = element_blank(),  # 
    panel.border = element_rect(color = "black", fill = NA),  # 
    axis.line = element_line(color = "black"),  # 
    aspect.ratio = 1,  # 
    plot.title = element_text(hjust = 0.5)  # 
  )
dev.off()

#
pdf("./pseudogenes/result/step3_deseq2/plot/pseudogenes_H_RA.pdf")
ggplot(pseudogene_diff_H_RA, aes(x = log2FoldChange, y = -log10(padj), color = significance)) +
  geom_point(alpha = 0.6, size = 1.5) +
  scale_color_manual(values = c("Stable" = "gray", 
                                "Up" = "red", 
                                "Down" = "blue")) +  # 
  labs(title = "HC-vs-RA volcano plot",
       x = "log2FoldChange",
       y = "-log10(padj)") +
  geom_vline(xintercept = c(-1, 1), lty = 4, col = "black", lwd = 0.5) +         # 
  geom_hline(yintercept = -log10(0.05), lty = 4, col = "black", lwd = 0.5) +     # 
  theme_minimal() +
  theme(
    text = element_text(size = 14),
    panel.grid = element_blank(),  # 
    panel.border = element_rect(color = "black", fill = NA),  # 
    axis.line = element_line(color = "black"),  # 
    aspect.ratio = 1,  # 
    plot.title = element_text(hjust = 0.5)  # 
  )
dev.off()
#
pdf("./pseudogenes/result/step3_deseq2/plot/pseudogenes_RA_eRA.pdf")
ggplot(pseudogene_diff_RA_eRA, aes(x = log2FoldChange, y = -log10(padj), color = significance)) +
  geom_point(alpha = 0.6, size = 1.5) +
  scale_color_manual(values = c("Stable" = "gray", 
                                "Up" = "red", 
                                "Down" = "blue")) +  # 
  labs(title = "RA-vs-eRA volcano plot",
       x = "log2FoldChange",
       y = "-log10(padj)") +
  geom_vline(xintercept = c(-1, 1), lty = 4, col = "black", lwd = 0.5) +         # 
  geom_hline(yintercept = -log10(0.05), lty = 4, col = "black", lwd = 0.5) +     # 
  theme_minimal() +
  theme(
    text = element_text(size = 14),
    panel.grid = element_blank(),  # 
    panel.border = element_rect(color = "black", fill = NA),  # 
    axis.line = element_line(color = "black"),  # 
    aspect.ratio = 1,  #  
    plot.title = element_text(hjust = 0.5)  # 
  )
dev.off()


###
need_1 <- pseudogene_diff_H_eRA %>% 
  dplyr::filter(significance !="Stable") %>% 
  dplyr::mutate(group = "H_eRA") %>% 
  dplyr::select(group,everything()) # 
write.table(need_1,"./pseudogenes/result/step3_deseq2/only_fold_change/HC_eRA.txt",quote = FALSE,sep = "\t")

need_2 <- pseudogene_diff_H_RA %>% 
  dplyr::filter(significance !="Stable") %>% 
  dplyr::mutate(group = "H_RA") %>% 
  dplyr::select(group,everything()) # 
write.table(need_2,"./pseudogenes/result/step3_deseq2/only_fold_change/HC_RA.txt",quote = FALSE,sep = "\t")

need_3 <- pseudogene_diff_RA_eRA %>% 
  dplyr::filter(significance !="Stable") %>% 
  dplyr::mutate(group = "RA_eRA") %>% 
  dplyr::select(group,everything()) # 
write.table(need_3,"./pseudogenes/result/step3_deseq2/only_fold_change/RA_eRA.txt",quote = FALSE,sep = "\t")

all <- rbind(need_1,need_2,need_3)
all_2 <- all %>% unique()
write.table(all,"./pseudogenes/result/step3_deseq2/only_fold_change/all_pseudogene_group.txt",quote = FALSE,sep = "\t")

# all <- all %>% 
#   dplyr::select(group,pseudogene_id) %>% unique() # 

length(unique(all$pseudogene_id)) # 
table(all_2$significance)


library(Biostrings)
library(dplyr)
library(stringr)
library(readr)

setwd("./pseudogenes")

tx_genecode <- readDNAStringSet("./pseudogenes/gencode.v49.transcripts.fa")
tx_genecode

# Build the FASTA header information table:
check <- names(tx_genecode) %>% as.data.frame()
colnames(check) <- "full_header"

check$ENST <- sub("\\|.*$", "", check$full_header)
head(check)
# Step 2: Extract the gene_id-to-transcript_id mapping from the pseudogene GTF
gtf_file <- "annotation/gencode.v49.pseudogene.clean.gtf"

gtf <- read.delim(
  gtf_file,
  header = FALSE,
  comment.char = "#",
  sep = "\t",
  stringsAsFactors = FALSE,
  quote = ""
)

colnames(gtf) <- c(
  "chr", "source", "feature", "start", "end",
  "score", "strand", "frame", "attribute"
)

pseudo_tx <- gtf %>%
  filter(feature == "transcript") %>%
  mutate(
    gene_id = str_match(attribute, 'gene_id "([^"]+)"')[, 2],
    transcript_id = str_match(attribute, 'transcript_id "([^"]+)"')[, 2],
    gene_name = str_match(attribute, 'gene_name "([^"]+)"')[, 2],
    gene_type = str_match(attribute, 'gene_type "([^"]+)"')[, 2],
    transcript_type = str_match(attribute, 'transcript_type "([^"]+)"')[, 2]
  ) %>%
  select(gene_id, transcript_id, gene_name, gene_type, transcript_type) %>%
  distinct()

head(pseudo_tx)
dim(pseudo_tx)

table(pseudo_tx$transcript_id %in% check$ENST)
###

all_pseudogene_group <- read.delim("~./pseudogenes/result/step3_deseq2/only_fold_change/all_pseudogene_group.txt")
all_pseudo_gene <- unique(all_pseudogene_group$pseudogene_id)
length(all_pseudo_gene)
head(all_pseudo_gene)
# Convert all differentially expressed pseudogene ENSG IDs to ENST IDs
all_de_pseudo_tx <- pseudo_tx %>%
  filter(gene_id %in% all_pseudo_gene)

dim(all_de_pseudo_tx)
head(all_de_pseudo_tx)

table(all_pseudo_gene %in% pseudo_tx$gene_id)

#
matched_check_all <- check %>%
  filter(ENST %in% all_de_pseudo_tx$transcript_id)

dim(matched_check_all)
head(matched_check_all)
# Extract full header
full_genecode_name_all <- matched_check_all$full_header %>% as.character()

all_de_pseudo_fasta <- tx_genecode[full_genecode_name_all]

all_de_pseudo_fasta

name_map_all <- all_de_pseudo_tx %>%
  select(gene_id, transcript_id, gene_name, gene_type) %>%
  distinct()

matched_check_all2 <- matched_check_all %>%
  left_join(name_map_all, by = c("ENST" = "transcript_id"))

names(all_de_pseudo_fasta) <- paste(
  matched_check_all2$ENST,
  matched_check_all2$gene_id,
  matched_check_all2$gene_name,
  matched_check_all2$gene_type,
  sep = "|"
)
name_map_all <- all_de_pseudo_tx %>%
  select(gene_id, transcript_id, gene_name, gene_type) %>%
  distinct()

matched_check_all2 <- matched_check_all %>%
  left_join(name_map_all, by = c("ENST" = "transcript_id"))

names(all_de_pseudo_fasta) <- paste(
  matched_check_all2$ENST,
  matched_check_all2$gene_id,
  matched_check_all2$gene_name,
  matched_check_all2$gene_type,
  sep = "|"
)
#
outdir <- "./pseudogenes/result/step4_seq/DE_pseudogene_fasta"
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

writeXStringSet(
  all_de_pseudo_fasta,
  filepath = file.path(outdir, "all_DE_pseudogene.transcripts.fasta"),
  format = "fasta"
)
#
my_group <- c("H_eRA", "H_RA", "RA_eRA")

outdir <- "./pseudogenes/result/step4_seq/DE_pseudogene_fasta/by_group"
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

deg_file <- "./pseudogenes/result/step3_deseq2/only_fold_change/all_pseudogene_group.txt"

de_pseudo_all <- read.delim(
  deg_file,
  header = TRUE,
  sep = "\t",
  check.names = FALSE,
  stringsAsFactors = FALSE
)

dim(de_pseudo_all)
head(de_pseudo_all)
colnames(de_pseudo_all)


table(de_pseudo_all$group)
table(de_pseudo_all$significance)
head(de_pseudo_all$pseudogene_id)
# 
# 
my_group <- c("H_eRA", "H_RA", "RA_eRA")

outdir <- "./pseudogenes/result/step4_seq/DE_pseudogene_fasta/by_group"
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

for (now_group in my_group) {
  
  message("Processing group: ", now_group)
  
## 1) Differentially expressed pseudogene gene_ids (ENSG) in the current group
  diff_pseudo_gene <- de_pseudo_all %>%
    dplyr::filter(group == now_group) %>%
    dplyr::pull(pseudogene_id) %>%
    unique()
  
## 2) Convert ENSG IDs to ENST IDs
  diff_pseudo_tx <- pseudo_tx %>%
    dplyr::filter(gene_id %in% diff_pseudo_gene)
  
## If no exact match is found, match after removing the version number
  if (nrow(diff_pseudo_tx) == 0) {
    tmp_gene_nover <- sub("\\..*$", "", diff_pseudo_gene)
    
    diff_pseudo_tx <- pseudo_tx %>%
      dplyr::mutate(gene_id_nover = sub("\\..*$", "", gene_id)) %>%
      dplyr::filter(gene_id_nover %in% tmp_gene_nover)
  }
  
## 3) Match these ENST IDs in the FASTA check table
  matched_check <- check %>%
    dplyr::filter(ENST %in% diff_pseudo_tx$transcript_id)
  
## 4) Retrieve the actual sequence names (full headers) from the FASTA file
  full_genecode_name <- matched_check$full_header %>% as.character()
  
## 5) Extract sequences from the XStringSet
  diff_pseudo_fasta <- tx_genecode[full_genecode_name]
  
## 6) Modify the FASTA headers for easier tracking
  name_map <- diff_pseudo_tx %>%
    dplyr::select(gene_id, transcript_id, gene_name, gene_type) %>%
    dplyr::distinct()
  
  matched_check2 <- matched_check %>%
    dplyr::left_join(name_map, by = c("ENST" = "transcript_id"))
  
  names(diff_pseudo_fasta) <- paste(
    matched_check2$ENST,
    matched_check2$gene_id,
    matched_check2$gene_name,
    matched_check2$gene_type,
    sep = "|"
  )
  
  ## 7) 
  out_fa <- file.path(outdir, paste0(now_group, ".pseudogene.transcripts.fasta"))
  writeXStringSet(diff_pseudo_fasta, filepath = out_fa, format = "fasta")
  
  message("  genes: ", length(diff_pseudo_gene))
  message("  transcripts: ", nrow(diff_pseudo_tx))
  message("  fasta sequences: ", length(diff_pseudo_fasta))
}

#  annotation table
de_pseudo_tx_group_map <- de_pseudo_all %>%
  dplyr::select(
    pseudogene_id,
    pseudogene_id,
    group,
    log2FoldChange,
    padj,
    significance
  ) %>%
  dplyr::left_join(
    pseudo_tx %>%
      dplyr::select(
        gene_id,
        transcript_id,
        gene_name,
        gene_type,
        transcript_type
      ),
    by = c("pseudogene_id" = "gene_id")
  )

readr::write_tsv(
  de_pseudo_tx_group_map,
  "./pseudogenes/result/step4_seq/DE_pseudogene_gene_transcript_group_map.tsv"
)

