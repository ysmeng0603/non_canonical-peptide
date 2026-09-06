# Load packages ----------------------------------------------------------------
library(dplyr)
library(tidyverse)
library(tidyr)
library(DESeq2)

####################################################################################################
# Read the count matrix --------------------------------------------------------
# Use the transcript-level count matrix
lncRNA_raw_counts <- read.delim("./lncRNA_known/StringTie/result/step_7_featureCounts/lncRNA_transcript_counts.txt", comment.char="#")
#counts_summary <- read.delim("./lncRNA_known/StringTie/result/step_7_featureCounts/lncRNA_raw_counts.txt.summary")
#
summary(colSums(lncRNA_raw_counts[,7:ncol(lncRNA_raw_counts)]))
# Extract SRR IDs
new_colnames <- colnames(lncRNA_raw_counts)
new_colnames

new_colnames[7:length(new_colnames)] <- gsub(
  ".*(SRR[0-9]+).*",
  "\\1",
  new_colnames[7:length(new_colnames)]
)
new_colnames

colnames(lncRNA_raw_counts) <- new_colnames
colnames(lncRNA_raw_counts)

# Differential expression analysis --------------------------------------------
# Construct the matrix
count_matrix <- lncRNA_raw_counts[, 7:ncol(lncRNA_raw_counts)]  # 
rownames(count_matrix) <- lncRNA_raw_counts$Geneid
# Prepare group information
SampleGroup2 <- read.delim("~/A_info/SampleGroup2.txt") # 
SampleGroup2 <- SampleGroup2 %>% 
  dplyr::filter(Sample %in% colnames(count_matrix)) # 
# 
SampleGroup2$Group <- as.factor(SampleGroup2$Group)
all(colnames(count_matrix) == SampleGroup2$Sample) # 
# 
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
res_df$lncRNA_id <- rownames(res_df) # 


res_df <- na.omit(res_df) # 
#
res_df$significance <- ifelse(res_df$padj < 0.05 & abs(res_df$log2FoldChange) > 1, 
                              ifelse(res_df$log2FoldChange > 0, "Up", "Down"), 
                              "Stable")
write.table(res_df,"./lncRNA_known/StringTie/result/step_8_DEseq2/lncRNA_diff_H_eRA.txt",quote = FALSE,sep = "\t")
# Hc-RA -------------------------------------------------------------------
res <- results(
  dds,
  contrast = c("Group", "RA", "healthy"),
  independentFiltering = FALSE   # 
)

res <- res[order(res$padj), ]
summary(res)
res_df <- as.data.frame(res)
res_df$lncRNA_id <- rownames(res_df) # 
res_df <- na.omit(res_df) # 
# 
res_df$significance <- ifelse(res_df$padj < 0.05 & abs(res_df$log2FoldChange) > 1, 
                              ifelse(res_df$log2FoldChange > 0, "Up", "Down"), 
                              "Stable")
write.table(res_df,"./lncRNA_known/StringTie/result/step_8_DEseq2/lncRNA_diff_H_RA.txt",quote = FALSE,sep = "\t")
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
res_df$lncRNA_id <- rownames(res_df) # 


res_df <- na.omit(res_df) # 
# 
res_df$significance <- ifelse(res_df$padj < 0.05 & abs(res_df$log2FoldChange) > 1, 
                              ifelse(res_df$log2FoldChange > 0, "Up", "Down"), 
                              "Stable")
write.table(res_df,"./lncRNA_known/StringTie/result/step_8_DEseq2/lncRNA_diff_RA_eRA.txt",quote = FALSE,sep = "\t")
#####
#  ------------------------------------------------------------------
setwd("./lncRNA_known/StringTie/result/step_8_DEseq2")
#
lncRNA_diff_H_eRA <- read.delim("./lncRNA_known/StringTie/result/step_8_DEseq2/lncRNA_diff_H_eRA.txt")
lncRNA_diff_H_RA <- read.delim("./lncRNA_known/StringTie/result/step_8_DEseq2/lncRNA_diff_H_RA.txt")
lncRNA_diff_RA_eRA <- read.delim("./lncRNA_known/StringTie/result/step_8_DEseq2/lncRNA_diff_RA_eRA.txt")
#
head(lncRNA_diff_H_eRA)
table(lncRNA_diff_H_eRA$significance)
# Down Stable     Up 
# 4501  10286   1619  
table(lncRNA_diff_H_RA$significance)
# Down Stable     Up 
# 3502  10642   2262  

#
#pdf("./lncRNA_known/StringTie/result/step_8_DEseq2/plot/lnc_H_eRA.pdf")
ggplot(lncRNA_diff_H_eRA, aes(x = log2FoldChange, y = -log10(padj), color = significance)) +
  geom_point(alpha = 0.6, size = 1.5) +
  scale_color_manual(values = c("Stable" = "gray", 
                                "Up" = "red", 
                                "Down" = "blue")) +  # 
  labs(title = "healthy-vs-eRA volcano plot",
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
#pdf("./lncRNA_known/StringTie/result/step_8_DEseq2/plot/lnc_H_RA.pdf")
ggplot(lncRNA_diff_H_RA, aes(x = log2FoldChange, y = -log10(padj), color = significance)) +
  geom_point(alpha = 0.6, size = 1.5) +
  scale_color_manual(values = c("Stable" = "gray", 
                                "Up" = "red", 
                                "Down" = "blue")) +  # 
  labs(title = "healthy-vs-RA volcano plot",
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
#pdf("./lncRNA_known/StringTie/result/step_8_DEseq2/plot/lnc_RA_eRA.pdf")
ggplot(lncRNA_diff_RA_eRA, aes(x = log2FoldChange, y = -log10(padj), color = significance)) +
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
need_1 <- lncRNA_diff_H_eRA %>% 
  dplyr::filter(significance !="Stable") %>% 
  dplyr::mutate(group = "H_eRA") %>% 
  dplyr::select(group,everything()) # 6120
#write.table(need_1,"./lncRNA_known/StringTie/result/step_8_DEseq2/only_fold_change/HC_eRA.txt",quote = FALSE,sep = "\t")

need_2 <- lncRNA_diff_H_RA %>% 
  dplyr::filter(significance !="Stable") %>% 
  dplyr::mutate(group = "H_RA") %>% 
  dplyr::select(group,everything()) # 5764
#write.table(need_2,"./lncRNA_known/StringTie/result/step_8_DEseq2/only_fold_change/HC_RA.txt",quote = FALSE,sep = "\t")

need_3 <- lncRNA_diff_RA_eRA %>% 
  dplyr::filter(significance !="Stable") %>% 
  dplyr::mutate(group = "RA_eRA") %>% 
  dplyr::select(group,everything()) # 2016
#write.table(need_3,"./lncRNA_known/StringTie/result/step_8_DEseq2/only_fold_change/RA_eRA.txt",quote = FALSE,sep = "\t")

all <- rbind(need_1,need_2,need_3)
all_2 <- all %>% unique()
#write.table(all,"./lncRNA_known/StringTie/result/step_8_DEseq2/only_fold_change/all_lncRNA_group.txt",quote = FALSE,sep = "\t")

# all <- all %>% 
#   dplyr::select(group,lncRNA_id) %>% unique() # 

length(unique(all$lncRNA_id)) # 
table(all_2$significance)

# TPM ------------------------------------------------------------------

# Calculate TPM for absolute lncRNA expression ----------------------------
library(tidyverse)
library(org.Hs.eg.db)
library(AnnotationDbi)
library(data.table)
fc <- fread("./lncRNA_known/StringTie/result/step_7_featureCounts/lncRNA_transcript_counts.txt")
# 
head(fc[, 1:8])
# Step 2：
# raw count matrix
count_mat <- as.matrix(fc[, 7:ncol(fc)])
rownames(count_mat) <- fc$Geneid

# transcript length (bp)
tx_len <- fc$Length
names(tx_len) <- fc$Geneid
# Step 3：TPM
calcTPM <- function(counts, lengths) {
  lengths_kb <- lengths / 1000
  rpk <- counts / lengths_kb  # 
  scaling_factor <- colSums(rpk) # 
  tpm <- sweep(rpk, 2, scaling_factor, "/") * 1e6 # 
  return(tpm)
}
# Step 4：
tpm_mat <- calcTPM(count_mat, tx_len)

# 
# Step 5：
summary(colSums(tpm_mat))
colnames(tpm_mat)
colnames(count_matrix)
colnames(tpm_mat) <- str_extract(colnames(tpm_mat), "SRR\\d+")

write.table(
  tpm_mat,
  file = "./lncRNA_known/lncRNA_TPM_matrix.txt",
  sep = "\t",
  quote = FALSE
)
# 
# Step 1：
tpm_df <- as.data.frame(tpm_mat)
tpm_df$lncRNA_id <- rownames(tpm_df)
#
SampleGroup2 <- read.delim("~/A_info/SampleGroup2.txt")
#
tpm_long <- tpm_df %>%
  pivot_longer(
    cols = -lncRNA_id,
    names_to = "Sample",
    values_to = "TPM"
  ) %>%
  left_join(SampleGroup2, by = "Sample")
# 
range(tpm_long$TPM)
# test <- tpm_long %>% 
#   dplyr::filter(TPM ==0)

# Step 2：log2(TPM + 1) 0 
tpm_long$logTPM <- log2(tpm_long$TPM + 1)
# Step 3：
library(RColorBrewer)
# 
# 
my_comparisons <- list(c("healthy", "eRA"),
                       c("healthy", "RA"),
                       c("eRA", "RA"))
#
tpm_long$Group <- factor(x = tpm_long$Group,levels = c("healthy","eRA","RA"))

sample_level <- tpm_long %>%
  group_by(Sample, Group) %>%
  summarise(
    #median_logTPM = median(logTPM, na.rm = TRUE),
    ave_logTPM = mean(logTPM),
    .groups = "drop"
  )

# 
#my_colors <- c("#add4a5","#faca88","#c9b8da")
my_colors <- c("#8A919799","#FED43999", "#709AE199")
#绘图
#pdf("./part2/boxplot/paper/papaer_CPM_box.pdf",width = 8,height = 6)
p <- ggplot(sample_level,aes(x=Group,y=ave_logTPM,fill = Group))+
  theme_classic()+
  #labs(title = "Comparison of CPM Across Disease Types")+
  geom_boxplot(color = "black",
               size = 0.75,# 
               fatten = 1 )+ # 
  scale_fill_manual(values = brewer.pal(3, "Pastel2"))+
  geom_jitter(position = position_jitter(0.2),color = "black",size = 1)+
  coord_cartesian(ylim = c(0, 4))+
  ylab("mean log2(TPM + 1)")+# 
  # 
  stat_compare_means(comparisons = my_comparisons,method = "wilcox.test",label = "p.signif",label.y = c(0.5,0.8,1.1))+
  scale_fill_manual(values = my_colors)+
  theme(
    #plot.title = element_text(hjust = 0.5, size = 22),  # 
    axis.title = element_text(size = 20), # 
    axis.text = element_text(size = 15),   # 
    panel.grid = element_blank(),                       # 
    #axis.line = element_line(color = "black"),  # 
    axis.line = element_blank(), # 
    panel.border = element_rect(color = "black", fill = NA,linewidth = 1),
    # ===== 图例 =====
    legend.title = element_text(size = 15),     
    legend.text = element_text(size = 13),      
    #legend.key.size = unit(0.4, "cm"),        
    # 
  )
p

dev.off()

colSums(featurecount[,7:ncol(featurecount)]) %>% as.data.frame() %>% View()

# Extract sequences -----------------------------------------------------------
# Step 10: Extract sequences of differentially expressed lncRNAs ----------------
lncRNA_genecode <- readDNAStringSet("./reference/gencode.v49.lncRNA_transcripts.fa")
lncRNA_genecode
check <- names(lncRNA_genecode) %>% as.data.frame() 
check$ENST <- sub("\\|.*$", "", check$.)
table(all$lncRNA_id %in% check$ENST) # all T 
# 
my_group <- c("H_eRA","H_RA","RA_eRA")
outdir <- "./lncRNA_known/StringTie/result/step_10_FC_fasta/"
# 
for (now_group in my_group) {
  ## 1) 
  diff_lncRNA <- all %>%
    dplyr::filter(group == now_group) %>%
    dplyr::pull(lncRNA_id) %>%
    unique()
  
  ## 2) 
  matched_check <- check %>%
    dplyr::filter(ENST %in% diff_lncRNA)
  
  ## 3) 
  full_genecode_name <- matched_check[[1]] %>% as.character()
  
  # ： matched_check$full_header
  ## 4) # Extract sequences from the XStringSet
  diff_lncRNA_fasta <- lncRNA_genecode[full_genecode_name]
  
  ## 5) # Export FASTA
  out_fa <- paste0(outdir, now_group, ".fasta")
  writeXStringSet(diff_lncRNA_fasta, filepath = out_fa, format = "fasta")
}

# 
test <- readDNAStringSet("./lncRNA_known/StringTie/result/step_10_FC_fasta/RA_eRA.fasta")
test

