#################################################################################
##############################   # Calculate CPM   #######################################
#################################################################################
.libPaths( c( "/mnt/data/hr1912/softs/R/x86_64-pc-linux-gnu-library/4.0" , .libPaths() ) )
library(ggpubr)
# known_circRNA <- read.delim("./circRNA/part2/known_circRNA.txt",header = TRUE)
# unknown_circRNA <- read.delim("./circRNA/part2/unknown_circRNA.txt",header = TRUE)

circFinder_allSRR_comb <- read.csv("./circRNA/venn/venn_need_data1/z_filterBSJ_allSRR/circFinder_allSRR_comb.txt", sep="",header = TRUE)
circFinder_allSRR_comb <- circFinder_allSRR_comb %>% mutate(comb = str_c(Chr,Start,End,chain,sep = "&"),.after = 4)

select_known_circRNA <- read.csv("./circRNA/venn/ResultData/Bed/select_known_circRNA.txt", sep="",header = TRUE)
remaining_unknown_circRNA <- read.csv("./circRNA/venn/ResultData/Bed/remaining_unknown_circRNA.txt", sep="",header = TRUE)

select_known_circRNA <- select_known_circRNA %>% mutate(comb = str_c(Chr,Start,End,chain,sep = "&"))
remaining_unknown_circRNA <- remaining_unknown_circRNA %>% mutate(comb = str_c(Chr,Start,End,chain,sep = "&")

known_circRNA <- circFinder_allSRR_comb[circFinder_allSRR_comb$comb %in% select_known_circRNA$comb,]
unknown_circRNA <- circFinder_allSRR_comb[circFinder_allSRR_comb$comb %in% remaining_unknown_circRNA$comb,]

#  Calculate the number of samples in which each circRNA is detected
known_circRNA[is.na(known_circRNA)] <- 0
unknown_circRNA[is.na(unknown_circRNA)] <- 0
# 
known_PerCirc_zero <- apply(known_circRNA[,6:ncol(known_circRNA)], 1,function(x) sum(x == 0))
known_circRNA <- known_circRNA %>% mutate(zero_count = known_PerCirc_zero,.after = 5)
known_circRNA <- known_circRNA %>% mutate(InSRR_num = 180 - zero_count,.after = 5)

unknown_PerCirc_zero <- apply(unknown_circRNA[,6:ncol(unknown_circRNA)], 1,function(x) sum(x == 0))
unknown_circRNA <- unknown_circRNA %>% mutate(zero_count = unknown_PerCirc_zero,.after = 5)
unknown_circRNA <- unknown_circRNA %>% mutate(InSRR_num = 180 - zero_count,.after = 5)
# Save the known and unknown circRNAs along with the number of samples in which each circRNA is detected
write.table(known_circRNA,"./circRNA/part2/Data/known_circRNA.txt",quote = FALSE,row.names = FALSE,sep = "\t")
write.table(unknown_circRNA,"./circRNA/part2/Data/unknown_circRNA.txt",quote = FALSE,row.names = FALSE,sep = "\t")
#
table(known_circRNA$InSRR_num)
table(unknown_circRNA$InSRR_num)

#
known_circRNA <- read.delim("./circRNA/part2/Data/known_circRNA.txt",header = TRUE)
unknown_circRNA <- read.delim("./circRNA/part2/Data/unknown_circRNA.txt",header = TRUE)

known_circRNA <- known_circRNA %>% mutate(group = "known") %>% select(group,everything())
unknown_circRNA <- unknown_circRNA %>% mutate(group = "unknown") %>% select(group,everything())

ToTable <- rbind(known_circRNA,unknown_circRNA)
ToTable <- ToTable %>% mutate(state1_count = NA,state1_zero = NA,state2_count = NA,state2_zero = NA,state3_count = NA,state3_zero = NA) %>% select(c(1:6,state1_count,state1_zero,state2_count,state2_zero,state3_count,state3_zero,everything()))

#定义3组范围
SampleGroup <- read.delim("./circRNA/a_info/SampleGroup.txt",header = TRUE)
state1 <- ToTable[,colnames(ToTable) %in% SampleGroup[1:28,1]] 
state2 <- ToTable[,colnames(ToTable) %in% SampleGroup[29:85,1]]
state3 <- ToTable[,colnames(ToTable) %in% SampleGroup[86:180,1]]

s1_zero <- apply(state1,1,function(x) sum(x == 0))
ToTable <- ToTable %>% mutate(state1_zero = s1_zero)
ToTable <- ToTable %>% mutate(state1_count = 28 - state1_zero)

s2_zero <- apply(state2,1,function(x) sum(x == 0))
ToTable <- ToTable %>% mutate(state2_zero = s2_zero)
ToTable <- ToTable %>% mutate(state2_count = 57 - state2_zero)

s3_zero <- apply(state3,1,function(x) sum(x ==0))
ToTable <- ToTable %>% mutate(state3_zero = s3_zero)
ToTable <- ToTable %>% mutate(state3_count = 95 - state3_zero)
#
library(scales)
ToTable <- ToTable %>% mutate(state1_percent = state1_count/28,.after = 8) %>% mutate(state1_percent = percent(state1_percent))
ToTable <- ToTable %>% mutate(state2_percent = state2_count/57,.after = 11) %>% mutate(state2_percent = percent(state2_percent))
ToTable <- ToTable %>% mutate(state3_percent = state3_count/95,.after = 14) %>% mutate(state3_percent = percent(state3_percent))
#
write.table(ToTable,"./circRNA/part2/ToTable.txt",quote = FALSE,row.names = FALSE,sep = "\t")

# Filtering criterion: retain circRNAs detected in >10% of samples in at least one of the three groups
# Here, a threshold of 15% was used
test <- ToTable %>%
  filter(
    as.numeric(gsub("%","",state1_percent)) > 15 |
    as.numeric(gsub("%","",state2_percent)) > 15 |
    as.numeric(gsub("%","",state3_percent)) > 15
  )

write.table(test,"./circRNA/part2/filter_ToTable.txt",quote = FALSE,row.names = FALSE,sep = "\t")
# Read the BAM files to obtain the total number of mapped reads
original_bam2 <- read_excel("./circRNA/part2/Data/circRNA_finder_bam/original_bam2.xlsx")

filter_ToTable <- read.delim("./circRNA/part2/filter_ToTable.txt",header = TRUE)
all_reads <- original_bam2$`all reads`

test <- filter_ToTable
df1 <- test[,c(1:17)]
df2 <- test[,c(18:ncol(test))]

df2_normli <- sweep(df2,2,all_reads,FUN = function(x,all_reads) x*1000000/all_reads)

test <- cbind(df1,df2_normli)
#
for_plot <- test[,c(1,6,18:ncol(test))]
write.table(for_plot,"./circRNA/part2/boxplot/for_plot.txt",quote = FALSE,row.names = FALSE,sep = "\t")
# 
for_plot_known <- for_plot %>% filter(group == "known")
for_plot_unknown <- for_plot %>% filter(group =="unknown")
########## 
#
test <- for_plot_known[,3:ncol(for_plot_known)]
mean_known <- apply(test,2,mean)
SampleGroup <- read.delim("./circRNA/a_info/SampleGroup.txt",header = TRUE)
plot_sample_mean <- SampleGroup
plot_sample_mean$Sample_mean <- mean_known
plot_sample_mean$ISknown <- "known"
#
test <- for_plot_unknown[,3:ncol(for_plot_unknown)]
mean_unknown <- apply(test,2,mean)
unknown_plot_sample_mean <- SampleGroup
unknown_plot_sample_mean$Sample_mean <- mean_unknown
unknown_plot_sample_mean$ISknown <- "unknown"

# 
PlotData <- rbind(plot_sample_mean,unknown_plot_sample_mean)
write.table(PlotData,"./circRNA/part2/boxplot/mean/PlotData.txt",quote = FALSE,row.names = FALSE,sep = "\t")

# 
test <- for_plot_known[,3:ncol(for_plot_known)]
sum_known <- apply(test,2,sum)
SampleGroup <- read.delim("./circRNA/a_info/SampleGroup.txt",header = TRUE)
plot_sample_sum <- SampleGroup
#
plot_sample_sum$Sample_sum <- sum_known
plot_sample_sum$ISknown <- "known"

# 
test <- for_plot_unknown[,3:ncol(for_plot_unknown)]
sum_unknown <- apply(test,2,sum)
unknown_plot_sample_sum <- SampleGroup
unknown_plot_sample_sum$Sample_sum <- sum_unknown
unknown_plot_sample_sum$ISknown <- "unknown"

PlotData_sum <- rbind(plot_sample_sum,unknown_plot_sample_sum)
write.table(PlotData_sum,"./circRNA/part2/boxplot/sum/PlotData_sum.txt",quote = FALSE,row.names = FALSE,sep = "\t")

######## 
#
# plotsum1 <- PlotData_sum %>% filter(ISknown =="known")
#
plot_sample_sum <- PlotData_sum %>% dplyr::filter(ISknown == "known")
my_comparisons <- list(c("Normal", "Rheumatoid arthritis(early)"),
                       c("Normal", "Rheumatoid arthritis(established)"),
                       c("Rheumatoid arthritis(early)", "Rheumatoid arthritis(established)"))

tiff("./circRNA/part2/boxplot/all_reads_known_sum.tiff", width = 2800, height = 2400, res = 300)
p <- ggplot(plot_sample_sum,aes(x=Group,y=Sample_sum,fill = Group))+
  theme_classic()+
  #labs(title = "Comparison of CPM Across Disease Types")+
  geom_boxplot(color = "black",size = 0.5)+
  scale_fill_manual(values = brewer.pal(3, "Pastel1"))+
  geom_jitter(position = position_jitter(0.2),color = "black",size = 1)+
  stat_compare_means(comparisons = my_comparisons,method = "wilcox.test",label = "p.signif",label.y = c(500,650,700))
p

dev.off()
#
#plotsum2 <- PlotData_sum %>% filter(ISknown =="unknown")
#
tiff("./circRNA/part2/boxplot/all_reads_unknown_sum.tiff", width = 2800, height = 2400, res = 300)
p <- ggplot(unknown_plot_sample_sum,aes(x=Group,y=Sample_sum,fill = Group))+
  theme_classic()+
  #labs(title = "Comparison of CPM Across Disease Types")+
  geom_boxplot(color = "black",size = 0.5,outlier.color = "red")+
  scale_fill_manual(values = brewer.pal(3, "Pastel2"))+
  geom_jitter(position = position_jitter(0.2),color = "black",size = 1)+
  coord_cartesian(ylim = c(0, 2.5)) +  #  
  stat_compare_means(comparisons = my_comparisons,method = "wilcox.test",label = "p.signif",label.y = c(1.5,2.3,1.8))
p

dev.off()

# # Assess whether the number of unmapped reads is associated with circRNA abundance ----------------------------------------

original_bam2 <- read.delim("./circRNA/part2/lm/original_bam2.txt")
SampleGroup2 <- read.delim("./circRNA/part2/DEseq2/need_data/SampleGroup2.txt",header = TRUE)
over_zero <- apply(filter_ToTable[18:ncol(filter_ToTable)],2,function(x) sum(x>0))
all_state1_BSJs <- colSums(filter_ToTable[,colnames(filter_ToTable) %in% SampleGroup2[SampleGroup2$Group=="healthy",1]])
all_state2_BSJs <- colSums(filter_ToTable[,colnames(filter_ToTable) %in% SampleGroup2[SampleGroup2$Group=="eRA",1]])
all_state3_BSJs <- colSums(filter_ToTable[,colnames(filter_ToTable) %in% SampleGroup2[SampleGroup2$Group=="RA",1]])


original_bam2$filter_circs <- over_zero
original_bam2$all_BSJs <- all_BSJs

write.table(original_bam2,"./circRNA/part2/lm/original_bam2.txt",quote = FALSE,row.names = FALSE,sep = "\t")

# original_bam2$all_state1_BSJs <- all_state1_BSJs
# original_bam2$all_state2_BSJs <- all_state2_BSJs
# original_bam2$all_state3_BSJs <- all_state3_BSJs
# y=BSJ
# x= unmapped reads

original_bam2 <- read.delim("./circRNA/part2/lm/original_bam2.txt",header = TRUE)
filter_ToTable <- read.delim("./circRNA/part2/filter_ToTable.txt",header = TRUE)
#--------------------------
y <- original_bam2$all_BSJs
x <- original_bam2$unmapped_reads
fit <- lm(y ~ x)

#tiff("./circRNA/part2/lm/picture/all.tiff",width = 2800, height = 2400, res = 300)
p <- ggplot(data = data.frame(x = x, y = y, disease = original_bam2$disease), aes(x = x, y = y, color = disease)) +
  geom_point() +  # 
  stat_smooth(method = "lm", col = "blue") +  # 
  stat_cor(method = "spearman") +  # 
  labs(title = "Scatter Plot with Linear Fit",
       x = "Unmap reads",
       y = "Junction reads") +
  theme_minimal()
#dev.off()
p

# ------------------------------------------------------------

pp <- ggplot(data = data.frame(x = x, y = y), aes(x = x, y = y)) +
  geom_point() +  # 
  stat_smooth(method = "lm", col = "blue") +  # 
  stat_cor(method = "spearman", aes(label = paste(..r.label.., ..p.label.., sep = "~`,`~")),
           label.x = 3, label.y = 15) +  # 
  labs(title = "Scatter Plot with Linear Fit",
       x = "Unmap reads",
       y = "Junction reads") +
  theme_minimal()
pp
#dev.off()

pp <- ggplot(data = data.frame(x = x, y = y), aes(x = x, y = y)) +
  geom_point() +  # 
  stat_smooth(method = "lm", col = "blue") +  # 
  stat_cor(method = "spearman", aes(label = paste(..r.label.., ..p.label.., sep = "~`,`~")),
           label.x = 3, label.y = 15) +  # 
       x = "Unmap reads",
       y = "Junction reads") +
  theme_minimal() +
  theme(
    panel.grid = element_blank(),               # 
    axis.line = element_line(color = "black"),  # 
    panel.border = element_rect(color = "black", fill = NA)  # 
  )

pp

pp <- ggplot(data = data.frame(x = x, y = y), aes(x = x, y = y)) +
  geom_point() +  # 
  stat_smooth(method = "lm", col = "blue") +  # 
  stat_cor(method = "spearman", aes(label = paste(..r.label.., ..p.label.., sep = "~`,`~")),
           label.x = 3, label.y = -15) +  # 
  labs(title = "Scatter Plot with Linear Fit",
       x = "Unmap reads",
       y = "Junction reads") +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 22),  # 
    axis.title = element_text(size = 16),               #   
    panel.grid = element_blank(),                       # 
    #axis.line = element_line(color = "black"),          # 
    panel.border = element_rect(color = "black", fill = NA,linewidth = 1)  # 
  )

pp


plot_grid(p,p1,p2,p3)
plot_grid(pp,NULL,p,NULL,
          nrow = 2)


library(ggpubr)

# 
my_comparisons <- list(c("Normal", "Rheumatoid arthritis(early)"),
                       c("Normal", "Rheumatoid arthritis(established)"),
                       c("Rheumatoid arthritis(early)", "Rheumatoid arthritis(established)"))

tiff("./circRNA/part2/boxplot/output_image.tiff", width = 2800, height = 2600, res = 300)
ggplot(plot_long2, aes(state, mean_CPM, fill = state)) +
  geom_boxplot() +
  geom_jitter(width = 0.2, alpha = 0.5) +  # 
  scale_fill_npg() +
  theme_classic() +
  facet_wrap(~ group) +  # 
  labs(x = 'group_state', y = 'CPM_count', 
       title = "Comparison of CPM Across Disease Types") +
  theme(text = element_text(size = 20)) +
  stat_compare_means(comparisons = my_comparisons, method = "wilcox.test", 
                     label = "p.signif", label.y = c(0, 0, 25))  # 


dev.off()

# 
# 
my_comparisons <- list(c("Normal", "Rheumatoid arthritis(early)"),
                       c("Normal", "Rheumatoid arthritis(established)"),
                       c("Rheumatoid arthritis(early)", "Rheumatoid arthritis(established)"))
                       
tiff("./circRNA/part2/boxplot/output_image.tiff", width = 2800, height = 2600, res = 300)
 ggplot(plot_long, aes(state, CPM, fill = state)) +
  geom_boxplot() +
  geom_jitter(width = 0.2, alpha = 0.5) +  # 
  scale_fill_npg() +
  theme_classic() +
  facet_wrap(~ group) +  # 
  labs(x = 'group_state', y = 'CPM_count', 
       title = "Comparison of CPM Across Disease Types") +
  theme(text = element_text(size = 20)) +
  stat_compare_means(comparisons = my_comparisons, method = "wilcox.test", 
                     label = "p.signif", label.y = c(0, 0, 25))  # 
dev.off()

###################
library(ggpubr)
my_comparisons <- list(c("Normal", "Rheumatoid arthritis(early)"),
                       c("Normal", "Rheumatoid arthritis(established)"),
                       c("Rheumatoid arthritis(early)", "Rheumatoid arthritis(established)"))
#   
tiff("./circRNA/part2/boxplot/output_image.tiff", width = 2800, height = 2600, res = 300)

# 
p <- ggplot(plot_long, aes(x = state, y = CPM, fill = state)) +
  geom_violin(trim = FALSE) +
  geom_boxplot(width = 0.1, position = position_dodge(0.9)) +  # 
  labs(x = "circRNA Group", y = "CPM") +
  theme_minimal() +
  facet_wrap(~ group)  # 

# 
p + stat_compare_means(aes(label = ..p.signif..), method = "wilcox.test")
# 
dev.off()
# 
circFinder_allSRR_comb[is.na(circFinder_allSRR_comb)] <- 0
to_use <- rbind(known_circRNA,unknown_circRNA)
# write.table(to_use,"./circRNA/part2/AfterTools_CPM.txt",quote = FALSE,row.names = FALSE,sep = "\t")
to_use[is.na(to_use)] <- 0
#  
CF_SRR_BSJ <- colSums(to_use[,6:ncol(to_use)])
#------------------------------------------------------------------------------------------------------
known_circRNA[is.na(known_circRNA)] <- 0
# 
df1 <- known_circRNA[,6:ncol(known_circRNA)]
#   
df1_normalized <- sweep(df1, 2, CF_SRR_BSJ, FUN = function(x, CF_SRR_BSJ) x * 1000000 / CF_SRR_BSJ)
CPM_known <- known_circRNA[,1:5]
CPM_known <- cbind(CPM_known,df1_normalized)

#------------------------------------------------------------------
unknown_circRNA[is.na(unknown_circRNA)] <- 0
# 
df2 <- unknown_circRNA[,6:ncol(unknown_circRNA)]
# #
# col_sums2 <- colSums(df2) 
df2_normalized <- sweep(df2,2, CF_SRR_BSJ, FUN = function(x, CF_SRR_BSJ) x * 1000000 / CF_SRR_BSJ)

CPM_unknown <- unknown_circRNA[,1:5]
CPM_unknown <- cbind(CPM_unknown,df2_normalized)

#--------------------------------
write.table(CPM_known,"./circRNA/part2/CPM_known.txt",quote = FALSE,row.names = FALSE,sep = "\t")
write.table(CPM_unknown,"./circRNA/part2/CPM_unknown.txt",quote = FALSE,row.names = FALSE,sep = "\t")

#  paper lm -------------------------------------------------------
original_bam2 <- read.delim("./circRNA/part2/lm/original_bam2.txt",header = TRUE)
#filter_ToTable <- read.delim("./circRNA/part2/filter_ToTable.txt",header = TRUE)
original_bam2 <- original_bam2 %>% 
  dplyr::mutate(new_Unmap_reads = unmapped_reads/1e7) %>% 
  dplyr::mutate(newe_BSJ = all_BSJs/1e4)
  
#--------------------------
y <- original_bam2$newe_BSJ
x <- original_bam2$new_Unmap_reads
fit <- lm(y ~ x)
write.table(original_bam2,"./circRNA/part2/lm/paper/paper_lm.txt",quote = FALSE,row.names = FALSE,sep = "\t")
#
pdf("./circRNA/part2/lm/paper/paper_lm.pdf",width = 7,height = 6)
pp <- ggplot(data = data.frame(x = x, y = y), aes(x = x, y = y)) +
  geom_point(size = 3) +  # 
  stat_smooth(method = "lm", col = "blue") +  # 
  stat_cor(method = "spearman", aes(label = paste(..r.label.., ..p.label.., sep = "~`,`~")),
           label.x = 6, label.y = 5, # 
          size = 6) +  # 
  labs(x = "Unmap reads (10⁷)",
       y = "Junction reads (10⁴)") +  = "Scatter Plot with Linear Fit",
  theme_minimal() +
  theme(
    #plot.title = element_text(hjust = 0.5, size = 22),  # 
    axis.title = element_text(size = 20),# 
    axis.text = element_text(size = 15),
    panel.grid = element_blank(),                   # 
    #axis.line = element_line(color = "black"),          # 
    panel.border = element_rect(color = "black", fill = NA,linewidth = 1)) + # 
  coord_cartesian(ylim = c(0, 6)) 

pp

dev.off()





















