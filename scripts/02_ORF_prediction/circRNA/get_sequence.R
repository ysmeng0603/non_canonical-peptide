# 用来取出circRNA序列

# 准备文件 --------------------------------------------------------------------
# 导入circbank数据库序列
# !!!! 这里必须用 read DNA ,否则读取会出现报错，取出的序列不全，
circbank_sequence <- readDNAStringSet("/mnt/data/ysmeng/circRNA/part2/ceRNA/circRNA/sequence_file/database_sequence/hg19_circbase_seq.fa")

unique_all_database_info <- read.delim("/mnt/data/ysmeng/circRNA/a_info/unique_all_database_info.txt",header = TRUE)
circAtlas_human_sequence_v3 <- read.table("/mnt/data/ysmeng/circRNA_peptide2/GetSequence/Three_databases/circAtlas_human_sequence_v3", quote="\"", comment.char="",header = FALSE)
colnames(circAtlas_human_sequence_v3) <- c("DB_name","sequence")
circbank_circbase <- read.delim("/mnt/data/ysmeng/circRNA/venn/venn_need_data2/three_database/circbank/circbank_circbase.txt")
# 先做H_eRA -----------------------------------------------------------------
# 该表用来筛选出差异且对应host gene不差异的circRNA
linear_circ_H_eRA <- read.delim("/mnt/data/ysmeng/circRNA/2025_DEseq2/z_HostGene/linear_circ_H_eRA.txt",header = TRUE)
test <- linear_circ_H_eRA %>% dplyr::filter(mRNA_significance == "Stable" & (significance %in% c("Up","Down"))) # 875

table(test$significance)# down 68,up 807
table(test$mRNA_significance) # 875

# 1.1 H_eRA_down ----------------------------------------------------------
H_eRA_down <- test %>% dplyr::filter(significance == "Down")
table(H_eRA_down$circRNA %in% circbank_circbase$circRNA) # F:1,TRUE:67
H_eRA_down <- H_eRA_down %>% dplyr::select(circRNA,significance)
H_eRA_down <- left_join(H_eRA_down,circbank_circbase,by = "circRNA")
# ori_H_eRA_down <- read.delim("/mnt/data/ysmeng/circRNA_peptide2/part1_coding_prob/CPAT/H_eRA_down.txt")

# 这组中跟原来的有一个不一样
table(ori_H_eRA_down$circRNA %in% H_eRA_down$circRNA)
table(H_eRA_down$circRNA %in% ori_H_eRA_down$circRNA)
# 
# 先查看circbase名称是否是字符
for (i in 1:nrow(H_eRA_down)) {
  if (!is.na(H_eRA_down[i, 4])) {
    circRNA_circbaseID  <- H_eRA_down[i,4]
    match_sequence <- circbank_sequence[grepl(circRNA_circbaseID,names(circbank_sequence))]
    writeXStringSet(match_sequence,filepath = paste0("/mnt/data/ysmeng/circRNA_peptide3/get_sequence/sequence_file/H_eRA_down/fasta/",
                                                     circRNA_circbaseID,".fasta"),format = "fasta")
    H_eRA_down[i,9] <- "done"
  } else {
    H_eRA_down[i,9] <- "not_found"
    
  }
  write.table(H_eRA_down,"/mnt/data/ysmeng/circRNA_peptide3/A_info/H_eRA_down.txt",quote = FALSE,row.names = FALSE,sep = "\t")
}
# 从circAtlas3中取出序列
colnames(H_eRA_down)[9] <- "state"
test <- H_eRA_down %>% filter(state == "not_found") %>% dplyr::select(circRNA,state)
test <- left_join(test,unique_all_database_info,by = "circRNA") 
test <- test %>% dplyr::filter(Database =="circAtlas3")
out <- circAtlas_human_sequence_v3 %>% dplyr::filter(DB_name == test[1,4])

output_file <- "/mnt/data/ysmeng/circRNA_peptide3/get_sequence/sequence_file/H_eRA_down/supple/hsa-GARNL3_0003.fasta"
fileConn <- file(output_file, "w")
writeLines(paste0(">",out[1,1]),fileConn)
writeLines(out[1,2],fileConn)
close(fileConn)

# 1.2 H_eRA_up ----------------------------------------------------------
test <- linear_circ_H_eRA %>% dplyr::filter(mRNA_significance == "Stable" & (significance %in% c("Up","Down"))) # 875
H_eRA_up <- test %>% dplyr::filter(significance == "Up")
length(unique(H_eRA_up$circRNA))

table(H_eRA_up$circRNA %in% circbank_circbase$circRNA) # F:10,TRUE:797
H_eRA_up <- H_eRA_up %>% dplyr::select(circRNA,significance)
H_eRA_up <- left_join(H_eRA_up,circbank_circbase,by = "circRNA",relationship = "many-to-many")
# ori_H_eRA_up <- read.delim("/mnt/data/ysmeng/circRNA_peptide2/part1_coding_prob/CPAT/H_eRA_up.txt")

# 
# 先查看circbase名称是否是字符
for (i in 1:nrow(H_eRA_up)) {
  if (!is.na(H_eRA_up[i, 4])) {
    circRNA_circbaseID  <- H_eRA_up[i,4]
    match_sequence <- circbank_sequence[grepl(circRNA_circbaseID,names(circbank_sequence))]
    writeXStringSet(match_sequence,filepath = paste0("/mnt/data/ysmeng/circRNA_peptide3/get_sequence/sequence_file/H_eRA_up/fasta/",
                                                     circRNA_circbaseID,".fasta"),format = "fasta")
    H_eRA_up[i,9] <- "done"
  } else {
    H_eRA_up[i,9] <- "not_found"
    
  }
  write.table(H_eRA_up,"/mnt/data/ysmeng/circRNA_peptide3/A_info/H_eRA_up.txt",quote = FALSE,row.names = FALSE,sep = "\t")
}
# 从circAtlas3中取出序列
colnames(H_eRA_up)[9] <- "state"
test <- H_eRA_up %>% filter(state == "not_found") %>% dplyr::select(circRNA,state)
test <- left_join(test,unique_all_database_info,by = "circRNA") 
#test <- test %>% dplyr::filter(Database =="circAtlas3")

# out <- circAtlas_human_sequence_v3 %>% dplyr::filter(DB_name == test[1,4])
# 
# output_file <- "/mnt/data/ysmeng/circRNA_peptide3/get_sequence/sequence_file/H_eRA_up/supple/hsa-GARNL3_0003.fasta"
# fileConn <- file(output_file, "w")
# writeLines(paste0(">",out[1,1]),fileConn)
# writeLines(out[1,2],fileConn)
# close(fileConn)
# 

# 手动下载剩余的两条补充
# chr3&121635734&121637256&-  HSA_CIRCpedia_43147

circbank_sequence <- readDNAStringSet("/mnt/data/ysmeng/circRNA/part2/ceRNA/circRNA/sequence_file/database_sequence/hg19_circbase_seq.fa")
unique_all_database_info <- read.delim("/mnt/data/ysmeng/circRNA/a_info/unique_all_database_info.txt",header = TRUE)
circAtlas_human_sequence_v3 <- read.table("/mnt/data/ysmeng/circRNA_peptide2/GetSequence/Three_databases/circAtlas_human_sequence_v3", quote="\"", comment.char="",header = FALSE)
colnames(circAtlas_human_sequence_v3) <- c("DB_name","sequence")
circbank_circbase <- read.delim("/mnt/data/ysmeng/circRNA/venn/venn_need_data2/three_database/circbank/circbank_circbase.txt")
# 先做H_RA -----------------------------------------------------------------
# 该表用来筛选出差异且对应host gene不差异的circRNA
linear_circ_H_RA <- read.delim("/mnt/data/ysmeng/circRNA/2025_DEseq2/z_HostGene/linear_circ_H_RA.txt",header = TRUE)
test <- linear_circ_H_RA %>% dplyr::filter(mRNA_significance == "Stable" & (significance %in% c("Up","Down"))) # 1608

table(test$significance)# down 37,up 1571
table(test$mRNA_significance) # 1608

# 2.1 H_RA_down ----------------------------------------------------------
H_RA_down <- test %>% dplyr::filter(significance == "Down")
table(H_RA_down$circRNA %in% circbank_circbase$circRNA) # F:1,TRUE:36
H_RA_down <- H_RA_down %>% dplyr::select(circRNA,significance)
H_RA_down <- left_join(H_RA_down,circbank_circbase,by = "circRNA")


# 先查看circbase名称是否是字符
for (i in 1:nrow(H_RA_down)) {
  if (!is.na(H_RA_down[i, 4])) {
    circRNA_circbaseID  <- H_RA_down[i,4]
    match_sequence <- circbank_sequence[grepl(circRNA_circbaseID,names(circbank_sequence))]
    writeXStringSet(match_sequence,filepath = paste0("/mnt/data/ysmeng/circRNA_peptide3/get_sequence/sequence_file/H_RA_down/fasta/",
                                                     circRNA_circbaseID,".fasta"),format = "fasta")
    H_RA_down[i,9] <- "done"
  } else {
    H_RA_down[i,9] <- "not_found"
    
  }
  write.table(H_RA_down,"/mnt/data/ysmeng/circRNA_peptide3/A_info/H_RA_down.txt",quote = FALSE,row.names = FALSE,sep = "\t")
}
# 从circAtlas3中取出序列
colnames(H_RA_down)[9] <- "state"
test <- H_RA_down %>% filter(state == "not_found") %>% dplyr::select(circRNA,state)
test <- left_join(test,unique_all_database_info,by = "circRNA") 
test <- test %>% dplyr::filter(Database =="circAtlas3")
out <- circAtlas_human_sequence_v3 %>% dplyr::filter(DB_name == test[1,4])

output_file <- "/mnt/data/ysmeng/circRNA_peptide3/get_sequence/sequence_file/H_RA_down/supple/hsa-GARNL3_0003.fasta"
fileConn <- file(output_file, "w")
writeLines(paste0(">",out[1,1]),fileConn)
writeLines(out[1,2],fileConn)
close(fileConn)

# 2.2 H_RA_up ----------------------------------------------------------
test <- linear_circ_H_RA %>% dplyr::filter(mRNA_significance == "Stable" & (significance %in% c("Up","Down"))) # 875
H_RA_up <- test %>% dplyr::filter(significance == "Up")
length(unique(H_RA_up$circRNA)) #1568

table(H_RA_up$circRNA %in% circbank_circbase$circRNA) # F:30,TRUE:1541
H_RA_up <- H_RA_up %>% dplyr::select(circRNA,significance)
H_RA_up <- left_join(H_RA_up,circbank_circbase,by = "circRNA",relationship = "many-to-many")
# ori_H_RA_up <- read.delim("/mnt/data/ysmeng/circRNA_peptide2/part1_coding_prob/CPAT/H_RA_up.txt")

# 
# 先查看circbase名称是否是字符
for (i in 1:nrow(H_RA_up)) {
  if (!is.na(H_RA_up[i, 4])) {
    circRNA_circbaseID  <- H_RA_up[i,4]
    match_sequence <- circbank_sequence[grepl(circRNA_circbaseID,names(circbank_sequence))]
    writeXStringSet(match_sequence,filepath = paste0("/mnt/data/ysmeng/circRNA_peptide3/get_sequence/sequence_file/H_RA_up/fasta/",
                                                     circRNA_circbaseID,".fasta"),format = "fasta")
    H_RA_up[i,9] <- "done"
  } else {
    H_RA_up[i,9] <- "not_found"
    
  }
  write.table(H_RA_up,"/mnt/data/ysmeng/circRNA_peptide3/A_info/H_RA_up.txt",quote = FALSE,row.names = FALSE,sep = "\t")
}
# 从circAtlas3中取出序列
colnames(H_RA_up)[9] <- "state"
test <- H_RA_up %>% filter(state == "not_found") %>% dplyr::select(circRNA,state)
test <- left_join(test,unique_all_database_info,by = "circRNA") 
#test <- test %>% dplyr::filter(Database =="circAtlas3")

length(unique(test$circRNA))#30

# out <- circAtlas_human_sequence_v3 %>% dplyr::filter(DB_name == test[1,4])
# 
# output_file <- "/mnt/data/ysmeng/circRNA_peptide3/get_sequence/sequence_file/H_RA_up/supple/hsa-GARNL3_0003.fasta"
# fileConn <- file(output_file, "w")
# writeLines(paste0(">",out[1,1]),fileConn)
# writeLines(out[1,2],fileConn)
# close(fileConn)
# 

# 手动下载剩余的两条补充
# chr3&121635734&121637256&-  HSA_CIRCpedia_43147


# 手动补
#
out <- circAtlas_human_sequence_v3 %>% dplyr::filter(DB_name == "hsa-ANKRD26_0004")
output_file <- "/mnt/data/ysmeng/circRNA_peptide3/get_sequence/sequence_file/H_RA_up/supple/supple2/hsa-ANKRD26_0004.fasta"
fileConn <- file(output_file, "w")
writeLines(paste0(">",out[1,1]),fileConn)
writeLines(out[1,2],fileConn)
close(fileConn)

#
out <- circAtlas_human_sequence_v3 %>% dplyr::filter(DB_name == "hsa-BLOC1S5-TXNDC5_0001")
output_file <- "/mnt/data/ysmeng/circRNA_peptide3/get_sequence/sequence_file/H_RA_up/supple/supple2/hsa-BLOC1S5-TXNDC5_0001.fasta"
fileConn <- file(output_file, "w")
writeLines(paste0(">",out[1,1]),fileConn)
writeLines(out[1,2],fileConn)
close(fileConn)

#
out <- circAtlas_human_sequence_v3 %>% dplyr::filter(DB_name == "hsa-CENPH_0002")
output_file <- "/mnt/data/ysmeng/circRNA_peptide3/get_sequence/sequence_file/H_RA_up/supple/supple2/hsa-CENPH_0002.fasta"
fileConn <- file(output_file, "w")
writeLines(paste0(">",out[1,1]),fileConn)
writeLines(out[1,2],fileConn)
close(fileConn)

#
out <- circAtlas_human_sequence_v3 %>% dplyr::filter(DB_name == "hsa-HCLS1_0002")
output_file <- "/mnt/data/ysmeng/circRNA_peptide3/get_sequence/sequence_file/H_RA_up/supple/supple2/hsa-HCLS1_0002.fasta"
fileConn <- file(output_file, "w")
writeLines(paste0(">",out[1,1]),fileConn)
writeLines(out[1,2],fileConn)
close(fileConn)


#
out <- circAtlas_human_sequence_v3 %>% dplyr::filter(DB_name == "hsa-PARP4_0006")
output_file <- "/mnt/data/ysmeng/circRNA_peptide3/get_sequence/sequence_file/H_RA_up/supple/supple2/hsa-PARP4_0006.fasta"
fileConn <- file(output_file, "w")
writeLines(paste0(">",out[1,1]),fileConn)
writeLines(out[1,2],fileConn)
close(fileConn)
#
out <- circAtlas_human_sequence_v3 %>% dplyr::filter(DB_name == "hsa-PPM1B_0010")
output_file <- "/mnt/data/ysmeng/circRNA_peptide3/get_sequence/sequence_file/H_RA_up/supple/supple2/hsa-PPM1B_0010.fasta"
fileConn <- file(output_file, "w")
writeLines(paste0(">",out[1,1]),fileConn)
writeLines(out[1,2],fileConn)
close(fileConn)

#

out <- circAtlas_human_sequence_v3 %>% dplyr::filter(DB_name == "hsa-PSMC6_0006")
output_file <- "/mnt/data/ysmeng/circRNA_peptide3/get_sequence/sequence_file/H_RA_up/supple/supple2/hsa-PSMC6_0006.fasta"
fileConn <- file(output_file, "w")
writeLines(paste0(">",out[1,1]),fileConn)
writeLines(out[1,2],fileConn)
close(fileConn)

# 3.1 RA_eRA_down ----------------------------------------------------------
# 导入circbank数据库序列
# !!!! 这里必须用 read DNA ,否则读取会出现报错，取出的序列不全，
circbank_sequence <- readDNAStringSet("/mnt/data/ysmeng/circRNA/part2/ceRNA/circRNA/sequence_file/database_sequence/hg19_circbase_seq.fa")
unique_all_database_info <- read.delim("/mnt/data/ysmeng/circRNA/a_info/unique_all_database_info.txt",header = TRUE)
circAtlas_human_sequence_v3 <- read.table("/mnt/data/ysmeng/circRNA_peptide2/GetSequence/Three_databases/circAtlas_human_sequence_v3", quote="\"", comment.char="",header = FALSE)
colnames(circAtlas_human_sequence_v3) <- c("DB_name","sequence")
circbank_circbase <- read.delim("/mnt/data/ysmeng/circRNA/venn/venn_need_data2/three_database/circbank/circbank_circbase.txt")
# 先做H_eRA -----------------------------------------------------------------
# 该表用来筛选出差异且对应host gene不差异的circRNA
linear_circ_RA_eRA <- read.delim("/mnt/data/ysmeng/circRNA/2025_DEseq2/z_HostGene/linear_circ_RA_eRA.txt",header = TRUE)
test <- linear_circ_RA_eRA %>% dplyr::filter(mRNA_significance == "Stable" & (significance %in% c("Up","Down"))) 

table(test$significance)# down 46
table(test$mRNA_significance) 
RA_eRA_down <- test %>% dplyr::filter(significance == "Down")
table(RA_eRA_down$circRNA %in% circbank_circbase$circRNA) # TRUE:46
RA_eRA_down <- RA_eRA_down %>% dplyr::select(circRNA,significance)
RA_eRA_down <- left_join(RA_eRA_down,circbank_circbase,by = "circRNA")


# 先查看circbase名称是否是字符
for (i in 1:nrow(RA_eRA_down)) {
  if (!is.na(RA_eRA_down[i, 4])) {
    circRNA_circbaseID  <- RA_eRA_down[i,4]
    match_sequence <- circbank_sequence[grepl(circRNA_circbaseID,names(circbank_sequence))]
    writeXStringSet(match_sequence,filepath = paste0("/mnt/data/ysmeng/circRNA_peptide3/get_sequence/sequence_file/RA_eRA_down/fasta/",
                                                     circRNA_circbaseID,".fasta"),format = "fasta")
    RA_eRA_down[i,9] <- "done"
  } else {
    RA_eRA_down[i,9] <- "not_found"
    
  }
  write.table(RA_eRA_down,"/mnt/data/ysmeng/circRNA_peptide3/A_info/RA_eRA_down.txt",quote = FALSE,row.names = FALSE,sep = "\t")
}




