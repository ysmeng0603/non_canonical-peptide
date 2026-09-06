
# The following files were generated from the four tools
# All samples were combined together
# "Chr"        "Start"      "End"        "chain"      "SRR4785812"

###------------------------------ circexplore2 ---------------------------------------------------------------------------------
circexplore2_info <- read.table("./circRNA/venn/venn_need_data1/circexplore2_info.txt", quote="\"", comment.char="",header = TRUE)
circexplore2_info <- circexplore2_info %>% mutate(file_name = sub("^\\./([^_]+)_.*","\\1",V9))
circexplore2_info <- circexplore2_info %>% mutate(file_path_to_use = sub("^\\./","",V9))
circexplore2_info <-  circexplore2_info %>% mutate(file_path = str_c("./circRNA/venn/venn_need_data1/circexplorer2/",file_path_to_use))
circexplore2_info <- circexplore2_info[,c(10:12)]
circexplore2_info <- circexplore2_info[order(circexplore2_info$file_name),]
rownames(circexplore2_info) <-c(1:nrow(circexplore2_info))
write.table(circexplore2_info,"./circRNA/venn/venn_need_data1/circexplore2_info.txt",row.names = FALSE,quote = FALSE)


# circexplore2_info <- circexplore2_info %>% mutate_all(~str_replace_all(.,"circexplore","circexplorer"))
circexplore2_info <- read.csv("./circRNA/venn/venn_need_data1/circexplore2_info.txt", sep="",header = TRUE)
circexplore2_info <-  circexplore2_info %>% mutate(file_path = str_c("./circRNA/venn/venn_need_data1/circexplorer2/",file_path_to_use))

circexplore2 <- read.delim(circexplore2_info[1,3], header=FALSE)
colnames(circexplore2) <- c("Chr","Start","End","chain",circexplore2_info[1,1])
for(i in 2:nrow(circexplore2_info)){
  df <- read.delim(circexplore2_info[i,3], header=FALSE)
  colnames(df) <- c("Chr","Start","End","chain",circexplore2_info[i,1])
  circexplore2 <- full_join(circexplore2,df,by =c("Chr","Start","End","chain"))
  
}
write.table(circexplore2,"./circRNA/venn/venn_need_data1/z_filterBSJ_allSRR/circexplore2_allSRR_comb.txt",row.names = FALSE,quote = FALSE)

###------------------------------ circRNAfinder ---------------------------------------------------------------------------------
circFinder_info <- read.table("./circRNA/venn/venn_need_data1/circFinder_info.txt", quote="\"", comment.char="",header = TRUE)
circFinder_info <- circFinder_info %>% mutate(file_name = sub("^\\./([^_]+)_.*","\\1",V9))
circFinder_info <- circFinder_info %>% mutate(file_path_to_use = sub("^\\./","",V9))
circFinder_info <-  circFinder_info %>% mutate(file_path = str_c("./circRNA/venn/venn_need_data1/circRNA_finder/",file_path_to_use))
circFinder_info <- circFinder_info[,c(10:12)]

circFinder_info <- circFinder_info[order(circFinder_info$file_name),]
rownames(circFinder_info) <-c(1:nrow(circFinder_info))
write.table(circFinder_info,"./circRNA/venn/venn_need_data1/circFinder_info.txt",row.names = FALSE,quote = FALSE)

circFinder_info <- read.csv("./circRNA/venn/venn_need_data1/circFinder_info.txt", sep="",header = TRUE)
circFinder <- read.delim(circFinder_info[1,3], header=FALSE)
circFinder <- circFinder[,c(1,2,3,5,4)]
colnames(circFinder) <- c("Chr","Start","End","chain",circFinder_info[1,1])

for(i in 2:nrow(circFinder_info)){
  df <- read.delim(circFinder_info[i,3], header=FALSE)
  df <- df[,c(1,2,3,5,4)]
  colnames(df) <- c("Chr","Start","End","chain",circFinder_info[i,1])
  circFinder <- full_join(circFinder,df,by =c("Chr","Start","End","chain")) #relationship = "many-to-many"
  
}
 
write.table(circFinder_me_comb,"./circRNA/venn/venn_need_data1/z_filterBSJ_allSRR/circFinder_allSRR_comb.txt")

###------------------------------ ciri2 ---------------------------------------------------------------------------------
# ciri2
ciri2_info <- read.table("./circRNA/venn/venn_need_data1/ciri2_info.txt", quote="\"", comment.char="",header = TRUE)
ciri2_info <- ciri2_info %>% mutate(file_name = sub("^\\./([^_]+)_.*","\\1",V9))
ciri2_info <- ciri2_info %>% mutate(file_path_to_use = sub("^\\./","",V9))
ciri2_info <-  ciri2_info %>% mutate(file_path = str_c("./circRNA/venn/venn_need_data1/ciri2/",file_path_to_use))
ciri2_info <- ciri2_info[,c(10:12)]

ciri2_info <- ciri2_info[order(ciri2_info$file_name),]
rownames(ciri2_info) <-c(1:nrow(ciri2_info))
write.table(ciri2_info,"./circRNA/venn/venn_need_data1/ciri2_info.txt",row.names = FALSE,quote = FALSE)

ciri2_info <- read.csv("./circRNA/venn/venn_need_data1/ciri2_info.txt", sep="",header = TRUE)
ciri2 <- read.delim(ciri2_info[1,3], header=FALSE)
ciri2 <- ciri2[,c(1,2,3,5,4)]
colnames(ciri2) <- c("Chr","Start","End","chain",ciri2_info[1,1])

for(i in 2:nrow(ciri2_info)){
  df <- read.delim(ciri2_info[i,3], header=FALSE)
  df <- df[,c(1,2,3,5,4)]
  colnames(df) <- c("Chr","Start","End","chain",ciri2_info[i,1])
  ciri2 <- full_join(ciri2,df,by =c("Chr","Start","End","chain")) #relationship = "many-to-many"
  
}

write.table(ciri2,"./circRNA/venn/venn_need_data1/z_filterBSJ_allSRR/ciri2_allSRR_comb.txt")

###------------------------------ find_circ ---------------------------------------------------------------------------------
# find_circ
find_circ_info <- read.table("./circRNA/venn/venn_need_data1/find_circ_info.txt", quote="\"", comment.char="",header = TRUE)
find_circ_info <- find_circ_info %>% mutate(file_name = sub("^\\./([^_]+)_.*","\\1",V9))
find_circ_info <- find_circ_info %>% mutate(file_path_to_use = sub("^\\./","",V9))
find_circ_info <-  find_circ_info %>% mutate(file_path = str_c("./circRNA/venn/venn_need_data1/find_circ/",file_path_to_use))
find_circ_info <- find_circ_info[,c(10:12)]

find_circ_info <- find_circ_info[order(find_circ_info$file_name),]
rownames(find_circ_info) <-c(1:nrow(find_circ_info))
write.table(find_circ_info,"./circRNA/venn/venn_need_data1/find_circ_info.txt",row.names = FALSE,quote = FALSE)

find_circ_info <- read.csv("./circRNA/venn/venn_need_data1/find_circ_info.txt", sep="",header = TRUE)
# find_circ_info <- find_circ_info %>% mutate_all(~str_replace_all(.,".txt",".bed"))
# find_circ_info <- find_circ_info %>% mutate_all(~str_replace_all(.,"_find_circ",""))


find_circ <- read.delim(find_circ_info[1,3], header=FALSE)
find_circ <- find_circ[,c(1,2,3,6,5)]
colnames(find_circ) <- c("Chr","Start","End","chain",find_circ_info[1,1])

for(i in 2:nrow(find_circ_info)){
  df <- read.delim(find_circ_info[i,3], header=FALSE)
  df <- df[,c(1,2,3,6,5)]
  colnames(df) <- c("Chr","Start","End","chain",find_circ_info[i,1])
  find_circ <- full_join(find_circ,df,by =c("Chr","Start","End","chain"),relationship = "many-to-many") 
  
}

write.table(find_circ,"./circRNA/venn/venn_need_data1/z_filterBSJ_allSRR/find_circ_allSRR_comb.txt")

##########################################################################################################################
# prepare the data required for the second Venn diagram and compare it against three databases
# Extract circRNAs from the tool with the largest number of predicted circRNAs
# circRNA_finder,that overlap with the results from other tools,
# thereby retaining circRNAs supported by at least two tools
# For circRNA_finder, keep the first four columns, including the "chain" column
circexplorer2 <- read.csv("./circRNA/venn/venn_need_data1/z_filterBSJ_allSRR/circexplorer2_allSRR_comb.txt", sep="")
circRNA_finder <- read.csv("./circRNA/venn/venn_need_data1/z_filterBSJ_allSRR/circFinder_allSRR_comb.txt", sep="")
ciri2 <- read.csv("./circRNA/venn/venn_need_data1/z_filterBSJ_allSRR/ciri2_allSRR_comb.txt", sep="")
find_circ <- read.csv("./circRNA/venn/venn_need_data1/z_filterBSJ_allSRR/find_circ_allSRR_comb.txt", sep="")

circRNA_finder_region <- paste(circRNA_finder$Chr,
                               circRNA_finder$Start, 
                               circRNA_finder$End,
                               circRNA_finder$chain,sep = '&')

find_circ_region  <- paste(find_circ$Chr, 
                           find_circ$Start,
                           find_circ$End,
                           find_circ$chain,sep = '&')

ciri2_region <- paste(ciri2$Chr, 
                      ciri2$Start, 
                      ciri2$End, 
                      ciri2$chain,sep = '&')

circexplorer2_region <- paste(circexplorer2$Chr,
                              circexplorer2$Start,
                              circexplorer2$End, 
                              circexplorer2$chain,sep = '&')


other3_comb <- union(circexplorer2_region,union(find_circ_region,ciri2_region)) #
circRNA_finder_in_other3 <- intersect(circRNA_finder_region,other3_comb) 
circRNA_finder_in_other3 <- as.data.frame(circRNA_finder_in_other3) #  

to_venn2_circRNA <- circRNA_finder_in_other3
to_venn2_circRNA <- to_venn2_circRNA %>% separate(col = circRNA_finder_in_other3, into = c("Chr", "Start","End","chain"), sep = "&")
write.table(to_venn2_circRNA,"./circRNA/venn/venn_need_data2/to_venn2_circRNA.txt",row.names = FALSE,quote = FALSE)

CIRCpedia_v2 <- read.csv("./circRNA/venn/venn_need_data2/three_database/CIRCpedia_v2.csv", header=FALSE)
# 
CIRCpedia_sample <- read.csv("./circRNA/venn/venn_need_data2/three_database/CIRCpedia_sample.txt")
colnames(CIRCpedia_v2) <- colnames(CIRCpedia_sample)
remove(CIRCpedia_sample)


circbase <- read.delim2("./circRNA/venn/venn_need_data2/three_database/circbase.txt")
circAtlas3 <- read.delim("./circRNA/venn/venn_need_data2/three_database/circAtlas3.txt")
# For the circAtlas3 database, subtract 1 from the start position
circAtlas3$Start <- circAtlas3$Start-1
write.table(circAtlas3,"./circRNA/venn/venn_need_data2/three_database/circAtlas3_match_start.txt",quote = FALSE,row.names = FALSE,sep = "\t")
circBank_circrna_annotation <- read.delim2("./circRNA/venn/venn_need_data2/three_database/circbank/circBank_circrna_annotation.txt")

# The circBank database uses hg19 coordinates; convert them to hg38 using liftOver

circbank_to_hg38 <- circBank_circrna_annotation
# colnames(circbank_to_hg38) <- "region"
circbank_to_hg38 <- circbank_to_hg38 %>%
  mutate(
    #SNP = NA,
    Chr = sub(":.*", "", position),               
    Start = as.numeric(sub(".*:(.*)-.*", "\\1", position)), # 
    Stop = as.numeric(sub(".*-", "", position))      # 
  )
circbank_to_hg38 <- circbank_to_hg38 %>% select(Chr,Start,Stop,everything())
circbank_to_hg38 <- circbank_to_hg38[,-c(4:6,9)]#
# 
circbank_to_hg38 <- circbank_to_hg38 %>% mutate(index = 1:nrow(circbank_to_hg38),.after = 3)
circbank_to_hg38 <- circbank_to_hg38 %>% mutate(index2 = ".",.after = 4)
# 
# 
circbank_to_hg38 <- circbank_to_hg38[,-8]
# 
circbank_to_hg38 <- circbank_to_hg38[,-8]
#  "Chr","Start","Stop","index","index2","strand","splicedSeqLength"
write.table(circbank_to_hg38,"./circRNA/venn/venn_need_data2/three_database/circbank/toLiftOver_index.bed",row.names = FALSE,col.names = FALSE,quote = FALSE)
#   liftOver
# linux
# liftOver ./circRNA/toLiftOver_index.bed ./liftOver/hg19Tohg38/chain/hg19ToHg38.over.chain.gz ./circRNA/map3.bed ./circRNA/unmap3.bed
circbank_liftover <- read.csv("./circRNA/venn/venn_need_data2/three_database/circbank/circbank_liftover.txt", sep="",header = TRUE)
colnames(circbank_liftover)[4] <- "index"
  
circBank_circrna_annotation$index <- 1:nrow(circBank_circrna_annotation)
circbank_test <- left_join(circBank_circrna_annotation,circbank_liftover,by = "index")
colnames(circbank_test)[3] <- "hg19_position"
circbank_test <- circbank_test %>% select(circBankID,circbaseID,hg19_position,V1,V2,V3,V8,everything())
colnames(circbank_test)[4:7] <- c("Chr","Start","End","IsMap")
colnames(circbank_test)[14:16] <- c("bed1","bed2","SpliceLength")
write.table(circbank_test,"./circRNA/venn/venn_need_data2/three_database/circbank/circbank_hg38.txt",row.names = FALSE,quote = FALSE,sep = "\t")
# 
circbank_test <- circbank_test[circbank_test$V8 == "map",]

######
# 
# 
test <- CIRCpedia_v2 %>% tidyr::separate(col = Location,into = c("chr","position1"),sep = ":")
test <- test %>% separate(col = position1,into = c("start","end"),sep = "-") 
test <- test %>% select(chr,start,end,everything())
write.table(test,"./circRNA/venn/venn_need_data2/three_database/CIRCpedia2.txt",row.names = FALSE,quote = FALSE)

# Define a circRNA as known if it is present in any of the databases
All_3database <- union(circAtlas3_region,union(CIRCpedia2_region,circbank_region)) # 
to_select <- intersect(circRNA_region,All_3database)

to_select <- to_select %>% as.data.frame() 

circRNA <- as.data.frame(circRNA_region)
colnames(circRNA) <- "."
# 
remaining <- anti_join(circRNA, to_select, by = ".") 

to_select <- to_select %>% separate(col = ".",into = c("Chr","Start","End"),sep = "&")
remaining <- remaining %>% separate(col = ".",into = c("Chr","Start","End"),sep = "&")
write.table(to_select,"./circRNA/venn/ResultData/Bed/select_known_circRNA.txt",row.names = FALSE,quote = FALSE)
write.table(remaining,"./circRNA/venn/ResultData/Bed/remaining_unknown_circRNA.txt",row.names = FALSE,quote = FALSE)
