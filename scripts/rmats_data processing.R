#
setwd("./AS_rMATS_jcast")

GSE89408_bam_path <- read.delim("~/A_info/GSE89408_bam_path.txt", header=FALSE)
SampleGroup2 <- read.delim("~/A_info/SampleGroup2.txt")

colnames(GSE89408_bam_path)[2] <- colnames(SampleGroup2)[1]
#colnames(SampleGroup2)
SampleGroup2 <- SampleGroup2 %>% 
  dplyr::filter(Sample %in% GSE89408_bam_path$Sample)
GSE89408_bam_path <- left_join(GSE89408_bam_path,SampleGroup2,by = "Sample")

#write.table(GSE89408_bam_path,"./A_info/bam_SRR_group.txt",quote = FALSE,row.names = FALSE,sep = "\t")

# 
bam_SRR_group <- read.delim("~/A_info/bam_SRR_group.txt")

Group1 <- bam_SRR_group %>% 
  dplyr::filter(Group == "healthy") %>%   # 
  dplyr::pull(V1) %>% 
  paste(collapse = ",") %>% as.data.frame()

write.table(Group1,"./AS_rMATS_jcast/rMATS/need_file/HC_bamPath.txt",
            quote = FALSE,row.names = FALSE,col.names = FALSE,sep = "\t")

Group2 <- bam_SRR_group %>% 
  dplyr::filter(Group == "eRA")%>%   # 
  dplyr::pull(V1) %>% 
  paste(collapse = ",") %>% as.data.frame()

write.table(Group2,"./AS_rMATS_jcast/rMATS/need_file/eRA_bamPath.txt",
            quote = FALSE,row.names = FALSE,col.names = FALSE,sep = "\t")

Group3 <- bam_SRR_group %>% 
  dplyr::filter(Group == "RA")%>% 
  dplyr::pull(V1) %>% 
  paste(collapse = ",") %>% as.data.frame()  # 94
write.table(Group3,"./AS_rMATS_jcast/rMATS/need_file/RA_bamPath.txt",
            quote = FALSE,row.names = FALSE,col.names = FALSE,sep = "\t")

# rmats result ---------------------------------------------------------------

library(data.table)
library(dplyr)
library(ggplot2)

dir.create("../rmats_plots", showWarnings = FALSE)

files <- c(
  SE   = "SE.MATS.JC.txt",
  MXE  = "MXE.MATS.JC.txt",
  A3SS = "A3SS.MATS.JC.txt",
  A5SS = "A5SS.MATS.JC.txt",
  RI   = "RI.MATS.JC.txt"
)
stopifnot(all(file.exists(files)))


# A2.  --------------------------------------------------
read_rmats_core <- function(path, type){
  cols <- names(fread(path, nrows = 1))
  id_col <- if ("ID" %in% cols) "ID" else cols[1]   #
  
  need <- c(id_col, "FDR", "IncLevelDifference",
            "IncLevel1", "IncLevel2",
            "IJC_SAMPLE_1","SJC_SAMPLE_1","IJC_SAMPLE_2","SJC_SAMPLE_2")
  # IncLevel1/2，IncLevelDifference（ΔPSI）
  need <- need[need %in% cols]
  
  x <- fread(path, select = need)
  setnames(x, id_col, "event_id")
  x[, event_type := type]
  x
}

rmats_list <- mapply(read_rmats_core, path = files, type = names(files), SIMPLIFY = FALSE)
rmats_all <- rbindlist(rmats_list, use.names = TRUE, fill = TRUE) # 

dim(rmats_all)
table(rmats_all$event_type)
# 
rmats_all <- rmats_all %>%
  dplyr::mutate(
    delta_psi = IncLevelDifference,
    pass_q = (FDR < 0.05),
    abs_dpsi = abs(delta_psi)
  )

table(rmats_all$pass_q) # 
summary(rmats_all$abs_dpsi)

rmats_all <- rmats_all %>% 
  dplyr::filter(pass_q) 

write.table(rmats_all,"./AS_rMATS_jcast/rMATS/result/rmats_all.txt",quote = FALSE,row.names = FALSE,sep = "\t")

# Step 1：---------------------------
df_sig <- rmats_all %>%
  dplyr::mutate(direction = ifelse(delta_psi >= 0,
                            "Higher in HC (delta PSI>0)",
                            "Higher in eRA (delta PSI<0)"))

# ：
tab_ie <- df_sig %>%
  dplyr::count(event_type, direction, name = "n")
# 
tab_sum <- df_sig %>%
  dplyr::count(event_type, name = "total")

#
total_all <- nrow(df_sig)
# 
tab_sum <- tab_sum %>%
  dplyr::mutate(pct = total / total_all)

tab_ie <- tab_ie %>%
  left_join(tab_sum, by = "event_type") %>%
  dplyr::mutate(label = sprintf("%d / %.2f", n, 100*pct))  # 

tab_ie
tab_sum

# Step 2： --------------------------------------------------------
type_order <- rev(c("SE","MXE","A5SS","A3SS","RI"))  
tab_ie$event_type <- factor(tab_ie$event_type, levels = type_order)
tab_sum$event_type <- factor(tab_sum$event_type, levels = type_order)

# Step 3： --------------------------------------------------

library(ggplot2)
#pdf("./AS_rMATS_jcast/rMATS/result/rmats_plots/barplots/plot.pdf",width = 12,height = 8)
p_hstack <- ggplot(tab_ie, aes(x = n, y = event_type, fill = direction)) +
  geom_col(width = 0.7) +
  
  # 
  geom_text(
    data = tab_sum,
    aes(x = total, y = event_type, label = total),
    inherit.aes = FALSE,
    hjust = -0.15,
    size = 4.5,                
    fontface = "bold",
    colour = "black"
  ) +
  
  theme_bw() +
  
  labs(
    title = "HC-eRA splicing events by type (FDR < 0.05)",
    x = "Number of events",
    y = NULL,
    fill = NULL                
  ) +
  
  # 
  scale_fill_manual(
    values = c("Higher in HC (delta PSI>0)" = "#c7cfb7",   # 
               "Higher in eRA (delta PSI<0)" = "#739559")   # 
  ) +
  
  coord_cartesian(clip = "off") +   # 
  
  theme(
    # 
    plot.title = element_text(hjust = 0.5, 
                              size = 16, 
                              face = "bold",
                              colour = "black"),
    
    # 
    axis.title.x = element_text(size = 14, 
                                face = "bold",
                                colour = "black"),
    
    # 
    axis.text = element_text(size = 12, 
                             colour = "black"),
    axis.text.y = element_text(size = 13),   # 
    
    # 
    legend.title = element_blank(),          # 
    legend.text = element_text(size = 12),
    legend.position = "top",                 # 
    
    # 
    plot.margin = margin(5.5, 50, 5.5, 5.5)  # 
  )

p_hstack
# dev.off()


# scatter plot ---------------------------------------------------------------------
# Step S1： ---------------------------------
mean_csv_num <- function(s){
  if (is.na(s) || s == "") return(NA_real_)
  vals <- as.numeric(strsplit(s, ",")[[1]])
  mean(vals, na.rm = TRUE)
}
# Step S2： --------------
df_scatter <- rmats_all %>%
  dplyr::mutate(
    psi_HC  = vapply(IncLevel1, mean_csv_num, numeric(1)), # PSI，IncLevel1/IncLevel2，
    psi_eRA = vapply(IncLevel2, mean_csv_num, numeric(1)),
    direction = ifelse(delta_psi >= 0,
                       "Higher in HC",
                       "Higher in eRA"),
    high_effect = abs_dpsi >= 0.1,
    size_map = ifelse(high_effect, abs_dpsi, 0.02)
  ) %>%
  dplyr::filter(!is.na(psi_HC) & !is.na(psi_eRA))
### 
nrow(df_scatter)
summary(df_scatter$psi_HC)
summary(df_scatter$psi_eRA)
table(df_scatter$high_effect)


# 散点图 ---------------------------------------------------------------------
#pdf("./AS_rMATS_jcast/rMATS/result/rmats_plots/dot_plot/HC_eRA.pdf",width = 8,height = 6)
p_scatter <- ggplot(df_scatter, aes(x = psi_HC, y = psi_eRA)) +
  geom_point(aes(color = direction, size = abs_dpsi, alpha = high_effect)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  scale_alpha_manual(values = c(`FALSE` = 0.10, `TRUE` = 0.70), guide = "none") +
  scale_size_continuous(name = "|ΔPSI|", range = c(0.2, 3.5)) +
  
  # 
  scale_color_manual(values = c("Higher in HC" = "#c7cfb7", "Higher in eRA" = "#739559"),  
                     name = "Direction") +  # 
  
  # 
  theme_bw(base_size = 14) +  #  
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16),  #
    plot.subtitle = element_text(hjust = 0.5, size = 12),              # 
    axis.title = element_text(size = 14, face = "bold"),               # 
    axis.text = element_text(size = 12),                              # 
    legend.title = element_text(size = 12, face = "bold"),            # 
    legend.text = element_text(size = 11)                             # 
  ) +
  
  labs(
    title = "PSI scatter (HC vs eRA, FDR < 0.05)",
    # subtitle = "Point size = |ΔPSI|, opacity highlights |ΔPSI| >= 0.1",
    x = "Mean PSI in HC (IncLevel1)",
    y = "Mean PSI in eRA (IncLevel2)",
    color = "Direction"   #   
  )

p_scatter
#dev.off()
table(df_scatter$)
####################

# Junction  ------------------------------------------------------------
# Step J0： --------------------------------------------
rmats_list_all <- mapply(read_rmats_core, path = files, type = names(files), SIMPLIFY = FALSE)
#  rmats_full
rmats_full <- rbindlist(rmats_list_all, use.names = TRUE, fill = TRUE)

rmats_full <- rmats_full %>%
  dplyr::mutate(
    delta_psi = IncLevelDifference,
    abs_dpsi = abs(delta_psi),
    sig = (FDR < 0.05)
  )
table(rmats_full$sig)
# Step J1： -----------------------------
sum_csv_int <- function(s){
  if (is.na(s) || s == "") return(NA_integer_)
  vals <- as.integer(strsplit(s, ",")[[1]])
  sum(vals, na.rm = TRUE)
}

rmats_full$SJC1_sum <- vapply(rmats_full$SJC_SAMPLE_1, sum_csv_int, integer(1))
rmats_full$SJC2_sum <- vapply(rmats_full$SJC_SAMPLE_2, sum_csv_int, integer(1))
#

rmats_full <- rmats_full %>%
  dplyr::mutate(
    SJC_total = SJC1_sum + SJC2_sum,
    SJC_min_group = pmin(SJC1_sum, SJC2_sum),
    logSJC = log1p(SJC_total),
    logSJCmin = log1p(SJC_min_group)
  )

summary(rmats_full$SJC_total)
summary(rmats_full$SJC_min_group)
# 
rmats_full <- rmats_full %>%
  dplyr::mutate(sig = factor(sig, levels = c("TRUE", "FALSE")))  # 


library(ggplot2)

# Step J2：
# pdf("./AS_rMATS_jcast/rMATS/result/rmats_plots/density/total_HC_eRA.pdf",width = 6,height = 5)
p_sjc_total <- ggplot(rmats_full, aes(x = logSJC, fill = sig)) +
  geom_density(alpha = 0.5) +
  theme_bw() +
  labs(
    title = "skipping junctions (SJC_total)",
    #subtitle = "Significant (FDR<0.05) vs non-significant events",
    x = "log(1 + SJC_total)", 
    y = "Density", 
    fill = "FDR<0.05"
  ) +
  theme(
    # 
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
    # 
    plot.subtitle = element_text(hjust = 0.5, size = 13),
    
    # 
    axis.title = element_text(size = 14, face = "bold"),
    
    # 
    axis.text = element_text(size = 12),
    
    # 
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 11),
    legend.position = "top"   #  
  )

p_sjc_total
# dev.off()

# Step J3：
# pdf("./AS_rMATS_jcast/rMATS/result/rmats_plots/density/min_HC_eRA.pdf",width = 6,height = 5)
p_sjc_min <- ggplot(rmats_full, aes(x = logSJCmin, fill = sig)) +
  geom_density(alpha = 0.5) +
  theme_bw() +
  labs(
    title = "skipping junctions (min SJC across groups)",
    #subtitle = "min(SJC_group1, SJC_group2) reduces one-sided low-coverage events",
    x = "log(1 + min(SJC_group1, SJC_group2))", 
    y = "Density", 
    fill = "FDR<0.05"
  ) +
  theme(
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5, size = 13),
    axis.title = element_text(size = 14, face = "bold"),
    axis.text = element_text(size = 12),
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 11),
    legend.position = "top"
  )

p_sjc_min
# dev.off()
#############################################################################
#  ----------------------------------------------------------------------
suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(ggplot2)
})
# ---------- 1)  ----------
mean_csv_num <- function(s){
  if (is.na(s) || s == "") return(NA_real_)
  vals <- as.numeric(strsplit(s, ",")[[1]])
  mean(vals, na.rm = TRUE)
}

sum_csv_int <- function(s){
  if (is.na(s) || s == "") return(NA_integer_)
  vals <- as.integer(strsplit(s, ",")[[1]])
  sum(vals, na.rm = TRUE)
}

read_rmats_core <- function(path, type){
  cols <- names(fread(path, nrows = 1))
  id_col <- if ("ID" %in% cols) "ID" else cols[1]
  
  need <- c(id_col, "FDR", "IncLevelDifference",
            "IncLevel1", "IncLevel2",
            "IJC_SAMPLE_1","SJC_SAMPLE_1","IJC_SAMPLE_2","SJC_SAMPLE_2")
  need <- need[need %in% cols]
  
  x <- fread(path, select = need)
  setnames(x, id_col, "event_id")
  x[, event_type := type]
  x
}
###
# ---------- 2)  ----------
# 
col_g1 <- "#c7cfb7"
col_g2 <- "#739559"

run_rmats_one <- function(rmats_dir,
                          group1 = "Group1",
                          group2 = "Group2",
                          out_root,
                          fdr_cut = 0.05,
                          dpsi_cut = 0.1,
                          draw_total_density = TRUE){
  
  stopifnot(dir.exists(rmats_dir))
  contrast <- basename(rmats_dir)
  
  # 
  out_bar  <- file.path(out_root, "barplots")
  out_dot  <- file.path(out_root, "dot_plot")
  out_den  <- file.path(out_root, "density")
  dir.create(out_bar, recursive = TRUE, showWarnings = FALSE)
  dir.create(out_dot, recursive = TRUE, showWarnings = FALSE)
  dir.create(out_den, recursive = TRUE, showWarnings = FALSE)
  
  # 
  files <- c(
    SE   = file.path(rmats_dir, "SE.MATS.JC.txt"),
    MXE  = file.path(rmats_dir, "MXE.MATS.JC.txt"),
    A3SS = file.path(rmats_dir, "A3SS.MATS.JC.txt"),
    A5SS = file.path(rmats_dir, "A5SS.MATS.JC.txt"),
    RI   = file.path(rmats_dir, "RI.MATS.JC.txt")
  )
  if (!all(file.exists(files))) {
    missing <- names(files)[!file.exists(files)]
    stop("Missing files in ", rmats_dir, ": ", paste(missing, collapse = ", "))
  }
  
  # ---- A)  
  rmats_full <- data.table::rbindlist(
    mapply(read_rmats_core, path = files, type = names(files), SIMPLIFY = FALSE),
    use.names = TRUE, fill = TRUE
  )
  
  rmats_full <- tibble::as_tibble(rmats_full)
  
  rmats_full <- rmats_full %>%
    dplyr::mutate(
      FDR = as.numeric(FDR),
      delta_psi = IncLevelDifference,
      abs_dpsi  = abs(delta_psi),
      sig = (!is.na(FDR) & FDR < fdr_cut & abs_dpsi > 0.1)
    )
  
  # ---- B) 
  rmats_sig <- rmats_full %>%
    dplyr::filter(sig) %>%
    dplyr::mutate(
      pass_q = TRUE,
      direction = ifelse(delta_psi >= 0, "Higher in group1", "Higher in group2")
      
    )
  tab_ie  <- dplyr::count(rmats_sig, event_type, direction, name = "n")
  tab_sum <- dplyr::count(rmats_sig, event_type, name = "total")
  
  total_all <- nrow(rmats_sig)
  tab_sum <- tab_sum %>%
    dplyr::mutate(pct = total / total_all)
  
  # 
  type_order <- rev(c("SE","MXE","A5SS","A3SS","RI"))
  tab_ie$event_type  <- factor(tab_ie$event_type,  levels = type_order)
  tab_sum$event_type <- factor(tab_sum$event_type, levels = type_order)
  
  p_bar <- ggplot(tab_ie, aes(x = n, y = event_type, fill = direction)) +
    geom_col(width = 0.7) +
    geom_text(
      data = tab_sum,
      aes(x = total, y = event_type, label = total),
      inherit.aes = FALSE,
      hjust = -0.15, size = 6, fontface = "bold", colour = "black"
    ) +
    theme_bw() +
    labs(
      title = paste0(contrast, " (FDR < ", fdr_cut, ")"),
      x = "Number of events", y = NULL, fill = NULL
    ) +
    scale_fill_manual(
      values = c("Higher in group1" = col_g1,
                 "Higher in group2" = col_g2),
      labels = c("Higher in group1" = paste0("Higher in ", group1, " (ΔPSI>0.1)"),
                 "Higher in group2" = paste0("Higher in ", group2, " (ΔPSI<0.1)"))
    ) +
    coord_cartesian(clip = "off") +
    theme(
      plot.title = element_text(hjust = 0.5, size = 22, face = "bold"),
      axis.title.x = element_text(size = 20, face = "bold"),
      axis.text = element_text(size = 18, colour = "black"),
      axis.text.y = element_text(size = 18),
      legend.text = element_text(size = 17),
      legend.position = "top",
      plot.margin = margin(5.5, 50, 5.5, 5.5)
    )
  
  # pdf(file.path(out_bar, paste0(contrast, ".barplot.pdf")), width = 7.5, height = 5)
  print(p_bar)
  # dev.off()
  
  # ============ ：PSI scatter（FDR<cut） ============
  df_scatter <- rmats_sig %>%
    dplyr::mutate(
      psi_g1 = vapply(IncLevel1, mean_csv_num, numeric(1)),
      psi_g2 = vapply(IncLevel2, mean_csv_num, numeric(1)),
      direction2 = ifelse(delta_psi >= 0, "Higher in group1", "Higher in group2"),
      high_effect = abs_dpsi >= dpsi_cut
    ) %>%
    dplyr::filter(!is.na(psi_g1) & !is.na(psi_g2))
  
  p_scatter <- ggplot(df_scatter, aes(x = psi_g1, y = psi_g2)) +
    geom_point(aes(color = direction2, size = abs_dpsi, alpha = high_effect)) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
    scale_alpha_manual(values = c(`FALSE` = 0.05, `TRUE` = 0.70), guide = "none") +
    scale_size_continuous(name = "|ΔPSI|", range = c(0.2, 3.5)) +
    scale_color_manual(
      values = c("Higher in group1" = col_g1,
                 "Higher in group2" = col_g2),
      labels = c("Higher in group1" = paste0("Higher in ", group1),
                 "Higher in group2" = paste0("Higher in ", group2)),
      name = "Direction"
    ) +
    theme_bw(base_size = 17) +
    labs(
      title = paste0(contrast, ", (FDR < ", fdr_cut, ")"),
      #subtitle = paste0("Point size = |ΔPSI|, opacity highlights |ΔPSI| >= ", dpsi_cut),
      x = paste0("Mean PSI in ", group1, " (IncLevel1)"),
      y = paste0("Mean PSI in ", group2, " (IncLevel2)")
    ) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 22),
      #plot.subtitle = element_text(hjust = 0.5, size = 17),
      axis.title = element_text(size = 19, face = "bold"),
      axis.text = element_text(size = 18),
      legend.title = element_text(size = 19, face = "bold"),
      legend.position = "right",
      legend.text = element_text(size = 18)
    )
  
  # pdf(file.path(out_dot, paste0(contrast, ".scatter.pdf")), width = 7, height = 5)
  print(p_scatter)
  # dev.off()
  
  
  # ============ junction density  ============
  rmats_full2 <- rmats_full %>%
    dplyr::mutate(
      SJC1_sum = vapply(SJC_SAMPLE_1, sum_csv_int, integer(1)),
      SJC2_sum = vapply(SJC_SAMPLE_2, sum_csv_int, integer(1)),
      SJC_total = SJC1_sum + SJC2_sum,
      SJC_min_group = pmin(SJC1_sum, SJC2_sum),
      logSJC = log1p(SJC_total),
      logSJCmin = log1p(SJC_min_group),
      sig = factor(sig, levels = c(TRUE, FALSE))
    )
  
  # 
  p_min <- ggplot(rmats_full2, aes(x = logSJCmin, fill = sig)) +
    geom_density(alpha = 0.5) +
    theme_bw() +
    labs(
      title = paste0(contrast),
      #subtitle = "log(1 + min(SJC_group1, SJC_group2))",
      x = "log(1 + min(SJC_group1, SJC_group2))",
      y = "Density",
      fill = paste0("FDR<", fdr_cut)
    ) +
    theme(
      plot.title = element_text(hjust = 0.5, size = 21, face = "bold"),
      #plot.subtitle = element_text(hjust = 0.5, size = 17),
      axis.title = element_text(size = 14, face = "bold"),
      axis.text = element_text(size = 16),
      legend.position = "right",
      legend.title = element_text(size = 13, face = "bold"),
      legend.text = element_text(size = 11)
    )
  
  # pdf(file.path(out_den, paste0(contrast, ".density.minSJC.pdf")), width = 5, height = 3.5)
  print(p_min)
  # dev.off()
  
  #
  if (isTRUE(draw_total_density)) {
    p_total <- ggplot(rmats_full2, aes(x = logSJC, fill = sig)) +
      geom_density(alpha = 0.5) +
      theme_bw() +
      labs(
        title = paste0("RNA support of skipping junctions (", contrast, ")"),
        #subtitle = "log(1 + SJC_total)",
        x = "log(1 + SJC_total)",
        y = "Density",
        fill = paste0("FDR<", fdr_cut)
      ) +
      theme(
        plot.title = element_text(hjust = 0.5, size = 15, face = "bold"),
        plot.subtitle = element_text(hjust = 0.5, size = 11),
        axis.title = element_text(size = 12, face = "bold"),
        axis.text = element_text(size = 11),
        legend.position = "top"
      )
    
    # pdf(file.path(out_den, paste0(contrast, ".density.SJCtotal.pdf")), width = 6, height = 5)
    print(p_total)
    # dev.off()
  }
  
  list(
    contrast = contrast,
    n_all = nrow(rmats_full),
    n_sig = nrow(rmats_sig),
    by_type_sig = dplyr::count(rmats_sig, event_type, name = "n_sig") %>%
      dplyr::arrange(dplyr::desc(n_sig)),
    dpsi_summary_sig = summary(rmats_sig$abs_dpsi)
  )
}

# ---------- 3)  ----------
contrasts <- list(
  list(dir = "./AS_rMATS_jcast/rMATS/result/hc_vs_eRA",
       g1 = "HC", g2 = "eRA"),
  list(dir = "./AS_rMATS_jcast/rMATS/result/hc_vs_RA",
       g1 = "HC", g2 = "RA"),
  list(dir = "./AS_rMATS_jcast/rMATS/result/RA_vs_eRA",
       g1 = "RA", g2 = "eRA")
)

out_root <- "./AS_rMATS_jcast/rMATS/result/rmats_plots"

stats <- lapply(contrasts, function(x){
  run_rmats_one(
    rmats_dir = x$dir,
    group1 = x$g1,
    group2 = x$g2,
    out_root = out_root,
    fdr_cut = 0.05,
    dpsi_cut = 0.1,
    draw_total_density = TRUE   # 
  )
})

for (s in stats) {
  cat("\n====================\n")
  cat("Contrast:", s$contrast, "\n")
  cat("N(all) :", s$n_all, "\n")
  cat("N(sig) :", s$n_sig, "\n")
  print(s$by_type_sig)
  print(s$dpsi_summary_sig)
}


# summary -----------------------------------------------------------
suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(tibble)
})

read_rmats_core <- function(path, type){
  cols <- names(fread(path, nrows = 1))
  id_col <- if ("ID" %in% cols) "ID" else cols[1]
  
  need <- c(id_col, "FDR", "IncLevelDifference",
            "IncLevel1", "IncLevel2",
            "IJC_SAMPLE_1","SJC_SAMPLE_1","IJC_SAMPLE_2","SJC_SAMPLE_2")
  need <- need[need %in% cols]
  
  x <- fread(path, select = need)
  setnames(x, id_col, "event_id")
  x[, event_type := type]
  x
}

run_rmats_one <- function(rmats_dir,
                          group1, group2,
                          fdr_cut = 0.05,
                          dpsi_cut = 0.1){
  
  stopifnot(dir.exists(rmats_dir))
  contrast <- basename(rmats_dir)
  
  files <- c(
    SE   = file.path(rmats_dir, "SE.MATS.JC.txt"),
    MXE  = file.path(rmats_dir, "MXE.MATS.JC.txt"),
    A3SS = file.path(rmats_dir, "A3SS.MATS.JC.txt"),
    A5SS = file.path(rmats_dir, "A5SS.MATS.JC.txt"),
    RI   = file.path(rmats_dir, "RI.MATS.JC.txt")
  )
  if (!all(file.exists(files))) {
    miss <- names(files)[!file.exists(files)]
    stop("Missing files in ", rmats_dir, ": ", paste(miss, collapse = ", "))
  }
  
  rmats_full_dt <- data.table::rbindlist(
    mapply(read_rmats_core, path = files, type = names(files), SIMPLIFY = FALSE),
    use.names = TRUE, fill = TRUE
  )
  
  rmats_full <- as_tibble(rmats_full_dt) %>%
    dplyr::mutate(
      contrast = contrast,
      group = paste0(group1, "_vs_", group2),
      FDR = as.numeric(FDR),
      delta_psi = as.numeric(IncLevelDifference),
      abs_dpsi = abs(delta_psi),
      sig = !is.na(FDR) & (FDR < fdr_cut) & (abs_dpsi > dpsi_cut)
    )
  
  # 
  rmats_sig <- rmats_full %>%
    dplyr::filter(sig) %>%
    dplyr::select(contrast, group, event_type, event_id, FDR, delta_psi, abs_dpsi,
           IncLevel1, IncLevel2, IJC_SAMPLE_1, SJC_SAMPLE_1, IJC_SAMPLE_2, SJC_SAMPLE_2)
  
  return(rmats_sig)
}

# -----------------  contrasts -----------------
contrasts <- list(
  list(dir = "./AS_rMATS_jcast/rMATS/result/hc_vs_eRA",
       g1 = "HC", g2 = "eRA"),
  list(dir = "./AS_rMATS_jcast/rMATS/result/hc_vs_RA",
       g1 = "HC", g2 = "RA"),
  list(dir = "./AS_rMATS_jcast/rMATS/result/RA_vs_eRA",
       g1 = "RA", g2 = "eRA")
)

# 1) bind
all_events_sig <- bind_rows(lapply(contrasts, \(x){
  run_rmats_one(
    rmats_dir = x$dir,
    group1 = x$g1,
    group2 = x$g2,
    fdr_cut = 0.05,
    dpsi_cut = 0.1
  )
}))

# 2) summary
summary_table <- all_events_sig %>%
  dplyr::count(contrast, group, event_type, name = "sig_count") %>%
  dplyr::group_by(contrast, group) %>%
  dplyr::mutate(
    total_sig = sum(sig_count),
    proportion = sig_count / total_sig
  ) %>%
  ungroup() %>%
  dplyr::arrange(contrast, desc(sig_count))



write.table(
  all_events_sig,
  "./AS_rMATS_jcast/rMATS/result/all_events_sig.txt",
  quote = FALSE,
  row.names = FALSE,
  sep = "\t"
)

# 
write.table(
  summary_table,
  "./AS_rMATS_jcast/rMATS/result/summary_table.txt",
  quote = FALSE,
  row.names = FALSE,
  sep = "\t"
)







