# ============================================================
# CPAT and IRESfinder analysis for circRNA coding potential
# ============================================================

library(dplyr)
library(stringr)
library(ggplot2)
library(ggvenn)
library(patchwork)

# ------------------------------------------------------------
# Global settings
# ------------------------------------------------------------

groups <- c(
  "H_eRA_down",
  "H_eRA_up",
  "H_RA_down",
  "H_RA_up",
  "RA_eRA_down"
)

paper_groups <- c("HC_eRA", "HC_RA", "RA_eRA")

cpat_cutoff <- 0.364

cpat_merge_dir <- "./circRNA_peptide3/part1_coding_prob/CPAT/A_merge"
cpat_result_dir <- "./circRNA_peptide3/part1_coding_prob/CPAT/Z_result"
cpat_plot_dir <- "./circRNA_peptide3/part1_coding_prob/CPAT/z_plot"

ires_result_dir <- "./circRNA_peptide3/part1_coding_prob/IRES_finder/Result"
ires_filter_dir <- "./circRNA_peptide3/part1_coding_prob/IRES_finder/z_Filter_result"
ires_plot_dir <- "./circRNA_peptide3/part1_coding_prob/IRES_finder/z_plot"

venn_plot_dir <- "./circRNA_peptide3/part1_coding_prob/A_venn/plot"
venn_data_dir <- "./circRNA_peptide3/part1_coding_prob/A_venn/data"
paper_dir <- "./circRNA_peptide3/part1_coding_prob/A_paper"


# ============================================================
# 1. Filter CPAT results
# ============================================================

read_cpat_merge <- function(group) {
  input_file <- file.path(cpat_merge_dir, paste0(group, ".txt"))

  df <- read.delim(input_file, header = FALSE)
  colnames(df) <- as.vector(df[1, ])

  df %>%
    dplyr::filter(seq_ID != "seq_ID")
}

filter_cpat <- function(group) {
  df <- read_cpat_merge(group) %>%
    dplyr::filter(Coding_prob >= cpat_cutoff)

  output_file <- file.path(cpat_result_dir, paste0("CPAT_", group, ".txt"))

  write.table(
    df,
    output_file,
    quote = FALSE,
    row.names = FALSE,
    sep = "\t"
  )

  invisible(df)
}

for (group in groups) {
  filter_cpat(group)
}


# ============================================================
# 2. Filter IRESfinder results
# ============================================================

read_ires_result <- function(group) {
  input_file <- file.path(ires_result_dir, paste0(group, "_IRES.txt"))
  read.delim(input_file)
}

filter_ires <- function(group) {
  df <- read_ires_result(group) %>%
    dplyr::filter(Index == "IRES")

  output_file <- file.path(ires_filter_dir, paste0("IRES_", group, ".txt"))

  write.table(
    df,
    output_file,
    quote = FALSE,
    row.names = FALSE,
    sep = "\t"
  )

  invisible(df)
}

for (group in groups) {
  filter_ires(group)
}


# ============================================================
# 3. Plot CPAT-positive proportions for individual groups
# ============================================================

plot_cpat_group <- function(group) {
  df <- read_cpat_merge(group)

  all_count <- nrow(df)
  positive_count <- sum(df$Coding_prob >= cpat_cutoff)

  plot_data <- data.frame(
    coding_prob = c("TRUE", "FALSE"),
    count = c(positive_count, all_count - positive_count),
    percentage = c(
      positive_count / all_count * 100,
      (all_count - positive_count) / all_count * 100
    )
  )

  p <- ggplot(plot_data, aes(x = coding_prob, y = percentage)) +
    geom_col(width = 0.5, fill = c("#ECA8A9", "#74AED4")) +
    geom_text(
      aes(label = paste0(round(percentage, 1), "%")),
      vjust = -0.5,
      size = 3.5
    ) +
    scale_y_continuous(breaks = seq(0, 100, 20), limits = c(0, 100)) +
    scale_x_discrete(limits = c("TRUE", "FALSE")) +
    ggtitle(paste("CPAT", group)) +
    theme(
      plot.title = element_text(hjust = 0.5),
      panel.grid = element_blank(),
      panel.background = element_rect(color = "black", fill = "transparent"),
      axis.line = element_line(linewidth = 0.0),
      aspect.ratio = 1
    )

  pdf(file.path(cpat_plot_dir, paste0("CPAT_", group, ".pdf")))
  print(p)
  dev.off()
}

for (group in groups) {
  plot_cpat_group(group)
}


# ============================================================
# 4. Plot IRES-positive proportions for individual groups
# ============================================================

plot_ires_group <- function(group) {
  df <- read_ires_result(group)

  all_count <- nrow(df)
  positive_count <- sum(df$Index == "IRES")

  plot_data <- data.frame(
    IRES = c("TRUE", "FALSE"),
    count = c(positive_count, all_count - positive_count),
    percentage = c(
      positive_count / all_count * 100,
      (all_count - positive_count) / all_count * 100
    )
  )

  p <- ggplot(plot_data, aes(x = IRES, y = percentage)) +
    geom_col(width = 0.5, fill = c("#ECA8A9", "#74AED4")) +
    geom_text(
      aes(label = paste0(round(percentage, 1), "%")),
      vjust = -0.5,
      size = 3.5
    ) +
    scale_y_continuous(breaks = seq(0, 100, 20), limits = c(0, 100)) +
    scale_x_discrete(limits = c("TRUE", "FALSE")) +
    ggtitle(paste("IRES", group)) +
    theme(
      plot.title = element_text(hjust = 0.5),
      panel.grid = element_blank(),
      panel.background = element_rect(color = "black", fill = "transparent"),
      axis.line = element_line(linewidth = 0.0),
      aspect.ratio = 1
    )

  pdf(file.path(ires_plot_dir, paste0("IRES_", group, ".pdf")))
  print(p)
  dev.off()
}

for (group in groups) {
  plot_ires_group(group)
}


# ============================================================
# 5. Compare CPAT and IRESfinder results for each group
# ============================================================

read_filtered_cpat <- function(group) {
  read.delim(file.path(cpat_result_dir, paste0("CPAT_", group, ".txt"))) %>%
    dplyr::mutate(
      seq_ID = str_replace_all(seq_ID, "HSA", "hsa"),
      seq_ID = str_replace_all(seq_ID, "CIRC", "circ")
    )
}

read_filtered_ires <- function(group) {
  read.delim(file.path(ires_filter_dir, paste0("IRES_", group, ".txt")))
}

plot_group_venn <- function(group) {
  cpat_df <- read_filtered_cpat(group)
  ires_df <- read_filtered_ires(group)

  venn_input <- list(
    CPAT = cpat_df$seq_ID,
    IRESfinder = ires_df$ID
  )

  p <- ggvenn(
    venn_input,
    c("CPAT", "IRESfinder"),
    fill_color = c("#113E7C", "#95D3D2"),
    stroke_alpha = 1,
    stroke_color = "transparent",
    text_size = 5
  ) +
    labs(title = group)

  pdf(file.path(venn_plot_dir, paste0("venn_", group, ".pdf")))
  print(p)
  dev.off()

  overlap <- intersect(cpat_df$seq_ID, ires_df$ID)

  write.table(
    overlap,
    file.path(venn_data_dir, paste0(group, ".txt")),
    quote = FALSE,
    row.names = FALSE,
    col.names = FALSE,
    sep = "\t"
  )

  invisible(overlap)
}

for (group in groups) {
  plot_group_venn(group)
}

# Four-panel Venn figure for the HC comparisons
draw_individual_venn <- function(group) {
  cpat_df <- read_filtered_cpat(group)
  ires_df <- read_filtered_ires(group)

  venn_input <- list(
    CPAT = cpat_df$seq_ID,
    IRESfinder = ires_df$ID
  )

  ggvenn(
    venn_input,
    c("CPAT", "IRESfinder"),
    set_name_size = 6,
    fill_color = c("#113E7C", "#95D3D2"),
    stroke_alpha = 1,
    stroke_color = "transparent",
    text_size = 5
  ) +
    labs(title = group) +
    theme(
      plot.title = element_text(
        face = "bold",
        size = 19,
        hjust = 0.5,
        vjust = 1
      ),
      plot.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "mm")
    ) +
    scale_y_continuous(expand = expansion(mult = .15))
}

hc_groups <- c("H_eRA_down", "H_eRA_up", "H_RA_down", "H_RA_up")
hc_venn_plots <- lapply(hc_groups, draw_individual_venn)

pdf(
  file.path("./circRNA_peptide3/part1_coding_prob/A_paper", "venn.pdf"),
  width = 10,
  height = 10
)
print(hc_venn_plots[[1]] + hc_venn_plots[[2]] + hc_venn_plots[[3]] + hc_venn_plots[[4]])
dev.off()


# ============================================================
# 6. Prepare combined groups for publication figures
# ============================================================

read_cpat_for_plot <- function(group) {
  read.delim(file.path(cpat_merge_dir, paste0(group, ".txt")))
}

cpat_raw <- lapply(groups, read_cpat_for_plot)
names(cpat_raw) <- groups

ires_raw <- lapply(groups, read_ires_result)
names(ires_raw) <- groups

cpat_combined <- list(
  HC_eRA = rbind(cpat_raw[["H_eRA_down"]], cpat_raw[["H_eRA_up"]]),
  HC_RA = rbind(cpat_raw[["H_RA_down"]], cpat_raw[["H_RA_up"]]),
  RA_eRA = cpat_raw[["RA_eRA_down"]]
)

ires_combined <- list(
  HC_eRA = rbind(ires_raw[["H_eRA_down"]], ires_raw[["H_eRA_up"]]),
  HC_RA = rbind(ires_raw[["H_RA_down"]], ires_raw[["H_RA_up"]]),
  RA_eRA = ires_raw[["RA_eRA_down"]]
)


# ============================================================
# 7. Plot combined CPAT results
# ============================================================

draw_cpat_plot <- function(need_data, need_group) {
  test <- need_data %>%
    dplyr::filter(seq_ID != "seq_ID")

  all_count <- nrow(test)

  positive <- test %>%
    dplyr::filter(Coding_prob >= cpat_cutoff)

  positive_count <- nrow(positive)

  df <- data.frame(
    coding_prob = c("TRUE", "FALSE"),
    count = c(positive_count, all_count - positive_count),
    percentage = c(
      positive_count / all_count * 100,
      (all_count - positive_count) / all_count * 100
    )
  )

  ggplot(df, aes(x = coding_prob, y = percentage)) +
    geom_col(width = 0.5, fill = c("#ECA8A9", "#74AED4")) +
    geom_text(
      aes(label = paste0(round(percentage, 1), "%")),
      vjust = -0.5,
      size = 5
    ) +
    scale_y_continuous(breaks = seq(0, 100, 20), limits = c(0, 100)) +
    scale_x_discrete(limits = c("TRUE", "FALSE")) +
    ggtitle(need_group) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 20),
      axis.title = element_text(size = 16),
      axis.text = element_text(size = 15),
      panel.grid = element_blank(),
      panel.background = element_rect(color = "black", fill = "transparent"),
      axis.line = element_line(linewidth = 0.0),
      aspect.ratio = 1
    )
}

pic_HC_eRA <- draw_cpat_plot(cpat_combined[["HC_eRA"]], paper_groups[1])
pic_HC_RA <- draw_cpat_plot(cpat_combined[["HC_RA"]], paper_groups[2])
pic_RA_eRA <- draw_cpat_plot(cpat_combined[["RA_eRA"]], paper_groups[3])

pdf(
  file.path(paper_dir, "new_CPAT_all.pdf"),
  width = 8,
  height = 8
)
print(pic_HC_eRA | pic_HC_RA | pic_RA_eRA)
dev.off()


# ============================================================
# 8. Plot combined IRESfinder results
# ============================================================

draw_ires_plot <- function(need_data, need_group) {
  positive <- need_data %>%
    dplyr::filter(Index == "IRES")

  all_count <- nrow(need_data)
  positive_count <- nrow(positive)

  df <- data.frame(
    IRES = c("TRUE", "FALSE"),
    count = c(positive_count, all_count - positive_count),
    percentage = c(
      positive_count / all_count * 100,
      (all_count - positive_count) / all_count * 100
    )
  )

  ggplot(df, aes(x = IRES, y = percentage)) +
    geom_col(width = 0.5, fill = c("#8DD3C7", "#D4B483")) +
    geom_text(
      aes(label = paste0(round(percentage, 1), "%")),
      vjust = -0.5,
      size = 5
    ) +
    scale_y_continuous(breaks = seq(0, 100, 20), limits = c(0, 100)) +
    scale_x_discrete(limits = c("TRUE", "FALSE")) +
    ggtitle(need_group) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 20),
      axis.title = element_text(size = 16),
      axis.text = element_text(size = 15),
      panel.grid = element_blank(),
      panel.background = element_rect(color = "black", fill = "transparent"),
      axis.line = element_line(linewidth = 0.0),
      aspect.ratio = 1
    )
}

pic_HC_eRA <- draw_ires_plot(ires_combined[["HC_eRA"]], paper_groups[1])
pic_HC_RA <- draw_ires_plot(ires_combined[["HC_RA"]], paper_groups[2])
pic_RA_eRA <- draw_ires_plot(ires_combined[["RA_eRA"]], paper_groups[3])

pdf(
  file.path(paper_dir, "new_IRES_all.pdf"),
  width = 8,
  height = 8
)
print(pic_HC_eRA + pic_HC_RA + pic_RA_eRA)
dev.off()


# ============================================================
# 9. Plot combined CPAT-IRESfinder Venn diagrams
# ============================================================

cpat_filtered <- lapply(groups, read_filtered_cpat)
names(cpat_filtered) <- groups

ires_filtered <- lapply(groups, read_filtered_ires)
names(ires_filtered) <- groups

cpat_venn_combined <- list(
  HC_eRA = bind_rows(cpat_filtered[["H_eRA_down"]], cpat_filtered[["H_eRA_up"]]),
  HC_RA = bind_rows(cpat_filtered[["H_RA_down"]], cpat_filtered[["H_RA_up"]]),
  RA_eRA = cpat_filtered[["RA_eRA_down"]]
)

ires_venn_combined <- list(
  HC_eRA = bind_rows(ires_filtered[["H_eRA_down"]], ires_filtered[["H_eRA_up"]]),
  HC_RA = bind_rows(ires_filtered[["H_RA_down"]], ires_filtered[["H_RA_up"]]),
  RA_eRA = ires_filtered[["RA_eRA_down"]]
)

make_venn_input <- function(cpat_df, ires_df) {
  list(
    CPAT = unique(cpat_df$seq_ID),
    IRESfinder = unique(ires_df$ID)
  )
}

venn_HC_eRA <- make_venn_input(
  cpat_venn_combined[["HC_eRA"]],
  ires_venn_combined[["HC_eRA"]]
)

venn_HC_RA <- make_venn_input(
  cpat_venn_combined[["HC_RA"]],
  ires_venn_combined[["HC_RA"]]
)

venn_RA_eRA <- make_venn_input(
  cpat_venn_combined[["RA_eRA"]],
  ires_venn_combined[["RA_eRA"]]
)

draw_venn_plot <- function(my_file, need_group) {
  ggvenn(
    my_file,
    c("CPAT", "IRESfinder"),
    set_name_size = 6,
    fill_color = c("#113E7C", "#95D3D2"),
    stroke_alpha = 1,
    stroke_color = "transparent",
    text_size = 5
  ) +
    labs(title = need_group) +
    theme(
      plot.title = element_text(
        face = "bold",
        size = 19,
        hjust = 0.5,
        vjust = 1
      ),
      plot.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "mm")
    ) +
    scale_y_continuous(expand = expansion(mult = .15))
}

pic_HC_eRA <- draw_venn_plot(venn_HC_eRA, paper_groups[1])
pic_HC_RA <- draw_venn_plot(venn_HC_RA, paper_groups[2])
pic_RA_eRA <- draw_venn_plot(venn_RA_eRA, paper_groups[3])

p_all <- pic_HC_eRA + pic_HC_RA + pic_RA_eRA +
  plot_layout(ncol = 3)

pdf(
  file.path(paper_dir, "new_venn.pdf"),
  width = 10,
  height = 4
)
print(p_all)
dev.off()


# ============================================================
# 10. Summarize intersections in the combined groups
# ============================================================

intersect_HC_eRA <- intersect(
  venn_HC_eRA$CPAT,
  venn_HC_eRA$IRESfinder
)

intersect_HC_RA <- intersect(
  venn_HC_RA$CPAT,
  venn_HC_RA$IRESfinder
)

intersect_RA_eRA <- intersect(
  venn_RA_eRA$CPAT,
  venn_RA_eRA$IRESfinder
)

intersect_count_df <- data.frame(
  group = c("HC_eRA", "HC_RA", "RA_eRA"),
  intersect_count = c(
    length(intersect_HC_eRA),
    length(intersect_HC_RA),
    length(intersect_RA_eRA)
  )
)

print(intersect_count_df)

all_intersect_id <- unique(c(
  intersect_HC_eRA,
  intersect_HC_RA,
  intersect_RA_eRA
))

total_intersect_count <- length(all_intersect_id)
print(total_intersect_count)
