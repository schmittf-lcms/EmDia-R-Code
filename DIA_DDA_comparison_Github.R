# Load required packages
library(dplyr)
library(ggplot2)

# Define base directory for DDA and DIA comparison (Mac vs Windows)
if (.Platform$OS.type == "unix") {
  distinctive_annotations <- file.path(
    "Filepath_macOS",
    "DIA_DDA_comparison.csv"
  )
} else {
  distinctive_annotations <- file.path(
    "Filepath_windows",
    "DIA_DDA_comparison.csv"
  )
}

# ---- Read CSV without headers ----
raw <- read.csv(distinctive_annotations, header = FALSE, stringsAsFactors = FALSE)

# Extract metadata
method_row <- as.character(unlist(raw[1, ]))       # Row 1 = DIA / DDA
strategytype_row <- as.character(unlist(raw[2, ])) # Row 2 = targeted / untargeted

# Actual data starts at row 3
df <- raw[-c(1, 2), ]
colnames(df) <- paste0(strategytype_row, "_", seq_along(strategytype_row)) # temporary names

# Function to extract reproducible metabolites (≥5 detections) in a single column
get_reproducible <- function(col_data) {
  clean <- col_data[col_data != "" & !is.na(col_data)]
  tab <- table(clean)
  tab[tab >= 5]  # return a named vector with metabolite and its count
}

# Get reproducible metabolites per strategy
reproducible_lists <- lapply(df, get_reproducible)

# ---- Build summary table ----
summary_df <- data.frame(
  Strategy = seq_along(reproducible_lists),
  Method = method_row,
  StrategyType = strategytype_row,
  Reproducible_Count = sapply(reproducible_lists, length),
  stringsAsFactors = FALSE
)

# ---- Build detailed metabolite table ----
detailed_df <- do.call(rbind, lapply(seq_along(reproducible_lists), function(i) {
  if (length(reproducible_lists[[i]]) == 0) return(NULL)
  data.frame(
    Strategy = i,
    Method = method_row[i],
    StrategyType = strategytype_row[i],
    Metabolite = names(reproducible_lists[[i]]),
    Count = as.numeric(reproducible_lists[[i]]),
    stringsAsFactors = FALSE
  )
}))

# Define output directory = same as input directory
output_dir <- dirname(distinctive_annotations)

# Save results
summary_file <- file.path(output_dir, "reproducibility_summary.csv")
detailed_file <- file.path(output_dir, "reproducible_metabolites.csv")

write.csv(summary_df, summary_file, row.names = FALSE)
write.csv(detailed_df, detailed_file, row.names = FALSE)

cat("\n✅ Analysis complete!\n")
cat("Summary saved to: ", summary_file, "\n")
cat("Detailed metabolite list saved to: ", detailed_file, "\n")

# ---- Visualization ----
# Set dodge width for separation between DDA and DIA
dodge_width <- 0.6

p <- ggplot(summary_df, aes(x = StrategyType, y = Reproducible_Count, fill = Method)) +
  geom_bar(
    stat = "identity",
    position = position_dodge(width = dodge_width),
    color = "black",  # black outline
    width = 0.5
  ) +
  scale_fill_manual(values = c("DDA" = "#1f78b4", "DIA" = "#d95f02")) +
  scale_y_continuous(
    expand = expansion(mult = c(0, 0.05)),  # start y at 0
    breaks = seq(0, max(summary_df$Reproducible_Count) + 50, by = 50)  # ticks every 50
  ) +
  theme_classic(base_size = 4) +  # classic theme includes axes with ticks
  theme(
    axis.text.x = element_text(size = 5, color = "black"),
    axis.text.y = element_text(size = 5, color = "black"),
    axis.title = element_text(size = 5, color = "black"),
    plot.title = element_text(size = 6, hjust = 0.5),
    legend.title = element_text(size = 5),
    legend.text = element_text(size = 5),
    legend.position = "right",
    legend.key.size = unit(0.15, "cm")
  ) +
  labs(
    title = "Replicate Annotations",
    x = "Annotation Strategy",
    y = "Number of Annotations",
    fill = "Acquisition"
  )

# Save the plot
plot_file <- file.path(output_dir, "reproducibility_barplot_paper.png")
ggsave(plot_file, p, width = 1.6, height = 1.2, dpi = 300)

cat("✅ Paper-ready grouped barplot with proper axes saved to: ", plot_file, "\n")




###
### Distinctive Annotations
###

# Load required packages
library(dplyr)
library(ggplot2)

# ---- Define base directory for DDA and DIA comparison ----
if (.Platform$OS.type == "unix") {
  distinctive_annotations <- file.path(
    "Filepath_macos",
    "DIA_DDA_comparison.csv"
  )
} else {
  distinctive_annotations <- file.path(
    "Filepath_windows",
    "DIA_DDA_comparison.csv"
  )
}

# ---- Read CSV without headers ----
raw <- read.csv(distinctive_annotations, header = FALSE, stringsAsFactors = FALSE)

# Extract metadata
method_row <- as.character(unlist(raw[1, ]))       # Row 1 = DIA / DDA
strategytype_row <- as.character(unlist(raw[2, ])) # Row 2 = targeted / untargeted

# Actual data starts at row 3
df <- raw[-c(1,2), ]
colnames(df) <- paste0(strategytype_row, "_", seq_along(strategytype_row)) # temporary names

# ---- Function to extract reproducible metabolites for a given threshold ----
get_reproducible <- function(col_data, threshold = 5) {
  clean <- col_data[col_data != "" & !is.na(col_data)]
  tab <- table(clean)
  tab[tab >= threshold]
}

# ---- Thresholds to compare ----
thresholds <- 2:5

# ---- Build summary for each threshold ----
summary_list <- lapply(thresholds, function(thresh) {
  counts <- sapply(df, function(col) length(get_reproducible(col, threshold = thresh)))
  data.frame(
    Method = method_row,
    StrategyType = strategytype_row,
    Threshold = paste0("≥", thresh),
    Reproducible_Count = counts,
    stringsAsFactors = FALSE
  )
})

summary_df <- do.call(rbind, summary_list)

# ---- Save summary as CSV ----
summary_file <- file.path(dirname(distinctive_annotations), "distinctive_annotations_summary.csv")
write.csv(summary_df, summary_file, row.names = FALSE)

# ---- Paper-ready grouped barplot ----
dodge_width <- 0.6

p <- ggplot(summary_df, aes(x = Threshold, y = Reproducible_Count, fill = Method)) +
  geom_bar(
    stat = "identity",
    position = position_dodge(width = dodge_width),
    color = "black",
    width = 0.5
  ) +
  facet_wrap(~StrategyType, nrow = 1) +  # separate panels for targeted / untargeted
  scale_fill_manual(values = c("DDA" = "#1f78b4", "DIA" = "#d95f02")) +
  scale_y_continuous(
    expand = expansion(mult = c(0, 0.05)),
    breaks = seq(0, max(summary_df$Reproducible_Count) + 50, by = 50)
  ) +
  theme_classic(base_size = 4) +
  theme(
    axis.text.x = element_text(size = 8, color = "black"),
    axis.text.y = element_text(size = 8, color = "black"),
    axis.title = element_text(size = 8, color = "black"),
    strip.text = element_text(size = 8, color = "black"),   # facet labels
    plot.title = element_text(size = 10, hjust = 0.5),
    legend.title = element_text(size = 8),
    legend.text = element_text(size = 8),
    legend.position = "right",
    legend.key.size = unit(0.15, "cm")
  ) +
  labs(
    title = "Matched annotations across technical replicate measurements of SRM 1950 extracts",
    x = "Minimum Occurrences out of 5",
    y = "Number of Annotations",
    fill = "Acquisition"
  )

# ---- Save plot ----
plot_file <- file.path(dirname(distinctive_annotations), "full_annotations_comparison.png")
ggsave(plot_file, p, width = 7, height = 4.2, dpi = 300)

cat("✅ Comparison plot saved to:", plot_file, "\n")
cat("✅ Summary table saved to:", summary_file, "\n")





### 2 out of 5 only


# ---- Filter untargeted data for threshold ≥2 ----
untargeted_2of5_df <- subset(summary_df, StrategyType == "untargeted" & Threshold == "≥2")

# Safety check: any data?
if(nrow(untargeted_2of5_df) == 0) {
  stop("⚠️ No untargeted metabolites meet the threshold of ≥2 occurrences.")
}

# ---- Plot 2-of-5 untargeted only ----
p_2of5_untargeted <- ggplot(untargeted_2of5_df, aes(x = Threshold, y = Reproducible_Count, fill = Method)) +
  geom_bar(
    stat = "identity",
    position = position_dodge(width = dodge_width),
    color = "black",
    width = 0.5
  ) +
  scale_fill_manual(values = c("DDA" = "#1f78b4", "DIA" = "#d95f02")) +
  scale_y_continuous(
    expand = expansion(mult = c(0, 0.05)),
    breaks = seq(0, max(untargeted_2of5_df$Reproducible_Count, na.rm = TRUE) + 50, by = 50)
  ) +
  theme_classic(base_size = 4) +
  theme(
    axis.text.x = element_text(size = 5, color = "black"),
    axis.text.y = element_text(size = 5, color = "black"),
    axis.title = element_text(size = 5, color = "black"),
    plot.title = element_text(size = 6, hjust = 0.5),
    legend.title = element_text(size = 5),
    legend.text = element_text(size = 5),
    legend.position = "right",
    legend.key.size = unit(0.15, "cm")
  ) +
  labs(
    title = "Distinctive Annotations",
    x = "Minimum Occurrences",
    y = "Number of Annotations",
    fill = "Acquisition"
  )

# ---- Save plot ----
plot_file_2of5 <- file.path(dirname(distinctive_annotations), "untargeted_2of5_barplot.png")
ggsave(plot_file_2of5, p_2of5_untargeted, width = 1.6, height = 1.2, dpi = 300)

cat("✅ 2-of-5 untargeted barplot saved to:", plot_file_2of5, "\n")



###
### Overlap of DDA and DIA Untargeted (2 out of 5) — Clean, Axis-Free ggplot Venn
###

cat("\n🔍 Generating axis-free Venn diagram for DDA and DIA untargeted (2-of-5 distinctive metabolites)...\n")

# Load required package
if (!requireNamespace("ggvenn", quietly = TRUE)) install.packages("ggvenn")
library(ggvenn)

# ---- Function to extract metabolites detected ≥2 times ----
get_reproducible_2of5 <- function(col_data, threshold = 2) {
  clean <- col_data[col_data != "" & !is.na(col_data)]
  tab <- table(clean)
  names(tab[tab >= threshold])
}

# ---- Extract DDA and DIA untargeted metabolite lists ----
dda_idx <- which(method_row == "DDA" & strategytype_row == "untargeted")
dia_idx <- which(method_row == "DIA" & strategytype_row == "untargeted")

if (length(dda_idx) == 0 | length(dia_idx) == 0) {
  stop("⚠️ Could not find both DDA and DIA untargeted datasets in the input file.")
}

dda_2of5 <- get_reproducible_2of5(df[[dda_idx]], threshold = 2)
dia_2of5 <- get_reproducible_2of5(df[[dia_idx]], threshold = 2)

# ---- Prepare data for ggvenn ----
venn_data <- list(
  `DDA Untargeted` = dda_2of5,
  `DIA Untargeted` = dia_2of5
)

# ---- Create publication-style Venn diagram ----
venn_plot <- ggvenn(
  data = venn_data,
  fill_color = c("#1f78b4", "#d95f02"),
  fill_alpha = 0.5,
  stroke_size = 0.5,
  set_name_size = 4,
  text_size = 4,
  show_percentage = FALSE
) +
  ggtitle("Overlap of DDA and DIA Untargeted Annotations (2 out of 5)") +
  theme_void(base_size = 10) +  # removes all axes, grids, and background
  theme(
    plot.title = element_text(size = 11, face = "bold", hjust = 0.5),
    legend.position = "bottom",
    legend.title = element_blank(),
    legend.text = element_text(size = 9),
    plot.margin = margin(10, 10, 10, 10)
  ) +
  scale_fill_manual(
    values = c("#1f78b4", "#d95f02"),
    name = "Acquisition Method"
  )

# ---- Save high-resolution image ----
venn_plot_file <- file.path(dirname(distinctive_annotations), "dda_dia_untargeted_overlap_2of5_venn_clean.png")
ggsave(venn_plot_file, venn_plot, width = 4.5, height = 4, dpi = 300)

cat("✅ Axis-free clean Venn diagram saved to:", venn_plot_file, "\n")
cat("🔹 DDA unique:", length(setdiff(dda_2of5, dia_2of5)), "\n")
cat("🔹 DIA unique:", length(setdiff(dia_2of5, dda_2of5)), "\n")
cat("🔹 Overlap:", length(intersect(dda_2of5, dia_2of5)), "\n")

















###
### Number of annotations across DDA, DIA and processing mode with SD error bars
###
# Define base directory for DDA and DIA comparison of metabolite annotations per sample (Mac vs Windows)
if (.Platform$OS.type == "unix") {
  per_sample_file <- file.path(
    "Filepath_macos",
    "DIA_DDA_per_sample_annotation.csv"
  )
} else {
  per_sample_file <- file.path(
    "Filepath_windows",
    "DIA_DDA_per_sample_annotation.csv"
  )
}


library(dplyr)
library(tidyr)
library(ggplot2)

# ---- Read CSV without headers ----
raw <- read.csv(per_sample_file, header = FALSE, stringsAsFactors = FALSE)

# ---- Extract metadata ----
method_row <- as.character(unlist(raw[1, ]))      # DDA/DIA
targeted_row <- as.character(unlist(raw[2, ]))    # targeted/untargeted
meta_row <- as.character(unlist(raw[3, ]))        # Sample_name / Metabolite_name

# ---- Actual data ----
df <- raw[-c(1,2,3), ]
colnames(df) <- paste(method_row, meta_row, sep = "_")

# ---- Identify sample and metabolite columns ----
sample_cols <- grep("_Sample_name$", colnames(df))
metabolite_cols <- grep("_Metabolite_name$", colnames(df))

# ---- Convert to long format ----
long_df <- data.frame(
  Sample = unlist(lapply(sample_cols, function(i) df[, i])),
  Metabolite = unlist(lapply(metabolite_cols, function(i) df[, i])),
  Method = rep(method_row[sample_cols], each = nrow(df)),
  Targeted = rep(targeted_row[sample_cols], each = nrow(df))
)

# Remove empty/NA metabolites
long_df <- long_df[!is.na(long_df$Metabolite) & long_df$Metabolite != "", ]

# ---- Count metabolites per sample ----
counts_df <- long_df %>%
  group_by(Sample, Method, Targeted) %>%
  summarise(Count = n(), .groups = "drop")

# ---- Compute mean ± SD per Method × Targeted ----
stats_df <- counts_df %>%
  group_by(Method, Targeted) %>%
  summarise(
    Mean = mean(Count),
    SD = sd(Count),
    .groups = "drop"
  )

# ---- Plot 1: Mean ± SD (original) ----
dodge_width <- 0.7

p <- ggplot(stats_df, aes(x = Targeted, y = Mean, fill = Method)) +
  geom_bar(
    stat = "identity",
    position = position_dodge(width = dodge_width),
    color = "black",
    width = 0.6
  ) +
  geom_errorbar(
    aes(ymin = Mean - SD, ymax = Mean + SD),
    width = 0.2,
    position = position_dodge(width = dodge_width)
  ) +
  scale_fill_manual(values = c("DDA" = "#1f78b4", "DIA" = "#d95f02")) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
  theme_classic(base_size = 10) +
  theme(
    axis.text.x = element_text(size = 8, color = "black"),
    axis.text.y = element_text(size = 8, color = "black"),
    axis.title = element_text(size = 8, color = "black"),
    plot.title = element_text(size = 10, hjust = 0.5, face = "bold"),
    legend.title = element_text(size = 8),
    legend.text = element_text(size = 8),
    legend.position = "right",
    legend.key.size = unit(0.3, "cm")
  ) +
  labs(
    title = "Average annotations across technical replicate measurements of SRM 1950 extracts",
    x = "Annotation Strategy",
    y = "Average No. of Annotations ± SD",
    fill = "Acquisition Method"
  )

# ---- Save original plot ----
plot_file <- file.path(dirname(per_sample_file), "per_sample_metabolite_counts.png")
ggsave(plot_file, p, width = 7, height = 4.2, dpi = 300)
cat("✅ Plot saved to:", plot_file, "\n")

# ---- Plot 2: Add SD labels above bars ----
p_labels <- p +
  geom_text(
    aes(
      y = Mean + SD,  # Position labels above the top of error bars
      label = paste0("SD=", round(SD, 1))
    ),
    position = position_dodge(width = dodge_width),
    vjust = -0.6,  # Slightly above the error bar
    size = 3
  ) +
  labs(
    title = "Average Metabolite Annotations per Sample"
  )


# ---- Save labeled plot ----
plot_file_sd <- file.path(dirname(per_sample_file), "per_sample_metabolite_counts_with_SD_labels.png")
ggsave(plot_file_sd, p_labels, width = 7, height = 4.2, dpi = 300)
cat("✅ Plot with SD labels saved to:", plot_file_sd, "\n")


# ---- Save detailed and summary data to CSV ----
counts_csv <- file.path(dirname(per_sample_file), "per_sample_metabolite_counts_raw.csv")
write.csv(counts_df, counts_csv, row.names = FALSE)
cat("✅ Raw per-sample counts saved to:", counts_csv, "\n")

stats_csv <- file.path(dirname(per_sample_file), "per_sample_metabolite_counts_summary.csv")
write.csv(stats_df, stats_csv, row.names = FALSE)
cat("✅ Summary statistics saved to:", stats_csv, "\n")
