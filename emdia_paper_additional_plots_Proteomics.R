
#R code for the Wiley Proteomics paper "Metabolic Profiling of the EmDia Cohort by LC-MS Reveals Empagliflozin-intake Associated Regulation of 1,5-anhydroglucitol and Urate".

#####
##### Raw data import and data rearrangement of LINEAR REGRESSION FILES
#####
{
library("ggplot2") 
library("dplyr")
library("titanic")
library("cowplot")
library("missForest")
library("matrixStats")
library("tidyverse")
library("KODAMA")
library("vegan")
library("MASS")
library("reshape2")
library("knitr")
library("ggfortify")
library("ggdendroplot")
library("pheatmap")
library("S4Vectors")
library("SummarizedExperiment")
library("pmp")
library("ggplot2")
library("reshape2")
library("gridExtra")
library("mdatools")
library("ropls")
library("mixOmics")
library("ggrepel")
library("pheatmap")
library("RColorBrewer")
library("scales")
library("colorRamp2")
library("circlize")
library("grid")
library(ggplot2)
library(readr)
library(stringr)
}

###
####Annotation targeted vs untargeted barcharts
###

# Define base directory for the HMDB Metabolite list
if (.Platform$OS.type == "unix") {
  distinctive_annotations <- file.path("filepatch_macOS", "targeted_untargeted_distinctive.csv")
} else {
  distinctive_annotations <- file.path("Filepath_windows", "targeted_untargeted_distinctive.csv")
}

# Import data
distinctive_annotations <- read_csv(distinctive_annotations)

# Reshape to long format for ggplot2
annotations_long <- pivot_longer(distinctive_annotations, everything(), names_to = "Method", values_to = "Count")

# Plot
ggplot(annotations_long, aes(x = Method, y = Count, fill = Method)) +
  geom_bar(stat = "identity", width = 0.5) +
  scale_fill_manual(values = c("#d95f02", "#1f78b4")) + # #e7298a for pink
  labs(
    title = "Distinctive annotations: DIA vs DDA untargeted",
    x = "",
    y = "Number of Annotations"
  ) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 250)) + # Ensure y-axis starts at 0
  theme_minimal() +
  theme(
    text = element_text(family = "Arial", size = 40),
    axis.text = element_text(size = 40, color = "black"),
    axis.title = element_text(size = 40, color = "black"),
    plot.title = element_blank(),
    legend.position = "none",
    axis.ticks = element_line(size = 2),  # Make ticks more prominent
    axis.line = element_line(size = 2),  # Ensure axis lines are visible
    axis.text.x = element_text(size = 40, color = "black"),  # x-axis labels
    axis.text.y = element_text(size = 40, color = "black"),  # y-axis labels
    axis.ticks.length = grid::unit(0.2, "inches"),  # Increase tick length
    panel.grid.major = element_blank(),  # Remove major grid lines
    panel.grid.minor = element_blank(),   # Remove minor grid lines
    plot.margin = margin(20, 20, 20, 20)  # Add space to the right (increase right margin)
  )


###USE this for grouping by untargeted:

# Reshape and split method names
annotations_long <- distinctive_annotations %>%
  pivot_longer(everything(), names_to = "Method", values_to = "Count") %>%
  separate(Method, into = c("Acquisition", "Targeting"), sep = " ")

ggplot(annotations_long, aes(x = Targeting, y = Count, fill = Acquisition)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.7), width = 0.5) +
  scale_fill_manual(values = c("DDA" = "#1f78b4", "DIA" = "#d95f02")) +
  labs(
    title = "Distinctive annotations",
    x = "Annotation strategy",
    y = "Number of Annotations",
    fill = "Acquisition"
  ) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 250)) +
  theme_minimal() +
  theme(
    text = element_text(family = "Arial", size = 40),
    axis.text = element_text(size = 40, color = "black"),
    axis.title = element_text(size = 40, color = "black"),
    plot.title = element_blank(),
    legend.position = "right",
    legend.title = element_text(size = 36),
    legend.text = element_text(size = 36),
    axis.ticks = element_line(size = 2),
    axis.line = element_line(size = 2),
    axis.text.x = element_text(size = 40, color = "black"),
    axis.text.y = element_text(size = 40, color = "black"),
    axis.ticks.length = grid::unit(0.2, "inches"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    plot.margin = margin(20, 20, 20, 20)
  )


###
#### Replicate annotations
###

# Define file path
if (.Platform$OS.type == "unix") {
  replicate_annotations <- file.path("filepath_macOS", "targeted_untargeted_replicate.csv")
} else {
  replicate_annotations <- file.path("filepath_windows", "targeted_untargeted_replicate.csv")
}

# Import data
replicate_annotations <- read_csv(replicate_annotations)

# Reshape and split method names
annotations_long <- replicate_annotations %>%
  pivot_longer(everything(), names_to = "Method", values_to = "Count") %>%
  separate(Method, into = c("Acquisition", "Targeting"), sep = " ")

# Define publication-friendly palette for Acquisition method
#pub_colors <- c("DDA" = "#1f78b4", "DIA" = "#e7298a") choose for blue and pink
pub_colors <- c("DDA" = "#1f78b4", "DIA" = "#d95f02")

# Plot: Group by Targeting, color by Acquisition
# Plot: Group by Targeting, color by Acquisition
ggplot(annotations_long, aes(x = Targeting, y = Count, fill = Acquisition)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.7), width = 0.6) +
  scale_fill_manual(values = pub_colors) +
  labs(
    title = "Replicate annotations",
    x = "Annotation strategy",
    y = "Number of Annotations",
    fill = "Acquisition"
  ) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, NA)) + # Ensure y-axis starts at 0
  theme_minimal() +
  theme(
    text = element_text(family = "Arial", size = 40),
    axis.text = element_text(size = 40, color = "black"),
    axis.title = element_text(size = 40, color = "black"),
    plot.title = element_blank(),
    legend.position = "right",
    axis.ticks = element_line(size = 2),  # Make ticks more prominent
    axis.line = element_line(size = 2),  # Ensure axis lines are visible
    axis.text.x = element_text(size = 40, color = "black", hjust = 0.5),  # x-axis labels
    axis.text.y = element_text(size = 40, color = "black"),  # y-axis labels
    axis.ticks.length = grid::unit(0.2, "inches"),  # Increase tick length
    panel.grid.major = element_blank(),  # Remove major grid lines
    panel.grid.minor = element_blank(),   # Remove minor grid lines
    plot.margin = margin(20, 20, 20, 20)  # Add space to the right (increase right margin)
  )



###
### Elastic net regression ORGAN DAMAGE
###

# Define base directory for the HMDB Metabolite list. Change csv to desired clinical parameter, such as glucose.csv, hba1c.csv etc.)
if (.Platform$OS.type == "unix") {
  elastic_net <- file.path("filepath_macOS", "hba1c.csv")
} else {
  elastic_net <- file.path("Filepath_windows", "hba1c.csv")
}

# Import data
elastic_net_data <- read_csv(elastic_net)

# Clean column names just in case
colnames(elastic_net_data) <- str_replace_all(colnames(elastic_net_data), " ", "")
colnames(elastic_net_data) <- str_replace_all(colnames(elastic_net_data), "\\[.*\\]", "SD")

# Keep only necessary columns
elastic_net_data <- elastic_net_data[, 1:3]
names(elastic_net_data) <- c("Variable", "direction", "Lambda_ratio")

# Convert direction into a factor for coloring
elastic_net_data$direction <- factor(elastic_net_data$direction, levels = c("+", "-"))

# Create a new column for Lambda_ratio with direction considered
elastic_net_data$Lambda_direction <- ifelse(elastic_net_data$direction == "-", 
                                            -elastic_net_data$Lambda_ratio, 
                                            elastic_net_data$Lambda_ratio)

# Plot horizontal bar chart with Lambda_ratio values based on direction
ggplot(elastic_net_data, aes(x = Lambda_direction, y = reorder(Variable, Lambda_direction), fill = direction)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = c("+" = "#1f78b4", "-" = "#d95f02")) +  # More neutral colors for publication
  labs(
    title = "Elastic Net Results for Metabolites explaining patient glucose levels at V1",
    x = "Lambda Ratio",
    y = "Metabolite",
    fill = "Direction"
  ) +
  theme_minimal() +
  theme(
    text = element_text(family = "Arial", size = 20),
    axis.text.y = element_text(size = 20),
    axis.title = element_text(size = 20),
    plot.title = element_text(size = 20, hjust = 0.5),
    legend.title = element_text(size = 20),
    legend.text = element_text(size = 20)
  )

#Plot only the top 20:

top20 <- elastic_net_data[order(-abs(elastic_net_data$Lambda_ratio)), ][1:20, ]

ggplot(top20, aes(x = Lambda_direction, y = reorder(Variable, Lambda_direction), fill = direction)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = c("+" = "#1f78b4", "-" = "#d95f02")) +
  labs(
    title = "Top 20 Elastic Net Predictors for Patient HbA1c Levels at V1",
    x = "Lambda Ratio",
    y = "Metabolite",
    fill = "Direction"
  ) +
  theme_minimal() +
  theme(
    text = element_text(family = "Arial", size = 20),
    axis.text.y = element_text(size = 20),
    axis.title = element_text(size = 20),
    plot.title = element_text(size = 20, hjust = 0.5),
    legend.title = element_text(size = 20),
    legend.text = element_text(size = 20)
  )


###
### Elastic net regression DISEASE
###

# Define base directory for the HMDB Metabolite list
if (.Platform$OS.type == "unix") {
  elastic_net <- file.path("filepath_macOS", "Liver Fibrosis.csv")
} else {
  elastic_net <- file.path("filepath_windows", "FLD.csv")
}

# Import data
elastic_net_data <- read_csv(elastic_net)

# Clean column names just in case
colnames(elastic_net_data) <- str_replace_all(colnames(elastic_net_data), " ", "")
colnames(elastic_net_data) <- str_replace_all(colnames(elastic_net_data), "\\[.*\\]", "SD")

# Keep only necessary columns
elastic_net_data <- elastic_net_data[, 1:3]
names(elastic_net_data) <- c("Variable", "direction", "Lambda_ratio")

# Convert direction into a factor for coloring
elastic_net_data$direction <- factor(elastic_net_data$direction, levels = c("+", "-"))

# Create a new column for Lambda_ratio with direction considered
elastic_net_data$Lambda_direction <- ifelse(elastic_net_data$direction == "-", 
                                            -elastic_net_data$Lambda_ratio, 
                                            elastic_net_data$Lambda_ratio)

# Plot horizontal bar chart with Lambda_ratio values based on direction
ggplot(elastic_net_data, aes(x = Lambda_direction, y = reorder(Variable, Lambda_direction), fill = direction)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = c("+" = "#1f78b4", "-" = "#d95f02")) +
  labs(
    title = "Elastic Net Results for Metabolites explaining Liver Fibrosis at V1",
    x = "Lambda Ratio",
    y = "Metabolite",
    fill = "Direction"
  ) +
  theme_minimal() +
  theme(
    text = element_text(family = "Arial", size = 20),
    axis.text.y = element_text(size = 20),
    axis.title = element_text(size = 20),
    plot.title = element_text(size = 20, hjust = 0.5),
    legend.title = element_text(size = 20),
    legend.text = element_text(size = 20)
  )


#Plot only the top 20:


top20 <- elastic_net_data[order(-abs(elastic_net_data$Lambda_ratio)), ][1:20, ]

ggplot(top20, aes(x = Lambda_direction, y = reorder(Variable, Lambda_direction), fill = direction)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = c("+" = "#1f78b4", "-" = "#d95f02")) +
  labs(
    title = "Top 20 Elastic Net Predictors for Liver Fibrosis at V1",
    x = "Lambda Ratio",
    y = "Metabolite",
    fill = "Direction"
  ) +
  theme_minimal() +
  theme(
    text = element_text(family = "Arial", size = 20),
    axis.text.y = element_text(size = 20),
    axis.title = element_text(size = 20),
    plot.title = element_text(size = 20, hjust = 0.5),
    legend.title = element_text(size = 20),
    legend.text = element_text(size = 20)
  )


###
### Load emdia corrected data POSITIVE!!!!
###

# Define base directory for the HMDB Metabolite list
if (.Platform$OS.type == "unix") {
  emdia_data <- file.path("filepath_macOS", "emdia_corrected_med_data.csv")
} else {
  emdia_data <- file.path("filepath_windows", "emdia_corrected_med_data.csv")
}

emdia_data <- read.csv2(emdia_data, header = TRUE, dec = ".", sep = ",")
head(emdia_data)


###
### Load emdia corrected data NEGATIVE!!!!
###


# Define base directory for the HMDB Metabolite list
if (.Platform$OS.type == "unix") {
  emdia_data <- file.path("filepath_macOS", "emdia_neg_corrected_med_data.csv")
} else {
  emdia_data <- file.path("filepath_windows", "emdia_neg_corrected_med_data.csv")
}

emdia_data <- read.csv2(emdia_data, header = TRUE, dec = ".", sep = ",")
head(emdia_data)


###
###
###CVs and Violin Plots of corrected emdia data
###
###

#transpose emdia_data
emdia_data_transposed <- as.data.frame(t(emdia_data))
#set name of column
emdia_data_transposed <- setNames(data.frame(t(emdia_data[,-1])), emdia_data[,1])
emdia_data_transposed <- filter(emdia_data_transposed, xor(class == "QC", class == "Sample"))
emdia_data_transposed_CV <- emdia_data_transposed[-c(2,4:6)]

# Convert row names to a column named "SampleID"
emdia_data_transposed_CV <- emdia_data_transposed_CV %>%
  rownames_to_column(var = "SampleID")

# Ensure the data contains the necessary columns
if (!all(c("SampleID", "batch", "class") %in% colnames(emdia_data_transposed_CV))) {
  stop("Data must contain SampleID, Batch, and Class columns")
}

# Reshape the data to long format for easier manipulation
data_long <- pivot_longer(emdia_data_transposed_CV, cols = -c(SampleID, batch, class), names_to = "Metabolite", values_to = "Value")
data_long[, c(5)] <- sapply(data_long[, c(5)], as.numeric)

# Calculate the coefficient of variation (CV) for each metabolite per batch and class
cv_data <- data_long %>%
  group_by(batch, class, Metabolite) %>%
  summarise(
    Mean_Value = mean(Value, na.rm = TRUE),
    SD_Value = sd(Value, na.rm = TRUE),
    CV = (SD_Value / Mean_Value) * 100,
    .groups = "drop"
  )

# Filter out rows where Mean_Value is zero to avoid infinite CVs
cv_data <- cv_data %>%
  filter(Mean_Value != 0)


# Updated violin plot
ggplot(cv_data, aes(x = class, y = CV, fill = batch)) +
  geom_violin(trim = FALSE, alpha = 0.7, color = NA) +
  geom_boxplot(
    width = 0.12,
    outlier.shape = NA,
    position = position_dodge(width = 0.9),
    color = "black"
  ) +
  geom_hline(yintercept = 25, linetype = "dotted", color = "red", linewidth = 1) +
  annotate(
    "text", x = 1.5, y = 10,
    label = "25% threshold",
    color = "red", size = 5, hjust = 0
  ) +
  coord_cartesian(ylim = c(0, 200)) +  # <- Y-axis limit here
  scale_fill_brewer(palette = "Set2") +
  labs(
    title = "CV of Metabolites by Batch and Class (Positive Ion Mode)",
    subtitle = "Violin and boxplot (CV capped at 200%)",
    x = "Class",
    y = "Coefficient of Variation (%)",
    fill = "Batch"
  ) +
  theme_minimal(base_size = 16) +
  theme(
    axis.text.x = element_text(size = 16, angle = 45, hjust = 1),
    axis.text.y = element_text(size = 16),
    axis.title.x = element_text(size = 18, margin = margin(t = 10)),
    axis.title.y = element_text(size = 18, margin = margin(r = 10)),
    legend.title = element_text(size = 16),
    legend.text = element_text(size = 14),
    plot.title = element_text(size = 20, hjust = 0.5, face = "bold"),
    plot.subtitle = element_text(size = 16, hjust = 0.5, margin = margin(b = 10))
  )


#EMDIA statistics:
###
###
###PLs-DA of positive data (visit condition)
###
###



#
#All emdia timepoints (V1-3)
#


custom_colors <- c("#1f78b4", "#33a02c")  # Blue and green, but cleaner and more distinct

#transpose emdia_data
emdia_data_transposed <- as.data.frame(t(emdia_data))
#set name of column
emdia_data_transposed <- setNames(data.frame(t(emdia_data[,-1])), emdia_data[,1])
#filter transposed dataframe for "Emdia" and "Placebo" group
emdia_groups_transposed <- filter(emdia_data_transposed, xor(Group == "Empagliflozin", Group == "Placebo"))
#remove column containing sample information for calculations
emdia_groups_transposed_values <- emdia_groups_transposed[-c(1:6)]
#define values as numeric for calculations
emdia_groups_transposed_values <- as.data.frame(apply(emdia_groups_transposed_values, 2, as.numeric))
#define variable für PLS-DA in this case the group info of each sample
emdia_groupinfo <- as.character(as.matrix(emdia_groups_transposed[c(5)]))
#

custom_colors <- c("#1f78b4", "#d95f02")  # Blue and green, but cleaner and more distinct

### >>> PCA SECTION <<< ###
# Perform PCA
pca_result <- prcomp(emdia_groups_transposed_values, center = TRUE, scale. = TRUE)

# Create a data frame for ggplot
pca_df <- as.data.frame(pca_result$x)
pca_df$Group <- emdia_groupinfo

ggplot(pca_df, aes(x = PC1, y = PC2, fill = Group)) +
  geom_point(shape = 21, size = 5, color = "black", alpha = 0.8, stroke = 0.3) +  # round filled dots with faint border
  scale_fill_manual(
    values = custom_colors,
    guide = guide_legend(override.aes = list(shape = 22, size = 5))  # square color-filled legend keys
  ) +
  geom_hline(yintercept = 0, color = "grey", linewidth = 1) +  # center y = 0
  geom_vline(xintercept = 0, color = "grey", linewidth = 1) +  # center x = 0
  labs(title = "PCA of EmDia V1-3 [ESI+", x = "PC1", y = "PC2") +
  theme_minimal(base_size = 16) +
  theme(
    panel.grid = element_blank(),            
    axis.line = element_line(color = "black"),
    axis.ticks = element_line(color = "black"),
    axis.title = element_text(size = 30),
    axis.text.x = element_text(size = 30),
    axis.text.y = element_text(size = 30),
    plot.title = element_text(size = 30, face = "bold"),  # <-- Title size updated here
    legend.title = element_blank(),
    legend.text = element_text(size = 30)
  )
### <<< END OF PCA SECTION <<< ###


# Run PLS-DA
plsda_result <- plsda(emdia_groups_transposed_values, emdia_groupinfo, ncomp = 2)

# Extract scores
plsda_df <- as.data.frame(plsda_result$variates$X)
plsda_df$Group <- emdia_groupinfo

# Plot with same styling as PCA
ggplot(plsda_df, aes(x = comp1, y = comp2, fill = Group)) +
  geom_point(shape = 21, size = 5, color = "black", alpha = 0.8, stroke = 0.3) +
  scale_fill_manual(
    values = custom_colors,
    guide = guide_legend(override.aes = list(shape = 22, size = 5))
  ) +
  geom_hline(yintercept = 0, color = "grey", linewidth = 1) +
  geom_vline(xintercept = 0, color = "grey", linewidth = 1) +
  labs(title = "PLS-DA EmDia V1-3 [ESI+]", x = "PLS-DA comp 1", y = "PLS-DA comp 2") +
  theme_minimal(base_size = 16) +
  theme(
    panel.grid = element_blank(),            
    axis.line = element_line(color = "black"),
    axis.ticks = element_line(color = "black"),
    axis.title = element_text(size = 30),
    axis.text.x = element_text(size = 30),
    axis.text.y = element_text(size = 30),
    plot.title = element_text(size = 30, face = "bold"),  # <-- Title size updated here
    legend.title = element_blank(),
    legend.text = element_text(size = 30)
  )

### >>> VIP PLOT SECTION <<< ###

# Get VIP scores
vip_scores <- vip(plsda_result)

# Convert to dataframe
vip_df <- data.frame(Variable = rownames(vip_scores), VIP = vip_scores[,1])

# Get the loadings to estimate group association (positive = one group, negative = other)
loadings <- plsda_result$loadings$X[,1]

# Calculate average intensities for each variable across all samples
intensities <- colMeans(emdia_groups_transposed_values, na.rm = TRUE)

# Combine everything
vip_df$Loading <- loadings[match(vip_df$Variable, names(loadings))]
vip_df$Intensity <- intensities[match(vip_df$Variable, names(intensities))]

# Keep top N VIPs
top_n <- 10
vip_df_top <- vip_df[order(vip_df$VIP, decreasing = TRUE), ][1:top_n, ]

### <<< END OF VIP PLOT SECTION <<< ###

# Determine group association from loading values
vip_df_top$Group <- ifelse(vip_df_top$Loading > 0, "Empagliflozin", "Placebo")

# Plot VIP with color by group
ggplot(vip_df_top, aes(x = reorder(Variable, VIP), y = VIP, color = Group)) +
  geom_point(stat = "identity", size = 5, alpha = 0.9) +
  scale_color_manual(values = c("Empagliflozin" = "#1f78b4", "Placebo" = "#33a02c")) +
  coord_flip() +
  labs(
    title = "Top VIP Scores from PLS-DA (EmDia V1-3)",
    x = "Variables",
    y = "VIP Score",
    color = "Upregulated in "
  ) +
  theme_minimal(base_size = 16) +
  theme(
    panel.grid = element_blank(),
    axis.line = element_line(color = "black"),
    axis.ticks = element_line(color = "black"),
    axis.title = element_text(size = 16),
    legend.title = element_text(size = 14),
    legend.position = "right",
    plot.title = element_text(face = "bold")
  )



#PLSDA for sex dependency in V1-3
emdia_sexinfo <- as.character(as.matrix(emdia_groups_transposed[c(4)]))
plsda_sex_result <- plsda(emdia_groups_transposed_values, emdia_sexinfo, ncomp = 2)

# Extract scores from the PLSDA result object
plsda_df_sex <- as.data.frame(plsda_sex_result$variates$X)
plsda_df_sex$Sex <- emdia_sexinfo  # Add sex group info

# Plot with ggplot in the same style as PCA
ggplot(plsda_df_sex, aes(x = comp1, y = comp2, fill = Sex)) +
  geom_point(shape = 21, size = 5, color = "black", alpha = 0.8, stroke = 0.3) +
  scale_fill_manual(
    values = custom_colors,
    guide = guide_legend(override.aes = list(shape = 22, size = 5))
  ) +
  geom_hline(yintercept = 0, color = "grey", linewidth = 1) +
  geom_vline(xintercept = 0, color = "grey", linewidth = 1) +
  labs(title = "PLS-DA: Sex Dependency in EmDia V1-3", x = "PLS-DA comp 1", y = "PLS-DA comp 2") +
  theme_minimal(base_size = 16) +
  theme(
    panel.grid = element_blank(),
    axis.line = element_line(color = "black"),
    axis.ticks = element_line(color = "black"),
    axis.title = element_text(size = 16),
    legend.title = element_blank()
  )



#
#Timepoint V2+3
#

#transpose emdia_data
emdia_data_transposed <- as.data.frame(t(emdia_data))
#set name of column
emdia_data_transposed <- setNames(data.frame(t(emdia_data[,-1])), emdia_data[,1])
emdia_groups_transposed_V23 <- filter(emdia_data_transposed, xor(Visit == "V2", Visit == "V3"))
emdia_groups_transposed_V23_values <- emdia_groups_transposed_V23[-c(1:6)]
emdia_groups_transposed_V23_values <- as.data.frame(apply(emdia_groups_transposed_V23_values, 2, as.numeric))
emdia_groupinfo_V23 <- as.character(as.matrix(emdia_groups_transposed_V23[c(5)]))

### >>> PCA SECTION <<< ###
# Perform PCA
pca_result_V23 <- prcomp(emdia_groups_transposed_V23_values, center = TRUE, scale. = TRUE)

# Create a data frame for ggplot
pca_df_V23 <- as.data.frame(pca_result_V23$x)
pca_df_V23$Group <- emdia_groupinfo_V23

ggplot(pca_df_V23, aes(x = PC1, y = PC2, fill = Group)) +
  geom_point(shape = 21, size = 5, color = "black", alpha = 0.8, stroke = 0.3) +  # round filled dots with faint border
  scale_fill_manual(
    values = custom_colors,
    guide = guide_legend(override.aes = list(shape = 22, size = 5))  # square color-filled legend keys
  ) +
  geom_hline(yintercept = 0, color = "grey", linewidth = 1) +  # center y = 0
  geom_vline(xintercept = 0, color = "grey", linewidth = 1) +  # center x = 0
  labs(title = "PCA of EmDia V2-3 [ESI+]", x = "PC1", y = "PC2") +
  theme_minimal(base_size = 16) +
  theme(
    panel.grid = element_blank(),            
    axis.line = element_line(color = "black"),
    axis.ticks = element_line(color = "black"),
    axis.title = element_text(size = 30),
    axis.text.x = element_text(size = 30),
    axis.text.y = element_text(size = 30),
    plot.title = element_text(size = 30, face = "bold"),  # <-- Title size updated here
    legend.title = element_blank(),
    legend.text = element_text(size = 30)
  )


# Run PLS-DA
plsda_result_V23 <- plsda(emdia_groups_transposed_V23_values, emdia_groupinfo_V23, ncomp = 2)

# Extract scores
plsda_df_V23 <- as.data.frame(plsda_result_V23$variates$X)
plsda_df_V23$Group <- emdia_groupinfo_V23

# Plot with same styling as PCA
ggplot(plsda_df_V23, aes(x = comp1, y = comp2, fill = Group)) +
  geom_point(shape = 21, size = 5, color = "black", alpha = 0.8, stroke = 0.3) +
  scale_fill_manual(
    values = custom_colors,
    guide = guide_legend(override.aes = list(shape = 22, size = 5))
  ) +
  geom_hline(yintercept = 0, color = "grey", linewidth = 1) +
  geom_vline(xintercept = 0, color = "grey", linewidth = 1) +
  labs(title = "PLS-DA EmDia V2-3 [ESI+]", x = "PLS-DA comp 1", y = "PLS-DA comp 2") +
  theme_minimal(base_size = 16) +
  theme(
    panel.grid = element_blank(),            
    axis.line = element_line(color = "black"),
    axis.ticks = element_line(color = "black"),
    axis.title = element_text(size = 30),
    axis.text.x = element_text(size = 30),
    axis.text.y = element_text(size = 30),
    plot.title = element_text(size = 30, face = "bold"),  # <-- Title size updated here
    legend.title = element_blank(),
    legend.text = element_text(size = 30),
    legend.position = "bottom"
  )


#PLSDA for sex dependency in V1
emdia_sexinfo_V23 <- as.character(as.matrix(emdia_groups_transposed_V23[c(4)]))
plsda_sex_result_V23 <- plsda(emdia_groups_transposed_V23_values, emdia_sexinfo_V23, ncomp = 2)

# Extract scores from the PLSDA result object
plsda_df_sex_V23 <- as.data.frame(plsda_sex_result_V23$variates$X)
plsda_df_sex_V23$Sex <- emdia_sexinfo_V23  # Add sex group info

# Plot with ggplot in the same style as PCA
ggplot(plsda_df_sex_V23, aes(x = comp1, y = comp2, fill = Sex)) +
  geom_point(shape = 21, size = 5, color = "black", alpha = 0.8, stroke = 0.3) +
  scale_fill_manual(
    values = custom_colors,
    guide = guide_legend(override.aes = list(shape = 22, size = 5))
  ) +
  geom_hline(yintercept = 0, color = "grey", linewidth = 1) +
  geom_vline(xintercept = 0, color = "grey", linewidth = 1) +
  labs(title = "PLS-DA: Sex Dependency in EmDia V2+3", x = "PLS-DA comp 1", y = "PLS-DA comp 2") +
  theme_minimal(base_size = 16) +
  theme(
    panel.grid = element_blank(),
    axis.line = element_line(color = "black"),
    axis.ticks = element_line(color = "black"),
    axis.title = element_text(size = 16),
    legend.title = element_blank()
  )




#
#Timepoint V1
#



#transpose emdia_data
emdia_data_transposed <- as.data.frame(t(emdia_data))
#set name of column
emdia_data_transposed <- setNames(data.frame(t(emdia_data[,-1])), emdia_data[,1])
emdia_groups_transposed_V1 <- filter(emdia_data_transposed, (Visit == "V1"))
emdia_groups_transposed_V1_values <- emdia_groups_transposed_V1[-c(1:6)]
emdia_groups_transposed_V1_values <- as.data.frame(apply(emdia_groups_transposed_V1_values, 2, as.numeric))
emdia_groupinfo_V1 <- as.character(as.matrix(emdia_groups_transposed_V1[c(5)]))


### >>> PCA SECTION <<< ###
# Perform PCA
pca_result_V1 <- prcomp(emdia_groups_transposed_V1_values, center = TRUE, scale. = TRUE)

# Create a data frame for ggplot
pca_df_V1 <- as.data.frame(pca_result_V1$x)
pca_df_V1$Group <- emdia_groupinfo_V1

ggplot(pca_df_V1, aes(x = PC1, y = PC2, fill = Group)) +
  geom_point(shape = 21, size = 5, color = "black", alpha = 0.8, stroke = 0.3) +  # round filled dots with faint border
  scale_fill_manual(
    values = custom_colors,
    guide = guide_legend(override.aes = list(shape = 22, size = 5))  # square color-filled legend keys
  ) +
  geom_hline(yintercept = 0, color = "grey", linewidth = 1) +  # center y = 0
  geom_vline(xintercept = 0, color = "grey", linewidth = 1) +  # center x = 0
  labs(title = "PCA of EmDia V1", x = "PC1", y = "PC2") +
  theme_minimal(base_size = 16) +
  theme(
    panel.grid = element_blank(),            # remove background grid lines
    axis.line = element_line(color = "black"), # show axis lines with ticks
    axis.ticks = element_line(color = "black"),
    axis.title = element_text(size = 16),
    legend.title = element_blank()
  )


# Run PLS-DA
plsda_result_V1 <- plsda(emdia_groups_transposed_V1_values, emdia_groupinfo_V1, ncomp = 2)

# Extract scores
plsda_df_V1 <- as.data.frame(plsda_result_V1$variates$X)
plsda_df_V1$Group <- emdia_groupinfo_V1

# Plot with same styling as PCA
ggplot(plsda_df_V1, aes(x = comp1, y = comp2, fill = Group)) +
  geom_point(shape = 21, size = 5, color = "black", alpha = 0.8, stroke = 0.3) +
  scale_fill_manual(
    values = custom_colors,
    guide = guide_legend(override.aes = list(shape = 22, size = 5))
  ) +
  geom_hline(yintercept = 0, color = "grey", linewidth = 1) +
  geom_vline(xintercept = 0, color = "grey", linewidth = 1) +
  labs(title = "PLS-DA EmDia V1", x = "PLS-DA comp 1", y = "PLS-DA comp 2") +
  theme_minimal(base_size = 16) +
  theme(
    panel.grid = element_blank(),
    axis.line = element_line(color = "black"),
    axis.ticks = element_line(color = "black"),
    axis.title = element_text(size = 16),
    legend.title = element_blank()
  )


#PLSDA for sex dependency in V1
emdia_sexinfo_V1 <- as.character(as.matrix(emdia_groups_transposed_V1[c(4)]))
plsda_sex_result_V1 <- plsda(emdia_groups_transposed_V1_values, emdia_sexinfo_V1, ncomp = 2)

# Extract scores from the PLSDA result object
plsda_df_sex_V1 <- as.data.frame(plsda_sex_result_V1$variates$X)
plsda_df_sex_V1$Sex <- emdia_sexinfo_V1  # Add sex group info

# Plot with ggplot in the same style as PCA
ggplot(plsda_df_sex_V1, aes(x = comp1, y = comp2, fill = Sex)) +
  geom_point(shape = 21, size = 5, color = "black", alpha = 0.8, stroke = 0.3) +
  scale_fill_manual(
    values = custom_colors,
    guide = guide_legend(override.aes = list(shape = 22, size = 5))
  ) +
  geom_hline(yintercept = 0, color = "grey", linewidth = 1) +
  geom_vline(xintercept = 0, color = "grey", linewidth = 1) +
  labs(title = "PLS-DA: Sex Dependency in EmDia V1", x = "PLS-DA comp 1", y = "PLS-DA comp 2") +
  theme_minimal(base_size = 16) +
  theme(
    panel.grid = element_blank(),
    axis.line = element_line(color = "black"),
    axis.ticks = element_line(color = "black"),
    axis.title = element_text(size = 16),
    legend.title = element_blank()
  )


#
#Timepoint V2
#

#transpose emdia_data
emdia_data_transposed <- as.data.frame(t(emdia_data))
#set name of column
emdia_data_transposed <- setNames(data.frame(t(emdia_data[,-1])), emdia_data[,1])
emdia_groups_transposed_V2 <- filter(emdia_data_transposed, (Visit == "V2"))
emdia_groups_transposed_V2_values <- emdia_groups_transposed_V2[-c(1:6)]
emdia_groups_transposed_V2_values <- as.data.frame(apply(emdia_groups_transposed_V2_values, 2, as.numeric))
emdia_groupinfo_V2 <- as.character(as.matrix(emdia_groups_transposed_V2[c(5)]))


### >>> PCA SECTION <<< ###
# Perform PCA
pca_result_V2 <- prcomp(emdia_groups_transposed_V2_values, center = TRUE, scale. = TRUE)

# Create a data frame for ggplot
pca_df_V2 <- as.data.frame(pca_result_V2$x)
pca_df_V2$Group <- emdia_groupinfo_V2

ggplot(pca_df_V2, aes(x = PC1, y = PC2, fill = Group)) +
  geom_point(shape = 21, size = 5, color = "black", alpha = 0.8, stroke = 0.3) +  # round filled dots with faint border
  scale_fill_manual(
    values = custom_colors,
    guide = guide_legend(override.aes = list(shape = 22, size = 5))  # square color-filled legend keys
  ) +
  geom_hline(yintercept = 0, color = "grey", linewidth = 1) +  # center y = 0
  geom_vline(xintercept = 0, color = "grey", linewidth = 1) +  # center x = 0
  labs(title = "PCA of EmDia V2", x = "PC1", y = "PC2") +
  theme_minimal(base_size = 16) +
  theme(
    panel.grid = element_blank(),            # remove background grid lines
    axis.line = element_line(color = "black"), # show axis lines with ticks
    axis.ticks = element_line(color = "black"),
    axis.title = element_text(size = 16),
    legend.title = element_blank()
  )


# Run PLS-DA for Visit condition 2
plsda_result_V2 <- plsda(emdia_groups_transposed_V2_values, emdia_groupinfo_V2, ncomp = 2)

# Extract scores
plsda_df_V2 <- as.data.frame(plsda_result_V2$variates$X)
plsda_df_V2$Group <- emdia_groupinfo_V2

# Plot with same styling as PCA
ggplot(plsda_df_V2, aes(x = comp1, y = comp2, fill = Group)) +
  geom_point(shape = 21, size = 5, color = "black", alpha = 0.8, stroke = 0.3) +
  scale_fill_manual(
    values = custom_colors,
    guide = guide_legend(override.aes = list(shape = 22, size = 5))
  ) +
  geom_hline(yintercept = 0, color = "grey", linewidth = 1) +
  geom_vline(xintercept = 0, color = "grey", linewidth = 1) +
  labs(title = "PLS-DA EmDia V2", x = "PLS-DA comp 1", y = "PLS-DA comp 2") +
  theme_minimal(base_size = 16) +
  theme(
    panel.grid = element_blank(),
    axis.line = element_line(color = "black"),
    axis.ticks = element_line(color = "black"),
    axis.title = element_text(size = 16),
    legend.title = element_blank()
  )


#PLSDA for sex dependency in V2
emdia_sexinfo_V2 <- as.character(as.matrix(emdia_groups_transposed_V2[c(4)]))
plsda_sex_result_V2 <- plsda(emdia_groups_transposed_V2_values, emdia_sexinfo_V2, ncomp = 2)

# Extract scores from the PLSDA result object
plsda_df_sex_V2 <- as.data.frame(plsda_sex_result_V2$variates$X)
plsda_df_sex_V2$Sex <- emdia_sexinfo_V2  # Add sex group info

# Plot with ggplot in the same style as PCA
ggplot(plsda_df_sex_V2, aes(x = comp1, y = comp2, fill = Sex)) +
  geom_point(shape = 21, size = 5, color = "black", alpha = 0.8, stroke = 0.3) +
  scale_fill_manual(
    values = custom_colors,
    guide = guide_legend(override.aes = list(shape = 22, size = 5))
  ) +
  geom_hline(yintercept = 0, color = "grey", linewidth = 1) +
  geom_vline(xintercept = 0, color = "grey", linewidth = 1) +
  labs(title = "PLS-DA: Sex Dependency in EmDia V2", x = "PLS-DA comp 1", y = "PLS-DA comp 2") +
  theme_minimal(base_size = 16) +
  theme(
    panel.grid = element_blank(),
    axis.line = element_line(color = "black"),
    axis.ticks = element_line(color = "black"),
    axis.title = element_text(size = 16),
    legend.title = element_blank()
  )

#
#Timepoint V3
#

#transpose emdia_data
emdia_data_transposed <- as.data.frame(t(emdia_data))
#set name of column
emdia_data_transposed <- setNames(data.frame(t(emdia_data[,-1])), emdia_data[,1])
emdia_groups_transposed_V3 <- filter(emdia_data_transposed, (Visit == "V3"))
emdia_groups_transposed_V3_values <- emdia_groups_transposed_V3[-c(1:6)]
emdia_groups_transposed_V3_values <- as.data.frame(apply(emdia_groups_transposed_V3_values, 2, as.numeric))
emdia_groupinfo_V3 <- as.character(as.matrix(emdia_groups_transposed_V3[c(5)]))


### >>> PCA SECTION <<< ###
# Perform PCA
pca_result_V3 <- prcomp(emdia_groups_transposed_V3_values, center = TRUE, scale. = TRUE)

# Create a data frame for ggplot
pca_df_V3 <- as.data.frame(pca_result_V3$x)
pca_df_V3$Group <- emdia_groupinfo_V3

ggplot(pca_df_V3, aes(x = PC1, y = PC2, fill = Group)) +
  geom_point(shape = 21, size = 5, color = "black", alpha = 0.8, stroke = 0.3) +  # round filled dots with faint border
  scale_fill_manual(
    values = custom_colors,
    guide = guide_legend(override.aes = list(shape = 22, size = 5))  # square color-filled legend keys
  ) +
  geom_hline(yintercept = 0, color = "grey", linewidth = 1) +  # center y = 0
  geom_vline(xintercept = 0, color = "grey", linewidth = 1) +  # center x = 0
  labs(title = "PCA of EmDia V3", x = "PC1", y = "PC2") +
  theme_minimal(base_size = 16) +
  theme(
    panel.grid = element_blank(),            # remove background grid lines
    axis.line = element_line(color = "black"), # show axis lines with ticks
    axis.ticks = element_line(color = "black"),
    axis.title = element_text(size = 16),
    legend.title = element_blank()
  )


# Run PLS-DA
plsda_result_V3 <- plsda(emdia_groups_transposed_V3_values, emdia_groupinfo_V3, ncomp = 2)

# Extract scores
plsda_df_V3 <- as.data.frame(plsda_result_V3$variates$X)
plsda_df_V3$Group <- emdia_groupinfo_V3

# Plot with same styling as PCA
ggplot(plsda_df_V3, aes(x = comp1, y = comp2, fill = Group)) +
  geom_point(shape = 21, size = 5, color = "black", alpha = 0.8, stroke = 0.3) +
  scale_fill_manual(
    values = custom_colors,
    guide = guide_legend(override.aes = list(shape = 22, size = 5))
  ) +
  geom_hline(yintercept = 0, color = "grey", linewidth = 1) +
  geom_vline(xintercept = 0, color = "grey", linewidth = 1) +
  labs(title = "PLS-DA EmDia V3", x = "PLS-DA comp 1", y = "PLS-DA comp 2") +
  theme_minimal(base_size = 16) +
  theme(
    panel.grid = element_blank(),
    axis.line = element_line(color = "black"),
    axis.ticks = element_line(color = "black"),
    axis.title = element_text(size = 16),
    legend.title = element_blank()
  )


#PLSDA for sex dependency in V3
emdia_sexinfo_V3 <- as.character(as.matrix(emdia_groups_transposed_V3[c(4)]))
plsda_sex_result_V3 <- plsda(emdia_groups_transposed_V3_values, emdia_sexinfo_V3, ncomp = 2)

# Extract scores from the PLSDA result object
plsda_df_sex_V3 <- as.data.frame(plsda_sex_result_V3$variates$X)
plsda_df_sex_V3$Sex <- emdia_sexinfo_V3  # Add sex group info

# Plot with ggplot in the same style as PCA
ggplot(plsda_df_sex_V3, aes(x = comp1, y = comp2, fill = Sex)) +
  geom_point(shape = 21, size = 5, color = "black", alpha = 0.8, stroke = 0.3) +
  scale_fill_manual(
    values = custom_colors,
    guide = guide_legend(override.aes = list(shape = 22, size = 5))
  ) +
  geom_hline(yintercept = 0, color = "grey", linewidth = 1) +
  geom_vline(xintercept = 0, color = "grey", linewidth = 1) +
  labs(title = "PLS-DA: Sex Dependency in EmDia V3", x = "PLS-DA comp 1", y = "PLS-DA comp 2") +
  theme_minimal(base_size = 16) +
  theme(
    panel.grid = element_blank(),
    axis.line = element_line(color = "black"),
    axis.ticks = element_line(color = "black"),
    axis.title = element_text(size = 16),
    legend.title = element_blank()
  )













#Barchart
### Step 1: Statistical Analysis for Significant Metabolites (p < 0.01)
# Perform t-test and filter metabolites with p-value < 0.01


###
###Visit dependent barcharts
###

#transpose emdia_data
emdia_data_transposed <- as.data.frame(t(emdia_data))
#set name of column
emdia_data_transposed <- setNames(data.frame(t(emdia_data[,-1])), emdia_data[,1])
#filter transposed dataframe for "Emdia" and "Placebo" group
emdia_groups_transposed <- filter(emdia_data_transposed, xor(Group == "Empagliflozin", Group == "Placebo"))
#filter transposed dataframe for Visit 2 e.g.ONLY USE ONE VIST FILTER"
#emdia_groups_transposed_v <- filter(emdia_groups_transposed, Visit == "V2")
#filter transposed dataframe for Visit 2+3 e.g. ONLY USE ONE VIST FILTER"
emdia_groups_transposed_v <- filter(emdia_groups_transposed, xor(Visit == "V2", Visit == "V3"))
#optional Sex filter:
#emdia_groups_transposed_v <- filter(emdia_groups_transposed_v, Sex == "2")
#remove column containing sample information for calculations except group
emdia_groups_transposed_v_t_test <- emdia_groups_transposed_v[-c(1,2,3,4,6)]
emdia_groups_transposed_v_t_test <- as.data.frame(emdia_groups_transposed_v_t_test)
emdia_groups_transposed_v_t_test[, -1] <- sapply(emdia_groups_transposed_v_t_test[, -1], as.numeric)
#Define Group as factor
emdia_groups_transposed_v_t_test$Group <- as.factor(emdia_groups_transposed_v_t_test$Group)

# Perform t-tests
results_V <- data.frame(Metabolite = character(), p.value = numeric(), log2_fold_change = numeric())

for (metabolite in colnames(emdia_groups_transposed_v_t_test)[-1]) {
  t_test <- t.test(emdia_groups_transposed_v_t_test[[metabolite]] ~ emdia_groups_transposed_v_t_test$Group)
  
  mean_group1 <- mean(emdia_groups_transposed_v_t_test[[metabolite]][emdia_groups_transposed_v_t_test$Group == levels(emdia_groups_transposed_v_t_test$Group)[1]])
  mean_group2 <- mean(emdia_groups_transposed_v_t_test[[metabolite]][emdia_groups_transposed_v_t_test$Group == levels(emdia_groups_transposed_v_t_test$Group)[2]])
  log2_fc <- log2(mean_group2 / mean_group1)
  
  results_V <- rbind(results_V, data.frame(Metabolite = metabolite, p.value = t_test$p.value, log2_fold_change = log2_fc))
}

# Adjust p-values. FDR correction: The method adjusts the p-values to account for multiple comparisons, reducing the likelihood of false positives while maintaining power in detecting true signals.

results_V <- results_V %>%
  mutate(
    neg_log10_p = -log10(p.value),
    p.adj = p.adjust(p.value, method = "fdr"),
    significance = case_when(
      p.adj < 0.0001 ~ "****",
      p.adj < 0.001 ~ "***",
      p.adj < 0.01 ~ "**",
      p.adj < 0.05 ~ "*",
      TRUE ~ ""
    )
  )

# Filter significant metabolites
significant_metabolites_V <- results_V %>%
  filter(p.adj < 0.05) %>%
  pull(Metabolite)


custom_colors <- c("#1f78b4", "#d95f02")  # Blue and green, but cleaner and more distinct

#(normalizing by max intensity per metabolite)
if (length(significant_metabolites_V) > 0) {
  
  # Normalize to max within each metabolite group
  data_long_V <- emdia_groups_transposed_v_t_test %>%
    pivot_longer(cols = -Group, names_to = "Metabolite", values_to = "Intensity") %>%
    filter(Metabolite %in% significant_metabolites_V) %>%
    group_by(Metabolite) %>%
    mutate(Relative_Intensity = Intensity / max(Intensity) * 100) %>%
    left_join(results_V[, c("Metabolite", "significance")], by = "Metabolite") %>%
    ungroup()
  
  # Set fixed y-position for annotations
  annotation_y <- 105
  annotation_line_y <- 102
  
  # Prepare annotations for significance bars and stars
  annotation_df <- data_long_V %>%
    as_tibble() %>%
    distinct(Metabolite, significance) %>%
    mutate(x = as.numeric(factor(Metabolite)),
           xstart = x - 0.3,
           xend = x + 0.3,
           y = annotation_y,
           yline = annotation_line_y)
  
  # Plot
  # Plot
  ggplot(data_long_V, aes(x = Metabolite, y = Relative_Intensity, fill = Group)) +
    geom_bar(stat = "summary", fun = mean, position = position_dodge(width = 0.9), width = 0.8) +
    
    # Apply custom colors
    scale_fill_manual(values = custom_colors) +
    
    # Line above bars
    geom_segment(data = annotation_df,
                 aes(x = xstart, xend = xend, y = yline, yend = yline),
                 inherit.aes = FALSE) +
    
    # Stars
    geom_text(data = annotation_df,
              aes(x = x, y = y, label = significance),
              inherit.aes = FALSE, vjust = 0, size = 5) +
    
    
    labs(
      title = "Relative Intensity Barplot of Significant Metabolites (Visit V2-3)",
      subtitle = "P-value significance: **** p < 0.0001 | *** p < 0.001 | ** p < 0.01 | * p < 0.05",
      y = "Relative Intensity (%)",
      x = "Metabolite"
    ) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 10))
  
  
} else {
  print("No significant metabolites found (p.adj < 0.05).")
}


# Plot
# Plot
# Plot
ggplot(data_long_V, aes(x = Metabolite, y = Relative_Intensity, fill = Group)) +
  geom_boxplot(position = position_dodge(width = 0.8), width = 0.7, outlier.shape = NA) +
  
  # Apply custom colors
  scale_fill_manual(values = custom_colors) +
  
  # Line above boxes
  geom_segment(data = annotation_df,
               aes(x = xstart, xend = xend, y = yline, yend = yline),
               inherit.aes = FALSE, size = 1.5) +
  
  # Rectangle around deoxyhexose (adjusted height)
  geom_rect(aes(xmin = which(levels(factor(data_long_V$Metabolite)) == "DEOXYHEXOSE") - 0.5,
                xmax = which(levels(factor(data_long_V$Metabolite)) == "DEOXYHEXOSE") + 0.5,
                ymin = quantile(data_long_V$Relative_Intensity, 0.5) - 50,  # Adjust ymin
                ymax = quantile(data_long_V$Relative_Intensity, 0.6) + 50), # Adjust ymax
            inherit.aes = FALSE, color = "#D62728", linetype = "dashed", size = 1, fill = NA) +   # No fill, transparent
  
  # Rectangle around urate (adjusted height)
  geom_rect(aes(xmin = which(levels(factor(data_long_V$Metabolite)) == "URATE") - 0.5,
                xmax = which(levels(factor(data_long_V$Metabolite)) == "URATE") + 0.5,
                ymin = quantile(data_long_V$Relative_Intensity, 0.5) - 50,  # Adjust ymin
                ymax = quantile(data_long_V$Relative_Intensity, 0.6) + 50), # Adjust ymax
            inherit.aes = FALSE, color = "#D62728", linetype = "dashed", size = 1, fill = NA) +   # No fill, transparent
  
  # Stars
  geom_text(data = annotation_df,
            aes(x = x, y = y, label = significance),
            inherit.aes = FALSE, vjust = 0.5, size = 8, family = "Arial") +
  
  # Title + subtitle with significance legend
  labs(
    title = "Boxplot of Significant Metabolites (Visit V2-3) [ESI-]",
    subtitle = "P-value significance: **** p < 0.0001 | *** p < 0.001 | ** p < 0.01 | * p < 0.05",
    y = "Relative Intensity (%)",
    x = ""
  ) +
  
  theme_minimal(base_family = "Arial") +
  theme(
    plot.title = element_text(size = 22, face = "bold"),
    plot.subtitle = element_text(size = 20),
    axis.title = element_text(size = 20),
    axis.text.x = element_text(angle = 90, hjust = 1, size = 20),
    axis.text.y = element_text(size = 20),
    legend.title = element_text(size = 20),
    legend.text = element_text(size = 20)
  )



#log transofrom:

if (length(significant_metabolites_V) > 0) {
  
  # Apply log transformation to intensity values (log2)
  data_long_V <- emdia_groups_transposed_v_t_test %>%
    pivot_longer(cols = -Group, names_to = "Metabolite", values_to = "Intensity") %>%
    filter(Metabolite %in% significant_metabolites_V) %>%
    group_by(Metabolite) %>%
    mutate(Log_Intensity = log2(Intensity + 1)) %>%  # Apply log2 transformation (adding 1 to avoid log(0))
    left_join(results_V[, c("Metabolite", "significance")], by = "Metabolite") %>%
    ungroup()
  
  # Set fixed y-position for annotations
  annotation_y <- 27
  annotation_line_y <- 25
  
  # Prepare annotations for significance bars and stars
  annotation_df <- data_long_V %>%
    as_tibble() %>%
    distinct(Metabolite, significance) %>%
    mutate(x = as.numeric(factor(Metabolite)),
           xstart = x - 0.3,
           xend = x + 0.3,
           y = annotation_y,
           yline = annotation_line_y)
  
  # Plot
  ggplot(data_long_V, aes(x = Metabolite, y = Log_Intensity, fill = Group)) +
    geom_bar(stat = "summary", fun = mean, position = position_dodge(width = 0.9), width = 0.8) +
    
    # Line above bars
    geom_segment(data = annotation_df,
                 aes(x = xstart, xend = xend, y = yline, yend = yline),
                 inherit.aes = FALSE) +
    
    # Stars
    geom_text(data = annotation_df,
              aes(x = x, y = y, label = significance),
              inherit.aes = FALSE, vjust = 0, size = 5) +
    
    # Add star legend in top-right corner
    annotate("text", x = Inf, y = Inf, label = "*** p < 0.001\n** p < 0.01\n* p < 0.05",
             hjust = 1, vjust = 3, size = 4, fontface = "italic") +
    
    labs(title = "Log2 Transformed Relative Intensity Barplot of Significant Metabolites (Visit V2-3)",
         y = "Log2 Transformed Intensity", x = "Metabolite") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 10))
  
} else {
  print("No significant metabolites found (p.adj < 0.05).")
}


#Boxplot with Fold change
ggplot(results_V %>% filter(Metabolite %in% significant_metabolites_V), 
       aes(x = reorder(Metabolite, log2_fold_change), y = log2_fold_change, fill = significance)) +
  geom_bar(stat = "identity", show.legend = FALSE) +
  coord_flip() + # Flip the axes to make it easier to read
  labs(title = "Log2 Fold Change of Significant Metabolites (Visit V2-3)",
       x = "Metabolite", y = "Log2 Fold Change") +
  theme_minimal() +
  theme(axis.text.y = element_text(size = 10))


#Dotplot of Foldchange
ggplot(results_V %>% filter(Metabolite %in% significant_metabolites_V), 
       aes(x = reorder(Metabolite, log2_fold_change), y = log2_fold_change, size = neg_log10_p, color = significance)) +
  geom_point() +
  coord_flip() +
  labs(title = "Dotplot of Fold Change and Significance (Visit V2-3)",
       x = "Metabolite", y = "Log2 Fold Change") +
  scale_size_continuous(name = "Significance (-log10 p-value)") +
  theme_minimal() +
  theme(axis.text.y = element_text(size = 10))



#Metabolite Correlation analysis (VISIT dependent): 

# Define the metabolite you want to correlate with others
selected_metabolite <- "DEOXYHEXOSE"  # replace with your metabolite of interest

# Ensure the metabolite exists in the dataframe
if (!(selected_metabolite %in% colnames(emdia_groups_transposed_v_t_test))) {
  stop("Selected metabolite not found in the data.")
}

metabolite_data_only <- emdia_groups_transposed_v_t_test[, !(names(emdia_groups_transposed_v_t_test) %in% "Group")]

# Calculate correlations with the selected metabolite
correlation_results <- data.frame(
  Metabolite = character(),
  Correlation = numeric(),
  p.value = numeric(),
  stringsAsFactors = FALSE
)

for (met in colnames(metabolite_data_only)) {
  if (met != selected_metabolite) {
    test <- cor.test(
      metabolite_data_only[[selected_metabolite]],
      metabolite_data_only[[met]],
      method = "pearson"
    )
    correlation_results <- rbind(
      correlation_results,
      data.frame(
        Metabolite = met,
        Correlation = test$estimate,
        p.value = test$p.value
      )
    )
  }
}

# Adjust p-values
correlation_results <- correlation_results %>%
  mutate(
    p.adj = p.adjust(p.value, method = "fdr"),
    significance = case_when(
      p.adj < 0.0001 ~ "****",
      p.adj < 0.001 ~ "***",
      p.adj < 0.01 ~ "**",
      p.adj < 0.05 ~ "*",
      TRUE ~ ""
    )
  )


# Optional: Plot significant correlations with p-value legend
ggplot(filter(correlation_results, p.adj < 0.05), 
       aes(x = reorder(Metabolite, Correlation), y = Correlation)) +
  geom_bar(stat = "identity", fill = "#1f78b4") +
  geom_text(aes(label = significance), vjust = 0.8, size = 10) +
  coord_flip() +
  theme_minimal() +
  labs(
    title = paste("Correlations with", selected_metabolite),
    subtitle = "P-value significance: **** p < 0.0001 | *** p < 0.001 | ** p < 0.01 | * p < 0.05",
    y = "Pearson Correlation",
    x = "Metabolite"
  )+
  
  theme_minimal(base_family = "Arial") +
  theme(
  plot.title = element_text(size = 20, face = "bold"),
plot.subtitle = element_text(size = 20),
axis.title = element_text(size = 20),
axis.text = element_text(size = 20),
panel.grid.major.y = element_blank(),
panel.grid.minor = element_blank()
)

# Lollipop plot with more x-axis ticks
ggplot(filter(correlation_results, p.adj < 0.05), 
       aes(x = reorder(Metabolite, Correlation), y = Correlation)) +
  geom_segment(aes(xend = Metabolite, yend = 0), color = "skyblue", size = 1.2) +
  geom_point(color = "#1f78b4", size = 5) +
  geom_text(aes(label = significance), hjust = -0.4, size = 8) +
  coord_flip() +
  scale_y_continuous(
    breaks = seq(-1, 1, by = 0.1),  # You can use seq(-0.5, 0.5, by = 0.1) if you prefer
    limits = c(-0.5, 0.5)           # Adjust this range as needed
  ) +
  theme_minimal(base_family = "Arial") +
  labs(
    title = paste("Correlations with", selected_metabolite, "at V2-3 [ESI+]"),
    subtitle = "P-value significance: **** p < 0.0001 | *** p < 0.001 | ** p < 0.01 | * p < 0.05",
    y = "Pearson Correlation",
    x = "Metabolite"
  ) +
  theme(
    plot.title = element_text(size = 25, face = "bold"),
    plot.subtitle = element_text(size = 25),
    axis.title = element_text(size = 25),
    axis.text = element_text(size = 25),
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank()
  )





###
###Sex condition
###


#transpose emdia_data
emdia_data_transposed <- as.data.frame(t(emdia_data))
#set name of column
emdia_data_transposed <- setNames(data.frame(t(emdia_data[,-1])), emdia_data[,1])
#filter transposed dataframe for "Emdia" and "Placebo" group
emdia_groups_transposed_s <- filter(emdia_data_transposed, xor(Sex == "1", Sex == "2"))
#optionally filter transposed dataframe for Visit 1
#emdia_groups_transposed_s <- filter(emdia_groups_transposed_s, Visit == "V2")
#remove column containing sample information for calculations except group
emdia_groups_transposed_s_t_test <- emdia_groups_transposed_s[-c(1,2,3,5,6)]
emdia_groups_transposed_s_t_test <- as.data.frame(emdia_groups_transposed_s_t_test)
emdia_groups_transposed_s_t_test[, -1] <- sapply(emdia_groups_transposed_s_t_test[, -1], as.numeric)
#Define Group as factor
emdia_groups_transposed_s_t_test$Sex <- as.factor(emdia_groups_transposed_s_t_test$Sex)

# Perform t-tests
results_S <- data.frame(Metabolite = character(), p.value = numeric(), log2_fold_change = numeric())

for (metabolite in colnames(emdia_groups_transposed_s_t_test)[-1]) {
  t_test <- t.test(emdia_groups_transposed_s_t_test[[metabolite]] ~ emdia_groups_transposed_s_t_test$Sex)
  
  mean_group1 <- mean(emdia_groups_transposed_s_t_test[[metabolite]][emdia_groups_transposed_s_t_test$Sex == levels(emdia_groups_transposed_s_t_test$Sex)[1]])
  mean_group2 <- mean(emdia_groups_transposed_s_t_test[[metabolite]][emdia_groups_transposed_s_t_test$Sex == levels(emdia_groups_transposed_s_t_test$Sex)[2]])
  log2_fc <- log2(mean_group2 / mean_group1)
  
  results_S <- rbind(results_S, data.frame(Metabolite = metabolite, p.value = t_test$p.value, log2_fold_change = log2_fc))
}

# Adjust p-values. FDR correction: The method adjusts the p-values to account for multiple comparisons, reducing the likelihood of false positives while maintaining power in detecting true signals.

results_S <- results_S %>%
  mutate(
    neg_log10_p = -log10(p.value),
    p.adj = p.adjust(p.value, method = "fdr"),
    significance = case_when(
      p.adj < 0.0001 ~ "****",
      p.adj < 0.001 ~ "***",
      p.adj < 0.01 ~ "**",
      p.adj < 0.05 ~ "*",
      TRUE ~ ""
    )
  )

# Filter significant metabolites
significant_metabolites_S <- results_S %>%
  filter(p.adj < 0.05) %>%
  pull(Metabolite)


#(normalizing by max intensity per metabolite)
if (length(significant_metabolites_S) > 0) {
  
  # Normalize to max within each metabolite group
  data_long_S <- emdia_groups_transposed_s_t_test %>%
    pivot_longer(cols = -Sex, names_to = "Metabolite", values_to = "Intensity") %>%
    filter(Metabolite %in% significant_metabolites_S) %>%
    group_by(Metabolite) %>%
    mutate(Relative_Intensity = Intensity / max(Intensity) * 100) %>%
    left_join(results_S[, c("Metabolite", "significance")], by = "Metabolite") %>%
    ungroup()
  
  # Set fixed y-position for annotations
  annotation_y <- 105
  annotation_line_y <- 102
  
  # Prepare annotations for significance bars and stars
  annotation_df <- data_long_S %>%
    as_tibble() %>%
    distinct(Metabolite, significance) %>%
    mutate(x = as.numeric(factor(Metabolite)),
           xstart = x - 0.3,
           xend = x + 0.3,
           y = annotation_y,
           yline = annotation_line_y)
  
  # Plot
  ggplot(data_long_S, aes(x = Metabolite, y = Relative_Intensity, fill = Sex)) +
    geom_bar(stat = "summary", fun = mean, position = position_dodge(width = 0.9), width = 0.8) +
    
    # Line above bars
    geom_segment(data = annotation_df,
                 aes(x = xstart, xend = xend, y = yline, yend = yline),
                 inherit.aes = FALSE) +
    
    # Stars
    geom_text(data = annotation_df,
              aes(x = x, y = y, label = significance),
              inherit.aes = FALSE, vjust = 0, size = 5) +
    
    # Add star legend in top-right corner
    annotate("text", x = Inf, y = Inf, label = "**** p < 0.0001\n*** p < 0.001\n** p < 0.01\n* p < 0.05",
             hjust = 1, vjust = 3, size = 4, fontface = "italic") +
    
    labs(title = "Relative Intensity Barplot of Significant Metabolites Male vs Female",
         y = "Relative Intensity (%)", x = "Metabolite") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 10))
  
} else {
  print("No significant metabolites found (p.adj < 0.05).")
}


#log transform:

if (length(significant_metabolites_S) > 0) {
  
  # Apply log transformation to intensity values (log2)
  data_long_S <- emdia_groups_transposed_s_t_test %>%
    pivot_longer(cols = -Sex, names_to = "Metabolite", values_to = "Intensity") %>%
    filter(Metabolite %in% significant_metabolites_S) %>%
    group_by(Metabolite) %>%
    mutate(Log_Intensity = log2(Intensity + 1)) %>%  # Apply log2 transformation (adding 1 to avoid log(0))
    left_join(results_S[, c("Metabolite", "significance")], by = "Metabolite") %>%
    ungroup()
  
  # Set fixed y-position for annotations
  annotation_y <- 105
  annotation_line_y <- 102
  
  # Prepare annotations for significance bars and stars
  annotation_df <- data_long_S %>%
    as_tibble() %>%
    distinct(Metabolite, significance) %>%
    mutate(x = as.numeric(factor(Metabolite)),
           xstart = x - 0.3,
           xend = x + 0.3,
           y = annotation_y,
           yline = annotation_line_y)
  
  # Plot
  ggplot(data_long_S, aes(x = Metabolite, y = Log_Intensity, fill = Sex)) +
    geom_bar(stat = "summary", fun = mean, position = position_dodge(width = 0.9), width = 0.8) +
    
    # Line above bars
    geom_segment(data = annotation_df,
                 aes(x = xstart, xend = xend, y = yline, yend = yline),
                 inherit.aes = FALSE) +
    
    # Stars
    geom_text(data = annotation_df,
              aes(x = x, y = y, label = significance),
              inherit.aes = FALSE, vjust = 0, size = 5) +
    
    # Add star legend in top-right corner
    annotate("text", x = Inf, y = Inf, label = "**** p < 0.0001\n*** p < 0.001\n** p < 0.01\n* p < 0.05",
             hjust = 1, vjust = 3, size = 4, fontface = "italic") +
    
    labs(title = "Log2 Transformed Relative Intensity Barplot of Significant Metabolites Male vs Female",
         y = "Log2 Transformed Intensity", x = "Metabolite") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 10))
  
} else {
  print("No significant metabolites found (p.adj < 0.05).")
}


#Boxplot with Fold change
ggplot(results_S %>% filter(Metabolite %in% significant_metabolites_S), 
       aes(x = reorder(Metabolite, log2_fold_change), y = log2_fold_change, fill = significance)) +
  geom_bar(stat = "identity", show.legend = FALSE) +
  coord_flip() + # Flip the axes to make it easier to read
  labs(title = "Log2 Fold Change of Significant Metabolites Male vs Female",
       x = "Metabolite", y = "Log2 Fold Change") +
  theme_minimal() +
  theme(axis.text.y = element_text(size = 10))


#Dotplot of Foldchange
ggplot(results_S %>% filter(Metabolite %in% significant_metabolites_S), 
       aes(x = reorder(Metabolite, log2_fold_change), y = log2_fold_change, size = neg_log10_p, color = significance)) +
  geom_point() +
  coord_flip() +
  labs(title = "Dotplot of Fold Change and Significance Male vs Female",
       x = "Metabolite", y = "Log2 Fold Change") +
  scale_size_continuous(name = "Significance (-log10 p-value)") +
  theme_minimal() +
  theme(axis.text.y = element_text(size = 10))


###
###Sex dependent changes of empa treatment ANOVA:
###
### Sex-dependent changes in Empa treatment (V1 to V2)


# Extract metabolite data (excluding metadata columns)
emdia_data_t <- as.data.frame(t(emdia_data))
emdia_data_t <- setNames(data.frame(t(emdia_data[,-1])), emdia_data[,1])
emdia_data_t <- filter(emdia_data_t, class == "Sample")

# Prepare the metabolite data
metabolite_data <- emdia_data_t[, -(1:3)]  # removes unneeded batch info columns
metabolite_data$SampleID <- rownames(metabolite_data)

# Pivot to long format for metabolite data
metabolite_long <- metabolite_data[, -(1:3)] %>%
  pivot_longer(-SampleID, names_to = "Metabolite", values_to = "Intensity")

# Add metadata back from the original df
meta_info <- emdia_data_t[, c("Sex", "Visit", "Group")]
meta_info$SampleID <- rownames(meta_info)

# Join metadata with metabolite data
long_data <- left_join(metabolite_long, meta_info, by = "SampleID")

# Clean up factor levels and data types
long_data <- long_data %>%
  mutate(
    Sex = factor(Sex, levels = c("1", "2"), labels = c("Male", "Female")),
    Visit = factor(Visit),
    Group = factor(Group),
    Intensity = as.numeric(Intensity)
  )

# Create a new factor combining Sex, Visit, and Group
long_data <- long_data %>%
  mutate(Sex_Visit_Group = paste(Sex, Visit, Group, sep = "_"))

# Perform ANOVA for each metabolite
results_anova <- data.frame(Metabolite = character(), p.value = numeric(), stringsAsFactors = FALSE)

for (metabolite in unique(long_data$Metabolite)) {
  # Run ANOVA for each metabolite based on the combined Sex_Visit_Group factor
  model <- aov(Intensity ~ Sex_Visit_Group, data = subset(long_data, Metabolite == metabolite))
  
  # Get the ANOVA summary for the model
  anova_summary <- summary(model)
  
  # Extract the p-value for the main effect of Sex_Visit_Group
  p_value <- anova_summary[[1]]["Sex_Visit_Group", "Pr(>F)"]
  
  # Store the results
  results_anova <- rbind(results_anova, data.frame(Metabolite = metabolite, p.value = p_value))
}

# Adjust p-values using FDR correction (False Discovery Rate)
results_anova$p.adj <- p.adjust(results_anova$p.value, method = "fdr")

# Filter significant metabolites based on adjusted p-value
significant_metabolites <- results_anova %>%
  filter(p.adj < 0.05) %>%
  pull(Metabolite)

# Subset long_data to include only significant metabolites
long_data_significant <- long_data %>% filter(Metabolite %in% significant_metabolites)

if (length(significant_metabolites) > 0) {
  
  # Ensure levels exist before palette generation
  long_data_significant$Sex_Visit_Group <- factor(long_data_significant$Sex_Visit_Group)
  group_levels <- levels(long_data_significant$Sex_Visit_Group)
  n_colors <- length(group_levels)
  custom_colors_anova <- colorRampPalette(brewer.pal(8, "Set2"))(n_colors)
  names(custom_colors_anova) <- group_levels
  
  ggplot(long_data_significant, aes(x = Sex_Visit_Group, y = Intensity, fill = Sex_Visit_Group)) +
    geom_boxplot(outlier.shape = NA) +
    geom_jitter(aes(color = Sex_Visit_Group),
                alpha = 0.4,
                size = 1,
                position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.75)) +
    facet_wrap(~ Metabolite, scales = "free_y") +
    labs(title = "Faceted Boxplots of Absolute Intensities for Significant Metabolites",
         x = NULL, y = "Absolute Intensity") +  # set x label to NULL
    scale_fill_manual(values = custom_colors_anova) +
    scale_color_manual(values = custom_colors_anova) +
    theme_minimal() +
    theme(
      strip.text = element_text(face = "bold"),
      axis.text.x = element_blank(),      # remove x-axis text
      axis.title.x = element_blank(),     # remove x-axis title
      legend.title = element_blank(),
      legend.position = "right"
    )
} else {
  print("No significant metabolites found (p.adj < 0.05).")
}



###
###Sex dependent changes of empa treatment Barcharts:
###
### Sex-dependent changes in Empa treatment (V1 to V2 and V2 to V3)
###

library(dplyr)
library(tidyr)
library(ggplot2)

{

### 1. Initial Setup ----

# Extract metabolite data (excluding metadata columns)
emdia_data_t <- as.data.frame(t(emdia_data))
emdia_data_t <- setNames(data.frame(t(emdia_data[,-1])), emdia_data[,1])
emdia_data_t <- filter(emdia_data_t, class == "Sample")

metabolite_data <- emdia_data_t[, -(1:3)] #removes unnessecary batch info data
metabolite_data$SampleID <- rownames(metabolite_data)

# Pivot to long format
metabolite_long <- metabolite_data[, -(1:3)] %>% #removes cols Sex, Group and Visit for long format
  pivot_longer(-SampleID, names_to = "Metabolite", values_to = "Intensity")

# Add metadata back in from original df
meta_info <- emdia_data_t[, c("Sex", "Visit", "Group")]
meta_info$SampleID <- rownames(meta_info)

# Join metadata
long_data <- left_join(metabolite_long, meta_info, by = "SampleID")

# Clean up factor levels and data types
long_data <- long_data %>%
  mutate(
    Sex = factor(Sex, levels = c("1", "2"), labels = c("Male", "Female")),
    Visit = factor(Visit),
    Group = factor(Group),
    Intensity = as.numeric(Intensity)
  )

# Step 1: Calculate mean intensities and SD for each metabolite by group, sex, and visit
mean_sd_data <- long_data %>%
  group_by(Metabolite, Sex, Group, Visit) %>%
  summarise(
    mean_intensity = mean(Intensity, na.rm = TRUE),
    sd_intensity = sd(Intensity, na.rm = TRUE),
    .groups = "drop"
  )

# Step 2: Pivot the data to calculate log2 fold change between visits (V1 -> V2, V2 -> V3, etc.)
mean_sd_pivot <- mean_sd_data %>%
  pivot_wider(names_from = Visit, values_from = c(mean_intensity, sd_intensity)) %>%
  arrange(Metabolite, Sex, Group)

# Step 3: Calculate log2 fold change between visits for each metabolite, sex, and group
log2_fc_data <- mean_sd_pivot %>%
  mutate(
    log2_fc_V1_V2 = log2(mean_intensity_V2 / mean_intensity_V1),
    log2_fc_V1_V3 = log2(mean_intensity_V3 / mean_intensity_V1),
    log2_fc_V2_V3 = log2(mean_intensity_V3 / mean_intensity_V2)
  )

# Visualizing the log2 fold change for each metabolite, grouped by Sex and Visit
ggplot(log2_fc_data, aes(x = Metabolite, y = log2_fc_V1_V2, color = Sex, shape = Group)) +
  geom_point(size = 4, position = position_jitter(width = 0.2, height = 0)) +
  geom_line(aes(group = interaction(Sex, Group)), linetype = "dashed") +
  labs(
    title = "Log2 Fold Change from V1 to V2 for Metabolites",
    x = "Metabolite",
    y = "Log2 Fold Change (V1 to V2)",
    color = "Sex",
    shape = "Group"
  ) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))


# Step 1: Calculate p-values for each group, sex, and visit comparison
p_value_data <- long_data %>%
  filter(Visit %in% c("V1", "V2", "V3")) %>%
  group_by(Metabolite, Sex, Group) %>%
  summarise(
    p_value_V1_V2 = t.test(Intensity[Visit == "V1"], Intensity[Visit == "V2"])$p.value,
    p_value_V2_V3 = t.test(Intensity[Visit == "V2"], Intensity[Visit == "V3"])$p.value,
    p_value_V1_V3 = t.test(Intensity[Visit == "V1"], Intensity[Visit == "V3"])$p.value,
    .groups = "drop"
  )



# Step 2: Merge the p-values back into the log2 fold change data
log2_fc_data <- log2_fc_data %>%
  left_join(p_value_data, by = c("Metabolite", "Sex", "Group"))

# Step 3: Adjust p-values for multiple comparisons (using FDR)
log2_fc_data <- log2_fc_data %>%
  mutate(
    p_value_V1_V2_FDR = p.adjust(p_value_V1_V2, method = "fdr"),
    p_value_V1_V3_FDR = p.adjust(p_value_V1_V3, method = "fdr"),
    p_value_V2_V3_FDR = p.adjust(p_value_V2_V3, method = "fdr")
  )

# Adjust p-values using FDR correction and add significance levels
log2_fc_data <- log2_fc_data %>%
  mutate(
    # Compute -log10 p-value for visualization
    neg_log10_p_V1_V2 = -log10(p_value_V1_V2_FDR),
    neg_log10_p_V1_V3 = -log10(p_value_V1_V3_FDR),
    neg_log10_p_V2_V3 = -log10(p_value_V2_V3_FDR),
    
    # Add significance levels based on FDR-corrected p-values
    significance_V1_V2 = case_when(
      p_value_V1_V2_FDR < 0.0001 ~ "****",
      p_value_V1_V2_FDR < 0.001 ~ "***",
      p_value_V1_V2_FDR < 0.01 ~ "**",
      p_value_V1_V2_FDR < 0.05 ~ "*",
      TRUE ~ "No Significance"
    ),
    
    significance_V1_V3 = case_when(
      p_value_V1_V3_FDR < 0.0001 ~ "****",
      p_value_V1_V3_FDR < 0.001 ~ "***",
      p_value_V1_V3_FDR < 0.01 ~ "**",
      p_value_V1_V3_FDR < 0.05 ~ "*",
      TRUE ~ "No Significance"
    ),
    
    significance_V2_V3 = case_when(
      p_value_V2_V3_FDR < 0.0001 ~ "****",
      p_value_V2_V3_FDR < 0.001 ~ "***",
      p_value_V2_V3_FDR < 0.01 ~ "**",
      p_value_V2_V3_FDR < 0.05 ~ "*",
      TRUE ~ "No Significance"
    )
  )

# Identify metabolites with at least one significant result
significant_metabolites_V1_V2 <- log2_fc_data %>%
  filter(significance_V1_V2 != "No Significance") %>%
  pull(Metabolite) %>%
  unique()

significant_metabolites_V1_V3 <- log2_fc_data %>%
  filter(significance_V1_V3 != "No Significance") %>%
  pull(Metabolite) %>%
  unique()

significant_metabolites_V2_V3 <- log2_fc_data %>%
  filter(significance_V2_V3 != "No Significance") %>%
  pull(Metabolite) %>%
  unique()


# Keep all rows for those metabolites
significant_data_V1_V2 <- log2_fc_data %>%
  filter(Metabolite %in% significant_metabolites_V1_V2)

significant_data_V1_V3 <- log2_fc_data %>%
  filter(Metabolite %in% significant_metabolites_V1_V3)

significant_data_V2_V3 <- log2_fc_data %>%
  filter(Metabolite %in% significant_metabolites_V2_V3)

}
# Create cleaner labels for Sex & Group combinations
significant_data_V1_V2$SexGroup <- interaction(significant_data_V1_V2$Sex, significant_data_V1_V2$Group, sep = " ")

# Set a fixed star height
fixed_star_height_V1_V2 <- max(significant_data_V1_V2$log2_fc_V1_V2, na.rm = TRUE) + 0.5

ggplot(significant_data_V1_V2, aes(x = Metabolite, y = log2_fc_V1_V2, fill = SexGroup)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7) +
  
  # Add stars at fixed height
  geom_text(
    aes(y = fixed_star_height_V1_V2,
        label = ifelse(significance_V1_V2 != "No Significance", significance_V1_V2, "")),
    position = position_dodge(width = 0.8),
    vjust = 0,
    size = 3.5
  ) +
  
  scale_fill_brewer(palette = "Set2", name = "Sex & Group") +
  
  labs(
    title = "Significant Metabolites (V1 to V2)",
    subtitle = "P-value significance: **** p < 0.0001 | *** p < 0.001 | ** p < 0.01 | * p < 0.05",
    x = "Metabolite",
    y = "Log2 Fold Change (V1 to V2)",
    fill = "Sex & Group"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1),
    legend.position = "right",
    plot.subtitle = element_text(size = 10, face = "italic")
  )


significant_data_V1_V3$SexGroup <- interaction(significant_data_V1_V3$Sex, significant_data_V1_V3$Group, sep = " ")

fixed_star_height_V1_V3 <- max(significant_data_V1_V3$log2_fc_V1_V3, na.rm = TRUE) + 0.5

ggplot(significant_data_V1_V3, aes(x = Metabolite, y = log2_fc_V1_V3, fill = SexGroup)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7) +
  geom_text(
    aes(y = fixed_star_height_V1_V3,
        label = ifelse(significance_V1_V3 != "No Significance", significance_V1_V3, "")),
    position = position_dodge(width = 0.8),
    vjust = 0,
    size = 3.5
  ) +
  scale_fill_brewer(palette = "Set2", name = "Sex & Group") +
  labs(
    title = "Significant Metabolites (V1 to V3)",
    subtitle = "P-value significance: **** p < 0.0001 | *** p < 0.001 | ** p < 0.01 | * p < 0.05",
    x = "Metabolite",
    y = "Log2 Fold Change (V1 to V3)",
    fill = "Sex & Group"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1),
    legend.position = "right",
    plot.subtitle = element_text(size = 10, face = "italic")
  )



significant_data_V2_V3$SexGroup <- interaction(significant_data_V2_V3$Sex, significant_data_V2_V3$Group, sep = " ")

# Set a fixed star height
fixed_star_height_V2_V3 <- max(significant_data_V2_V3$log2_fc_V2_V3, na.rm = TRUE) + 0.5

ggplot(significant_data_V2_V3, aes(x = Metabolite, y = log2_fc_V2_V3, fill = SexGroup)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7) +
  
  # Add stars at fixed height
  geom_text(
    aes(y = fixed_star_height_V2_V3,
        label = ifelse(significance_V2_V3 != "No Significance", significance_V2_V3, "")),
    position = position_dodge(width = 0.8),
    vjust = 0,
    size = 3.5
  ) +
  
  scale_fill_brewer(palette = "Set2", name = "Sex & Group") +
  
  # Add p-value legend as caption
  labs(
    title = "Significant Metabolites (V2 to V3)",
    subtitle = "P-value significance: **** p < 0.0001 | *** p < 0.001 | ** p < 0.01 | * p < 0.05",
    x = "Metabolite",
    y = "Log2 Fold Change (V2 to V3)",
    fill = "Sex & Group"
  ) +
  
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1),
    legend.position = "right",
    plot.subtitle = element_text(size = 10, face = "italic")
  )





###
###Revision
###


{
###
###V2 Barplot Empagliflozin
###

#transpose emdia_data
emdia_data_transposed <- as.data.frame(t(emdia_data))
#set name of column
emdia_data_transposed <- setNames(data.frame(t(emdia_data[,-1])), emdia_data[,1])
#filter transposed dataframe for "Emdia" and "Placebo" group
emdia_groups_transposed <- filter(emdia_data_transposed, xor(Group == "Empa", Group == "Placebo"))
#filter transposed dataframe for Visit 1
emdia_groups_transposed_v2 <- filter(emdia_data_transposed, xor(Visit == "V2", Visit == "V3"))
#remove column containing sample information for calculations except group
emdia_groups_transposed_v2_t_test <- emdia_groups_transposed_v2[-c(1,2,3,4,6)]
emdia_groups_transposed_v2_t_test <- as.data.frame(emdia_groups_transposed_v2_t_test)
emdia_groups_transposed_v2_t_test[, c(2:136)] <- sapply(emdia_groups_transposed_v2_t_test[, c(2:136)], as.numeric)
#Define Group as factor
emdia_groups_transposed_v2_t_test$Group <- as.factor(emdia_groups_transposed_v2_t_test$Group)

#
# Volcanoplot of (log-transormed p-values vs. log2foldchange)
#

#Create new data.frame for metabolite+p-value+mean_diff
results_V2 <- data.frame(Metabolite = character(), p.value = numeric(), log2_fold_change = numeric(), stringsAsFactors = FALSE)

#Goes through data (Column 2 to n)
for (metabolite in colnames(emdia_groups_transposed_v2_t_test)[-1]) {
  #perform t-test
  t_test <- t.test(emdia_groups_transposed_v2_t_test[[metabolite]] ~ emdia_groups_transposed_v2_t_test$Group)
  
  # calculate log2(Fold Change)
  mean_group1 <- mean(emdia_groups_transposed_v2_t_test[[metabolite]][emdia_groups_transposed_v2_t_test$Group == levels(emdia_groups_transposed_v2_t_test$Group)[1]])
  mean_group2 <- mean(emdia_groups_transposed_v2_t_test[[metabolite]][emdia_groups_transposed_v2_t_test$Group == levels(emdia_groups_transposed_v2_t_test$Group)[2]])
  log2_fold_change <- log2(mean_group2 / mean_group1)
  
  #Save results in newly created data.frame
  results_V2 <- rbind(results_V2, data.frame(Metabolite = metabolite, p.value = t_test$p.value, log2_fold_change = log2_fold_change))
}

#Calculate neg log10 of p-values
results_V2 <- results_V2 %>%
  mutate(neg_log10_p = -log10(p.value))


#Visualization of the average peak area of a specific metabolite per group at V2.


# Define the specific metabolite you're interested in
selected_metabolite <- "EMPAGLIFLOZIN"  # Replace with the actual metabolite name

# Step 1: Reshape the data to long format
data_long_V2 <- emdia_groups_transposed_v2_t_test %>%
  pivot_longer(cols = -Group, names_to = "Metabolite", values_to = "Intensity")

# Step 2: Filter the data for the specific metabolite
filtered_data_V2 <- data_long_V2 %>%
  filter(Metabolite == selected_metabolite) %>%
  group_by(Group) %>%
  summarize(
    mean_intensity = mean(Intensity, na.rm = TRUE),
    sd_intensity = sd(Intensity, na.rm = TRUE),  # Standard deviation
    n = n(),  # Count of samples in each group
    se_intensity = sd_intensity / sqrt(n)  # Standard error
  )



###
###V3 Barplot Empagliflozin
###


###
###V3
###

#transpose emdia_data
emdia_data_transposed <- as.data.frame(t(emdia_data))
#set name of column
emdia_data_transposed <- setNames(data.frame(t(emdia_data[,-1])), emdia_data[,1])
#filter transposed dataframe for "Emdia" and "Placebo" group
emdia_groups_transposed <- filter(emdia_data_transposed, xor(Group == "Empa", Group == "Placebo"))
#filter transposed dataframe for Visit 1
emdia_groups_transposed_v3 <- filter(emdia_data_transposed, Visit == "V3")
#remove column containing sample information for calculations except group
emdia_groups_transposed_v3_t_test <- emdia_groups_transposed_v3[-c(1,2,3,4,6)]
emdia_groups_transposed_v3_t_test <- as.data.frame(emdia_groups_transposed_v3_t_test)
emdia_groups_transposed_v3_t_test[, c(2:136)] <- sapply(emdia_groups_transposed_v3_t_test[, c(2:136)], as.numeric)
#Define Group as factor
emdia_groups_transposed_v3_t_test$Group <- as.factor(emdia_groups_transposed_v3_t_test$Group)

#
# Volcanoplot of (log-transormed p-values vs. log2foldchange)
#

#Create new data.frame for metabolite+p-value+mean_diff
results_V3 <- data.frame(Metabolite = character(), p.value = numeric(), log2_fold_change = numeric(), stringsAsFactors = FALSE)

#Goes through data (Column 2 to n)
for (metabolite in colnames(emdia_groups_transposed_v3_t_test)[-1]) {
  #perform t-test
  t_test <- t.test(emdia_groups_transposed_v3_t_test[[metabolite]] ~ emdia_groups_transposed_v3_t_test$Group)
  
  # calculate log2(Fold Change)
  mean_group1 <- mean(emdia_groups_transposed_v3_t_test[[metabolite]][emdia_groups_transposed_v3_t_test$Group == levels(emdia_groups_transposed_v3_t_test$Group)[1]])
  mean_group2 <- mean(emdia_groups_transposed_v3_t_test[[metabolite]][emdia_groups_transposed_v3_t_test$Group == levels(emdia_groups_transposed_v3_t_test$Group)[2]])
  log2_fold_change <- log2(mean_group2 / mean_group1)
  
  #Save results in newly created data.frame
  results_V3 <- rbind(results_V3, data.frame(Metabolite = metabolite, p.value = t_test$p.value, log2_fold_change = log2_fold_change))
}

#Calculate neg log10 of p-values
results_V3 <- results_V3 %>%
  mutate(neg_log10_p = -log10(p.value))


#Visualization of the average peak area of a specific metabolite per group at V3.


# Define the specific metabolite you're interested in
selected_metabolite <- "EMPAGLIFLOZIN"  # Replace with the actual metabolite name

# Step 1: Reshape the data to long format
data_long_V3 <- emdia_groups_transposed_v3_t_test %>%
  pivot_longer(cols = -Group, names_to = "Metabolite", values_to = "Intensity")

# Step 2: Filter the data for the specific metabolite
filtered_data_V3 <- data_long_V3 %>%
  filter(Metabolite == selected_metabolite) %>%
  group_by(Group) %>%
  summarize(
    mean_intensity = mean(Intensity, na.rm = TRUE),
    sd_intensity = sd(Intensity, na.rm = TRUE),  # Standard deviation
    n = n(),  # Count of samples in each group
    se_intensity = sd_intensity / sqrt(n)  # Standard error
  )


# After each ggplot call, add ggsave:


# Define save path depending on OS
if (.Platform$OS.type == "unix") {
  save_path <- file.path("filepath_macOS")
} else {
  save_path <- file.path("filepath_windows", "R images")
}

# After each ggplot call, add ggsave:

# --- For V2 ---
p_V2 <- ggplot(filtered_data_V2, aes(x = Group, y = mean_intensity, fill = Group)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8),
           width = 0.7, color = "black", size = 0.5) +   # thinner outline
  geom_errorbar(aes(ymin = mean_intensity - se_intensity, ymax = mean_intensity + se_intensity),
                position = position_dodge(width = 0.8), width = 0.25, size = 0.5) +
  theme_minimal() +
  scale_fill_manual(values = c("#1a80bb", "#ea801c")) +
  labs(title = paste("Average Intensity of", selected_metabolite, "(V2)"),
       x = "Group", y = "Average Intensity") +
  theme(
    axis.text.x = element_text(size = 10),
    axis.text.y = element_text(size = 10),
    axis.title.x = element_text(size = 10),
    axis.title.y = element_text(size = 10),
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 10),
    legend.position = "right",
    plot.title = element_text(size = 12, hjust = 0.3)
  )

ggsave(filename = file.path(save_path, paste0("Barplot_", selected_metabolite, "_V2.png")),
       plot = p_V2, width = 3.5, height = 3.5, units = "in", dpi = 300)


# --- For V3 ---
p_V3 <- ggplot(filtered_data_V3, aes(x = Group, y = mean_intensity, fill = Group)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8),
           width = 0.7, color = "black", size = 0.5) +   # thinner outline
  geom_errorbar(aes(ymin = mean_intensity - se_intensity, ymax = mean_intensity + se_intensity),
                position = position_dodge(width = 0.8), width = 0.25, size = 0.5) +
  theme_minimal() +
  scale_fill_manual(values = c("#1a80bb", "#ea801c")) +
  labs(title = paste("Average Intensity of", selected_metabolite, "(V3)"),
       x = "Group", y = "Average Intensity") +
  theme(
    axis.text.x = element_text(size = 10),
    axis.text.y = element_text(size = 10),
    axis.title.x = element_text(size = 10),
    axis.title.y = element_text(size = 10),
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 10),
    legend.position = "right",
    plot.title = element_text(size = 12, hjust = 0.3)
  )

ggsave(filename = file.path(save_path, paste0("Barplot_", selected_metabolite, "_V3.png")),
       plot = p_V3, width = 3.5, height = 3.5, units = "in", dpi = 300)

  }
