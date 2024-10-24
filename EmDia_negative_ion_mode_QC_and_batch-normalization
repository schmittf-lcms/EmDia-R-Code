#####
##### Raw data import and data rearrangement
#####


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


#data import - MAC
emdia_neg <- "/Users/filepath_to_csv/Emdia MS1 RAW data neg.csv" 
emdia_neg <- read.csv2(emdia_neg, header = TRUE, dec = ".", sep = ";")
head(emdia_neg)

#data import - Windows
emdia_neg <- "C:/Users/filepath_to_csv/Emdia MS1 RAW data neg.csv" 
emdia_neg <- read.csv2(emdia_neg, header = TRUE, dec = ".", sep = ",")
head(emdia_neg)

###
###
###PCA of non-corrected data
###
###

emdia_transposed <- as.data.frame(t(emdia_neg))
#set name of column
emdia_transposed <- setNames(data.frame(t(emdia_neg[,-1])), emdia_neg[,1])
#filter for QC and Sample
emdia_transposed <- filter(emdia_transposed, xor(class == "QC", class == "Sample"))
#removes index, batch and class information
emdia_transposed_values <- emdia_transposed[-c(1:3)]
emdia_transposed_values <- as.data.frame(apply(emdia_transposed_values, 2, as.numeric))
#calculates pc
emdia_transposed_values_pca <- prcomp(emdia_transposed_values, scale = T)
#plots pca
autoplot(emdia_transposed_values_pca, data = emdia_transposed, colour = "batch", shape = "class", size = 4)+
  ggtitle("PCA of non-corrected EmDia (negative ion mode)")+
  theme(axis.text.x = element_text(size=12, angle = 45, vjust = 1, hjust=1), legend.position = "right") +
  theme(legend.text=element_text(size=20))+
  theme(panel.background = element_rect(fill = "white"),
        panel.grid.major = element_line(color = "grey", linetype = "dotted", linewidth = 0.75))+
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=20,face="bold"))+
  theme(title=element_text(size=20,face="bold"))


###
###
###CVs and Violin Plots of uncorrected GHS data
###
###

#transpose emdia_data
emdia_transposed <- as.data.frame(t(emdia_neg))
#set name of column
emdia_transposed <- setNames(data.frame(t(emdia_neg[,-1])), emdia_neg[,1])
emdia_transposed <- filter(emdia_transposed, xor(class == "QC", class == "Sample"))
emdia_transposed_CV <- emdia_transposed[-c(2)]

# Convert row names to a column named "SampleID"
emdia_transposed_CV <- emdia_transposed_CV %>%
  rownames_to_column(var = "SampleID")

# Ensure the data contains the necessary columns
if (!all(c("SampleID", "batch", "class") %in% colnames(emdia_transposed_CV))) {
  stop("Data must contain SampleID, Batch, and Class columns")
}

# Reshape the data to long format for easier manipulation
data_long <- pivot_longer(emdia_transposed_CV, cols = -c(SampleID, batch, class), names_to = "Metabolite", values_to = "Value")
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

# Convert batch column to a factor and ensure the levels are in the correct numeric order
cv_data$batch <- factor(cv_data$batch, levels = sort(as.numeric(unique(cv_data$batch))))

# Visualize the CVs using violin plots with boxplots and points for outliers
ggplot(cv_data, aes(x = class, y = CV, fill = batch)) +
  ylim(0,100)+
  geom_boxplot(width = 0.1, outlier.shape = NA, position = position_dodge(width = 0.9)) +
  geom_hline(yintercept = 6, linetype = "dotted", color = "red") +
  annotate("text", x = Inf, y = 5, label = "6%", vjust = -1, hjust = 1.05, color = "red", size = 5) +
  #geom_point(data = filter(cv_data, is_outlier), aes(x = class, y = CV, color = batch), 
  #position = position_jitter(width = 0.1), size = 1.5) +
  scale_color_manual(values = c("Treatment" = "red", "Placebo" = "blue")) +
  labs(
    title = "Coefficient of Variation of Metabolites by Batch and Class non-corrected EmDia data (negative ion mode)",
    x = "Class",
    y = "Coefficient of Variation (%)"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(size = 20, angle = 45, hjust = 1),
    axis.text.y = element_text(size = 20),
    axis.title.x = element_text(size = 20),
    axis.title.y = element_text(size = 20),
    legend.title = element_text(size = 20),
    legend.text = element_text(size = 20),
    plot.title = element_text(size = 20, hjust = 0.5)
  )

# Print the CV data
print(cv_data)


###
###
###Batch correction of skyline emdia data export (manual curation of index, class, batch in the csv rows required!!!)
###
####

#transform emdia data set
emdia_t <- t(emdia_neg) 
#set row and column names
emdia_t= setNames(data.frame(t(emdia_neg[,-1])), emdia_neg[,1])
#filter for QC and Samples
emdia_t_f <- filter(emdia_t, (class == "QC" | class == "Sample"))
#classifies the batch information for QC-RSC
batch <- emdia_t_f$batch
#classifies the class information for QC-RSC
class <- emdia_t_f$class
#classifies the AUC value information for QC-RSC
data <- emdia_t_f[-c(1:3)]
data <- as.data.frame(apply(data, 2, as.numeric))

sample_order <- c(1:nrow(emdia_t_f))
#data filtering for QC class, as sample normalization is performed in relation to averaged QC intensities
data_filtered <- filter_peaks_by_fraction(df=data, classes=class, method="QC",
                                          qc_label="QC", min_frac=0.8)
#data normalization using the QC-RSC algorrithm
corrected_data <- QCRSC(df=data_filtered, order=sample_order, batch=batch, 
                        classes=class, spar=0, minQC=4)
#creates plots the normalization on a per metabolite basis, indexes=c() (for all) can be changed to a specific metabolite e.g. indexes=c(1))
plots <- sbc_plot (df=data, corrected_df=corrected_data, classes=class, 
                   batch=batch, output=NULL, indexes=c(1))
#plots created plots individually
plots

#plots all wanted plots at once in one figure
autoplot(plots)+
  geom_point(size = 2)+
  theme(axis.text.x = element_text(size=12, angle = 45, vjust = 1, hjust=1), legend.position = "right") +
  theme(legend.text=element_text(size=20))+
  theme(panel.background = element_rect(fill = "white"),
        panel.grid.major = element_line(color = "grey", linetype = "dotted", linewidth = 0.75))+
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=20,face="bold"))+
  theme(title=element_text(size=20,face="bold"))

#replaces possible NAs with 0
corrected_data[is.na(corrected_data)] <- 0
corrected_data_t <- t(corrected_data)

#creates putput df format for csv
emdia_header <- emdia_t_f[c(1:3)]
emdia_corrected_data_t <- cbind(emdia_header, corrected_data_t)
emdia_corrected_data <- t(emdia_corrected_data_t)

write.table(emdia_corrected_data, 
            file = "C:/Users/filepath_to_csv/emdia_neg_corrected.csv", 
            sep = ",", 
            dec = ".", 
            row.names = TRUE, 
            col.names = TRUE)
###
###
###PCA of corrected data
###
###

emdia_corrected_data <- "C:/Users/filepath_to_csv/emdia_neg_corrected.csv" 
emdia_corrected_data <- read.csv2(emdia_corrected_data, header = TRUE, dec = ".", sep = ",")
emdia_corrected_data_transposed <- as.data.frame(t(emdia_corrected_data))
emdia_corrected_data_transposed <- filter(emdia_corrected_data_transposed, xor(class == "QC", class == "Sample"))
emdia_corrected_data_transposed_values <- emdia_corrected_data_transposed[-c(1:3)]
emdia_corrected_data_transposed_values <- as.data.frame(apply(emdia_corrected_data_transposed_values, 2, as.numeric))

pca_res <- prcomp(emdia_corrected_data_transposed_values, scale = T)

autoplot(pca_res, data = emdia_corrected_data_transposed, colour = "batch", shape = "class", size = 4)+
  ggtitle("PCA of batch-corrected EmDia data (negative ion mode)")+
  theme(axis.text.x = element_text(size=12, angle = 45, vjust = 1, hjust=1), legend.position = "right") +
  theme(legend.text=element_text(size=20))+
  theme(panel.background = element_rect(fill = "white"),
        panel.grid.major = element_line(color = "grey", linetype = "dotted", linewidth = 0.75))+
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=20,face="bold"))+
  theme(title=element_text(size=20,face="bold"))


###
###
###CVs and Violin Plots of batch-corrected EmDia data
###
###

emdia_corrected_data <- "C:/Users/filepath_to_csv/emdia_neg_corrected.csv" 
emdia_corrected_data <- read.csv2(emdia_corrected_data, header = TRUE, dec = ".", sep = ",")
emdia_corrected_data_transposed <- as.data.frame(t(emdia_corrected_data))
emdia_corrected_data_transposed <- filter(emdia_corrected_data_transposed, xor(class == "QC", class == "Sample"))
emdia_corrected_data_transposed_CV <- emdia_corrected_data_transposed[-c(2)]

# Convert row names to a column named "SampleID"
emdia_corrected_data_transposed_CV <- emdia_corrected_data_transposed_CV %>%
  rownames_to_column(var = "SampleID")

# Ensure the data contains the necessary columns
if (!all(c("SampleID", "batch", "class") %in% colnames(emdia_corrected_data_transposed_CV))) {
  stop("Data must contain SampleID, Batch, and Class columns")
}

# Reshape the data to long format for easier manipulation
data_long <- pivot_longer(emdia_corrected_data_transposed_CV, cols = -c(SampleID, batch, class), names_to = "Metabolite", values_to = "Value")
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

# Convert batch column to a factor and ensure the levels are in the correct numeric order
cv_data$batch <- factor(cv_data$batch, levels = sort(as.numeric(unique(cv_data$batch))))

# Visualize the CVs using violin plots with boxplots and points for outliers
ggplot(cv_data, aes(x = class, y = CV, fill = batch)) +
  ylim(0,100)+
  geom_boxplot(width = 0.1, outlier.shape = NA, position = position_dodge(width = 0.9)) +
  geom_hline(yintercept = 6, linetype = "dotted", color = "red") +
  annotate("text", x = Inf, y = 5, label = "6%", vjust = -1, hjust = 1.05, color = "red", size = 5) +
  #geom_point(data = filter(cv_data, is_outlier), aes(x = class, y = CV, color = batch), 
  #position = position_jitter(width = 0.1), size = 1.5) +
  scale_color_manual(values = c("Treatment" = "red", "Placebo" = "blue")) +
  labs(
    title = "Coefficient of Variation of Metabolites by Batch and Class batch-corrected EmDia data (negative ion mode)",
    x = "Class",
    y = "Coefficient of Variation (%)"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(size = 20, angle = 45, hjust = 1),
    axis.text.y = element_text(size = 20),
    axis.title.x = element_text(size = 20),
    axis.title.y = element_text(size = 20),
    legend.title = element_text(size = 20),
    legend.text = element_text(size = 20),
    plot.title = element_text(size = 20, hjust = 0.5)
  )

# Print the CV data
print(cv_data)



