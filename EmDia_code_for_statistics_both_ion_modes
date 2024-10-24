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
library("mdatools")
library("ropls")
library("mixOmics")
library("ggrepel")
library("pheatmap")
library("RColorBrewer")
library("scales")
library("ComplexHeatmap")
library("colorRamp2")
library("circlize")
library("grid")


#data import positive corrected data
emdia_data <- "C:/Users/filepath_to_csv/emdia_corrected_med_data.csv" 
emdia_data <- read.csv2(emdia_data, header = TRUE, dec = ".", sep = ",")
head(emdia_data)

#data import negative corrected data
emdia_data_neg <- "C:/Users/filepath_to_csv/emdia_neg_corrected_med_data.csv" 
emdia_data_neg <- read.csv2(emdia_data_neg, header = TRUE, dec = ".", sep = ",")
head(emdia_data_neg)


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


# Visualize the CVs using violin plots with boxplots and points for outliers
ggplot(cv_data, aes(x = class, y = CV, fill = batch)) +
  ylim(0,100)+
  geom_violin(trim = FALSE) +
  geom_boxplot(width = 0.1, outlier.shape = NA, position = position_dodge(width = 0.9)) +
  geom_hline(yintercept = 6, linetype = "dotted", color = "red") +
  annotate("text", x = Inf, y = 5, label = "6%", vjust = -1, hjust = 1.05, color = "red", size = 5) +
  #geom_point(data = filter(cv_data, is_outlier), aes(x = class, y = CV, color = batch), 
             #position = position_jitter(width = 0.1), size = 1.5) +
  scale_color_manual(values = c("Treatment" = "red", "Placebo" = "blue")) +
  labs(
    title = "Coefficient of Variation of Metabolites by Batch and Class (positive ion mode)",
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
###CVs of corrected negative ion mode data
###


#transpose emdia_data
emdia_data_neg_transposed <- as.data.frame(t(emdia_data_neg))
#set name of column
emdia_data_neg_transposed <- setNames(data.frame(t(emdia_data_neg[,-1])), emdia_data_neg[,1])
emdia_data_neg_transposed <- filter(emdia_data_neg_transposed, xor(class == "QC", class == "Sample"))
emdia_data_neg_transposed_CV <- emdia_data_neg_transposed[-c(2,4:6)]

# Convert row names to a column named "SampleID"
emdia_data_neg_transposed_CV <- emdia_data_neg_transposed_CV %>%
  rownames_to_column(var = "SampleID")

# Ensure the data contains the necessary columns
if (!all(c("SampleID", "batch", "class") %in% colnames(emdia_data_neg_transposed_CV))) {
  stop("Data must contain SampleID, Batch, and Class columns")
}

# Reshape the data to long format for easier manipulation
data_long <- pivot_longer(emdia_data_neg_transposed_CV, cols = -c(SampleID, batch, class), names_to = "Metabolite", values_to = "Value")
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

# Visualize the CVs using violin plots with boxplots and points for outliers
ggplot(cv_data, aes(x = class, y = CV, fill = batch)) +
  ylim(0,100)+
  geom_violin(trim = FALSE) +
  geom_boxplot(width = 0.1, outlier.shape = NA, position = position_dodge(width = 0.9)) +
  geom_hline(yintercept = 6, linetype = "dotted", color = "red") +
  annotate("text", x = Inf, y = 5, label = "6%", vjust = -1, hjust = 1.05, color = "red", size = 5) +
  #geom_point(data = filter(cv_data, is_outlier), aes(x = class, y = CV, color = batch), 
  #position = position_jitter(width = 0.1), size = 1.5) +
  scale_color_manual(values = c("Treatment" = "red", "Placebo" = "blue")) +
  labs(
    title = "Coefficient of Variation of Metabolites by Batch and Class (negative ion mode)",
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
###PCA of positive ion mode emdia data
####
###


#transpose emdia_data
emdia_data_transposed <- as.data.frame(t(emdia_data))
#set name of column
emdia_data_transposed <- setNames(data.frame(t(emdia_data[,-1])), emdia_data[,1])
#emdia_data_transposed <- filter(emdia_data_transposed, (Visit == "V3"))

#remove sample information columns
emdia_data_transposed_values <- emdia_data_transposed[-c(1:6)]
#set values as numeric for pca analaysis
emdia_data_transposed_values <- as.data.frame(apply(emdia_data_transposed_values, 2, as.numeric))

#calculate pca of transposed emdia data
pca_emdia <- prcomp(emdia_data_transposed_values, scale = T)
#plot the pca
autoplot(pca_emdia, data = emdia_data_transposed, colour = "batch", shape = "class")+
  ggtitle("PCA of batch-corrected data") +
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

###
###PCA of negative ion mode emdia data
####


#transpose emdia_data
emdia_data_neg_transposed <- as.data.frame(t(emdia_data_neg))
#set name of column
emdia_data_neg_transposed <- setNames(data.frame(t(emdia_data_neg[,-1])), emdia_data_neg[,1])

#emdia_data_neg_transposed <- filter(emdia_data_neg_transposed, (Visit == "V2"))

#remove sample information columns
emdia_data_neg_transposed_values <- emdia_data_neg_transposed[-c(1:6)]
#set values as numeric for pca analaysis
emdia_data_neg_transposed_values <- as.data.frame(apply(emdia_data_neg_transposed_values, 2, as.numeric))

#calculate pca of transposed emdia data
pca_emdia_neg <- prcomp(emdia_data_neg_transposed_values, scale = T)
#plot the pca
autoplot(pca_emdia_neg, data = emdia_data_neg_transposed, colour = "batch", shape = "class")+
  ggtitle("PCA of negative batch-corrected data") +
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



###
###
###PLS-DA of positive data (visit condition)
###
###

#
#All emdia timepoints (V1-3)
#

#transpose emdia_data
emdia_data_transposed <- as.data.frame(t(emdia_data))
#set name of column
emdia_data_transposed <- setNames(data.frame(t(emdia_data[,-1])), emdia_data[,1])
#filter transposed dataframe for "Emdia" and "Placebo" group
emdia_groups_transposed <- filter(emdia_data_transposed, xor(Group == "Empa", Group == "Placebo"))
#remove column containing sample information for calculations
emdia_groups_transposed_values <- emdia_groups_transposed[-c(1:6)]
#define values as numeric for calculations
emdia_groups_transposed_values <- as.data.frame(apply(emdia_groups_transposed_values, 2, as.numeric))
#define variable für PLS-DA in this case the group info of each sample
emdia_groupinfo <- as.character(as.matrix(emdia_groups_transposed[c(5)]))

#PLSDA with ropls
emdia.plsda <- opls(emdia_groups_transposed_values, emdia_groupinfo)

#PLSDA with mixomics
#calculate distance
plsda_result <- plsda(emdia_groups_transposed_values, emdia_groupinfo, ncomp = 2)
#plot PLS-DA
plotIndiv(plsda_result, ind.names = FALSE, legend=TRUE,
          comp=c(1,2), ellipse = TRUE, 
          title = 'PLS-DA EmDia V1-3',
          X.label = 'PLS-DA comp 1', Y.label = 'PLS-DA comp 2',
          size.title = 20, size.xlabel = 20, size.ylabel = 20, size.axis = 20, size.legend = 20,
          size.legend.title = 20)
#plot correlation circle plot
plotVar(plsda_result, title = "PLS Correlations - cutoff of 0.5", 
        cutoff = 0.5, legend = F, style = 'ggplot2')

#plot loading for comp 1
plotLoadings(plsda_result, comp = 1, method = 'mean', contrib = 'max')
#plot loading for comp 2
plotLoadings(plsda_result, comp = 2, method = 'mean', contrib = 'max')


#PLSDA for sex dependency in V1-3
emdia_sexinfo <- as.character(as.matrix(emdia_groups_transposed[c(4)]))

plsda_sex_result <- plsda(emdia_groups_transposed_values, emdia_sexinfo, ncomp = 2)
plotIndiv(plsda_sex_result, ind.names = FALSE, legend=TRUE,
          comp=c(1,2), ellipse = TRUE, 
          title = 'PLS-DA EmDia Influnce of Sex V1-3',
          X.label = 'PLS-DA comp 1', Y.label = 'PLS-DA comp 2',
size.title = 20, size.xlabel = 20, size.ylabel = 20, size.axis = 20, size.legend = 20,
size.legend.title = 20)

plotVar(plsda_sex_result)
plotLoadings(plsda_sex_result, comp = 1, method = 'mean', contrib = 'max')
plotLoadings(plsda_sex_result, comp = 2, method = 'mean', contrib = 'max')


#
#Timepoint V1
#


custom_colors <- c("#1a80bb", "#ea801c")

#transpose emdia_data
emdia_data_transposed <- as.data.frame(t(emdia_data))
#set name of column
emdia_data_transposed <- setNames(data.frame(t(emdia_data[,-1])), emdia_data[,1])
emdia_groups_transposed_V1 <- filter(emdia_data_transposed, (Visit == "V1"))
emdia_groups_transposed_V1_values <- emdia_groups_transposed_V1[-c(1:6)]
emdia_groups_transposed_V1_values <- as.data.frame(apply(emdia_groups_transposed_V1_values, 2, as.numeric))
emdia_groupinfo_V1 <- as.character(as.matrix(emdia_groups_transposed_V1[c(5)]))


#PLSDA with mixomics
plsda_result_V1 <- plsda(emdia_groups_transposed_V1_values, emdia_groupinfo_V1, ncomp = 2)
plotIndiv(plsda_result_V1, col = custom_colors, ind.names = FALSE, legend=TRUE,
          comp=c(1,2), ellipse = TRUE, 
          title = 'PLS-DA EmDia V1',
          X.label = 'PLS-DA comp 1', Y.label = 'PLS-DA comp 2',
          size.title = 30, size.xlabel = 30, size.ylabel = 30, size.axis = 30, size.legend = 30,
          size.legend.title = 30,
          pch = 16)
plotVar(plsda_result_V1, col = custom_colors)
plotLoadings(plsda_result_V1, comp = 1, method = 'mean', contrib = 'max')
plotLoadings(plsda_result_V1, comp = 2, method = 'mean', contrib = 'max')

#PLSDA for sex dependency in V1
emdia_sexinfo_V1 <- as.character(as.matrix(emdia_groups_transposed_V1[c(4)]))

plsda_sex_result_V1 <- plsda(emdia_groups_transposed_V1_values, emdia_sexinfo_V1, ncomp = 2)
plotIndiv(plsda_sex_result_V1, ind.names = FALSE, legend=TRUE,
          comp=c(1,2), ellipse = TRUE, 
          title = 'PLS-DA EmDia Influnce of Sex V1',
          X.label = 'PLS-DA comp 1', Y.label = 'PLS-DA comp 2',
          size.title = 20, size.xlabel = 20, size.ylabel = 20, size.axis = 20, size.legend = 20,
          size.legend.title = 20)
plotVar(plsda_sex_result_V1)
plotLoadings(plsda_sex_result_V1, comp = 1, method = 'mean', contrib = 'max')
plotLoadings(plsda_sex_result_V1, comp = 2, method = 'mean', contrib = 'max')


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


#PLSDA with mixomics
plsda_result_V2 <- plsda(emdia_groups_transposed_V2_values, emdia_groupinfo_V2, ncomp = 2)
plotIndiv(plsda_result_V2, col = custom_colors, ind.names = FALSE, legend=TRUE,
          comp=c(1,2), ellipse = TRUE, 
          title = 'PLS-DA EmDia V2',
          X.label = 'PLS-DA comp 1', Y.label = 'PLS-DA comp 2',
          size.title = 30, size.xlabel = 30, size.ylabel = 30, size.axis = 30, size.legend = 30,
          size.legend.title = 20,
          pch = 16)
plotVar(plsda_result_V2)
plotLoadings(plsda_result_V2, comp = 1, method = 'mean', contrib = 'max')
plotLoadings(plsda_result_V2, comp = 2, method = 'mean', contrib = 'max')

#PLSDA for sex dependency in V2
emdia_sexinfo_V2 <- as.character(as.matrix(emdia_groups_transposed_V2[c(4)]))

plsda_sex_result_V2 <- plsda(emdia_groups_transposed_V2_values, emdia_sexinfo_V2, ncomp = 2)
plotIndiv(plsda_sex_result_V2, ind.names = FALSE, legend=TRUE,
          comp=c(1,2), ellipse = TRUE, 
          title = 'PLS-DA EmDia Influnce of Sex V2',
          X.label = 'PLS-DA comp 1', Y.label = 'PLS-DA comp 2',
          size.title = 20, size.xlabel = 20, size.ylabel = 20, size.axis = 20, size.legend = 20,
          size.legend.title = 20)
plotVar(plsda_sex_result_V2)
plotLoadings(plsda_sex_result_V2, comp = 1, method = 'mean', contrib = 'max')
plotLoadings(plsda_sex_result_V2, comp = 2, method = 'mean', contrib = 'max')

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


#PLSDA with mixomics
plsda_result_V3 <- plsda(emdia_groups_transposed_V3_values, emdia_groupinfo_V3, ncomp = 2)
plotIndiv(plsda_result_V3, col = custom_colors, ind.names = FALSE, legend=TRUE,
          comp=c(1,2), ellipse = TRUE, 
          title = 'PLS-DA EmDia V3',
          X.label = 'PLS-DA comp 1', Y.label = 'PLS-DA comp 2',
          size.title = 30, size.xlabel = 30, size.ylabel = 30, size.axis = 30, size.legend = 30,
          size.legend.title = 20,
          pch = 16)
plotVar(plsda_result_V3)
plotLoadings(plsda_result_V3, comp = 1, method = 'mean', contrib = 'max')
plotLoadings(plsda_result_V3, comp = 2, method = 'mean', contrib = 'max')

#PLSDA for sex dependency in V3
emdia_sexinfo_V3 <- as.character(as.matrix(emdia_groups_transposed_V3[c(4)]))

plsda_sex_result_V3 <- plsda(emdia_groups_transposed_V3_values, emdia_sexinfo_V3, ncomp = 2)
plotIndiv(plsda_sex_result_V3, ind.names = FALSE, legend=TRUE,
          comp=c(1,2), ellipse = TRUE, 
          title = 'PLS-DA EmDia Influnce of Sex V3',
          X.label = 'PLS-DA comp 1', Y.label = 'PLS-DA comp 2',
          size.title = 20, size.xlabel = 20, size.ylabel = 20, size.axis = 20, size.legend = 20,
          size.legend.title = 20)
plotVar(plsda_sex_result_V3)
plotLoadings(plsda_sex_result_V3, comp = 1, method = 'mean', contrib = 'max')
plotLoadings(plsda_sex_result_V3, comp = 2, method = 'mean', contrib = 'max')


###
###
###PLs-DA of negative ion mode data (visit condition)
###
###

#
#All emdia timepoints (V1-3)
#

#transpose emdia_data
emdia_data_neg_transposed <- as.data.frame(t(emdia_data_neg))
#set name of column
emdia_data_neg_transposed <- setNames(data.frame(t(emdia_data_neg[,-1])), emdia_data_neg[,1])
emdia_neg_groups_transposed <- filter(emdia_data_neg_transposed, xor(Group == "Empa", Group == "Placebo"))
emdia_neg_groups_transposed_values <- emdia_neg_groups_transposed[-c(1:6)]
emdia_neg_groups_transposed_values <- as.data.frame(apply(emdia_neg_groups_transposed_values, 2, as.numeric))
emdia_neg_groupinfo <- as.character(as.matrix(emdia_neg_groups_transposed[c(5)]))


#PLSDA with mixomics
plsda_result_neg <- plsda(emdia_neg_groups_transposed_values, emdia_neg_groupinfo, ncomp = 2)
plotIndiv(plsda_result_neg, ind.names = FALSE, legend=TRUE,
          comp=c(1,2), ellipse = TRUE, 
          title = 'PLS-DA negative EmDia V1-3',
          X.label = 'PLS-DA comp 1', Y.label = 'PLS-DA comp 2',
          size.title = 20, size.xlabel = 20, size.ylabel = 20, size.axis = 20, size.legend = 20,
          size.legend.title = 20)
plotVar(plsda_result_neg)
plotLoadings(plsda_result_neg, comp = 1, method = 'mean', contrib = 'max')
plotLoadings(plsda_result_neg, comp = 2, method = 'mean', contrib = 'max')

#PLSDA for sex dependency in V1-3
emdia_neg_sexinfo <- as.character(as.matrix(emdia_neg_groups_transposed[c(4)]))

plsda_neg_sex_result <- plsda(emdia_neg_groups_transposed_values, emdia_neg_sexinfo, ncomp = 2)
plotIndiv(plsda_neg_sex_result, ind.names = FALSE, legend=TRUE,
          comp=c(1,2), ellipse = TRUE, 
          title = 'PLS-DA EmDia neg Influnce of Sex V1-3',
          X.label = 'PLS-DA comp 1', Y.label = 'PLS-DA comp 2',
          size.title = 20, size.xlabel = 20, size.ylabel = 20, size.axis = 20, size.legend = 20,
          size.legend.title = 20)
plotVar(plsda_neg_sex_result)
plotLoadings(plsda_neg_sex_result, comp = 1, method = 'mean', contrib = 'max')
plotLoadings(plsda_neg_sex_result, comp = 2, method = 'mean', contrib = 'max')

#
#Timepoint V1
#

#transpose emdia_data
emdia_data_neg_transposed <- as.data.frame(t(emdia_data_neg))
#set name of column
emdia_data_neg_transposed <- setNames(data.frame(t(emdia_data_neg[,-1])), emdia_data_neg[,1])
emdia_neg_groups_transposed_V1 <- filter(emdia_data_neg_transposed, (Visit == "V1"))
emdia_neg_groups_transposed_V1_values <- emdia_neg_groups_transposed_V1[-c(1:6)]
emdia_neg_groups_transposed_V1_values <- as.data.frame(apply(emdia_neg_groups_transposed_V1_values, 2, as.numeric))
emdia_neg_groupinfo_V1 <- as.character(as.matrix(emdia_neg_groups_transposed_V1[c(5)]))


#PLSDA with mixomics
plsda_result_neg_V1 <- plsda(emdia_neg_groups_transposed_V1_values, emdia_neg_groupinfo_V1, ncomp = 2)
plotIndiv(plsda_result_neg_V1, col = custom_colors, ind.names = FALSE, legend=TRUE,
          comp=c(1,2), ellipse = TRUE, 
          title = 'PLS-DA EmDia V1 negative',
          X.label = 'PLS-DA comp 1', Y.label = 'PLS-DA comp 2',
          size.title = 30, size.xlabel = 30, size.ylabel = 30, size.axis = 30, size.legend = 30,
          size.legend.title = 30,
          pch = 16)
plotVar(plsda_result_neg_V1)
plotLoadings(plsda_result_neg_V1, comp = 1, method = 'mean', contrib = 'max')
plotLoadings(plsda_result_neg_V1, comp = 2, method = 'mean', contrib = 'max')

#PLSDA for sex dependency in V1
emdia_neg_sexinfo_V1 <- as.character(as.matrix(emdia_neg_groups_transposed_V1[c(4)]))

plsda_neg_sex_result_V1 <- plsda(emdia_neg_groups_transposed_V1_values, emdia_neg_sexinfo_V1, ncomp = 2)
plotIndiv(plsda_neg_sex_result_V1, ind.names = FALSE, legend=TRUE,
          comp=c(1,2), ellipse = TRUE, 
          title = 'PLS-DA EmDia neg Influnce of Sex V1',
          X.label = 'PLS-DA comp 1', Y.label = 'PLS-DA comp 2',
          size.title = 20, size.xlabel = 20, size.ylabel = 20, size.axis = 20, size.legend = 20,
          size.legend.title = 20)
plotVar(plsda_neg_sex_result_V1)
plotLoadings(plsda_neg_sex_result_V1, comp = 1, method = 'mean', contrib = 'max')
plotLoadings(plsda_neg_sex_result_V1, comp = 2, method = 'mean', contrib = 'max')

#
#Timepoint V2
#

#transpose emdia_data
emdia_data_neg_transposed <- as.data.frame(t(emdia_data_neg))
#set name of column
emdia_data_neg_transposed <- setNames(data.frame(t(emdia_data_neg[,-1])), emdia_data_neg[,1])
emdia_neg_groups_transposed_V2 <- filter(emdia_data_neg_transposed, (Visit == "V2"))
emdia_neg_groups_transposed_V2_values <- emdia_neg_groups_transposed_V2[-c(1:6)]
emdia_neg_groups_transposed_V2_values <- as.data.frame(apply(emdia_neg_groups_transposed_V2_values, 2, as.numeric))
emdia_neg_groupinfo_V2 <- as.character(as.matrix(emdia_neg_groups_transposed_V2[c(5)]))


#PLSDA with mixomics
plsda_result_neg_V2 <- plsda(emdia_neg_groups_transposed_V2_values, emdia_neg_groupinfo_V2, ncomp = 2)
plotIndiv(plsda_result_neg_V2, col = custom_colors, ind.names = FALSE, legend=TRUE,
          comp=c(1,2), ellipse = TRUE, 
          title = 'PLS-DA EmDia V2 negative',
          X.label = 'PLS-DA comp 1', Y.label = 'PLS-DA comp 2',
          size.title = 30, size.xlabel = 30, size.ylabel = 30, size.axis = 30, size.legend = 30,
          size.legend.title = 30,
          pch = 16)
plotVar(plsda_result_neg_V2)
plotLoadings(plsda_result_neg_V2, comp = 1, method = 'mean', contrib = 'max')
plotLoadings(plsda_result_neg_V2, comp = 2, method = 'mean', contrib = 'max')

#PLSDA for sex dependency in V2
emdia_neg_sexinfo_V2 <- as.character(as.matrix(emdia_neg_groups_transposed_V2[c(4)]))

plsda_neg_sex_result_V2 <- plsda(emdia_neg_groups_transposed_V2_values, emdia_neg_sexinfo_V2, ncomp = 2)
plotIndiv(plsda_neg_sex_result_V2, ind.names = FALSE, legend=TRUE,
          comp=c(1,2), ellipse = TRUE, 
          title = 'PLS-DA EmDia neg Influnce of Sex V2',
          X.label = 'PLS-DA comp 1', Y.label = 'PLS-DA comp 2',
          size.title = 20, size.xlabel = 20, size.ylabel = 20, size.axis = 20, size.legend = 20,
          size.legend.title = 20)
plotVar(plsda_neg_sex_result_V2)
plotLoadings(plsda_neg_sex_result_V2, comp = 1, method = 'mean', contrib = 'max')
plotLoadings(plsda_neg_sex_result_V2, comp = 2, method = 'mean', contrib = 'max')

#
#Timepoint V3
#

#transpose emdia_data
emdia_data_neg_transposed <- as.data.frame(t(emdia_data_neg))
#set name of column
emdia_data_neg_transposed <- setNames(data.frame(t(emdia_data_neg[,-1])), emdia_data_neg[,1])
emdia_neg_groups_transposed_V3 <- filter(emdia_data_neg_transposed, (Visit == "V3"))
emdia_neg_groups_transposed_V3_values <- emdia_neg_groups_transposed_V3[-c(1:6)]
emdia_neg_groups_transposed_V3_values <- as.data.frame(apply(emdia_neg_groups_transposed_V3_values, 2, as.numeric))
emdia_neg_groupinfo_V3 <- as.character(as.matrix(emdia_neg_groups_transposed_V3[c(5)]))


#PLSDA with mixomics
plsda_result_neg_V3 <- plsda(emdia_neg_groups_transposed_V3_values, emdia_neg_groupinfo_V3, ncomp = 2)
plotIndiv(plsda_result_neg_V3, col = custom_colors, ind.names = FALSE, legend=TRUE,
          comp=c(1,2), ellipse = TRUE, 
          title = 'PLS-DA EmDia V3 negative',
          X.label = 'PLS-DA comp 1', Y.label = 'PLS-DA comp 2',
          size.title = 30, size.xlabel = 30, size.ylabel = 30, size.axis = 30, size.legend = 30,
          size.legend.title = 30,
          pch = 16)
plotVar(plsda_result_neg_V3)
plotLoadings(plsda_result_neg_V3, comp = 1, method = 'mean', contrib = 'max')
plotLoadings(plsda_result_neg_V3, comp = 2, method = 'mean', contrib = 'max')

#PLSDA for sex dependency in V3
emdia_neg_sexinfo_V3 <- as.character(as.matrix(emdia_neg_groups_transposed_V3[c(4)]))

plsda_neg_sex_result_V3 <- plsda(emdia_neg_groups_transposed_V3_values, emdia_neg_sexinfo_V3, ncomp = 2)
plotIndiv(plsda_neg_sex_result_V3, ind.names = FALSE, legend=TRUE,
          comp=c(1,2), ellipse = TRUE, 
          title = 'PLS-DA EmDia neg Influnce of Sex V3',
          X.label = 'PLS-DA comp 1', Y.label = 'PLS-DA comp 2',
          size.title = 20, size.xlabel = 20, size.ylabel = 20, size.axis = 20, size.legend = 20,
          size.legend.title = 20)
plotVar(plsda_neg_sex_result_V3)
plotLoadings(plsda_neg_sex_result_V3, comp = 1, method = 'mean', contrib = 'max')
plotLoadings(plsda_neg_sex_result_V3, comp = 2, method = 'mean', contrib = 'max')




###
###
###t-test and boxplots
###
###


###
###V1-3 positive ion mode data
###

#transpose emdia_data
emdia_data_transposed <- as.data.frame(t(emdia_data))
#set name of column
emdia_data_transposed <- setNames(data.frame(t(emdia_data[,-1])), emdia_data[,1])
#filter transposed dataframe for "Emdia" and "Placebo" group
emdia_groups_transposed <- filter(emdia_data_transposed, xor(Group == "Empa", Group == "Placebo"))
#remove column containing sample information for calculations except group
emdia_groups_transposed_t_test <- emdia_groups_transposed[-c(1,2,3,4,6)]
emdia_groups_transposed_t_test <- as.data.frame(emdia_groups_transposed_t_test)
emdia_groups_transposed_t_test[, c(2:132)] <- sapply(emdia_groups_transposed_t_test[, c(2:132)], as.numeric)
#Define Group as factor
emdia_groups_transposed_t_test$Group <- as.factor(emdia_groups_transposed_t_test$Group)


###Calculation with mean difference

#Create new data.frame for metabolite+p-value+mean_diff
results <- data.frame(Metabolite = character(), p.value = numeric(), mean_diff = numeric(), stringsAsFactors = FALSE)

#Goes through data (Column 2 to n)
for (metabolite in colnames(emdia_groups_transposed_t_test)[-1]) {
  #perform t-test
  t_test <- t.test(emdia_groups_transposed_t_test[[metabolite]] ~ emdia_groups_transposed_t_test$Group)
  
  #Calculate mean difference
  mean_diff <- mean(emdia_groups_transposed_t_test[[metabolite]][emdia_groups_transposed_t_test$Group == levels(emdia_groups_transposed_t_test$Group)[1]]) - 
    mean(emdia_groups_transposed_t_test[[metabolite]][emdia_groups_transposed_t_test$Group == levels(emdia_groups_transposed_t_test$Group)[2]])
  
  #save results in newly created data.frame
  results <- rbind(results, data.frame(Metabolite = metabolite, p.value = t_test$p.value, mean_diff = mean_diff))
}


#Save results in .csv file
write.csv(results, "t_test_results.csv", row.names = FALSE)

#
#Visulization of p-values mean difference
#

#Histogram of p-values
ggplot(results, aes(x = p.value)) +
  geom_histogram(binwidth = 0.05, fill = "blue", color = "black") +
  theme_minimal() +
  labs(title = "Histogram of p-values", x = "p-value", y = "frequency")

# Boxplot of metabolites per group for significant metabolites (p < 0.05)
significant_metabolites <- results %>% filter(p.value < 0.05) %>% pull(Metabolite)

if(length(significant_metabolites) > 0) {
  data_long <- emdia_groups_transposed_t_test %>%
    pivot_longer(cols = -Group, names_to = "Metabolite", values_to = "Intensity") %>%
    filter(Metabolite %in% significant_metabolites)
  
  ggplot(data_long, aes(x = Group, y = Intensity, fill = Group)) +
    geom_boxplot() +
    facet_wrap(~ Metabolite, scales = "free") +
    theme_minimal() +
    labs(title = "Boxplot of significant metabolites V1-3", x = "Group", y = "Intensity")+
  theme(
    axis.text.x = element_text(size = 10),
    axis.text.y = element_text(size = 10),
    axis.title.x = element_text(size = 10),
    axis.title.y = element_text(size = 10),
    legend.title = element_text(size = 20),
    legend.text = element_text(size = 20),
    legend.position = "right",
    plot.title = element_text(size = 20, hjust = 0.5))
} else {
  print("No significant metaboolites found (p < 0.05).")
}

#
# Volcanoplot of (log-transormed p-values vs. mean difference)
#
results <- results %>%
  mutate(neg_log10_p = -log10(p.value))

ggplot(results, aes(x = mean_diff, y = neg_log10_p)) +
  geom_point(alpha = 0.7) +
  theme_minimal() +
  labs(title = "Volcanoplot of t-test results EmDia V1-3", x = "mean difference", y = "-log10(p-Wert)") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red") +
  geom_text_repel(aes(label = ifelse(neg_log10_p > -log10(0.05), Metabolite, "")), 
                  max.overlaps = 10) +
  theme(text = element_text(size = 12))


#
# Volcanoplot of (log-transormed p-values vs. log2foldchange)
#

#Create new data.frame for metabolite+p-value+mean_diff
results <- data.frame(Metabolite = character(), p.value = numeric(), log2_fold_change = numeric(), stringsAsFactors = FALSE)

#Goes through data (Column 2 to n)
for (metabolite in colnames(emdia_groups_transposed_t_test)[-1]) {
  #perform t-test
  t_test <- t.test(emdia_groups_transposed_t_test[[metabolite]] ~ emdia_groups_transposed_t_test$Group)
  
  # calculate log2(Fold Change)
  mean_group1 <- mean(emdia_groups_transposed_t_test[[metabolite]][emdia_groups_transposed_t_test$Group == levels(emdia_groups_transposed_t_test$Group)[1]])
  mean_group2 <- mean(emdia_groups_transposed_t_test[[metabolite]][emdia_groups_transposed_t_test$Group == levels(emdia_groups_transposed_t_test$Group)[2]])
  log2_fold_change <- log2(mean_group2 / mean_group1)
  
  #Save results in newly created data.frame
  results <- rbind(results, data.frame(Metabolite = metabolite, p.value = t_test$p.value, log2_fold_change = log2_fold_change))
}

#Calculate neg log10 of p-values
results <- results %>%
  mutate(neg_log10_p = -log10(p.value))

#Volcanoplot of log2foldchange
ggplot(results, aes(x = log2_fold_change, y = neg_log10_p)) +
  geom_point(alpha = 1) +
  theme_classic() +
  labs(title = "Volcanoplot of EmDia V1-3", x = "log2(Fold Change)", y = "-log10(p-Wert)") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red") +
  geom_vline(xintercept=c(-0.6, 0.6), col="red") +
  geom_text_repel(aes(label = ifelse(neg_log10_p > -log10(0.05), Metabolite, "")), 
                  max.overlaps = 20) +
  theme(text = element_text(size = 20))

#
#Visulization of p-values V1-3
#

#Histogram of p-values
ggplot(results, aes(x = p.value)) +
  geom_histogram(binwidth = 0.05, fill = "blue", color = "black") +
  theme_minimal() +
  labs(title = "Histogram of p-values V1", x = "p-value", y = "frequency")

# Boxplot of metabolites per group for significant metabolites (p < 0.05)
significant_metabolites <- results %>% filter(p.value < 0.05) %>% pull(Metabolite)

if(length(significant_metabolites) > 0) {
  data_long <- emdia_groups_transposed_t_test %>%
    pivot_longer(cols = -Group, names_to = "Metabolite", values_to = "Intensity") %>%
    filter(Metabolite %in% significant_metabolites)
  
  ggplot(data_long, aes(x = Group, y = Intensity, fill = Group)) +
    geom_boxplot() +
    facet_wrap(~ Metabolite, scales = "free") +
    theme_minimal() +
    labs(title = "Boxplot of significant metabolites V1-3", x = "Group", y = "Intensity")+
    theme(
      axis.text.x = element_text(size = 10),
      axis.text.y = element_text(size = 10),
      axis.title.x = element_text(size = 10),
      axis.title.y = element_text(size = 10),
      legend.title = element_text(size = 20),
      legend.text = element_text(size = 20),
      legend.position = "right",
      plot.title = element_text(size = 20, hjust = 0.5))
} else {
  print("No significant metaboolites found (p < 0.05).")
}
  
###
###V1
###

#transpose emdia_data
emdia_data_transposed <- as.data.frame(t(emdia_data))
#set name of column
emdia_data_transposed <- setNames(data.frame(t(emdia_data[,-1])), emdia_data[,1])
#filter transposed dataframe for "Emdia" and "Placebo" group
emdia_groups_transposed <- filter(emdia_data_transposed, xor(Group == "Empa", Group == "Placebo"))
#filter transposed dataframe for Visit 1
emdia_groups_transposed_v1 <- filter(emdia_data_transposed, Visit == "V1")
#remove column containing sample information for calculations except group
emdia_groups_transposed_v1_t_test <- emdia_groups_transposed_v1[-c(1,2,3,4,6)]
emdia_groups_transposed_v1_t_test <- as.data.frame(emdia_groups_transposed_v1_t_test)
emdia_groups_transposed_v1_t_test[, c(2:132)] <- sapply(emdia_groups_transposed_v1_t_test[, c(2:132)], as.numeric)
#Define Group as factor
emdia_groups_transposed_v1_t_test$Group <- as.factor(emdia_groups_transposed_v1_t_test$Group)
  
#
# Volcanoplot of (log-transormed p-values vs. log2foldchange)
#

#Create new data.frame for metabolite+p-value+mean_diff
results_V1 <- data.frame(Metabolite = character(), p.value = numeric(), log2_fold_change = numeric(), stringsAsFactors = FALSE)

#Goes through data (Column 2 to n)
for (metabolite in colnames(emdia_groups_transposed_v1_t_test)[-1]) {
  #perform t-test
  t_test <- t.test(emdia_groups_transposed_v1_t_test[[metabolite]] ~ emdia_groups_transposed_v1_t_test$Group)
  
  # calculate log2(Fold Change)
  mean_group1 <- mean(emdia_groups_transposed_v1_t_test[[metabolite]][emdia_groups_transposed_v1_t_test$Group == levels(emdia_groups_transposed_v1_t_test$Group)[1]])
  mean_group2 <- mean(emdia_groups_transposed_v1_t_test[[metabolite]][emdia_groups_transposed_v1_t_test$Group == levels(emdia_groups_transposed_v1_t_test$Group)[2]])
  log2_fold_change <- log2(mean_group2 / mean_group1)
  
  #Save results in newly created data.frame
  results_V1 <- rbind(results_V1, data.frame(Metabolite = metabolite, p.value = t_test$p.value, log2_fold_change = log2_fold_change))
}

#Calculate neg log10 of p-values
results_V1 <- results_V1 %>%
  mutate(neg_log10_p = -log10(p.value))

#Volcanoplot of log2foldchange
ggplot(results_V1, aes(x = log2_fold_change, y = neg_log10_p)) +
  geom_point(alpha = 1) +
  theme_classic() +
  labs(title = "Volcanoplot of EmDia V1", x = "log2(Fold Change)", y = "-log10(p-Wert)") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red") +
  geom_vline(xintercept=c(-0.6, 0.6), col="red") +
  geom_text_repel(aes(label = ifelse(neg_log10_p > -log10(0.05), Metabolite, "")), 
                  max.overlaps = 20) +
  theme(text = element_text(size = 20))



###Visualization of average peak intensity of a specific metabolite.


# Define the specific metabolite you're interested in
selected_metabolite <- "EMPAGLIFLOZIN"  # Replace with the actual metabolite name

# Step 1: Reshape the data to long format
data_long_V1 <- emdia_groups_transposed_v1_t_test %>%
  pivot_longer(cols = -Group, names_to = "Metabolite", values_to = "Intensity")

# Step 2: Filter the data for the specific metabolite
filtered_data_V1 <- data_long_V1 %>%
  filter(Metabolite == selected_metabolite) %>%
  group_by(Group) %>%
  summarize(
    mean_intensity = mean(Intensity, na.rm = TRUE),
    sd_intensity = sd(Intensity, na.rm = TRUE),  # Standard deviation
    n = n(),  # Count of samples in each group
    se_intensity = sd_intensity / sqrt(n)  # Standard error
  )

# Step 3: Create a bar chart with error bars for the selected metabolite
ggplot(filtered_data_V1, aes(x = Group, y = mean_intensity, fill = Group)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.7) +  # Create the bars
  geom_errorbar(aes(ymin = mean_intensity - se_intensity, ymax = mean_intensity + se_intensity), 
                position = position_dodge(width = 0.8), width = 0.25) +  # Error bars
  theme_minimal() +
  labs(title = paste("Average Intensity of", selected_metabolite, "per Group (V1)"), 
       x = "Group", y = "Average Intensity") +
  theme(
    axis.text.x = element_text(size = 10),
    axis.text.y = element_text(size = 10),
    axis.title.x = element_text(size = 15),
    axis.title.y = element_text(size = 15),
    legend.title = element_text(size = 15),
    legend.text = element_text(size = 12),
    legend.position = "right",
    plot.title = element_text(size = 20, hjust = 0.5)
  )


#
#Visulization of p-values V1
#

#Histogram of p-values
ggplot(results_V1, aes(x = p.value)) +
  geom_histogram(binwidth = 0.05, fill = "blue", color = "black") +
  theme_minimal() +
  labs(title = "Histogram of p-values V1", x = "p-value", y = "frequency")

# Boxplot of metabolites per group for significant metabolites (p < 0.05)
significant_metabolites <- results_V1 %>% filter(p.value < 0.05) %>% pull(Metabolite)

if(length(significant_metabolites) > 0) {
  data_long_V1 <- emdia_groups_transposed_v1_t_test %>%
    pivot_longer(cols = -Group, names_to = "Metabolite", values_to = "Intensity") %>%
    filter(Metabolite %in% significant_metabolites)
  
  ggplot(data_long_V1, aes(x = Group, y = Intensity, fill = Group)) +
    geom_boxplot() +
    facet_wrap(~ Metabolite, scales = "free") +
    theme_minimal() +
    labs(title = "Boxplot of significant metabolites V1", x = "Group", y = "Intensity")+
    theme(
      axis.text.x = element_text(size = 10),
      axis.text.y = element_text(size = 10),
      axis.title.x = element_text(size = 10),
      axis.title.y = element_text(size = 10),
      legend.title = element_text(size = 20),
      legend.text = element_text(size = 20),
      legend.position = "right",
      plot.title = element_text(size = 20, hjust = 0.5))
} else {
  print("No significant metaboolites found (p < 0.05).")
}





###
###V2
###

#transpose emdia_data
emdia_data_transposed <- as.data.frame(t(emdia_data))
#set name of column
emdia_data_transposed <- setNames(data.frame(t(emdia_data[,-1])), emdia_data[,1])
#filter transposed dataframe for "Emdia" and "Placebo" group
emdia_groups_transposed <- filter(emdia_data_transposed, xor(Group == "Empa", Group == "Placebo"))
#filter transposed dataframe for Visit 1
emdia_groups_transposed_v2 <- filter(emdia_data_transposed, Visit == "V2")
#remove column containing sample information for calculations except group
emdia_groups_transposed_v2_t_test <- emdia_groups_transposed_v2[-c(1,2,3,4,6)]
emdia_groups_transposed_v2_t_test <- as.data.frame(emdia_groups_transposed_v2_t_test)
emdia_groups_transposed_v2_t_test[, c(2:132)] <- sapply(emdia_groups_transposed_v2_t_test[, c(2:132)], as.numeric)
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

#Volcanoplot of log2foldchange
ggplot(results_V2, aes(x = log2_fold_change, y = neg_log10_p)) +
  geom_point(alpha = 1) +
  theme_classic() +
  labs(title = "Volcanoplot of EmDia V2", x = "log2(Fold Change)", y = "-log10(p-Wert)") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red") +
  geom_vline(xintercept=c(-0.6, 0.6), col="red") +
  geom_text_repel(aes(label = ifelse(neg_log10_p > -log10(0.05), Metabolite, "")), 
                  max.overlaps = 20) +
  theme(text = element_text(size = 20))



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

# Step 3: Create a bar chart with error bars for the selected metabolite
ggplot(filtered_data_V2, aes(x = Group, y = mean_intensity, fill = Group)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.7, color = "black", size = 2) +  # Thicker outline for bars
  geom_errorbar(aes(ymin = mean_intensity - se_intensity, ymax = mean_intensity + se_intensity), 
                position = position_dodge(width = 0.8), width = 0.25, size = 2) +  # Thicker error bars
  theme_minimal() +
  scale_fill_manual(values = c("#1a80bb", "#ea801c")) +
  labs(title = paste("Average Intensity of", selected_metabolite, "(V2)"), 
       x = "Group", y = "Average Intensity") +
  theme(
    axis.text.x = element_text(size = 40),
    axis.text.y = element_text(size = 40),
    axis.title.x = element_text(size = 40),
    axis.title.y = element_text(size = 40),
    legend.title = element_text(size = 40),
    legend.text = element_text(size = 40),
    legend.position = "right",
    plot.title = element_text(size = 30, hjust = 0.3)
  )


#
#Visulization of p-values V2
#

#Histogram of p-values
ggplot(results_V2, aes(x = p.value)) +
  geom_histogram(binwidth = 0.05, fill = "blue", color = "black") +
  theme_minimal() +
  labs(title = "Histogram of p-values V2", x = "p-value", y = "frequency")

# Boxplot of metabolites per group for significant metabolites (p < 0.05)
significant_metabolites <- results_V2 %>% filter(p.value < 0.05) %>% pull(Metabolite)

if(length(significant_metabolites) > 0) {
  data_long_V2 <- emdia_groups_transposed_v2_t_test %>%
    pivot_longer(cols = -Group, names_to = "Metabolite", values_to = "Intensity") %>%
    filter(Metabolite %in% significant_metabolites)
  
  ggplot(data_long_V2, aes(x = Group, y = Intensity, fill = Group)) +
    geom_boxplot() +
    facet_wrap(~ Metabolite, scales = "free") +
    theme_minimal() +
    labs(title = "Boxplot of significant metabolites V2", x = "Group", y = "Intensity")+
    theme(
        axis.text.x = element_text(size = 40),
        axis.text.y = element_text(size = 40),
        axis.title.x = element_text(size = 40),
        axis.title.y = element_text(size = 40),
        legend.title = element_text(size = 40),
        legend.text = element_text(size = 40),
        legend.position = "right",
        plot.title = element_text(size = 20, hjust = 0.5))
} else {
  print("No significant metaboolites found (p < 0.05).")
}

###
# Single Boxplots of metabolites per group for significant metabolites (p < 0.05)
###

# Filter for significant metabolites
significant_metabolites <- results_V2 %>% filter(p.value < 0.05) %>% pull(Metabolite)

# Check if there are significant metabolites
if(length(significant_metabolites) > 0) {
  
  # Transform data to long format and filter for significant metabolites
  data_long_V2 <- emdia_groups_transposed_v2_t_test %>%
    pivot_longer(cols = -Group, names_to = "Metabolite", values_to = "Intensity") %>%
    filter(Metabolite %in% significant_metabolites)
  
  # Loop through each significant metabolite and create a boxplot
  for (metabolite in significant_metabolites) {
    
    # Filter data for the current metabolite
    metabolite_data <- data_long_V2 %>% filter(Metabolite == metabolite)
    
    # Create the boxplot with custom colors (purple and blue)
    p <- ggplot(metabolite_data, aes(x = Group, y = Intensity, fill = Group)) +
      geom_boxplot(lwd = 1.3) +
      scale_fill_manual(values = c("#1a80bb", "#ea801c")) +  # Custom color palette
      theme_minimal() +
      labs(title = paste("Boxplot of", metabolite),
           x = "Group", y = "Intensity") +
      theme(
        axis.text.x = element_text(size = 40),
        axis.text.y = element_text(size = 40),
        axis.title.x = element_text(size = 40),
        axis.title.y = element_text(size = 40),
        legend.title = element_text(size = 40),
        legend.text = element_text(size = 40),
        legend.position = "right",
        plot.title = element_text(size = 20, hjust = 0.5))
    
    # Print the plot
    print(p)
  }
  
} else {
  print("No significant metabolites found (p < 0.05).")
}



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
emdia_groups_transposed_v3_t_test[, c(2:132)] <- sapply(emdia_groups_transposed_v3_t_test[, c(2:132)], as.numeric)
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

#Volcanoplot of log2foldchange
ggplot(results_V3, aes(x = log2_fold_change, y = neg_log10_p)) +
  geom_point(alpha = 1) +
  theme_classic() +
  labs(title = "Volcanoplot of EmDia V3", x = "log2(Fold Change)", y = "-log10(p-Wert)") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red") +
  geom_vline(xintercept=c(-0.6, 0.6), col="red") +
  geom_text_repel(aes(label = ifelse(neg_log10_p > -log10(0.05), Metabolite, "")), 
                  max.overlaps = 20) +
  theme(text = element_text(size = 20))



#Visualization of the average peak area of a specific metabolite per group at V3.


# Define the specific metabolite you're interested in
selected_metabolite <- "EMPAGLIFLOZIN"  # Replace with the actual metabolite name

# Step 1: Reshape the data to long format
data_long_V3 <- emdia_groups_transposed_v3_t_test %>%
  pivot_longer(cols = -Group, names_to = "Metabolite", values_to = "Intensity")

# Step 2: Filter the data for the specific metabolite
filtered_data_V3 <- data_long_V2 %>%
  filter(Metabolite == selected_metabolite) %>%
  group_by(Group) %>%
  summarize(
    mean_intensity = mean(Intensity, na.rm = TRUE),
    sd_intensity = sd(Intensity, na.rm = TRUE),  # Standard deviation
    n = n(),  # Count of samples in each group
    se_intensity = sd_intensity / sqrt(n)  # Standard error
  )

# Step 3: Create a bar chart with error bars for the selected metabolite
ggplot(filtered_data_V3, aes(x = Group, y = mean_intensity, fill = Group)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.7, color = "black", size = 2) +  # Thicker outline for bars
  geom_errorbar(aes(ymin = mean_intensity - se_intensity, ymax = mean_intensity + se_intensity), 
                position = position_dodge(width = 0.8), width = 0.25, size = 2) +  # Thicker error bars
  theme_minimal() +
  scale_fill_manual(values = c("#1a80bb", "#ea801c")) +
  labs(title = paste("Average Intensity of", selected_metabolite, "(V3)"), 
       x = "Group", y = "Average Intensity") +
  theme(
    axis.text.x = element_text(size = 40),
    axis.text.y = element_text(size = 40),
    axis.title.x = element_text(size = 40),
    axis.title.y = element_text(size = 40),
    legend.title = element_text(size = 40),
    legend.text = element_text(size = 40),
    legend.position = "right",
    plot.title = element_text(size = 30, hjust = 0.3)
  )

#
#Visulization of p-values V3
#

#Histogram of p-values
ggplot(results_V3, aes(x = p.value)) +
  geom_histogram(binwidth = 0.05, fill = "blue", color = "black") +
  theme_minimal() +
  labs(title = "Histogram of p-values V3", x = "p-value", y = "frequency")

# Boxplot of metabolites per group for significant metabolites (p < 0.05)
significant_metabolites <- results_V3 %>% filter(p.value < 0.05) %>% pull(Metabolite)

if(length(significant_metabolites) > 0) {
  data_long_V3 <- emdia_groups_transposed_v3_t_test %>%
    pivot_longer(cols = -Group, names_to = "Metabolite", values_to = "Intensity") %>%
    filter(Metabolite %in% significant_metabolites)
  
  ggplot(data_long_V3, aes(x = Group, y = Intensity, fill = Group)) +
    geom_boxplot() +
    facet_wrap(~ Metabolite, scales = "free") +
    theme_minimal() +
    labs(title = "Boxplot of significant metabolites V3", x = "Group", y = "Intensity")+
    theme(
      axis.text.x = element_text(size = 10),
      axis.text.y = element_text(size = 10),
      axis.title.x = element_text(size = 10),
      axis.title.y = element_text(size = 10),
      legend.title = element_text(size = 20),
      legend.text = element_text(size = 20),
      legend.position = "right",
      plot.title = element_text(size = 20, hjust = 0.5))
} else {
  print("No significant metaboolites found (p < 0.05).")
}


###
###V1-3 negative
###

#transpose emdia_data
emdia_data_neg_transposed <- as.data.frame(t(emdia_data_neg))
#set name of column
emdia_data_neg_transposed <- setNames(data.frame(t(emdia_data_neg[,-1])), emdia_data_neg[,1])
#filter transposed dataframe for "Emdia" and "Placebo" group
emdia_data_neg_transposed <- filter(emdia_data_neg_transposed, xor(Group == "Empa", Group == "Placebo"))
#remove column containing sample information for calculations except group
emdia_data_neg_transposed_t_test <- emdia_data_neg_transposed[-c(1,2,3,4,6)]
emdia_data_neg_transposed_t_test <- as.data.frame(emdia_data_neg_transposed_t_test)
emdia_data_neg_transposed_t_test[, c(2:121)] <- sapply(emdia_data_neg_transposed_t_test[, c(2:121)], as.numeric)
#Define Group as factor
emdia_data_neg_transposed_t_test$Group <- as.factor(emdia_data_neg_transposed_t_test$Group)


#
# Volcanoplot of (log-transormed p-values vs. log2foldchange)
#

#Create new data.frame for metabolite+p-value+mean_diff
results <- data.frame(Metabolite = character(), p.value = numeric(), log2_fold_change = numeric(), stringsAsFactors = FALSE)

#Goes through data (Column 2 to n)
for (metabolite in colnames(emdia_data_neg_transposed_t_test)[-1]) {
  #perform t-test
  t_test <- t.test(emdia_data_neg_transposed_t_test[[metabolite]] ~ emdia_data_neg_transposed_t_test$Group)
  
  # calculate log2(Fold Change)
  mean_group1 <- mean(emdia_data_neg_transposed_t_test[[metabolite]][emdia_data_neg_transposed_t_test$Group == levels(emdia_data_neg_transposed_t_test$Group)[1]])
  mean_group2 <- mean(emdia_data_neg_transposed_t_test[[metabolite]][emdia_data_neg_transposed_t_test$Group == levels(emdia_data_neg_transposed_t_test$Group)[2]])
  log2_fold_change <- log2(mean_group2 / mean_group1)
  
  #Save results in newly created data.frame
  results <- rbind(results, data.frame(Metabolite = metabolite, p.value = t_test$p.value, log2_fold_change = log2_fold_change))
}

#Calculate neg log10 of p-values
results <- results %>%
  mutate(neg_log10_p = -log10(p.value))

#Volcanoplot of log2foldchange
ggplot(results, aes(x = log2_fold_change, y = neg_log10_p)) +
  geom_point(alpha = 1) +
  theme_classic() +
  labs(title = "Volcanoplot of EmDia V1-3 negative", x = "log2(Fold Change)", y = "-log10(p-Wert)") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red") +
  geom_vline(xintercept=c(-0.6, 0.6), col="red") +
  geom_text_repel(aes(label = ifelse(neg_log10_p > -log10(0.05), Metabolite, "")), 
                  max.overlaps = 20) +
  theme(text = element_text(size = 20))

#
#Visulization of p-values V1-3
#

#Histogram of p-values
ggplot(results, aes(x = p.value)) +
  geom_histogram(binwidth = 0.05, fill = "blue", color = "black") +
  theme_minimal() +
  labs(title = "Histogram of p-values V1-3 negative", x = "p-value", y = "frequency")

# Boxplot of metabolites per group for significant metabolites (p < 0.05)
significant_metabolites <- results %>% filter(p.value < 0.05) %>% pull(Metabolite)

if(length(significant_metabolites) > 0) {
  data_long <- emdia_data_neg_transposed_t_test %>%
    pivot_longer(cols = -Group, names_to = "Metabolite", values_to = "Intensity") %>%
    filter(Metabolite %in% significant_metabolites)
  
  ggplot(data_long, aes(x = Group, y = Intensity, fill = Group)) +
    geom_boxplot() +
    facet_wrap(~ Metabolite, scales = "free") +
    theme_minimal() +
    labs(title = "Boxplot of significant metabolites V1-3 negative", x = "Group", y = "Intensity")+
    theme(
      axis.text.x = element_text(size = 10),
      axis.text.y = element_text(size = 10),
      axis.title.x = element_text(size = 10),
      axis.title.y = element_text(size = 10),
      legend.title = element_text(size = 20),
      legend.text = element_text(size = 20),
      legend.position = "right",
      plot.title = element_text(size = 20, hjust = 0.5))
} else {
  print("No significant metaboolites found (p < 0.05).")
}


###
###V1
###

#transpose emdia_data
emdia_data_neg_transposed <- as.data.frame(t(emdia_data_neg))
#set name of column
emdia_data_neg_transposed <- setNames(data.frame(t(emdia_data_neg[,-1])), emdia_data_neg[,1])
#filter transposed dataframe for "Emdia" and "Placebo" group
emdia_data_neg_transposed <- filter(emdia_data_neg_transposed, xor(Group == "Empa", Group == "Placebo"))
#filter transposed dataframe for Visit 1
emdia_data_neg_transposed_v1 <- filter(emdia_data_neg_transposed, Visit == "V1")
#remove column containing sample information for calculations except group
emdia_data_neg_transposed_v1_t_test <- emdia_data_neg_transposed_v1[-c(1,2,3,4,6)]
emdia_data_neg_transposed_v1_t_test <- as.data.frame(emdia_data_neg_transposed_v1_t_test)
emdia_data_neg_transposed_v1_t_test[, c(2:121)] <- sapply(emdia_data_neg_transposed_v1_t_test[, c(2:121)], as.numeric)
#Define Group as factor
emdia_data_neg_transposed_v1_t_test$Group <- as.factor(emdia_data_neg_transposed_v1_t_test$Group)

#
# Volcanoplot of (log-transormed p-values vs. log2foldchange)
#

#Create new data.frame for metabolite+p-value+mean_diff
results_V1 <- data.frame(Metabolite = character(), p.value = numeric(), log2_fold_change = numeric(), stringsAsFactors = FALSE)

#Goes through data (Column 2 to n)
for (metabolite in colnames(emdia_data_neg_transposed_v1_t_test)[-1]) {
  #perform t-test
  t_test <- t.test(emdia_data_neg_transposed_v1_t_test[[metabolite]] ~ emdia_data_neg_transposed_v1_t_test$Group)
  
  # calculate log2(Fold Change)
  mean_group1 <- mean(emdia_data_neg_transposed_v1_t_test[[metabolite]][emdia_data_neg_transposed_v1_t_test$Group == levels(emdia_data_neg_transposed_v1_t_test$Group)[1]])
  mean_group2 <- mean(emdia_data_neg_transposed_v1_t_test[[metabolite]][emdia_data_neg_transposed_v1_t_test$Group == levels(emdia_data_neg_transposed_v1_t_test$Group)[2]])
  log2_fold_change <- log2(mean_group2 / mean_group1)
  
  #Save results in newly created data.frame
  results_V1 <- rbind(results_V1, data.frame(Metabolite = metabolite, p.value = t_test$p.value, log2_fold_change = log2_fold_change))
}

#Calculate neg log10 of p-values
results_V1 <- results_V1 %>%
  mutate(neg_log10_p = -log10(p.value))

#Volcanoplot of log2foldchange
ggplot(results_V1, aes(x = log2_fold_change, y = neg_log10_p)) +
  geom_point(alpha = 1) +
  theme_classic() +
  labs(title = "Volcanoplot of EmDia V1 negative", x = "log2(Fold Change)", y = "-log10(p-Wert)") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red") +
  geom_vline(xintercept=c(-0.6, 0.6), col="red") +
  geom_text_repel(aes(label = ifelse(neg_log10_p > -log10(0.05), Metabolite, "")), 
                  max.overlaps = 20) +
  theme(text = element_text(size = 20))


#
#Visulization of p-values V1
#

#Histogram of p-values
ggplot(results_V1, aes(x = p.value)) +
  geom_histogram(binwidth = 0.05, fill = "blue", color = "black") +
  theme_minimal() +
  labs(title = "Histogram of p-values V1", x = "p-value", y = "frequency")

# Boxplot of metabolites per group for significant metabolites (p < 0.05)
significant_metabolites <- results_V1 %>% filter(p.value < 0.05) %>% pull(Metabolite)

if(length(significant_metabolites) > 0) {
  data_long_V1 <- emdia_groups_transposed_v1_t_test %>%
    pivot_longer(cols = -Group, names_to = "Metabolite", values_to = "Intensity") %>%
    filter(Metabolite %in% significant_metabolites)
  
  ggplot(data_long_V1, aes(x = Group, y = Intensity, fill = Group)) +
    geom_boxplot() +
    facet_wrap(~ Metabolite, scales = "free") +
    theme_minimal() +
    labs(title = "Boxplot of significant metabolites V1 negative", x = "Group", y = "Intensity")+
    theme(
      axis.text.x = element_text(size = 10),
      axis.text.y = element_text(size = 10),
      axis.title.x = element_text(size = 10),
      axis.title.y = element_text(size = 10),
      legend.title = element_text(size = 20),
      legend.text = element_text(size = 20),
      legend.position = "right",
      plot.title = element_text(size = 20, hjust = 0.5))
} else {
  print("No significant metaboolites found (p < 0.05).")
}





###
###V2
###

#transpose emdia_data
emdia_data_neg_transposed <- as.data.frame(t(emdia_data_neg))
#set name of column
emdia_data_neg_transposed <- setNames(data.frame(t(emdia_data_neg[,-1])), emdia_data_neg[,1])
#filter transposed dataframe for "Emdia" and "Placebo" group
emdia_data_neg_transposed <- filter(emdia_data_neg_transposed, xor(Group == "Empa", Group == "Placebo"))
#filter transposed dataframe for Visit 1
emdia_data_neg_transposed_v2 <- filter(emdia_data_neg_transposed, Visit == "V2")
#remove column containing sample information for calculations except group
emdia_data_neg_transposed_v2_t_test <- emdia_data_neg_transposed_v2[-c(1,2,3,4,6)]
emdia_data_neg_transposed_v2_t_test <- as.data.frame(emdia_data_neg_transposed_v2_t_test)
emdia_data_neg_transposed_v2_t_test[, c(2:121)] <- sapply(emdia_data_neg_transposed_v2_t_test[, c(2:121)], as.numeric)
#Define Group as factor
emdia_data_neg_transposed_v2_t_test$Group <- as.factor(emdia_data_neg_transposed_v2_t_test$Group)

#
# Volcanoplot of (log-transormed p-values vs. log2foldchange)
#

#Create new data.frame for metabolite+p-value+mean_diff
results_V2 <- data.frame(Metabolite = character(), p.value = numeric(), log2_fold_change = numeric(), stringsAsFactors = FALSE)

#Goes through data (Column 2 to n)
for (metabolite in colnames(emdia_groups_transposed_v2_t_test)[-1]) {
  #perform t-test
  t_test <- t.test(emdia_groups_transposed_v2_t_test[[metabolite]] ~ emdia_data_neg_transposed_v2_t_test$Group)
  
  # calculate log2(Fold Change)
  mean_group1 <- mean(emdia_data_neg_transposed_v2_t_test[[metabolite]][emdia_data_neg_transposed_v2_t_test$Group == levels(emdia_data_neg_transposed_v2_t_test$Group)[1]])
  mean_group2 <- mean(emdia_data_neg_transposed_v2_t_test[[metabolite]][emdia_data_neg_transposed_v2_t_test$Group == levels(emdia_data_neg_transposed_v2_t_test$Group)[2]])
  log2_fold_change <- log2(mean_group2 / mean_group1)
  
  #Save results in newly created data.frame
  results_V2 <- rbind(results_V2, data.frame(Metabolite = metabolite, p.value = t_test$p.value, log2_fold_change = log2_fold_change))
}

#Calculate neg log10 of p-values
results_V2 <- results_V2 %>%
  mutate(neg_log10_p = -log10(p.value))

#Volcanoplot of log2foldchange
ggplot(results_V2, aes(x = log2_fold_change, y = neg_log10_p)) +
  geom_point(alpha = 1) +
  theme_classic() +
  labs(title = "Volcanoplot of EmDia V2 negative", x = "log2(Fold Change)", y = "-log10(p-Wert)") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red") +
  geom_vline(xintercept=c(-0.6, 0.6), col="red") +
  geom_text_repel(aes(label = ifelse(neg_log10_p > -log10(0.05), Metabolite, "")), 
                  max.overlaps = 20) +
  theme(text = element_text(size = 20))


#
#Visualization of p-values V2
#

#Histogram of p-values
ggplot(results_V2, aes(x = p.value)) +
  geom_histogram(binwidth = 0.05, fill = "blue", color = "black") +
  theme_minimal() +
  labs(title = "Histogram of p-values V1", x = "p-value", y = "frequency")

# Boxplot of metabolites per group for significant metabolites (p < 0.05)
significant_metabolites <- results_V2 %>% filter(p.value < 0.05) %>% pull(Metabolite)

if(length(significant_metabolites) > 0) {
  data_long_V2 <- emdia_data_neg_transposed_v2_t_test %>%
    pivot_longer(cols = -Group, names_to = "Metabolite", values_to = "Intensity") %>%
    filter(Metabolite %in% significant_metabolites)
  
  ggplot(data_long_V2, aes(x = Group, y = Intensity, fill = Group)) +
    geom_boxplot() +
    facet_wrap(~ Metabolite, scales = "free") +
    theme_minimal() +
    labs(title = "Boxplot of significant metabolites V2 negative", x = "Group", y = "Intensity")+
    theme(
      axis.text.x = element_text(size = 10),
      axis.text.y = element_text(size = 10),
      axis.title.x = element_text(size = 10),
      axis.title.y = element_text(size = 10),
      legend.title = element_text(size = 20),
      legend.text = element_text(size = 20),
      legend.position = "right",
      plot.title = element_text(size = 20, hjust = 0.5))
} else {
  print("No significant metaboolites found (p < 0.05).")
}


###
# Single Boxplots of metabolites per group for significant metabolites (p < 0.05) negative
###

# Filter for significant metabolites
significant_metabolites <- results_V2 %>% filter(p.value < 0.05) %>% pull(Metabolite)

# Check if there are significant metabolites
if(length(significant_metabolites) > 0) {
  
  # Transform data to long format and filter for significant metabolites
  data_long_V2 <- emdia_data_neg_transposed_v2_t_test %>%
    pivot_longer(cols = -Group, names_to = "Metabolite", values_to = "Intensity") %>%
    filter(Metabolite %in% significant_metabolites)
  
  # Loop through each significant metabolite and create a boxplot
  for (metabolite in significant_metabolites) {
    
    # Filter data for the current metabolite
    metabolite_data <- data_long_V2 %>% filter(Metabolite == metabolite)
    
    # Create the boxplot with custom colors (purple and blue)
    p <- ggplot(metabolite_data, aes(x = Group, y = Intensity, fill = Group)) +
      geom_boxplot(lwd = 1.3) +
      scale_fill_manual(values = c("#1a80bb", "#ea801c")) +  # Custom color palette
      theme_minimal() +
      labs(title = paste("Boxplot of", metabolite),
           x = "Group", y = "Intensity") +
      theme(
        axis.text.x = element_text(size = 40),
        axis.text.y = element_text(size = 40),
        axis.title.x = element_text(size = 40),
        axis.title.y = element_text(size = 40),
        legend.title = element_text(size = 40),
        legend.text = element_text(size = 40),
        legend.position = "right",
        plot.title = element_text(size = 20, hjust = 0.5))
    
    # Print the plot
    print(p)
  }
  
} else {
  print("No significant metabolites found (p < 0.05).")
}



###
###V3
###

#transpose emdia_data
emdia_data_neg_transposed <- as.data.frame(t(emdia_data_neg))
#set name of column
emdia_data_neg_transposed <- setNames(data.frame(t(emdia_data_neg[,-1])), emdia_data_neg[,1])
#filter transposed dataframe for "Emdia" and "Placebo" group
emdia_data_neg_transposed <- filter(emdia_data_neg_transposed, xor(Group == "Empa", Group == "Placebo"))
#filter transposed dataframe for Visit 1
emdia_data_neg_transposed_v3 <- filter(emdia_data_neg_transposed, Visit == "V2")
#remove column containing sample information for calculations except group
emdia_data_neg_transposed_v3_t_test <- emdia_data_neg_transposed_v3[-c(1,2,3,4,6)]
emdia_data_neg_transposed_v3_t_test <- as.data.frame(emdia_data_neg_transposed_v3_t_test)
emdia_data_neg_transposed_v3_t_test[, c(2:121)] <- sapply(emdia_data_neg_transposed_v3_t_test[, c(2:121)], as.numeric)
#Define Group as factor
emdia_data_neg_transposed_v3_t_test$Group <- as.factor(emdia_data_neg_transposed_v3_t_test$Group)

#
# Volcanoplot of (log-transormed p-values vs. log2foldchange)
#

#Create new data.frame for metabolite+p-value+mean_diff
results_V3 <- data.frame(Metabolite = character(), p.value = numeric(), log2_fold_change = numeric(), stringsAsFactors = FALSE)

#Goes through data (Column 2 to n)
for (metabolite in colnames(emdia_data_neg_transposed_v3_t_test)[-1]) {
  #perform t-test
  t_test <- t.test(emdia_data_neg_transposed_v3_t_test[[metabolite]] ~ emdia_data_neg_transposed_v3_t_test$Group)
  
  # calculate log2(Fold Change)
  mean_group1 <- mean(emdia_data_neg_transposed_v3_t_test[[metabolite]][emdia_data_neg_transposed_v3_t_test$Group == levels(emdia_data_neg_transposed_v3_t_test$Group)[1]])
  mean_group2 <- mean(emdia_data_neg_transposed_v3_t_test[[metabolite]][emdia_data_neg_transposed_v3_t_test$Group == levels(emdia_data_neg_transposed_v3_t_test$Group)[2]])
  log2_fold_change <- log2(mean_group2 / mean_group1)
  
  #Save results in newly created data.frame
  results_V3 <- rbind(results_V3, data.frame(Metabolite = metabolite, p.value = t_test$p.value, log2_fold_change = log2_fold_change))
}

#Calculate neg log10 of p-values
results_V3 <- results_V3 %>%
  mutate(neg_log10_p = -log10(p.value))

#Volcanoplot of log2foldchange
ggplot(results_V3, aes(x = log2_fold_change, y = neg_log10_p)) +
  geom_point(alpha = 1) +
  theme_classic() +
  labs(title = "Volcanoplot of EmDia V3 negative", x = "log2(Fold Change)", y = "-log10(p-Wert)") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red") +
  geom_vline(xintercept=c(-0.6, 0.6), col="red") +
  geom_text_repel(aes(label = ifelse(neg_log10_p > -log10(0.05), Metabolite, "")), 
                  max.overlaps = 20) +
  theme(text = element_text(size = 20))


#
#Visulization of p-values V3
#

#Histogram of p-values
ggplot(results_V3, aes(x = p.value)) +
  geom_histogram(binwidth = 0.05, fill = "blue", color = "black") +
  theme_minimal() +
  labs(title = "Histogram of p-values V3", x = "p-value", y = "frequency")

# Boxplot of metabolites per group for significant metabolites (p < 0.05)
significant_metabolites <- results_V3 %>% filter(p.value < 0.05) %>% pull(Metabolite)

if(length(significant_metabolites) > 0) {
  data_long_V3 <- emdia_data_neg_transposed_v3_t_test %>%
    pivot_longer(cols = -Group, names_to = "Metabolite", values_to = "Intensity") %>%
    filter(Metabolite %in% significant_metabolites)
  
  ggplot(data_long_V3, aes(x = Group, y = Intensity, fill = Group)) +
    geom_boxplot() +
    facet_wrap(~ Metabolite, scales = "free") +
    theme_minimal() +
    labs(title = "Boxplot of significant metabolites V3 negative", x = "Group", y = "Intensity")+
    theme(
      axis.text.x = element_text(size = 10),
      axis.text.y = element_text(size = 10),
      axis.title.x = element_text(size = 10),
      axis.title.y = element_text(size = 10),
      legend.title = element_text(size = 20),
      legend.text = element_text(size = 20),
      legend.position = "right",
      plot.title = element_text(size = 20, hjust = 0.5))
} else {
  print("No significant metaboolites found (p < 0.05).")
}



###
###
###Patternhunter - correlation of one specific metabolite against a set of metabolites
###
###



#transpose emdia_data
emdia_data_transposed <- as.data.frame(t(emdia_data))
#set name of column
emdia_data_transposed <- setNames(data.frame(t(emdia_data[,-1])), emdia_data[,1])
#filter transposed dataframe for "Emdia" and "Placebo" group
emdia_groups_transposed <- filter(emdia_data_transposed, xor(Group == "Empa", Group == "Placebo"))
emdia_groups_transposed <- filter(emdia_groups_transposed, Visit == "V2")


emdia_groups_transposed_m_correlation <- emdia_groups_transposed[-c(1,2,3,4,6)]


# Example data structure:
# Rows: Samples
# Columns: Metabolites (with an optional "Group" column indicating patient group)
# Remove the Group column if it exists
if ("Group" %in% colnames(emdia_groups_transposed_m_correlation)) {
  emdia_groups_transposed_m_correlation <- emdia_groups_transposed_m_correlation[, -which(names(emdia_groups_transposed_m_correlation) == "Group")]
  emdia_groups_transposed_m_correlation <- as.data.frame(apply(emdia_groups_transposed_m_correlation,2,as.numeric))
  }

# Transpose data so that metabolites are columns and samples are rows if needed
# data <- t(data)

# Normalize the data (optional, but can help with scaling issues)
#emdia_groups_transposed_m_correlation <- scale(emdia_groups_transposed_m_correlation)

# Select the specific metabolite of interest
metabolite_of_interest <- "DEOXYHEXOSE"  # Replace with the name of your specific metabolite

# Ensure the metabolite of interest is in the dataset
if (!(metabolite_of_interest %in% colnames(emdia_groups_transposed_m_correlation))) {
  stop("Metabolite of interest not found in the dataset")
}

# Calculate the correlation of the selected metabolite with all other metabolites
correlations <- sapply(emdia_groups_transposed_m_correlation, function(x) cor(emdia_groups_transposed_m_correlation[, metabolite_of_interest], x, use = "pairwise.complete.obs"))

# Remove the correlation of the metabolite with itself
correlations <- correlations[names(correlations) != metabolite_of_interest]

# Create a data frame for easy plotting
correlation_df <- data.frame(
  Metabolite = names(correlations),
  Correlation = correlations
)

# Select the top 20 metabolites by absolute correlation value
top_25 <- correlation_df %>%
  mutate(AbsCorrelation = abs(Correlation)) %>%
  top_n(10, AbsCorrelation) %>%
  arrange(desc(AbsCorrelation))

# Plot the top 20 correlations
# Plot the top 20 correlations with conditional coloring
ggplot(top_25, aes(x = reorder(Metabolite, Correlation), y = Correlation, fill = Correlation > 0)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = c("TRUE" = "#298c8c", "FALSE" = "#800074")) +
  coord_flip() +
  labs(
    title = paste("Top 10 Correlations of", metabolite_of_interest, "with other Metabolites"),
    x = "Metabolites",
    y = "Correlation"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(size = 20),
    axis.text.y = element_text(size = 20),
    axis.title.x = element_text(size = 20),
    axis.title.y = element_text(size = 20),
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 10),
    legend.position = "none",
    plot.title = element_text(size = 15, hjust = 0, vjust = 0))



# Print the top 20 correlation values
print(top_25)



###
###Pattern Hunter negative
###



#transpose emdia_data
emdia_data_neg_transposed <- as.data.frame(t(emdia_data_neg))
#set name of column
emdia_data_neg_transposed <- setNames(data.frame(t(emdia_data_neg[,-1])), emdia_data_neg[,1])
#filter transposed dataframe for "Emdia" and "Placebo" group
emdia_groups_neg_transposed <- filter(emdia_data_neg_transposed, xor(Group == "Empa", Group == "Placebo"))
emdia_groups_neg_transposed <- filter(emdia_groups_neg_transposed, Visit == "V2")


emdia_groups_neg_transposed_m_correlation <- emdia_groups_neg_transposed[-c(1,2,3,4,6)]


# Example data structure:
# Rows: Samples
# Columns: Metabolites (with an optional "Group" column indicating patient group)
# Remove the Group column if it exists
if ("Group" %in% colnames(emdia_groups_neg_transposed_m_correlation)) {
  emdia_groups_neg_transposed_m_correlation <- emdia_groups_neg_transposed_m_correlation[, -which(names(emdia_groups_neg_transposed_m_correlation) == "Group")]
  emdia_groups_neg_transposed_m_correlation <- as.data.frame(apply(emdia_groups_neg_transposed_m_correlation,2,as.numeric))
}

# Transpose data so that metabolites are columns and samples are rows if needed
# data <- t(data)

# Normalize the data (optional, but can help with scaling issues)
#emdia_groups_transposed_m_correlation <- scale(emdia_groups_transposed_m_correlation)

# Select the specific metabolite of interest
metabolite_of_interest <- "CITRATE"  # Replace with the name of your specific metabolite

# Ensure the metabolite of interest is in the dataset
if (!(metabolite_of_interest %in% colnames(emdia_groups_neg_transposed_m_correlation))) {
  stop("Metabolite of interest not found in the dataset")
}

# Calculate the correlation of the selected metabolite with all other metabolites
correlations <- sapply(emdia_groups_neg_transposed_m_correlation, function(x) cor(emdia_groups_neg_transposed_m_correlation[, metabolite_of_interest], x, use = "pairwise.complete.obs"))

# Remove the correlation of the metabolite with itself
correlations <- correlations[names(correlations) != metabolite_of_interest]

# Create a data frame for easy plotting
correlation_df <- data.frame(
  Metabolite = names(correlations),
  Correlation = correlations
)

# Select the top 20 metabolites by absolute correlation value
top_25 <- correlation_df %>%
  mutate(AbsCorrelation = abs(Correlation)) %>%
  top_n(10, AbsCorrelation) %>%
  arrange(desc(AbsCorrelation))

# Plot the top 20 correlations
# Plot the top 20 correlations with conditional coloring
ggplot(top_25, aes(x = reorder(Metabolite, Correlation), y = Correlation, fill = Correlation > 0)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = c("TRUE" = "lightgreen", "FALSE" = "lightblue")) +
  coord_flip() +
  labs(
    title = paste("Top 25 Correlations of", metabolite_of_interest, "with other Metabolites"),
    x = "Metabolites",
    y = "Correlation"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(size = 15),
    axis.text.y = element_text(size = 15),
    axis.title.x = element_text(size = 15),
    axis.title.y = element_text(size = 15),
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 10),
    legend.position = "none",
    plot.title = element_text(size = 15, hjust = 0, vjust = 0))

# Print the top 20 correlation values
print(top_25)


