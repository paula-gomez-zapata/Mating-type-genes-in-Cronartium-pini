##############################################
# Plot proportion of heterozygous SNPs per sample
#############################################

# Clear global environment
rm(list = ls())

# Define a user-specific library directory
setwd("C:/Users/pago0001/Documents/Personal Documents/Paula's notes Meetings Project Progress Follow-up/Report 31/vcf_biallelic_filtered_removed_include_variants_only_forPCA")
getwd()

# Install and load necessary libraries
#install.packages(c("ggplot2", "withr", "labeling", "farver", "gridExtra", "dplyr", "adegenet"))
#if (!require("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install("VariantAnnotation")
library(vcfR)
library(withr)
library(ggplot2)
library(labeling)
library(farver)
library(gridExtra)
library(dplyr)
library(adegenet)
library(reshape2)
library(ape)
library(ggtree)
library(pinfsc50)
library(VariantAnnotation)

#Read the vcf file
vcf <- read.vcfR("vcf_filtered_no_constants_snps_biallelic_18samples_morethan80percent_without_missdata.vcf.gz")
vcf

#here we use the function vcfR2genlight() to convert our vcfR object to a genlight object. 
#This makes our VCF data available to the analyses in adegenet.
vcf_gl<- vcfR2genlight(vcf)
vcf_gl

#Genotype count per sample

#get sample names
sample_names <- indNames(vcf_gl)
sample_names

scaffolds <- unique(chromosome(vcf_gl))

#extract genotype matrix
genotype<-t(as.matrix(vcf_gl))
genotype
#genotype[1,] 

#count the number of occurrences for first column where the genotype is equal to 2, 1 or 0.
#0 means the sample is homozygous to the reference: 0/0
#1 means the sample is heterozygous, carrying one copy of each of the reference and alternate allele: 0/1
#2 means homozygous alternate: 1/1
length(which(genotype[,1] == 2))
length(which(genotype[,1] == 1))
length(which(genotype[,1] == 0))


# Initialize vectors to store counts for each value (0, 1, 2)
count_data <- data.frame(Sample = sample_names, Count_0 = 0, Count_1 = 0, Count_2 = 0)
#######################################
# Iterate over each column (sample)
for (i in 1:ncol(genotype)) {
  # Get counts for 0, 1, 2 in the current column
  counts <- table(genotype[, i])
  
  # Create a data frame with all possible genotype counts
  count_df <- data.frame(Genotype = 0:2, Count = 0)
  
  # Update counts for existing genotype values
  count_df$Count[as.character(count_df$Genotype) %in% names(counts)] <- counts
  
  # Assign counts to count_data
  count_data[i, 2:4] <- count_df$Count
}

# Reshape data for ggplot
count_data_melt <- melt(count_data, id.vars = "Sample")

count_data_melt$Sample <- as.factor(count_data_melt$Sample)

head(count_data_melt)

# Define a darker color palette
dark_palette <- c("#1f78b4", "#33a02c", "#e31a1c", "#ff7f00", "#6a3d9a")

# Define legend labels
legend_labels <- c("Homozygous to ref", "Heterozygous", "Homozygous alt")
#0 means the sample is homozygous to the reference: 0/0
#1 means the sample is heterozygous, carrying one copy of each of the reference and alternate allele: 0/1
#2 means homozygous alternate: 1/1

# Plot using ggplot with stacked bars and modified legend
a <- ggplot(count_data_melt, aes(x = Sample, y = value, fill = variable)) +
  geom_bar(stat = "identity") +  # Use geom_bar for stacked bars
  labs(x = "Sample", y = "SNPs Count", fill = "Genotype") +
  scale_fill_manual(values = dark_palette, labels = legend_labels) +  # Modify legend labels
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10)) +  # Rotating x-axis labels
  scale_y_continuous(labels = scales::comma)

###################################################################################################

# Add a new column representing the sum of count_1 and count_2
count_data$snps <- count_data$Count_1 + count_data$Count_2
#Add a new column representing the proportion of count_1 (heterozygous) among the total of SNPs
count_data$prop_het <- count_data$Count_1 / count_data$snps * 100

count_data$Form <- c("Autoecious", "Autoecious",	"Autoecious",	"Autoecious",	"Heteroecious",	"Heteroecious",	"Autoecious",	"Heteroecious",	"Heteroecious",	"Heteroecious",	"Heteroecious",	"Autoecious",	"Heteroecious",	"Heteroecious",	"Autoecious",	"Autoecious",	"Autoecious",	"Autoecious")

# Here, I will plot the proportion of heterozygous SNPs per sample
b <- ggplot(count_data, aes(x = Sample, y = prop_het, fill = Form)) +
  geom_bar(stat = "identity") +
  labs(x = "Sample", y = "Proportion of Heterozygous SNPs (%)") +
  scale_fill_manual(values = c("Heteroecious" = "purple", "Autoecious" = "#ff7f00")) +  # Specify colors for each category
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10))
