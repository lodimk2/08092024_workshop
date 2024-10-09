# Install necessary packages
# ggplot2 is very useful/customizable for visualization!
install.packages("ggplot2")
# BiocManager allows users to install Bioconductor packages, which is where several Bioinformatics specific packages are posted and maintained 
install.packages("BiocManager")

BiocManager::install(c("GEOquery", "DESeq2"))

# Load libraries
library(GEOquery)
library(ggplot2)
library(DESeq2)
library(dplyr)

# Get the GEO dataset
gse <- getGEO("GSE152418", GSEMatrix = TRUE)

# Extract expression matrix and metadata
# The raw counts are located in supplementary files; for simplicity, we can use processed count data
gse_data <- gse[[1]]
metadata <- pData(gse_data)

# Load in counts file 
counts_file <- read.csv("GSE152418_raw_counts_GRCh38.p13_NCBI.tsv", sep = "\t", row.names = 1)

# Inspect metadata and counts
head(metadata)
head(counts)

unique(metadata$`disease state:ch1`)

# Remove what doesn't belong 
# Create a condition column

# Remove the convalescent column from count_file
counts <- counts_file[, metadata$`disease state:ch1` != "Convalescent"]
# Remove it from the metadata file 
metadata <- metadata %>% filter(!`disease state:ch1` == "Convalescent")


# Check the dimensions of counts and metadata to make sure they overlap

dim(counts)
dim(metadata)


metadata$condition <- metadata$`disease state:ch1`
# Check the condition column to identify patients and controls
unique(metadata$condition)


# PERFORM DESeq2 ANALYSIS
# First, let's check the DESeq2 Vignette to understand how to use DESeq2

condition <- factor(metadata$condition)  # Adjust the column name for disease status
coldata <- data.frame(row.names = colnames(counts), condition)
coldata$condition <- relevel(coldata$condition, ref = "COVID-19")
dds <- DESeqDataSetFromMatrix(countData = counts, colData = coldata, design = ~ condition)
dds <- DESeq(dds)
res <- results(dds)

# Make the res table a dataframe
results_table <- as.data.frame(res)

# We need to first sort by p value and then by log2FoldChange
results_table <- results_table[order(results_table$pvalue), ]
results_table <- results_table[order(results_table$log2FoldChange), ]

# Filter the non significant results
results_table <- results_table %>% filter(padj < 0.05)
results_table_filtered <- results_table %>%
    filter(abs(log2FoldChange) > 2)
# View the comparison that DESEQ2 made 
resultsNames(dds)

# How can we interpret negative log2fold change?

# Let's write a function to save the DE file 

# Verify the current working directory 
getwd()

# Set the working directory to 
setwd("/Users/lodimk2/Documents/0806_workshop/")
write.csv(results_table_filtered, "results.csv") 
