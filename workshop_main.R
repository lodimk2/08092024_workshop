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
counts_file <-  # Function to load in the counts file 

# Inspect metadata and counts
# What's the best function to get a quick overview of the data?
##
##

# How can we look at the current disease states in the dataset, to make the resulting comparison?

unique(metadata$`disease state:ch1`)

# Now, lets view the metadata and see what conditions are present in the dataset 
# Remove what doesn't belong 

## 
# Remove the convalescent column from count_file
counts <- counts_file[, metadata$`disease state:ch1` != ""]
# Remove it from the metadata file 
metadata <- metadata %>% #Which dplyr function can go here? (!`disease state:ch1` == "")


# Check the dimensions of counts and metadata to make sure they overlap. What needs to overlap?
##
##


# Create a condition column in the metadata
## 

# Sanity check to make sure the conditions are correct, what function can we use?
## 

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
# Order using dplyr
## 
## 

# Filter the non significant results
## Which function can go here?



results_table_filtered <- results_table %>%
    filter(abs(log2FoldChange) > 2)
# View the comparison that DESEQ2 made 
resultsNames(dds)

# DESIGN QUESTION: How can we interpret negative log2fold change?


# Verify the current working directory. How can we do this in RStudio?
getwd()

# Set the working directory to 
setwd("/Users/lodimk2/Documents/0806_workshop/")
write.csv(results_table_filtered, "results.csv") 

