#!/usr/bin/Rscript

# Purpose of script is to perform differential analysis on the genentic data in 
# joined_table. Using the program DeSeq2, deifferential expression analysis can
# be performed on a DeSeq object using a desired design. Log transformations
# can then be perfromed on these analysis outputs to be used in the generation
# of clustering plots such as PCA plots. 

# package info and download:
# https://www.bioconductor.org/packages/release/bioc/html/DESeq2.html

# load packages as needed
library("dplyr")
library("tidyr")
library("knitr")
library("ggplot2")
library("magrittr")
library("tibble")
library("DESeq2")
library("RColorBrewer")
library("ggfortify")
library("cluster")

# this package allows for the easy inclusion of literature citations in our Rmd
# more info here: https://github.com/crsh/citr
# and here:
# http://rmarkdown.rstudio.com/authoring_bibliographies_and_citations.html
library("citr")

# load data produced from analysis scripts using
# something like load("output/processed_data.Rdata")
load("output/final_compiled_counts/joined_count_data.RData")

# test that it loaded correctly before proceeding
stopifnot(exists("joined_table"))

# Round the read_abund_for_deseq column because it is required to make DESeq dataset
rounded <- joined_table %>%
  mutate_if(is.numeric, ~round(., 0))

# Make a counts table of read_abund_for_deseq from rounded dataset
# to read into DESeq dataset and use for PCA plot
cts <- rounded %>%
  mutate_if(is.numeric, ~round(., 1)) %>%
  select(sample, genename, read_abund_for_deseq) %>%
  arrange(sample) %>%
  spread(sample, read_abund_for_deseq) %>%
  column_to_rownames(var = "genename")

# Make a data column table of sample names and corresponding disease_state
# from rounded dataset to read into DESeq dataset
coldata <- rounded %>%
  select(sample, disease_state) %>%
  distinct(sample, .keep_all = TRUE) %>%
  arrange(sample)

# Use DESeq to read in the coldata and cts tables into a single DESeq dataset
# using the chosen design based on the disease_state to make graphs
dds <- DESeqDataSetFromMatrix(countData = cts,
                              colData = coldata,
                              design = ~ disease_state)

# Perform differential expression analysis
dds <- DESeq(dds)

# Generate results table which extracts the log2 fold changes, p values, and adjusted p values. 
res <- results(dds)

# Perform log transformation of DESeq dataset for clustering, heatmaps, etc. 
rld <- rlogTransformation(dds)

# Save the file so it can be read in by knit in the .Rmd
save(rld,
     file = "./output/deseq2_data/rld_deseq2_data.RData")
