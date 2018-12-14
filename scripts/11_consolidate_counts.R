#!/usr/bin/Rscript

# This purpose of this R script is to read in all of the transcript counts
# made by sailfish individually on each sample, aggregate those counts (which
# are at the transcript level) to the gene level, normalize them based on
# length and read depth, and then build a table with rows as genes and columns
# as samples, with each cell being the normalized count for that gene in that
# sample. Then, we want to join in the two other metadata files we have,
# from SRA and from the supplementatry table from the original manuscript that
# has information like patient age and smoking behavior. We have to do some
# cleaning along the way, and the gathered (melted) output file is rather big
# so we only include interesting columns, and write out the output in compressed
# csv format as well as binary RData format

# package info and download:
# http://bioconductor.org/packages/release/bioc/html/tximport.html

## the basic approach is based around the
## instructions here:
## http://bioconductor.org/packages/release/bioc/vignettes/tximport/inst/doc/tximport.html

# set to TRUE to get some extra output to help troubleshoot
debug <- FALSE

## install needed packages
## first thing get biocLite: source("https://bioconductor.org/biocLite.R")
## load the packages with: biocLite("tximport")
## load the packages with: biocLite("tximportData")

# Load packages
library("tximport")
library("readr")
library("tximportData")
library("tidyr")
library("dplyr")
library("tibble")
library("magrittr")

# Load the transcript to gene lookup table we constructed
# this needs to be tx ID, then gene ID
tx2gene <- read.table(paste0("./output/gene_name_lookup_files/",
                             "transcript_to_gene_lookup_2cols.txt"),
                      header = FALSE,
                      col.names =
                          c("TX-NAME", "GENE-NAME"))

# read in the SRA sample metadata table, downloaded from the SRA Run Selector
samples <- read.table("./data/metadata/SraRunTable.txt",
                      header = TRUE,
                      sep = "\t")

# construct vector of paths to each of the sailfish quant files to read in
# and add names so that the sample names get added to the constructed table
files <- file.path(".", "output", "sailfish_quants", samples$Run, "quant.sf")
names(files) <- samples$Run

# error out if the files don't exist
stopifnot(all(file.exists(files)))

# load in and summarize all of the sailfish files, join into a single
# table summarized by gene (instead of transcript) and scale counts
txi <- tximport(files,
                type = "sailfish",
                tx2gene = tx2gene,
                countsFromAbundance = "lengthScaledTPM")

# print out a peak at this dataset to make sure it's ok
if (debug) {
  txi$counts[1:4, 1:4]
  str(txi$counts)
}

# build the final melted/gathered long format table that includes all the
# counts for all the samples, by gene, joined to the two metadata tables
# we read in previously, filter out any where we don't have good metadata on
# the patients, since for this analysis we want to look at those variables
# and then finally select out the columns of interest, discarding the rest
# which are largely redundant and make the file size unacceptably large
final_table <- txi$counts %>%
    as.data.frame() %>%
    rownames_to_column() %>%
    rename(GeneName = rowname) %>%
    gather(key = "Sample",
           value = "counts_lengthScaledTPM",
           -GeneName) %>%
    left_join(samples, by = c("Sample" = "Run")) %>%
    mutate_if(is.factor, as.character) %>%
    select(GeneName,
           Sample,
           counts_lengthScaledTPM,
           disease_state,
           source_name,
           Organism,
           tissue)

names(final_table) %<>% tolower


# create raw abundance table (not normalized) as this is what DESEQ2 needs
txi_for_deseq <- tximport(files,
                          type = "sailfish",
                          tx2gene = tx2gene)

# turn it into a data frame
deseq_table <- txi_for_deseq$counts %>%
    as.data.frame() %>%
    rownames_to_column() %>%
    rename(GeneName = rowname) %>%
    gather(key = "Sample",
           value = "read_abund_for_deseq",
           -GeneName)

# join to earlier table based on GeneName and Sample
joined_table <- final_table %>%
    left_join(deseq_table, by = c("genename" = "GeneName", "sample" = "Sample"))

# write out the final table in zipped form to save space and make it fit into
# GitHub's file size limits
write.csv(joined_table,
      file = gzfile("./output/final_compiled_counts/joined_count_data.csv.zip"),
      row.names = FALSE)

# also save the tibble as an RData file, which will be much faster to read in
# later and doesn't require the computationally expensive decompression step
# this is where the analysis will begin, using dplyr and ggplot
save(joined_table,
     file = "./output/final_compiled_counts/joined_count_data.RData")

# to keep going with differential expression analysis, the next step would
# be to move this data over into DESeq2 or edgeR to calculate differential
# expression and etc.

# some references for that:
# https://www.biostars.org/p/143458/
# http://www.gettinggeneticsdone.com/2015/12/tutorial-rna-seq-differential.html
