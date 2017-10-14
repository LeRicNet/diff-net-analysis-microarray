###############
# Description #
###############

# AUTHOR: Eric Prince
# DATE: 2017-10-13
# DETAILS: This is a work-flow for using the dna package to investigate network relationships between RNA sequence data.
#         One pair of networks will look at the effect of shRNA knockdown and the other will look at the effect of drug treatment.

library(dna)
library(tidyverse)
library(magrittr)
library(reshape2)

#################
# Sample Layout #
#################

# In-vitro cells
# Sample ID   Treatment
#----------------------
# 1, 2, 3     shNULL
# 4, 5, 6     shTARGET
# 7, 8, 9     NoTx
# 10, 11, 12  Tx


############################
# Import & Format Seq Data #
############################

directoryList <- list.dirs("/Volumes/EP/__masked__/for")
directoryList <- directoryList[c(2:5, 7:14)] # removing unwanted directories; could do this in list.dirs..

# bring in all of the FPKM files 
fpkmList <- list()
for (i in 1:length(directoryList)) {
  file.path = list.files(path = directoryList[i], pattern = "genes.fpkm_tracking", full.names = TRUE)
  sampleID = sapply(strsplit(file.path, "for/"), '[', 2)
  sampleID = sapply(strsplit(as.character(sampleID), "_"), '[', 1)
  fpkmList[[paste0("sample", sampleID)]] <- read_tsv(file.path)
}

# collapse the list into one tibble
fpkm <- bind_rows(fpkmList, .id = "sample")

# filter OK reads, select relevant attributes, and arrange data to meet args(test.individual.genes) requirements
fpkm_spread <- fpkm %>%
  filter(FPKM_status == "OK") %>%
  select(sample, gene_short_name, FPKM) %>%
  dcast(gene_short_name ~ sample, mean) %>%
  t() %>%
  as_data_frame()

# turn first row into names, then remove row 1
as.character(fpkm_spread[1,]) -> names(fpkm_spread)
names(fpkm_spread)
fpkm_spread <- fpkm_spread[c(2:nrow(fpkm_spread)),]

# Convert to numeric matrix and change NA to 0
fpkm_spread <- apply(fpkm_spread, 2, as.numeric)


###################
# Quality Control #
###################

# Checking for NAs
if (any(is.na(fpkm_spread))) {
  fpkm_spread[is.na(fpkm_spread)] <- 0
} else { print("No <NA>s in Dataset; <NA> converted to 0") }

# Set row/column color scales and check heatmap for shNull/shTarget
sg1.idx <- c(1:6)
rc <- rainbow(nrow(fpkm_spread[sg1.idx,]), start = 0, end = 0.3)
cc <- rainbow(ncol(fpkm_spread[sg1.idx,]), start = 0, end = 0.3)
heatmap(fpkm_spread[sg1.idx,], Colv=NA, scale = "column", RowSideColors = rc, ColSideColors = cc)
# Heatmap for shNull/shTarget suggest that samples 1 and 4 have been incorrectly labeled.
null.idx <- c(2,3,4)
target.idx <- c(1, 5, 6)


# Set row/column color scales and check heatmap for NoTx/Tx
sg2.idx <- c(7:12)
rc <- rainbow(nrow(fpkm_spread[sg2.idx,]), start = 0, end = 0.3)
cc <- rainbow(ncol(fpkm_spread[sg2.idx,]), start = 0, end = 0.3)
heatmap(fpkm_spread[sg2.idx,], Colv=NA, scale = "column", RowSideColors = rc, ColSideColors = cc)
# Samples appear to be okay
notx.idx <- c(7:9)
tx.idx <- c(10:12)


#####################
# Diff Net Analysis #
#####################

# test individual genes (tig)
# The test.individual.genes calls are very long; estimated at approximately 45 hours each on a 2015 iMac with i5[?] and 32GB RAM
#                                                             ***************************
shNull = fpkm_spread[null.idx,]
shTarget = fpkm_spread[target.idx,]
tig.res.null_target <- test.individual.genes(X1 = shNull, 
                                 X2 = shSIRT2, 
                                 scores = "PLS", 
                                 distance = "abs", 
                                 num.permutations = 1000, 
                                 check.networks = TRUE)

NoTx = fpkm_spread[notx.idx,]
Tx = fpkm_spread[tx.idx,]
tig.res.notx_tx <- test.individual.genes(X1 = NoTx, 
                                         X2 = Tx, 
                                         scores = "PLS", 
                                         distance = "abs", 
                                         num.permutations = 1000, 
                                         check.networks = TRUE)

