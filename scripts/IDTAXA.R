#!/usr/bin/env Rscript

#usage: Rscript --vanilla IDTAXA.R /path-to-otu-seqs/ /path-to-training-set/
#alternative usage: R CMD BATCH IDTAXA.R (but need to specify the path directory in the script)

#define argumens from command line
args <- commandArgs(trailingOnly=TRUE)

DirPath <- args[1]
TrainingSet <- args[2]

#install DECIPHER which is designed to complement other tools that are part of BioConductor, a suite of packages for the R statistical programming language
#	install.packages("BiocManager")
#	BiocManager::install("DECIPHER")
	#get error: package ‘BiocManager’ is not available (for R version 3.4.4)
#this works for R version 3.4.4 in Ubuntu 18.0
#	source("https://bioconductor.org/biocLite.R")
#	biocLite("DECIPHER")

# load the DECIPHER library in R
library(DECIPHER)

# specify the path to the FASTA file (in quotes)
setwd(DirPath)
fas <- "all.otus_BLCA-unclassified.fasta"

# load the sequences from the file
seqs <- readDNAStringSet(fas)

# remove any gaps (if needed)
seqs <- RemoveGaps(seqs)

# load a training set object (trainingSet)
# see http://DECIPHER.codes/Downloads.html
load(TrainingSet)

# classify the sequences
ids <- IdTaxa(seqs, trainingSet,type="extended", strand="top", threshold=70, bootstraps=100, processors=88)

# look at the results
print(ids)
#plot(ids, trainingSet)

#write tables
output <- sapply(ids, function (id) { paste(id$rank, id$taxon, id$confidence, sep=";", collapse=";") })

write.table(output, file = "all.otus_BLCA-unclassified_taxonomy_IDTAXA_0.7_RDP-V16-modified_original.txt", sep="\t")

#quit and save the R objects for future view
quit(save="no")
