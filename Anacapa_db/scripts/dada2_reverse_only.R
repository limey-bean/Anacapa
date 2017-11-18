#!/usr/bin/env Rscript

##### command arguments for running script!
## the dada2 commands are sight modifications of those from https://benjjneb.github.io/dada2/tutorial.html
 args = commandArgs(trailingOnly=TRUE)
# 
 barC = args[1]  #barcode target
 odirpath = args[2]  #path to the fastq files
 barC_length = args[3] # expected seq length of the barcode.

## path to output
path = paste(odirpath,  "/unpaired_R/",barC,  sep='')
outpath=paste(odirpath,"/", barC, "/dada2_bowtie2/unpaired_R/", sep='')


############################################################################################Forward Reads reads

# Install packages that are not currently installed -----
# This chunk may need attention, putting it here as a starting point - gsk
 .cran_packages  <-  c("ggplot2", "plyr", "dplyr","seqRFLP", "reshape2", "tibble")
 .bioc_packages <- c("phyloseq", "genefilter", "impute", "Biostrings", "dada2")
 
 .inst <- .cran_packages %in% installed.packages()
 if (any(!.inst)) {
   install.packages(.cran_packages[!.inst], repos = "http://cran.rstudio.com/")
 }
# 
 .inst <- .bioc_packages %in% installed.packages()
 if (any(!.inst)) {
   source("http://bioconductor.org/biocLite.R")
   biocLite(.bioc_packages[!.inst])
 }
####################################################################################### process paired end reads

library("dada2")
cat(paste("dada2 package version:", packageVersion("dada2")))
library("seqRFLP")
library("plyr")
library("Biostrings")
library("reshape2")
library("dplyr")
library("tibble")
library("ggplot2")
####################
### Show where reads are kept
####################

list.files(path)

####################
### Filter and Trim
####################

# Sort ensures reverse/reverse reads are in same order
fnRo <- sort(list.files(path, pattern="_Paired_2_singletons.fastq"))
# Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fastq
sample.names <- sapply(strsplit(fnRo, "_Paired_2_singletons.fastq"), `[`, 1)
# Specify the full path to the fnRo and fnRs
fnRo <- file.path(path, fnRo)

# look at plots
plotQualityProfile(fnRo[1:2])

filt_path <- file.path(path, "filtered") # Place filtered files in filtered/ subdirectory
filtRos <- file.path(filt_path, paste0(sample.names, "_r_filto.fastq.gz"))


out <- filterAndTrim(fnRo, filtRos, minLen = 70,
                     maxN=0, maxEE=c(2), truncQ=0, rm.phix=TRUE,
                     compress=TRUE, multithread=TRUE) # On Windows set multithread=FALSE
head(out)

####################
### learn errors
####################

errRo <- learnErrors(filtRos, multithread=TRUE)
plotErrors(errRo, nominalQ=TRUE)

####################
### dereplicate
####################

derepFso <- derepFastq(filtRos, verbose=TRUE)
# Name the derep-class objects by the sample names
names(derepFso) <- sample.names

####################
### sample inference
####################

dadaRos <- dada(derepFso, err=errRo, multithread=TRUE)
dadaRos[[1]]



####################
### Make Sequence table
####################
seqtabRo <- makeSequenceTable(derepFso)
dim(seqtabRo)
table(nchar(getSequences(seqtabRo)))

####################
### Remove Chimeras
####################

seqtabRo.nochim <- removeBimeraDenovo(seqtabRo, method="consensus", multithread=TRUE, verbose=TRUE)
dim(seqtabRo.nochim)

sum(seqtabRo.nochim)/sum(seqtabRo)

####################
### Track Reads
####################

getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaRos, getN), rowSums(seqtabRo), rowSums(seqtabRo.nochim))
colnames(track) <- c("input", "filtered", "denoised", "tabled", "nonchim")
rownames(track) <- sample.names
head(track)

###########################################
### make output files for the merged reads
############################################

# make variables, and file paths for output
reverse_barC =  paste("reverse_", barC , sep='')
reversebarCseqnum = paste("reverse_", barC , "_seq_number", sep = '')
nochime_r_fname.fasta = paste(outpath,"/", "nochim_reverse",barC,".fasta", sep='')
nochime_r_fname.txt = paste(outpath,"/", "nochim_reverse",barC,".txt", sep='')

#modify table so that there are unique specific names for each read, and that also include sample realtive abundance.
makes.sense.seqtabRo.nochim <- t(seqtabRo.nochim)
makes.sense.seqtabRo.nochim <- cbind(seqeunces = rownames(makes.sense.seqtabRo.nochim), makes.sense.seqtabRo.nochim)
rownames(makes.sense.seqtabRo.nochim) <- NULL
makes.sense.seqtabRo.nochim <- as.data.frame(makes.sense.seqtabRo.nochim)
makes.sense.seqtabRo.nochim$seqnum <- 1:nrow(makes.sense.seqtabRo.nochim) 
namevector <- c(reverse_barC)
makes.sense.seqtabRo.nochim[ , namevector] <- reverse_barC
makes.sense.seqtabRo.nochim[[reversebarCseqnum]] <- paste(makes.sense.seqtabRo.nochim[[reverse_barC]] , makes.sense.seqtabRo.nochim$seqnum,sep="_")
makes.sense.seqtabRo.nochim$seqnum <- NULL
makes.sense.seqtabRo.nochim[[reverse_barC]] <- NULL

# final data transfermations, and export of reads as fasta files and accompanying data table.
nochim_reverse  <- makes.sense.seqtabRo.nochim[,c(which(colnames(makes.sense.seqtabRo.nochim)==reversebarCseqnum),which(colnames(makes.sense.seqtabRo.nochim)!=reversebarCseqnum))]
nochim_reverse_seq <- data.frame(nochim_reverse[[reversebarCseqnum]],nochim_reverse$seqeunces)
nochim_reverse_seq.fasta = dataframe2fas(nochim_reverse_seq, file= nochime_r_fname.fasta)
write.table(nochim_reverse, file = nochime_r_fname.txt, row.names=FALSE, sep="\t", quote=FALSE)
