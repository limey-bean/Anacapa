#!/usr/bin/env Rscript

##### command arguments for running script!
## the dada2 commands are sight modifications of those from https://benjjneb.github.io/dada2/tutorial.html
 args = commandArgs(trailingOnly=TRUE)
# 
 barC = args[1]  #barcode target
 odirpath = args[2]  #path to the fastq files
 barC_length = args[3] # expected seq length of the barcode.

## path to output
path = paste(odirpath,  "/unpaired_F/" ,barC, sep='')
outpath=paste(odirpath, "/", barC, "/dada2_bowtie2/unpaired_F", sep='')


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

# Sort ensures forward/reverse reads are in same order
fnFo <- sort(list.files(path, pattern="_Paired_1_singletons.fastq"))
# Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fastq
sample.names <- sapply(strsplit(fnFo, "_Paired_1_singletons.fastq"), `[`, 1)
# Specify the full path to the fnFo and fnRs
fnFo <- file.path(path, fnFo)

# look at plots
plotQualityProfile(fnFo[1:2])

filt_path <- file.path(path, "filtered") # Place filtered files in filtered/ subdirectory
filtFos <- file.path(filt_path, paste0(sample.names, "_F_filto.fastq.gz"))


out <- filterAndTrim(fnFo, filtFos, minLen = 70,
                     maxN=0, maxEE=c(2), truncQ=0, rm.phix=TRUE,
                     compress=TRUE, multithread=TRUE) # On Windows set multithread=FALSE
head(out)

####################
### learn errors
####################

errFo <- learnErrors(filtFos, multithread=TRUE)
plotErrors(errFo, nominalQ=TRUE)

####################
### dereplicate
####################

derepFso <- derepFastq(filtFos, verbose=TRUE)
# Name the derep-class objects by the sample names
names(derepFso) <- sample.names

####################
### sample inference
####################

dadaFos <- dada(derepFso, err=errFo, multithread=TRUE)
dadaFos[[1]]



####################
### Make Sequence table
####################
seqtabFo <- makeSequenceTable(derepFso)
dim(seqtabFo)
table(nchar(getSequences(seqtabFo)))

####################
### Remove Chimeras
####################

seqtabFo.nochim <- removeBimeraDenovo(seqtabFo, method="consensus", multithread=TRUE, verbose=TRUE)
dim(seqtabFo.nochim)

sum(seqtabFo.nochim)/sum(seqtabFo)

####################
### Track Reads
####################

getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFos, getN), rowSums(seqtabFo), rowSums(seqtabFo.nochim))
colnames(track) <- c("input", "filtered", "denoised", "tabled", "nonchim")
rownames(track) <- sample.names
head(track)

###########################################
### make output files for the merged reads
############################################

# make variables, and file paths for output
forward_barC =  paste("forward_", barC , sep='')
forwardbarCseqnum = paste("forward_", barC , "_seq_number", sep = '')
nochime_f_fname.fasta = paste(outpath,"/", "nochim_forward",barC,".fasta", sep='')
nochime_f_fname.txt = paste(outpath,"/", "nochim_forward",barC,".txt", sep='')

#modify table so that there are unique specific names for each read, and that also include sample realtive abundance.
makes.sense.seqtabFo.nochim <- t(seqtabFo.nochim)
makes.sense.seqtabFo.nochim <- cbind(sequences = rownames(makes.sense.seqtabFo.nochim), makes.sense.seqtabFo.nochim)
rownames(makes.sense.seqtabFo.nochim) <- NULL
makes.sense.seqtabFo.nochim <- as.data.frame(makes.sense.seqtabFo.nochim)
makes.sense.seqtabFo.nochim$seqnum <- 1:nrow(makes.sense.seqtabFo.nochim) 
namevector <- c(forward_barC)
makes.sense.seqtabFo.nochim[ , namevector] <- forward_barC
makes.sense.seqtabFo.nochim[[forwardbarCseqnum]] <- paste(makes.sense.seqtabFo.nochim[[forward_barC]] , makes.sense.seqtabFo.nochim$seqnum,sep="_")
makes.sense.seqtabFo.nochim$seqnum <- NULL
makes.sense.seqtabFo.nochim[[forward_barC]] <- NULL

# final data transfermations, and export of reads as fasta files and accompanying data table.
nochim_forward  <- makes.sense.seqtabFo.nochim[,c(which(colnames(makes.sense.seqtabFo.nochim)==forwardbarCseqnum),which(colnames(makes.sense.seqtabFo.nochim)!=forwardbarCseqnum))]
nochim_forward_seq <- data.frame(nochim_forward[[forwardbarCseqnum]],nochim_forward$sequences)
nochim_forward_seq.fasta = dataframe2fas(nochim_forward_seq, file= nochime_f_fname.fasta)
write.table(nochim_forward, file = nochime_f_fname.txt, row.names=FALSE, sep="\t", quote=FALSE)
