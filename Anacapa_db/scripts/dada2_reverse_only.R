#!/usr/bin/env Rscript

##### command arguments for running script!

 args = commandArgs(trailingOnly=TRUE)

 barC = arg[1]  #barcode target
 path = arg[2]  #path to the fastq files

# barC = "CO1"
# path <- paste("/Users/limeybean/Downloads/unpaired_2/",barC,sep="")


############################################################################################Forward Reads reads

library(dada2); packageVersion("dada2")
library("seqRFLP")

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

sum(seqtabRo.nochim)/sum(seqtab)

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
nochime_r_fname.fasta = paste(path,"/", "nochim_reverse",barC,".fasta", sep='')
nochime_r_fname.txt = paste(path,"/", "nochim_reverse",barC,".txt", sep='')

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

