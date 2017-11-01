#!/usr/bin/env Rscript

##### command arguments for running script ADD Soon!

args = commandArgs(trailingOnly=TRUE)

barC = args[1]  #barcode target
path = args[2]  #path to the fastq files
barC_length = args[3] # expected seq length of the barcode.

#barC = "CO1"
#path <- paste("/Users/limeybean/Downloads/paired/",barC,sep="")
#barC_length = "500"

####################################################################################### process paired end reads

library(dada2); packageVersion("dada2")
library("seqRFLP")
library(plyr)
library("Biostrings")
library(reshape2)
####################
### Show where reads are kept
####################

list.files(path)

####################
### Filter and Trim
####################

# Sort ensures forward/reverse reads are in same order
fnFs <- sort(list.files(path, pattern="_Paired_1_sorted.fastq"))
fnRs <- sort(list.files(path, pattern="_Paired_2_sorted.fastq"))
# Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fastq
sample.names <- sapply(strsplit(fnFs, "_Paired"), `[`, 1)
# Specify the full path to the fnFs and fnRs
fnFs <- file.path(path, fnFs)
fnRs <- file.path(path, fnRs)

# look at plots
plotQualityProfile(fnFs[1:2])
plotQualityProfile(fnRs[1:2])

filt_path <- file.path(path, "filtered") # Place filtered files in filtered/ subdirectory
filtFs <- file.path(filt_path, paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(filt_path, paste0(sample.names, "_R_filt.fastq.gz"))


out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, minLen = 70,
                     maxN=0, maxEE=c(2,2), truncQ=0, rm.phix=TRUE,
                     compress=TRUE, multithread=TRUE) # On Windows set multithread=FALSE
head(out)

####################
### learn errors
####################

errF <- learnErrors(filtFs, multithread=TRUE)
errR <- learnErrors(filtRs, multithread=TRUE)
plotErrors(errF, nominalQ=TRUE)
plotErrors(errR, nominalQ=TRUE)

####################
### dereplicate
####################

derepFs <- derepFastq(filtFs, verbose=TRUE)
derepRs <- derepFastq(filtRs, verbose=TRUE)
# Name the derep-class objects by the sample names
names(derepFs) <- sample.names
names(derepRs) <- sample.names

####################
### sample inference
####################

dadaFs <- dada(derepFs, err=errF, multithread=TRUE)
dadaRs <- dada(derepRs, err=errR, multithread=TRUE)
dadaFs[[1]]

####################
### Merge Pairs
####################

mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs,returnRejects=TRUE, verbose=TRUE,minOverlap = 20, maxMismatch = 2)
# Inspect the merger data.frame from the first sample
head(mergers[[1]])


####################
### Make Sequence table
####################
seqtab <- makeSequenceTable(mergers)
dim(seqtab)
table(nchar(getSequences(seqtab)))

####################
### Remove Chimeras
####################

seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
dim(seqtab.nochim)

sum(seqtab.nochim)/sum(seqtab)

####################
### Track Reads
####################

getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(mergers, getN), rowSums(seqtab), rowSums(seqtab.nochim))
colnames(track) <- c("input", "filtered", "denoised", "merged", "tabled", "nonchim")
rownames(track) <- sample.names
head(track)



###########################################
### make output files for the merged reads
############################################

mergedbarC =  paste("merged_", barC , sep='')
mergedbarCseqnum = paste("merged_", barC , "_seq_number", sep = '')
nochime_fname.fasta = paste(path,"/", "nochim_merged",barC,".fasta", sep='')
nochime_fname.txt = paste(path,"/", "nochim_merged",barC,".txt", sep='')


makes.sense.seqtab.nochim <- t(seqtab.nochim)
makes.sense.seqtab.nochim <- cbind(seqeunces = rownames(makes.sense.seqtab.nochim), makes.sense.seqtab.nochim)
rownames(makes.sense.seqtab.nochim) <- NULL
makes.sense.seqtab.nochim <- as.data.frame(makes.sense.seqtab.nochim)
makes.sense.seqtab.nochim$seqnum <- 1:nrow(makes.sense.seqtab.nochim) 
namevector <- c(mergedbarC)
makes.sense.seqtab.nochim[ , namevector] <- mergedbarC
makes.sense.seqtab.nochim[[mergedbarCseqnum]] <- paste(makes.sense.seqtab.nochim[[mergedbarC]] , makes.sense.seqtab.nochim$seqnum,sep="_")
makes.sense.seqtab.nochim$seqnum <- NULL
makes.sense.seqtab.nochim[[mergedbarC]] <- NULL

nochim_merged  <- makes.sense.seqtab.nochim[,c(which(colnames(makes.sense.seqtab.nochim)==mergedbarCseqnum),which(colnames(makes.sense.seqtab.nochim)!=mergedbarCseqnum))]
nochim_merged_seq <- data.frame(nochim_merged[[mergedbarCseqnum]],nochim_merged$seqeunces)
nochim_merged_seq.fasta = dataframe2fas(nochim_merged_seq, file= nochime_fname.fasta)
write.table(nochim_merged, file = nochime_fname.txt, row.names=FALSE, sep="\t", quote=FALSE)


########################################################################################### processed unmerged paired end reads
### We need to remove the reads that did not merge and process those that did not overlap.  


##############################################################
# make dataframe from merged reads to pull ot those that failed!
##############################################################
try <- ldply(mergers)  #merge all dataframes resulting from merger

pairedsum.table <- cbind(try$.id,try$forward,try$reverse, try$abundance, try$accept)
colnames(pairedsum.table) <- c("id","forward", "reverse","abundance","accept")
test <- as.data.frame(pairedsum.table)
namevector <- c("sequenceF","sequenceR")
test[ , namevector] <- "NA"
pairedsum.unmerged.table <- subset(test, accept=="FALSE", select=c(id,forward,reverse,abundance,sequenceF,sequenceR))

##############################################################
## grab the correct F and R reads from the dadaF and dadR files
##############################################################

mine_unmerged <- function(df, seqs_list, forward_or_reverse) {
  seq_obj <- seqs_list[[(df[1])]]
  
  if (forward_or_reverse == "forward") {
    seq <- seq_obj$sequence[as.numeric(df[2])]
  } else {
    seq <- seq_obj$sequence[as.numeric(df[3])]
  }
  
  return(seq)
}

pairedsum.unmerged.table$sequenceF <- apply(pairedsum.unmerged.table, 1, function(x) 
  mine_unmerged(x, seqs_list = dadaFs, forward_or_reverse = "forward"))

pairedsum.unmerged.table$sequenceR <- apply(pairedsum.unmerged.table, 1, function(x) 
  mine_unmerged(x, seqs_list = dadaRs, forward_or_reverse = "reverse"))


##############################################################
# add seqeunce length for forwards and reversed to unmerged dereplicated data
##############################################################


charmatchesF<-nchar(gsub("[a-z]","",pairedsum.unmerged.table$sequenceF))
charmatchesR<-nchar(gsub("[a-z]","",pairedsum.unmerged.table$sequenceR))
pairedsum.unmerged.table$lengthF<-charmatchesF
pairedsum.unmerged.table$lengthR<-charmatchesR
pairedsum.unmerged.table$totalseq<-(pairedsum.unmerged.table$lengthF + pairedsum.unmerged.table$lengthR)

# add expected length of amplicon
pairedsum.unmerged.table$expected_amplicon <- barC_length
pairedsum.unmerged.table$keep[pairedsum.unmerged.table$totalseq>=pairedsum.unmerged.table$expected_amplicon] <- FALSE
pairedsum.unmerged.table$keep[pairedsum.unmerged.table$totalseq<pairedsum.unmerged.table$expected_amplicon] <- TRUE

pairedsum.unmerged.table$sequenceRc <- sapply((sapply(sapply(pairedsum.unmerged.table$sequenceR, DNAString), reverseComplement)), toString)

pairedsum.unmerged.table$NNNNNN <- "AAAAAAAAAATTCTTAAAAAAAAAA"
pairedsum.unmerged.table$sequenceF_N_Rrc <- paste(pairedsum.unmerged.table$sequenceF,pairedsum.unmerged.table$NNNNNN,pairedsum.unmerged.table$sequenceRc,sep="")

pairedsum.unmerged.dada2 <- subset(pairedsum.unmerged.table, keep=="TRUE")
pairedsum.unmerged.dada2  <- pairedsum.unmerged.dada2[,c(which(colnames(pairedsum.unmerged.dada2)=="sequenceF_N_Rrc"),which(colnames(pairedsum.unmerged.dada2)!="sequenceF_N_Rrc"))]

pairedsum.unmerged.dada2$forward <- NULL
pairedsum.unmerged.dada2$reverse <- NULL
pairedsum.unmerged.dada2$sequenceF <- NULL
pairedsum.unmerged.dada2$sequenceR <- NULL
pairedsum.unmerged.dada2$sequenceRc <- NULL
pairedsum.unmerged.dada2$lengthF <- NULL
pairedsum.unmerged.dada2$lengthR <- NULL
pairedsum.unmerged.dada2$keep <- NULL
pairedsum.unmerged.dada2$totalseq <- NULL
pairedsum.unmerged.dada2$expected_amplicon <- NULL
pairedsum.unmerged.dada2$name <- NULL
pairedsum.unmerged.dada2$NNNNNN  <- NULL



#View(pairedsum.unmerged.dada2 %>% group_by(sequenceF_N_Rrc) %>% summarise(sum= n()))

# convert abundance to a numeric vector
pairedsum.unmerged.dada2$abundance <- as.integer(pairedsum.unmerged.dada2$abundance)

# Spread the dataframe, and sum up the abundances per id
umerged.seq.tab <- dcast(pairedsum.unmerged.dada2, sequenceF_N_Rrc ~ id, fun.aggregate = sum)

rownames(umerged.seq.tab) <- umerged.seq.tab$sequenceF_N_Rrc

umerged.seq.tab$sequenceF_N_Rrc <- NULL
unmerged.seq.tab <-t(as.vector(umerged.seq.tab))

rownames <- row.names(unmerged.seq.tab)
tt<- apply(unmerged.seq.tab, 2 , as.integer)
row.names(tt) <- rownames

##########################################
# Run bimera detection on unmerged reads -> discard bimeras
##########################################
unmerged.seq.tab.nochim <- removeBimeraDenovo(tt, method="consensus", multithread=TRUE, verbose=TRUE)


##########################################
# Reformat data for Bowtie2
##########################################

# get barcode name info for numbering / nameing sequences in final fasta files / tables
unmergedbarC =  paste("unmerged_", barC , sep='')
unmergedbarCseqnum = paste("unmerged_", barC , "_seq_number", sep = '')


#transpose table, make rownames into a column -> sequnces
unmerged.seq.tab.nochim <- t(unmerged.seq.tab.nochim)
unmerged.seq.tab.nochim <- cbind(sequences = rownames(unmerged.seq.tab.nochim), unmerged.seq.tab.nochim)
rownames(unmerged.seq.tab.nochim) <- c()
unmerged.seq.tab.nochim <- as.data.frame(unmerged.seq.tab.nochim)

## split sequence column into two -> sequences F and Rrc (reverse complement), split on dummy sequnces used to merge reads initially
unmerged.seq.tab.nochim$sequencesF <- sapply(strsplit(as.character(unmerged.seq.tab.nochim$sequences),'AAAAAAAAAATTCTTAAAAAAAAAA'), "[", 1)
unmerged.seq.tab.nochim$sequencesRrc <- sapply(strsplit(as.character(unmerged.seq.tab.nochim$sequences),'AAAAAAAAAATTCTTAAAAAAAAAA'), "[", 2)
unmerged.seq.tab.nochim$sequences <- NULL

# reverse complement the reverse read for Bowtie2...
unmerged.seq.tab.nochim$sequencesR <- sapply((sapply(sapply(unmerged.seq.tab.nochim$sequencesRrc, DNAString), reverseComplement)), toString)
unmerged.seq.tab.nochim$sequencesRrc <- NULL

#Order the columns
unmerged.seq.tab.nochim <- unmerged.seq.tab.nochim[,c(which(colnames(unmerged.seq.tab.nochim)=="sequencesR"),which(colnames(unmerged.seq.tab.nochim)!="sequencesR"))]
unmerged.seq.tab.nochim <- unmerged.seq.tab.nochim[,c(which(colnames(unmerged.seq.tab.nochim)=="sequencesF"),which(colnames(unmerged.seq.tab.nochim)!="sequencesF"))]

# add new first column that gives the barcode, the type of read (unmerged), and the read number
unmerged.seq.tab.nochim$seqnum <- 1:nrow(unmerged.seq.tab.nochim) 
namevector <- c(unmergedbarC)
unmerged.seq.tab.nochim[ , namevector] <- unmergedbarC
unmerged.seq.tab.nochim[[unmergedbarCseqnum]] <- paste(unmerged.seq.tab.nochim[[unmergedbarC]] , unmerged.seq.tab.nochim$seqnum,sep="_")
unmerged.seq.tab.nochim$seqnum <- NULL
unmerged.seq.tab.nochim[[unmergedbarC]] <- NULL
unmerged.seq.tab.nochim <- unmerged.seq.tab.nochim[,c(which(colnames(unmerged.seq.tab.nochim)==unmergedbarCseqnum),which(colnames(unmerged.seq.tab.nochim)!=unmergedbarCseqnum))]

########################################
### make output files
#########################################

# make file paths
nochime_unfnameF.fasta = paste(path,"/", "nochim_unmerged",barC,"F",".fasta", sep='')
nochime_unfnameR.fasta = paste(path,"/", "nochim_unmerged",barC,"R",".fasta", sep='')
nochime_unfname.txt = paste(path,"/", "nochim_unmerged",barC,".txt", sep='')

# make data frames for the soon to be made fasta files  with reads and read names
nochim_unmerged_seq_F <- data.frame(unmerged.seq.tab.nochim[[unmergedbarCseqnum]],unmerged.seq.tab.nochim$sequencesF)
nochim_unmerged_seq_R <- data.frame(unmerged.seq.tab.nochim[[unmergedbarCseqnum]],unmerged.seq.tab.nochim$sequencesR)

# export the fasta file ready dataframes
nochim_unmerged_seq_F.fasta = dataframe2fas(nochim_unmerged_seq_F, file= nochime_unfnameF.fasta)
nochim_unmerged_seq_R.fasta = dataframe2fas(nochim_unmerged_seq_R, file= nochime_unfnameR.fasta)

# write summary table 
write.table(unmerged.seq.tab.nochim, file = nochime_unfname.txt, row.names=FALSE, sep="\t", quote=FALSE)





