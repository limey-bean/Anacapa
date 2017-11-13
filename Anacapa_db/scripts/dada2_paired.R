#!/usr/bin/env Rscript

##### command arguments for running script ADD Soon!
## the dada2 commands come directly from https://benjjneb.github.io/dada2/tutorial.html

 args = commandArgs(trailingOnly=TRUE)
# 
 barC = args[1]  #barcode target
 odirpath = args[2]  #path to the fastq files
 barC_length = args[3] # expected seq length of the barcode.

## path to output
path = paste(odirpath, "/paired/" ,barC, sep='')
mergedoutpath=paste(odirpath, "/dada2_out/paired/merged/",barC, sep='')
unmergedoutpath=paste(odirpath, "/dada2_out/paired/unmerged/",barC, sep='')

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
fnFs <- sort(list.files(path, pattern="_Paired_1_sorted.fastq"))
fnRs <- sort(list.files(path, pattern="_Paired_2_sorted.fastq"))
# Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fastq
sample.names <- sapply(strsplit(fnFs, "_Paired"), `[`, 1)
# Specify the full path to the fnFs and fnRs
fnFs <- file.path(path, fnFs)
fnRs <- file.path(path, fnRs)

# save plots in a plots directory
dir.create(path = paste0(path, "/plots"))
plotQualityProfile(fnFs[1:2]) + ggsave(filename = paste0(path, "/plots/forward_qualities.pdf"))
plotQualityProfile(fnRs[1:2]) + ggsave(filename = paste0(path, "/plots/reverse_qualities.pdf"))

filt_path <- file.path(path, "filtered") # Place filtered files in filtered/ subdirectory
filtFs <- file.path(filt_path, paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(filt_path, paste0(sample.names, "_R_filt.fastq.gz"))


out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, minLen = 70,
                     maxN=0, maxEE=c(2,2), truncQ=0, rm.phix=TRUE,
                     compress=F, multithread=F) # On Windows set multithread=FALSE
head(out)

####################
### learn errors
####################

errF <- learnErrors(filtFs, multithread=TRUE)
errR <- learnErrors(filtRs, multithread=TRUE)
plotErrors(errF, nominalQ=TRUE) + ggsave(filename = paste0(path, "/plots/forward_error_profile.pdf"))
plotErrors(errR, nominalQ=TRUE) + ggsave(filename = paste0(path, "/plots/reverse_error_profile.pdf"))

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

cat("Distribution of merged sequence lengths:"); table(nchar(getSequences(seqtab)))

####################
### Remove Chimeras
####################

seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
dim(seqtab.nochim)

cat("Proportion of sequences kept after chimera removal: "); cat(sum(seqtab.nochim)/sum(seqtab))

####################
### Track Reads
####################

getN <- function(x) {sum(getUniques(x))}
track <- data.frame(input = out[,"reads.in"], filtered = out[,"reads.out"], 
                    denoised =sapply(dadaFs, getN), 
                    merged = sapply(mergers, getN), 
                    tabled = rowSums(seqtab),
                    nochim = rowSums(seqtab.nochim))
rownames(track) <- sample.names
track

###########################################
### make output files (FASTA and abundance matrix) for the merged reads
############################################

# Note to Emily- these variable names should have an explanation
# I.e. how will each be used?
mergedbarC =  paste("merged_", barC , sep='')
mergedbarCseqnum = paste("merged_", barC , "_seq_number", sep = '')
nochime_fname.fasta = paste(path,"/", "nochim_merged",barC,".fasta", sep='')
nochime_fname.txt = paste(path,"/", "nochim_merged",barC,".txt", sep='')

makes.sense.seqtab.nochim <- t(seqtab.nochim)
nochim_merged <- makes.sense.seqtab.nochim %>% data.frame %>%
       rownames_to_column %>% rename(sequence = rowname) %>% # make sequences into a column
       mutate(!!mergedbarCseqnum := paste0(mergedbarC,"_",row_number())) %>% # Make a new column w seq number
       select(!!mergedbarCseqnum,everything()) # reorder the columns

# Save this table 
write.table(nochim_merged, file = nochime_fname.txt, row.names=FALSE, sep="\t", quote=FALSE)

# Make a fasta file out of the sequences in this table
nochim_merged_seq <- nochim_merged %>% select(!!mergedbarCseqnum, sequence)

#
nochim_merged_seq.fasta = dataframe2fas(nochim_merged_seq, file= nochime_fname.fasta)



########################################################################################### 
### processed unmerged paired end reads
### We need to remove the reads that did not merge and process those that did not overlap.  


##############################################################
# make dataframe from merged reads to pull ot those that failed!
##############################################################
try <- ldply(mergers)  #merge all dataframes resulting from merger


pairedsum.unmerged.table <- try %>% select(id = .id, forward, reverse, abundance, accept) %>% data.frame %>% 
  filter(accept == FALSE) %>% select(-accept) # remove the ones that worked well during merge

##############################################################
## grab the correct F and R reads from the dadaF and dadR files
##############################################################

mine_unmerged <- function(df, seqs_list, forward_or_reverse) {
  seq_obj <- seqs_list[[(df[1])]]
  
  # Make sure nothing silly is happening
  if (!(forward_or_reverse %in% c("forward","reverse"))){
    stop("forward_or_reverse must be specified as 'forward' or 'reverse'")
  }
  
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


pairedsum.unmerged.table$lengthF <- nchar(gsub("[a-z]","",pairedsum.unmerged.table$sequenceF))
pairedsum.unmerged.table$lengthR <- nchar(gsub("[a-z]","",pairedsum.unmerged.table$sequenceR))
pairedsum.unmerged.table$totalseq <- pairedsum.unmerged.table$lengthF + pairedsum.unmerged.table$lengthR

# add expected length of amplicon
pairedsum.unmerged.table$keep[pairedsum.unmerged.table$totalseq>=barC_length] <- FALSE
pairedsum.unmerged.table$keep[pairedsum.unmerged.table$totalseq<barC_length] <- TRUE

# Get the reverse complement of the R sequence
pairedsum.unmerged.table$sequenceRc <- sapply(sapply(sapply(pairedsum.unmerged.table$sequenceR, DNAString), reverseComplement), toString)

# pairedsum.unmerged.table$NNNNNN <- "AAAAAAAAAATTCTTAAAAAAAAAA"
pairedsum.unmerged.table$sequenceF_N_Rrc <- paste(pairedsum.unmerged.table$sequenceF,"AAAAAAAAAATTCTTAAAAAAAAAA",pairedsum.unmerged.table$sequenceRc,sep="")


pairedsum.unmerged.dada2 <- pairedsum.unmerged.table %>% filter(keep == TRUE) %>% select(sequenceF_N_Rrc, id, abundance)


#View(pairedsum.unmerged.dada2 %>% group_by(sequenceF_N_Rrc) %>% summarise(sum= n()))

# Spread the dataframe, and sum up the abundances per id
unmerged.seq.tab <- dcast(pairedsum.unmerged.dada2, sequenceF_N_Rrc ~ id, fun.aggregate = sum) %>% 
  data.frame %>% column_to_rownames( "sequenceF_N_Rrc") %>% t()


##########################################
# Run bimera detection on unmerged reads -> discard bimeras
##########################################
unmerged.seq.tab.nochim <- removeBimeraDenovo(unmerged.seq.tab, method="consensus", multithread=TRUE, verbose=TRUE)


##########################################
# Reformat data for Bowtie2
##########################################

# get barcode name info for numbering / nameing sequences in final fasta files / tables
unmergedbarC =  paste("unmerged_", barC , sep='')
unmergedbarCseqnum = paste("unmerged_", barC , "_seq_number", sep = '')


#transpose table, make rownames into a column -> sequnces
unmerged.seq.tab.nochim <- unmerged.seq.tab.nochim %>% t %>% data.frame %>% rownames_to_column("sequences")


## split sequence column into two -> sequences F and Rrc (reverse complement), split on dummy sequnces used to merge reads initially
unmerged.seq.tab.nochim$sequencesF <- sapply(strsplit(as.character(unmerged.seq.tab.nochim$sequences),'AAAAAAAAAATTCTTAAAAAAAAAA'), "[", 1)
unmerged.seq.tab.nochim$sequencesRrc <- sapply(strsplit(as.character(unmerged.seq.tab.nochim$sequences),'AAAAAAAAAATTCTTAAAAAAAAAA'), "[", 2)
unmerged.seq.tab.nochim$sequences <- NULL

# reverse complement the reverse read for Bowtie2...
unmerged.seq.tab.nochim$sequencesR <- sapply((sapply(sapply(unmerged.seq.tab.nochim$sequencesRrc, DNAString), reverseComplement)), toString)
unmerged.seq.tab.nochim$sequencesRrc <- NULL

# Reorder the columns (seqF, seqR, all the samples)

unmerged.seq.tab.nochim <- unmerged.seq.tab.nochim %>% select(sequencesF, sequencesR, everything()) %>% 
  mutate(!!unmergedbarCseqnum := paste0(unmergedbarC, "_", row_number())) %>% select(!!unmergedbarCseqnum, everything())


# ########################################
# ### make output files
# #########################################
# 
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
