#!/usr/bin/env Rscript
# Interpret the input variables -----

##### command arguments for running script ADD Soon!
## the dada2 commands come directly from https://benjjneb.github.io/dada2/tutorial.html


args = commandArgs(trailingOnly=TRUE)
# check that there's four arguments
if (length(args) != 5) {
  stop("please make sure there are 5 arguments: the barcode name, the path to fastq files, the expected seq length of the barcode, and the type of reads being processed, minimum ASV abundance")
}


barC = args[1]  #barcode target
odirpath = args[2]  #path to the fastq files
barC_length = args[3] # expected seq length of the barcode.
paired_or_not = args[4] # type of reads- should be "paired", "forward", or "reverse
min_asv_abundance = as.numeric(args[5]) # minimum number of times an ASV needs to appear to be kept in output files

# confirm that the user has specified paired_or_not properly
if (!(paired_or_not %in% c("paired", "forward", "reverse"))) {
  cat("Please specify sequence type as 'paired', 'forward', or 'reverse'")
  quit()
}

## path to output

if (paired_or_not == "paired") {
  path = paste(odirpath,"/", barC, "/", barC, "_sort_by_read_type/paired",  sep='')
  mergedoutpath=paste(odirpath,"/", barC, "/", barC, "dada2_out", sep='')
  unmergedoutpath=paste(odirpath,"/", barC, "/", barC, "dada2_out", sep='')
} else if(paired_or_not == "forward") {
  path = paste(odirpath,"/", barC, "/", barC, "_sort_by_read_type/unpaired_F", sep='')
  outpath=paste(odirpath,"/", barC, "/", barC, "dada2_out", sep='')
} else {
  path = paste(odirpath,"/", barC, "/", barC, "_sort_by_read_type/unpaired_R", sep='')
  outpath=paste(odirpath,"/", barC, "/", barC, "dada2_out", sep='')

}

# Confirm that the user has Write access to the path
if (file.access(path, mode = 2) != 0) {
  stop("Please make sure that you have write access to the supplied path.")
}

# Manage packages -----

#1. Download packages from CRAN
.cran_packages  <-  c("ggplot2", "plyr", "dplyr","seqRFLP", "reshape2", "tibble", "devtools", "Matrix", "mgcv")
.inst <- .cran_packages %in% installed.packages()
if (any(!.inst)) {
  install.packages(.cran_packages[!.inst], repos = "http://cran.rstudio.com/")
}

# 2. Download packages from biocLite
.bioc_packages <- c("phyloseq", "genefilter", "impute", "Biostrings")
.inst <- .bioc_packages %in% installed.packages()
if (any(!.inst)) {
  source("http://bioconductor.org/biocLite.R")
  biocLite(.bioc_packages[!.inst])
}

.dada_version = "1.6.0"
.dada_version_gh = "v1.6"
if("dada2" %in% installed.packages()){
  if(packageVersion("dada2") == .dada_version) {
    cat("congrats, right version of dada2")
  } else {
    devtools::install_github("benjjneb/dada2", ref=.dada_version_gh)
  }
}

if(!("dada2" %in% installed.packages())){
  # if the user doesn't have dada2 installed, install version 1.6 from github
  devtools::install_github("benjjneb/dada2", ref=.dada_version_gh)
}

library("dada2")
cat(paste("dada2 package version:", packageVersion("dada2")))
if(packageVersion("dada2") != '1.6.0') {
  stop("Please make sure you have dada version ", .dada_version, " installed")
}

library("seqRFLP")
library("plyr")
library("Biostrings")
library("reshape2")
library("dplyr")
library("tibble")
library("ggplot2")

# Set up paths to files ----------

list.files(path)

if(paired_or_not == "paired") {
  fnFs <- sort(list.files(path, pattern="_Paired_1_pairs_R1.fastq"))
  fnRs <- sort(list.files(path, pattern="_Paired_2_pairs_R2.fastq"))
  all_sample_names <- sapply(strsplit(fnFs, "_Paired"), `[`, 1)
  fnFs <- file.path(path, fnFs)
  fnRs <- file.path(path, fnRs)

} else {
  fnFs <- sort(list.files(path, pattern="_Paired_\\d_singles.fastq"))
  all_sample_names <- sapply(strsplit(fnFs, "_Paired_\\d_singles.fastq"), `[`, 1)
  fnFs <- file.path(path, fnFs)
}

# Make plots of sequence quality ------
#dir.create(path = paste0(path, "/plots"))

# Make a vector of two randomly sampled files with non-zero amounts of sequence data
#files_for_qual_plot <- sample(fnFs[file.size(fnFs) > 0], 2, replace = F)
#plotQualityProfile(files_for_qual_plot) + ggsave(filename = paste0(path, "/plots/forward_qualities.pdf"))
#if(paired_or_not == "paired") {
#  files_for_qual_plot_R <- sample(fnRs[file.size(fnRs) > 0], 2, replace = F)
#  plotQualityProfile(files_for_qual_plot_R) + ggsave(filename = paste0(path, "/plots/reverse_qualities.pdf"))
#}


# Make the path to which filtered sequences should be outputted ---------
filt_path <- file.path(path, "filtered")

if(paired_or_not == "paired") {
  filtered_seqs_name <- file.path(filt_path, paste0(all_sample_names, "_F_filt.fastq.gz"))
  filtered_seqs_name_R <- file.path(filt_path, paste0(all_sample_names, "_R_filt.fastq.gz"))
} else if (paired_or_not == "forward") {
  filtered_seqs_name <- file.path(filt_path, paste0(all_sample_names, "_F_filt.fastq.gz"))
} else {
  filtered_seqs_name <- file.path(filt_path, paste0(all_sample_names, "_R_filt.fastq.gz"))
}

# Run the filtering step ------
if(paired_or_not == "paired") {
  filtered_seqs <- filterAndTrim(fnFs, filtered_seqs_name, fnRs, filtered_seqs_name_R, minLen = 10,
                                 maxN=0, maxEE=c(2,2), truncQ=0, rm.phix=TRUE,
                                 compress=F,matchIDs=TRUE, multithread=TRUE) # On Windows set multithread=FALSE
} else {
  filtered_seqs <- filterAndTrim(fnFs, filtered_seqs_name, minLen = 10,
                                 maxN=0, maxEE=c(2), truncQ=0, rm.phix=TRUE,
                                 compress=TRUE, multithread=TRUE) # On Windows set multithread=FALSE

}
head(filtered_seqs)


# Check for cases where the filtering step left zero sequences in the output file----

if (paired_or_not == "paired") {
  exists <- file.exists(filtered_seqs_name) & file.exists(filtered_seqs_name_R)
  filtered_seqs_name <- filtered_seqs_name[exists]
  filtered_seqs_name_R <- filtered_seqs_name_R[exists]
  filtered_sample_names <- sapply(strsplit(basename(filtered_seqs_name), "_F_filt.fastq.gz"), `[`, 1)
} else if (paired_or_not == "forward") {
  exists <- file.exists(filtered_seqs_name)
  filtered_seqs_name <- filtered_seqs_name[exists]
  filtered_sample_names <- sapply(strsplit(basename(filtered_seqs_name), "_F_filt.fastq.gz"), `[`, 1)
} else {
  exists <- file.exists(filtered_seqs_name)
  filtered_seqs_name <- filtered_seqs_name[exists]
  filtered_sample_names <- sapply(strsplit(basename(filtered_seqs_name), "_R_filt.fastq.gz"), `[`, 1)
}

# Learn errors and save plots---------
if (paired_or_not == "paired"){
  error_profile <- learnErrors(filtered_seqs_name, multithread=TRUE)
  error_profile_R <- learnErrors(filtered_seqs_name_R, multithread=TRUE)
 # plotErrors(error_profile, nominalQ=TRUE) + ggsave(filename = paste0(path, "/plots/forward_error_profile.pdf"))
 # plotErrors(error_profile_R, nominalQ=TRUE) + ggsave(filename = paste0(path, "/plots/reverse_error_profile.pdf"))
} else {
  error_profile <- learnErrors(filtered_seqs_name, multithread=TRUE)
 #plotErrors(error_profile, nominalQ=TRUE) + ggsave(filename = paste0(path, "/plots/singleton_error_profile.pdf"))
}

# Dereplicate sequences -----------
if (paired_or_not == "paired"){
  derep_seqs <- derepFastq(filtered_seqs_name, verbose=TRUE)
  derep_seqs_R <- derepFastq(filtered_seqs_name_R, verbose=TRUE)
  names(derep_seqs) <- filtered_sample_names
  names(derep_seqs_R) <- filtered_sample_names
} else {
  derep_seqs <- derepFastq(filtered_seqs_name, verbose=TRUE)
  names(derep_seqs) <- filtered_sample_names
}

# Run dada on the dereplicated sequences ------
if (paired_or_not == "paired"){
  dada_output <- dada(derep_seqs, err=error_profile, multithread=TRUE)
  dada_output_R <- dada(derep_seqs_R, err=error_profile_R, multithread=TRUE)
} else {
  dada_output <- dada(derep_seqs, err=error_profile, multithread=TRUE)
}

# Merge F and R if paired, and make sequence table -----
if (paired_or_not == "paired") {
  mergers <- mergePairs(dada_output, derep_seqs, dada_output_R, derep_seqs_R,
                        returnRejects=TRUE, verbose=TRUE,minOverlap = 20, maxMismatch = 2)
  # Inspect the merger data.frame from the first sample
  head(mergers[[1]])

  seqtab <- makeSequenceTable(mergers)

} else {
  seqtab <- makeSequenceTable(derep_seqs)
}

# Remove chimeras ----
seqtab_nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)

# Track the status of the sequences through the steps ----
getN <- function(x) {sum(getUniques(x))}

rownames(filtered_seqs) <- all_sample_names

filtered_seqs <- as.data.frame(filtered_seqs) %>% rownames_to_column("name")

filtered_seqs_nonzeros <- filtered_seqs %>% filter(name %in% filtered_sample_names)
functionally_useless <- filtered_seqs %>% filter(!(name %in% filtered_sample_names)) %>% rename(input = reads.in, filtered = reads.out)
# functionally_useless <- data.frame(input = filtered_seqs[-(which(rownames(filtered_seqs) %in% filtered_sample_names)),1],
#                                    filtered = filtered_seqs[-(which(rownames(filtered_seqs) %in% filtered_sample_names)),2])

if (paired_or_not == "paired") {
  track <- data.frame(name = filtered_seqs_nonzeros[,"name"],
                      input = filtered_seqs_nonzeros[,"reads.in"],
                      filtered = filtered_seqs_nonzeros[,"reads.out"],
                      denoised = sapply(dada_output, getN),
                      merged = sapply(mergers, getN),
                      tabled = rowSums(seqtab),
                      nochim = rowSums(seqtab_nochim))
} else {
  track <- data.frame(name = filtered_seqs_nonzeros[,"name"],
                      input = filtered_seqs_nonzeros[,"reads.in"],
                      filtered = filtered_seqs_nonzeros[,"reads.out"],
                      denoised = sapply(dada_output, getN),
                      tabled = rowSums(seqtab),
                      nochim = rowSums(seqtab_nochim))
}
track <- plyr::rbind.fill(track, functionally_useless)

head(track)

# write.csv(track, ...)

# Make output FASTA files and abundance tables for the processed reads------

if (paired_or_not == "paired") {
  working_barcode_name =  paste("merged_", barC , sep='')
  working_barcode_seqnum = paste("merged_", barC , "_seq_number", sep = '')
  nochim_fname.fasta = paste(mergedoutpath,"/", "nochim_merged",barC,".fasta", sep='')
  nochim_fname.txt = paste(mergedoutpath,"/", "nochim_merged",barC,".txt", sep='')
} else if(paired_or_not == "forward") {
  working_barcode_name =  paste("forward_", barC , sep='')
  working_barcode_seqnum = paste("forward_", barC , "_seq_number", sep = '')
  nochim_fname.fasta = paste(outpath,"/", "nochim_forward",barC,".fasta", sep='')
  nochim_fname.txt = paste(outpath,"/", "nochim_forward",barC,".txt", sep='')
} else {
  working_barcode_name =  paste("reverse_", barC , sep='')
  working_barcode_seqnum = paste("reverse_", barC , "_seq_number", sep = '')
  nochim_fname.fasta = paste(outpath,"/", "nochim_reverse",barC,".fasta", sep='')
  nochim_fname.txt = paste(outpath,"/", "nochim_reverse",barC,".txt", sep='')
}

seqtab_nochim <- t(seqtab_nochim)
nochim_merged <- seqtab_nochim %>% data.frame %>%
  rownames_to_column %>% rename(sequence = rowname) %>% # make sequences into a column
  mutate(!!working_barcode_seqnum := paste0(working_barcode_name,"_",row_number())) %>% # Make a new column w seq number
  select(!!working_barcode_seqnum,everything()) # reorder the columns

# Filter out ASVs with fewer occurrences than the threshold
nochim_merged$total <- rowSums(nochim_merged[,3:ncol(nochim_merged)])
nochim_merged <- nochim_merged %>% filter(total > min_asv_abundance) %>% select(-total)

# Save this table
if(paired_or_not == "paired") {
  dir.create(mergedoutpath, recursive = T)
  write.table(nochim_merged, file = nochim_fname.txt, row.names=FALSE, sep="\t", quote=FALSE)
} else {
  dir.create(outpath, recursive = T)
  write.table(nochim_merged, file = nochim_fname.txt, row.names=FALSE, sep="\t", quote=FALSE)
}

# Make a fasta file out of the sequences in this table
nochim_merged_seq <- nochim_merged %>% select(!!working_barcode_seqnum, sequence)
nochim_merged_seq.fasta = dataframe2fas(nochim_merged_seq, file= nochim_fname.fasta)

# If working with unpaired reads, this is the end! --------------
if(paired_or_not != "paired"){
  cat("Done with analyzing your unpaired reads!\n\n")
  quit()
}

# If working with paired data, there's more to be done ----------
try <- ldply(mergers)  #merge all dataframes resulting from merger


pairedsum_unmerged_table <- try %>% select(id = .id, forward, reverse, abundance, accept) %>% data.frame %>%
  filter(accept == FALSE) %>% select(-accept) # remove the ones that worked well during merge

# grab the correct F and R reads from the dadaF and dadaR files ------


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

pairedsum_unmerged_table$sequenceF <- apply(pairedsum_unmerged_table, 1, function(x)
  mine_unmerged(x, seqs_list = dada_output, forward_or_reverse = "forward"))

pairedsum_unmerged_table$sequenceR <- apply(pairedsum_unmerged_table, 1, function(x)
  mine_unmerged(x, seqs_list = dada_output_R, forward_or_reverse = "reverse"))



# add seqeunce length for forwards and reversed to unmerged dereplicated data ------



pairedsum_unmerged_table$lengthF <- nchar(gsub("[a-z]","",pairedsum_unmerged_table$sequenceF))
pairedsum_unmerged_table$lengthR <- nchar(gsub("[a-z]","",pairedsum_unmerged_table$sequenceR))
pairedsum_unmerged_table$totalseq <- pairedsum_unmerged_table$lengthF + pairedsum_unmerged_table$lengthR

# add expected length of amplicon
pairedsum_unmerged_table$keep[pairedsum_unmerged_table$totalseq>=barC_length] <- FALSE
pairedsum_unmerged_table$keep[pairedsum_unmerged_table$totalseq<barC_length] <- TRUE

# Get the reverse complement of the R sequence
pairedsum_unmerged_table$sequenceRc <- sapply(sapply(sapply(pairedsum_unmerged_table$sequenceR, DNAString), reverseComplement), toString)

# pairedsum_unmerged_table$NNNNNN <- "AAAAAAAAAATTCTTAAAAAAAAAA"
pairedsum_unmerged_table$sequenceF_N_Rrc <- paste(pairedsum_unmerged_table$sequenceF,"AAAAAAAAAATTCTTAAAAAAAAAA",pairedsum_unmerged_table$sequenceRc,sep="")


pairedsum_unmerged_dada2 <- pairedsum_unmerged_table %>% filter(keep == TRUE) %>% select(sequenceF_N_Rrc, id, abundance)



if (nrow(pairedsum_unmerged_dada2) == 0) {
  cat("Done with analyzing your paired reads! None of the paired-but-unmerged reads were kept.\n\n")
  quit()
}

# Spread the dataframe, and sum up the abundances per id
unmerged_seqtab <- dcast(pairedsum_unmerged_dada2, sequenceF_N_Rrc ~ id, fun.aggregate = sum) %>%
  data.frame %>% column_to_rownames( "sequenceF_N_Rrc") %>% t()



# Run bimera detection on unmerged reads and discard bimeras -----

unmerged_seqtab_nochim <- removeBimeraDenovo(unmerged_seqtab, method="consensus", multithread=TRUE, verbose=TRUE)



# Reformat data for Bowtie2 ---------


# get barcode name info for numbering / nameing sequences in final fasta files / tables
unmergedbarC =  paste("unmerged_", barC , sep='')
unmergedbarCseqnum = paste("unmerged_", barC , "_seq_number", sep = '')


#transpose table, make rownames into a column -> sequnces
unmerged_seqtab_nochim <- unmerged_seqtab_nochim %>% t %>% data.frame %>% rownames_to_column("sequences")


## split sequence column into two -> sequences F and Rrc (reverse complement), split on dummy sequnces used to merge reads initially
unmerged_seqtab_nochim$sequencesF <- sapply(strsplit(as.character(unmerged_seqtab_nochim$sequences),'AAAAAAAAAATTCTTAAAAAAAAAA'), "[", 1)
unmerged_seqtab_nochim$sequencesRrc <- sapply(strsplit(as.character(unmerged_seqtab_nochim$sequences),'AAAAAAAAAATTCTTAAAAAAAAAA'), "[", 2)
unmerged_seqtab_nochim$sequences <- NULL

# reverse complement the reverse read for Bowtie2...
unmerged_seqtab_nochim$sequencesR <- sapply((sapply(sapply(unmerged_seqtab_nochim$sequencesRrc, DNAString), reverseComplement)), toString)
unmerged_seqtab_nochim$sequencesRrc <- NULL

# Reorder the columns (seqF, seqR, all the samples)

unmerged_seqtab_nochim <- unmerged_seqtab_nochim %>% select(sequencesF, sequencesR, everything()) %>%
  mutate(!!unmergedbarCseqnum := paste0(unmergedbarC, "_", row_number())) %>% select(!!unmergedbarCseqnum, everything())

# filter out ASVs with fewer occurrences than the total
unmerged_seqtab_nochim$total <- rowSums(unmerged_seqtab_nochim[,4:ncol(unmerged_seqtab_nochim)])
unmerged_seqtab_nochim <- unmerged_seqtab_nochim %>% filter(total > min_asv_abundance) %>% select(-total)


# ########################################
# ### make output files
# #########################################
#
# make file paths

nochim_unfnameF.fasta = paste(unmergedoutpath,"/", "nochim_unmerged",barC,"F",".fasta", sep='')
nochim_unfnameR.fasta = paste(unmergedoutpath,"/", "nochim_unmerged",barC,"R",".fasta", sep='')
nochim_unfname.txt = paste(unmergedoutpath,"/", "nochim_unmerged",barC,".txt", sep='')

# make data frames for the soon to be made fasta files  with reads and read names
nochim_unmerged_seq_F <- data.frame(unmerged_seqtab_nochim[[unmergedbarCseqnum]],unmerged_seqtab_nochim$sequencesF)
nochim_unmerged_seq_R <- data.frame(unmerged_seqtab_nochim[[unmergedbarCseqnum]],unmerged_seqtab_nochim$sequencesR)

# export the fasta file ready dataframes
nochim_unmerged_seq_F.fasta = dataframe2fas(nochim_unmerged_seq_F, file= nochim_unfnameF.fasta)
nochim_unmerged_seq_R.fasta = dataframe2fas(nochim_unmerged_seq_R, file= nochim_unfnameR.fasta)

# write summary table
write.table(unmerged_seqtab_nochim, file = nochim_unfname.txt, row.names=FALSE, sep="\t", quote=FALSE)
cat("Done with analzing your paired reads!\n\n")
quit()
