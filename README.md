# Anacapa Toolkit

### Anacapa last updated 1-10-2018

#### Written by Emily Curd (eecurd@g.ucla.edu), Jesse Gomer (jessegomer@gmail.com), Zack Gold (zjgold@ucla.edu), Gaurav Kandlikar (gkandlikar@ucla.edu), Baochen Shi (biosbc@gmail.com), Rachel Turba (rturba@ucla.edu) and Rachel Meyer (rsmeyer@ucla.edu).
#### Developed at UCLA for the University of California Conservation Consortium's CALeDNA Program

## Introduction
Anacapa Island's name is derived from the Chumash __Ennepah__ or __Anyapakh__ which means "mirage island" (Gudde and Bright, 2004). Much like Anacapa Island, using environmental DNA (eDNA) to uncover Biodiversity seems like an illusion on the horizon. Anacapa toolkit processes eDNA reads and assigns taxonomy using existing software or modifications to existing software.

Anacapa is an automated metabarcoding read processing toolkit. It is a modular toolkit that is is designed to analyze multiple samples and metabarcodes simultaneously. Anacapa accomplishes this in four steps: 1) building reference libraries using CRUX: Creating Reference libraries Using eXisting tools, 2) running QC and assigning Amplicon Sequence Variants (ASV) using Dada2, 3) assigning taxonomy using bowtie2 and a bowtie2 specific Bayesian Least Common Ancestor (BLCA) and 4) Generating ecological diversity summary statistics for a set of samples.  This toolkit does not require that paired reads overlap, or that both reads in a pair pass QC.  Taxonomy results are generated for all read types and the user can decide which read types they wish to retain.

#### Step 1: CRUX: Creating Reference libraries Using eXisting tools
For details on building reference libraries using CRUX, please refer to the following: https://github.com/limey-bean/CRUX_Creating-Reference-libraries-Using-eXisting-tools. The output of CRUX consists of two reference libraries. Each Library contains unique metabarcode specific reads that correspond NCBI accession version numbers. Libraries consist of fasta and taxonomy files and a bowtie2 index library.  The libraries are unfiltered or filtered.  Unfiltered libraries contain every dereplicated read found during the BLAST searches. The filtered library contains only reads with robust taxonomic assignments.  In CRUX robust means that any read with the following in their taxonomic path: 'uncultured', 'environmental', 'sample', or 'NA;NA;NA;NA' is excluded from the library. Prebuilt CRUX reference libraries (CO1-Leray, 12S-MiFish, 16S-V4, 18S-EMP, Fungal ITS (FITS), Plant ITS2 (PITS), and PPM cytochrome Oxidase [add refs and primer sequence] can be found at [link to dryad].

#### Step 2: Running QC and assigning Amplicon Sequence Variants (ASV) using dada2
This step of the toolkit simultaneously processes raw fastq reads for single or multiple samples with single or multiple metabarcode targets generated on Illumina HiSeq and MiSeq machines. It is not required that all samples contain reads for each metabarcode. This script trims nextera and truseq adapters (**cutadapt**; Martin 2011), removes low quality reads, and sorts reads by metabarcode primer sequence (**cutadapt**). Clean reads are processed with dada2 to generate ASV for each metabarcode. Dada2 denoises, dereplicates, merges (where possible), and removes chimeric sequences from the data set. The input is raw Illumina metabarcode sequence data [\*.fastq.gz] reads and outputs are site frequency tables of ASV and species count data for multiple samples and metabarcodes. Successful implementation of this step requires: 1) raw illumina sequencing data and 2) a set of fasta formatted forward and reverse fasta format files that include the metabarcoding primers used to generate sequence data.

#### Step 3: Assigning taxonomy using bowtie2 and a bowtie2 specific Bayesian Least Common Ancestor (BLCA) algorithm
The Anacapa toolkit determines the best taxonomic hits for an ASV using bowtie2.  Because Anacapa considers paired merged, paired unmerged, and unpaired sequencing reads, a fast and flexible read aligner is required. Bowtie2 handles all possible kinds of data resulting from sequencing runs, and it uses CRUX databases as reference libraries. All reads are globally aligned against the CRUX database. Any reads that fail to align are then aligned locally. The best hits (the top 100 bowtie 2 returns) are then processed with BLCA to assign taxonomy.  The bowtie2 BLCA algorithm was adapted from (add blast BLCA). Add explanation of BLCA.... [Jesse is this accurate?  more text?]


#### Step 4: Generating ecological diversity summary statistics for a set of samples

Zack and Gaurav do you wanna take a stab at this?

## add this stuff to the above text
The workflow: Anacapa takes raw Illumina fastq format reads and preprocesses them to assess file corruption (**md5sum**) and uncompresses (**gunzip**) and then renames the files for readability  readable.  Reads are next trimmed using **cutadapt** (Martin 2011) to remove sequencing adapters from the 5' ends and sequnging adapters and primers from the 3' end of reads.  **Fastx-toolkit** (Gordon and Hannon, 2010) is then used to processed for quality control. Read are retained if they have a Q ≥ 35 and are at lease 100bp after adapter and 3' primer trimming. **Cutadapt** is next used to sort reads by primer, and to trim additional basepairs from the end of read to increase quality going into **dada2**. Prior to running **dada2** there is a paired end read checking step.  Reads that have pairs that passed qc are run separately from unpaired F or unpaired R reads.  **dada2** is then used to denoise, dereplicate, merge (where possible), and remove chimeric sequences from the data set.   **dada2** processed reads are then assigned taxonomy using **Bowtie2** (Langmead and Salzberg, 2012).
*	elaborate on dada2 and the bowtie2 process and **R** (Team, RC, 2000)


## Anacapa relies on many programs and databases to run properly.
**__First Download the Anacapa_db folder.__**
* This folder contains:
	* Two folders
		* adapters_and_PrimAdapt_rc
		* scripts
	* Two files:
		* forward_primers.txt
		* reverse_primers.txt
		* __The two files are examples of how to format the primer forward and reverse input files.  It is VERY IMPORTANT that you modify these files or make new files to reflect your data set!__

* adapters_and_PrimAdapt_rc folder contains:
	* the forward and reverse nextera adapters
	* the forward and reverse trueseq adapters.


* scripts folder contains:
	* anacapa_bowtie2_blca.sh
	* anacapa_config.sh
	* anacapa_format_primers_cutadapt.py
	* anacapa_QC_dada2.sh
	* anacapa_vars_nextV.sh
	* append_blca_to_summary.py
	* append_bowtie_to_summary.py
	* blca_from_bowtie.py
	* check_paired.py
	* dada2_unified_script.R
	* group_alignments_to_db.py
	* group_alignments_to_files_p_mod.py
	* group_alignments_to_files.py
	* merge_asv.py
	* reformat_summary_for_r.py
	* run_bowtie2_blca.sh
	* run_dada2.sh
	* sqlite_storage.py
	* summarize_bowtie2_hits_full_taxonomy.py
	* summarize_bowtie2_hits.py
	* and a directory for downstream-analyses



**__Programs__** [Jesse can I get you to verify these?]
To run Anacap, you need verify that the full path to each of the following programs is correctly indicated in the anacapa_config.sh file.  

1. cutadapt: http://cutadapt.readthedocs.io/en/stable/index.html

2. fastxtoolkit (add versions)

3. anaconda/python2-4.2
	* make sure biopython is installed

4. R 3.4.2

5. Bowtie2: http://bowtie-bio.sourceforge.net/bowtie2/index.shtml
	* Bowtie2 does not need to be installed in the crux_release_V1_db folder, however you will need to verify that the Crux_config.sh is modified for you computing environment.

6. muscle: https://www.drive5.com/muscle/downloads.htm


**__CRUX Databases to download__**
Download taxonomy reference libraries from this google drive folder: https://drive.google.com/drive/folders/0BycoA83WF7aNOEFFV2Z6bC1GM1E?usp=sharing

Users can also make their own libraries using CRUX.  Sliva and greengees libraries can easily be converted to Anacapa compatible libraries.  [make documentation for this.]

**_Hoffman users running the QC dada2 need to do the following before dada2 will run_**
```
qrsh
module load module load R/3.4.2
module load gcc/6.3.0
R # will open R in terminal to go back to bach quit() I think...
source("https://bioconductor.org/biocLite.R")
biocLite(suppressUpdates = FALSE)
```
if given this option "would you like to use a personal library" say "y"
This is the reason that it is not possible to install this in the R script, because it requires user input.
```
biocLite("ShortRead", suppressUpdates = FALSE)
biocLite("devtools")
```
this bit could take a very long time so no worries....


## How to run the QC / dada2 step:
```
sh ~/Anacapa_db/scripts/anacapa_QC_dada2.sh -i <input_dir> -o <out_dir> -d <database_directory> -u <hoffman_account_user_name> -f <fasta file of forward primers> -r <fasta file of reverse primers> -a <adapter type (nextera or truseq)>  -t <illumina run type HiSeq or MiSeq> -l (add flag (-l), no text required, if running locally)
```

## How to run the bowtie2 blca step:

```
# sh ~/Anacapa_db/scripts/anacapa_bowtie2_blca.sh -o <out_dir_for_anacapa_QC_run> -d <database_directory> -u <hoffman_account_user_name> -l (add flag (-l), no text required, if running locally)
```



This is a preliminary attempt a documentation...  To run the latest version of Anacapa using the 6 primer sets that are commonly used for CALeDNA projects, you will need to download the following:
* Download the Anacapa_db folder
* Download taxonomy reference libraries (updated 10-16-2017) from this google drive folder: https://drive.google.com/drive/folders/0BycoA83WF7aNOEFFV2Z6bC1GM1E?usp=sharing

Once you download these folders, transfer the reference library folders to the Anacapa_db folder.  You should have the following subdirectories in the Anacapa_db folder:
* 12S  
* 16S  
* 18S  
* adapters_and_PrimAdapt_rc  
* CO1  
* FITS  
* PITS  
* scripts

The script to run Anacapa is in the scripts directory.  It is called: anacapa_release_V1.sh

To run the script you need to run the following command:

sh ~/Anacapa_db/scripts/anacapa_release_V1.sh -i <input_dir> -o <out_dir> -d <database_directory> -u <hoffman_account_user_name> -f <fasta file of forward primers> -r <fasta file of reverse primers> -a <adapter type ("nextera" or "truseq")>  -t <illumina run type HiSeq or MiSeq>

### More to come, and it might be a bit buggy.
* if you choose to take this on...  Good Luck!



## References
Caporaso, J.G., Kuczynski, J., Stombaugh, J., Bittinger, K., Bushman, F.D., Costello, E.K., Fierer, N., Peña, A.G., Goodrich, J.K., Gordon, J.I. and Huttley, G.A., 2010. QIIME allows analysis of high-throughput community sequencing data. Nature methods, 7(5), pp.335-336.

Edgar, R.C., 2010. Search and clustering orders of magnitude faster than BLAST. Bioinformatics, 26(19), pp.2460-2461.

Gordon, A. and Hannon, G.J., 2010. Fastx-toolkit. FASTQ/A short-reads preprocessing tools (unpublished) http://hannonlab. cshl. edu/fastx_toolkit, 5.

Gudde, Erwin; William Bright (2004). California Place Names (Fourth ed.). University of California Press. p. 12. ISBN 0-520-24217-3.

Langmead, B. and Salzberg, S.L., 2012. Fast gapped-read alignment with Bowtie 2. Nature methods, 9(4), pp.357-359.

Martin, M., 2011. Cutadapt removes adapter sequences from high-throughput sequencing reads. EMBnet. journal, 17(1), pp.pp-10.

Team, R.C., 2000. R language definition. Vienna, Austria: R foundation for statistical computing.

Zhang, J., Kobert, K., Flouri, T. and Stamatakis, A., 2013. PEAR: a fast and accurate Illumina Paired-End reAd mergeR. Bioinformatics, 30(5), pp.614-620.
