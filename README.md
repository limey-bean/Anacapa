# Anacapa Toolkit

### Anacapa last updated 1-09-2017, Big Bug in run_dada2_bowtie2.  Fixing today

#### Written by Emily Curd (eecurd@g.ucla.edu), Jesse Gomer (jessegomer@gmail.com), Zack Gold (zjgold@ucla.edu), Gaurav Kandlikar (gkandlikar@ucla.edu), Baochen Shi (biosbc@gmail.com), Rachel Turba (rturba@ucla.edu ) and Rachel Meyer (rsmeyer@ucla.edu).
#### Developed at UCLA for the University of California Conservation Consortium's CALeDNA Program

## Introduction
Anacapa Island's name is derived from the Chumash __Ennepah__ or __Anyapakh__ which means "mirage island" (Gudde and Bright, 2004). Much like Anacapa Island, using environmental DNA (eDNA) to uncover Biodiversity seems like an illusion on the horizon. Anacapa toolkit processes eDNA reads and assigns taxonomy using existing software or modifications to existing software.

Anacapa is an automated metabarcoding read processing toolkit. It is designed to analyze multiple samples and metabarcodes simultaneously. Anacapa accomplishes this in four steps: 1) building reference libraries using CRUX: Creating Reference libraries Using eXisting tools, 2) running QC and assigning Amplicon Sequence Variants (ASV) using Dada2, 3) assigning taxonomy using bowtie2 and a bowtie2 specific Bayesian Least Common Ancestor (BLCA) and 4) Generating ecological diversity summary statistics for a set of samples.  

For details on building reference libraries using CRUX, please refer to the following: https://github.com/limey-bean/CRUX_Creating-Reference-libraries-Using-eXisting-tools

It processes raw fastq reads generated on Illumina HiSeq and MiSeq machines. This pipeline does not require that the raw paired reads overlap, or that both reads in a pair pass qc.  Taxonomy results are generated for all read types, however taxonomy is reported for overlapping (assembled) and non overlapping (unassembled) paired end reads, and also single end reads (discarded forward or discarded reverse) where one of the pairs fails qc. The input is raw Illumina metabarcode sequence data and outputs are species count data for multiple samples and metabarcodes. Successful implementation of Anacap requires: 1) raw illumina sequencing data, 2) a set of fasta formatted forward and reverse fasta format files that include the metabarcoding primers used to generate sequence data, 3) reference libraries that correspond to the metabarcodes of interest made using CRUX, and 4) the dependencies indicated below.

The workflow: Anacapa takes raw Illumina fastq format reads and preprocesses them to assess file corruption (**md5sum**) and uncompresses (**gunzip**) and then renames the files for readability  readable.  Reads are next trimmed using **cutadapt** (Martin 2011) to remove sequencing adapters from the 5' ends and sequnging adapters and primers from the 3' end of reads.  **Fastx-toolkit** (Gordon and Hannon, 2010) is then used to processed for quality control. Read are retained if they have a Q ≥ 35 and are at lease 100bp after adapter and 3' primer trimming. **Cutadapt** is next used to sort reads by primer, and to trim additional basepairs from the end of read to increase quality going into **dada2**. Prior to running **dada2** there is a paired end read checking step.  Reads that have pairs that passed qc are run separately from unpaired F or unpaired R reads.  **dada2** is then used to denoise, dereplicate, merge (where possible), and remove chimeric sequences from the data set.   **dada2** processed reads are then assigned taxonomy using **Bowtie2** (Langmead and Salzberg, 2012).


*	elaborate on dada2 and the bowtie2 process and **R** (Team, RC, 2000)

(Zack wants to add a paragraph about ESV's VS OTUs. We sort of ignore both really, and go for defined groups like species, genus, family, etc... Is that a shortcomming of the method?)  


## Anacapa relies on many programs and databases to run properly.
**__First Download the Anacapa_db folder.__**
* This folder contains:
	* Two folders
		* adapters_and_PrimAdapt_rc
		* scripts
	* Two files:
		* forward_primers.txt
		* reverse_primers.txt

* adapters_and_PrimAdapt_rc folder contains:
	* the forward and reverse nextera adapters
	* add forward and reverse trueseq (add) adapters.

* scripts folder contains:
	* anacapa_config.sh
	* anacapa_format_primers_cutadapt.py
	* anacapa_release_V1_dada2_plus_bowtie2.sh
	* anacapa_vars.sh
	* check_paired.pl
	* dada2_paired.R
	* dada2_unpaired_F.R
	* dada2_unpaired_R.R
	* run_dada2_bowtie2_paired.sh
	* run_dada2_bowtie2_unpaired_F.sh
	* run_dada2_bowtie2_unpaired_R.sh
	* pending -> bowtie2summary script

* The two files are examples of how to format the primer forward and reverse input files.  It is VERY IMPORTANT that you modify these files to reflect your data set!

**__Programs__**
To run Anacap, you need verify that the full path to each of the following programs is correctly indicated in the anacapa_config.sh file.  

1. cutadapt: http://cutadapt.readthedocs.io/en/stable/index.html

2. fastxtoolkit (add versions)

3. Python (add versions)

4. Perl (add versions)

5. R (add versions)

3. Bowtie2: http://bowtie-bio.sourceforge.net/bowtie2/index.shtml
	* Bowtie2 does not need to be installed in the crux_release_V1_db folder, however you will need to verify that the Crux_config.sh is modified for you computing environment.


**__Databases to download__**

## Running this sucker...

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
