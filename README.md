# Anacapa

### Anacapa_release_V1		10-17-2017
#### Written by Emily Curd (eecurd@g.ucla.edu), Jesse Gomer (jessegomer@gmail.com), Baochen Shi (biosbc@gmail.com), Zack Gold (zjgold@ucla.edu), Gaurav Kandlikar (gkandlikar@ucla.edu), Rachel Turba (rturba@ucla.edu ) and Rachel Meyer (rsmeyer@ucla.edu). 
#### Developed at UCLA for the University of California Conservation Consortium's CALeDNA Program

## Introduction
Anacapa Island's name is derived from the Chumash __Ennepah__ or __Anyapakh__ which means "mirage island" (Gudde and Bright, 2004). Much like Anacapa Island, using environmental DNA (eDNA) to uncover Biodiversity seems like an issusion on the horizon. This pipeline... (something amazing)  

Anacapa is an automated metabarcoding read processing pipeline.  It is designed to analyze multiple samples and metabarcodes simultaneously. It processed fastq reads generated on Illumina HiSeq and MiSeq machines. This pipleline does not require that the raw paired reads overlap, or that both reads in a pair pass qc.  Taxonomy results are generated for all read types, however taxomony is reported for overlapping (assembled) and non overlapping (unassembled) paired end reads, and also single end reads (discarded forward or discarded reverse) where one of the pairs fails qc. The input is raw Illumina metabarcode sequence data and outputs are species count data for multiple samples and metabarcodes. Sucessful implementaion of Anacap requires: 1) raw illumina data, 2) a set of fasta formatted forward and reverse fasta format files that indlude the metabarcodes used to generate seqeunce data, 3) reference libraries that correspond to the metabrcodes of interest made using CRUX, and 4) the dependencies indicated below. 

The workflow: Anacapa takes raw Illumina fastq format reads and preprocesses them to assess file corruption (md5sum) and uncompresses (gunzip) and renames files and reads within files.  Reads are next processed for quality control using cutadapt, where read are retained if they have a Q ≥ 30 and are at lease 100bp after adapter and 3' primer trimming.  PEAR (Zhang et al., 2013) is used to merge paired reads.  Metabarcode reads are sorted by 5' primers and fastq files are converted to fasta files with cutadapt (Martin, 2011). Some reads (unassembled reverse reads) are reverse complemented using Fastx-toolkit (Gordon and Hannon, 2010), prior to sorting. Sorted reads are assigned taxonomy using Bowtie2, and Qiime is used to merge and summarize results tables. 

Bowtie2 taxonomic assignment: Reference libraries for many metabarcode primers suffer from several problems: low coverage for many taxa, incomplete seqeunce coverage of amplicon region, species have low variability and are not distinct from other species with in a genus or family, etc. Anacapa processes reads iteratively and takes into account 1) the read overhang between the sample read and the reference sequences, 2) the percent identity of the sample read to a reference, and 3) the number of equal best hits between a sample read and a given reference library. Anacapa sorts reads iteratively using bowtie2. Reads are initially sorted by read overhang (the default is bins with reads containing <= 25 bp, 25 < reads >= 50, 50 < reads >= 75, and 75 < reads >= 100). Within each overhang bin, Anacapa uses bowtie2 to identify reads with single hits (99% identity or better) to a reference library clustered at 99%, the remaining reads are then sorted using a reference library clustered at 97% and the single best hits (97% identity or better) are retained, the remaining reads are sorted with 95% reference libraries, and so on for 90, 85, and 80% reference libraries. Sample reads that have a forward and reverse read are sorted into bins based on the read (in the pair) with the largest overhand or lowest percent single hit to a reference database. Reads are summarized for each bin / percent reference library combination, and then summarized after merging all percent reference library results within bins, and then summarized for all data from every bin / percent reference library combination.  

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
	* anacapa_release_V1.sh
	* anacapa_vars.sh
	* check_paired.pl
	* group_alignments_to_files.py *- not currently used (Emily github knowledge problem)*
	* group_alignments_to_files_p_mod.py
	* pick_open_otus_and_summ.sh *- not currently used (for qiime processing)*
	* run_bowtie2_make_3_Sfolders.sh
	* summarize_bowtie2_hits.py *- not currently used (Emily github knowledge problem)*
	* summarize_bowtie2_hits_full_taxonomy.py

* The two files are examples of how to format the primer forward and reverse input files.  It is VERY IMPORTANT that you modify these files to reflect your data set!

**__Programs__**
To run Anacap, you need verify that the full path to each of the following programs is correctly indicated in the anacapa_config.sh file.  

1. cutadapt: http://cutadapt.readthedocs.io/en/stable/index.html

2. PEAR: 

3. fastxtoolkit

4. Perl

3. Qiime 1: http://qiime.org/index.html
	* Installation information can be found here: http://qiime.org/install/install.html
	* We will transition to Qiime 2 by January 01, 2018. 
	
4. Bowtie2: http://bowtie-bio.sourceforge.net/bowtie2/index.shtml
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

sh ~/Anacapa_db/scripts/anacapa_release_V1.sh -i <input_dir> -o <out_dir> -d <database_directory> -u <hoffman_account_user_name> -f <fasta file of forward primers> -r <fasta file of reverse primers> -a <adapter type ("nextera" or "truseq")>
 
### More to come, and it might be a bit buggy.
* if you choose to take this on...  Good Luck!



## References
Caporaso, J.G., Kuczynski, J., Stombaugh, J., Bittinger, K., Bushman, F.D., Costello, E.K., Fierer, N., Peña, A.G., Goodrich, J.K., Gordon, J.I. and Huttley, G.A., 2010. QIIME allows analysis of high-throughput community sequencing data. Nature methods, 7(5), pp.335-336.

Edgar, R.C., 2010. Search and clustering orders of magnitude faster than BLAST. Bioinformatics, 26(19), pp.2460-2461.

Gordon, A. and Hannon, G.J., 2010. Fastx-toolkit. FASTQ/A short-reads preprocessing tools (unpublished) http://hannonlab. cshl. edu/fastx_toolkit, 5.

Gudde, Erwin; William Bright (2004). California Place Names (Fourth ed.). University of California Press. p. 12. ISBN 0-520-24217-3.

Martin, M., 2011. Cutadapt removes adapter sequences from high-throughput sequencing reads. EMBnet. journal, 17(1), pp.pp-10.

Zhang, J., Kobert, K., Flouri, T. and Stamatakis, A., 2013. PEAR: a fast and accurate Illumina Paired-End reAd mergeR. Bioinformatics, 30(5), pp.614-620.
