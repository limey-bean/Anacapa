# Anacapa Toolkit

### Anacapa last updated 1-10-2018

#### Written by Emily Curd (eecurd@g.ucla.edu), Jesse Gomer (jessegomer@gmail.com), Zack Gold (zjgold@ucla.edu), Gaurav Kandlikar (gkandlikar@ucla.edu), Baochen Shi (biosbc@gmail.com), Rachel Turba (rturba@ucla.edu) and Rachel Meyer (rsmeyer@ucla.edu).
#### Developed at UCLA for the University of California Conservation Consortium's CALeDNA Program

## Introduction
Anacapa is an iconic island off the Southern Californian coast that has significant cultural and biodiversity importance. Anacapa Island's name is actually derived from the Chumash __Ennepah__ or __Anyapakh__ which means "mirage island" (Gudde and Bright, 2004). Much like the name Anacapa, using environmental DNA (eDNA) to uncover biodiversity seems like an illusion on the horizon. Anacapa toolkit processes eDNA reads and assigns taxonomy using existing software or modifications to existing software.

Anacapa is an automated metabarcoding read processing toolkit. This modular toolkit is designed to analyze multiple samples and metabarcodes simultaneously from Ilumina sequencing platforms. Anacapa accomplishes this in four steps: 1) building reference libraries using CRUX: Creating Reference libraries Using eXisting tools, 2) running quality control (QC) and assigning Amplicon Sequence Variants (ASV) using Dada2, 3) assigning taxonomy using bowtie2 and a bowtie2 specific Bayesian Least Common Ancestor (BLCA) and 4) Running exploratory data analysis and generating ecological diversity summary statistics for a set of samples. A significant advantage of the Anacapa toolkit is that it does not require that paired reads overlap, or that both reads in a pair pass QC.  Taxonomy results are generated for all read types and the user can decide which read types they wish to retain for downstream analysis.

#### Step 1: CRUX: Creating Reference libraries Using eXisting tools
For full details on building reference libraries using CRUX, please refer to the following: https://github.com/limey-bean/CRUX_Creating-Reference-libraries-Using-eXisting-tools.

This first part of the toolkit generates reference libraries needed for taxonomic assignment using CRUX.  The output of CRUX consists of two reference libraries, either unfiltered or filtered. Unfiltered libraries contain every dereplicated read found during the BLAST searches. The filtered library contains only reads with robust taxonomic assignments. Specifically we refer to robust taxonomic assignments as any reads that do not have the following in their taxonomic path: 'uncultured', 'environmental', 'sample', or 'NA;NA;NA;NA'. Prebuilt CRUX reference libraries (12S-MiFish, CO1, 16S-EMP, 18S V4, 18S V8-9, 18S-EMP, Fungal ITS (FITS), Plant ITS2 (PITS), and PPM cytochrome Oxidase [see Table 1] can be found at [link to dryad]. Each library contains unique metabarcode specific reads that correspond to NCBI accession version numbers. Libraries consist of fasta,  taxonomy files, and a bowtie2 index library.

<p align="center">
<img src="/figures-and-tables-for-the-Github/Table_1.png">
</p>

We acknowledge that users may wish to substitute their own reference libraries or add additional samples to a CRUX generated reference library to improve taxonomic assignments. Please refer to the CRUX page above for instructions to create a CRUX formatted reference library for use in this toolkit.

#### Step 2: Running QC and assigning Amplicon Sequence Variants (ASV) using dada2
This next step of the toolkit aims to conduct standard sequence QC and then generate amplicon sequence variants (ASV) from Ilumina data using dada2. ASVs are a novel solution to identifying biologically informative unique sequences in metabarcoding samples that replaces the operational taxonomic unit (OTU) framework. Unlike OTUs which cluster sequences using an arbitrary sequence similarity (ex 97%), ASVs are unique sequence reads determined using Bayesian probabilities of known sequencing error. These unique sequences can be as little as 2 bp different, providing improved taxonomic resolution and an increase in observed diversity. Please see [insert citation to Dada2 and deblur] for further discussion.

An strong advantage of the Ancapa toolkit is that is can simultaneously processes raw fastq reads for samples with single or multiple metabarcode targets generated on Illumina HiSeq and MiSeq machines. It is also not required that all samples contain reads for each metabarcode, thus allowing users to combine multiple projects or targets on the same sequencing run while only running the pipeline once.

Anacapa takes raw Illumina fastq format reads and preprocesses them to assess file corruption (**md5sum**) and uncompresses (**gunzip**) and then renames the files for readability  readable. The QC portion of this script trims nextera and truseq adapters (**cutadapt**; Martin 2011), removes low quality reads **Fastx-toolkit**, and sorts reads by metabarcode primer sequence (**cutadapt**). Reads are trimmed using **cutadapt** (Martin 2011) to remove sequencing adapters from the 5' ends and sequencing adapters and primers from the 3' end of reads.  **Fastx-toolkit** (Gordon and Hannon, 2010) is then used to processed for quality control. Read are retained if they have a Q ≥ 35 and are at lease 100bp after adapter and 3' primer trimming. **cutadapt** is next used to sort reads by primer, and to trim additional basepairs from the end of read to increase quality going into **dada2**. Prior to running **dada2** a custom python script sorts reads into unpaired F, unpaired R and unmerged read files.  The files are passed separately into **dada2*** where they are denoised, dereplicated, merged (where possible), and  chimeric sequences removed from the data set.  


The input is raw Illumina metabarcode sequence data [\*.fastq.gz] reads and outputs are site frequency tables of ASVs and species count data for multiple samples and metabarcodes (ASV table). Successful implementation of this step requires: 1) raw illumina sequencing data and 2) a set of fasta formatted forward and reverse fasta format files that include the metabarcoding primers used to generate sequence data.

#### Step 3: Assigning taxonomy using bowtie2 and a bowtie2 specific Bayesian Least Common Ancestor (BLCA) algorithm
This next module of the pipeline assigns taxonomy to ASVs.

The Anacapa toolkit determines the best taxonomic hits for an ASV using bowtie2 (Langmead and Salzberg, 2012). Anacapa considers paired merged, paired unmerged, and unpaired sequencing reads, and thus a fast and flexible read aligner, such as bowtie2, is required to handle all four read types. This script uses bowtie2, CRUX reference libraries, and BLCA to then assign taxonomy. All reads are first globally aligned against the CRUX database using bowtie2. Any reads that fail to align are then aligned locally. The best hits (the top 100 bowtie2 returns) are then processed with BLCA script to assign taxonomy. The bowtie2 BLCA algorithm was adapted from https://github.com/qunfengdong/BLCA. BLCA uses pairwise sequence alignment to calculate sequence similarity between query sequences and reference library hits. Taxonomy is assigned based on the lowest common ancestor of multiple reference library hits for each query sequence. The reliability of each taxonomic assignment is then evaluated through bootstrap confidence scores [BLCA citation].

Successful implementation of this script requires an ASV table (site frequency table) with the following columns (ASV number, Sequence, Samples). The output is multiple ASV tables for each metabarcode and both global and local alignments.


#### Step 4: Generating ecological diversity summary statistics for a set of samples

The last step of the Anacapa Pipeline conducts exploratory data analysis to provide a first pass look at sequencing depth, taxonomic assignments, and generated data tables. This analysis is not meant for publication, but solely as a first stab at visualization of your data. This is helpful in identifying potential glaring errors or contamination, and identifying patterns worth investigating further through more robust analysis. We highly encourage data exploration before further analysis as different parameters within the Anacapa pipeline may produce differences in downstream results and these parameters will vary by project, stringency of taxonomic assignment, and users opinions. The exploratory_analysis.R script uses a variety of R packages, relying heavily on __phyloseq__, __vegan__, and __ggplot2__. See below for full list of R package dependencies and scripts. The output from reformat_summary_for_R.py is an ASV site frequency table with assigned taxonomy and is the input used for the exploratory_analysis.R script. In addition, the user supplies an input metadata table that only requires the first column be sample names. Users can include any type of metadata including categorical, continuous, and discete variables. The first step of the R script is to convert the input files into a Phyloseq class object. We then generate bar plots looking at total number of observed classes and relative abundance of each class. We then generate rarefaction curves, alpha diversity boxplots to observe total number of taxa and Shannon diversity, and alpha diversity statistics. In addition, we calculate jaccard and bray-curtis distance matrices and conduct NMDS ordination plots, network map, heat maps, and ward-linkage maps. Each of the above analyses are repeated with different grouping for each metadata column. In addition we conduct two betadiversity statistical tests, pairwise adonis and betadisp from the vegan package. Again each analyses is repeated across groupings of each metadata column.


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
	* version muscle3.8.31
	* **__muscle must be installed within the anacapa_db folder__**


**__CRUX Databases to download__**
Download taxonomy reference libraries from this google drive folder: https://drive.google.com/drive/folders/0BycoA83WF7aNOEFFV2Z6bC1GM1E?usp=sharing

Users can also make their own libraries using CRUX.  Sliva and greengeens libraries can easily be converted to Anacapa compatible libraries.  [make documentation for this.]

**_Hoffman users running the QC dada2 need to do the following before dada2 will run_**
```
qrsh
module load R/3.4.2
module load gcc/6.3.0
R
source("https://bioconductor.org/biocLite.R")
biocLite(suppressUpdates = FALSE)
```
if given this option "would you like to use a personal library" say "y"
This is the reason that it is not possible to install this in the R script, because it requires user input.
```
biocLite("ShortRead", suppressUpdates = FALSE)
biocLite("devtools")
```
you may also need to install MASS, mgcv, and rpart
```
install.packages("MASS")
install.packages("mgcv")
install.packages("rpart")
```


this bit could take a very long time so no worries....

also make sure that biopython is installed.

```
module load anconda
pip install biopython --user
```
"user" not your user name

## How to run the QC / dada2 step:
```
sh ~/Anacapa_db/anacapa_QC_dada2.sh -i <input_dir> -o <out_dir> -d <database_directory> -u <hoffman_account_user_name> -f <fasta file of forward primers> -r <fasta file of reverse primers> -a <adapter type (nextera or truseq)>  -t <illumina run type HiSeq or MiSeq>
```

## How to run the bowtie2 blca step:

```
# ~/Anacapa_db/anacapa_bowtie2_blca.sh -o <out_dir_for_anacapa_QC_run> -d <database_directory> -u <hoffman_account_user_name>
```

! describe parameter b !

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
