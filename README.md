# Anacapa Toolkit

### Anacapa last updated 3-21-2018

#### Written by Emily Curd (eecurd@g.ucla.edu), Zack Gold (zjgold@ucla.edu), Gaurav Kandlikar (gkandlikar@ucla.edu), Jesse Gomer (jessegomer@gmail.com), , Baochen Shi (biosbc@gmail.com), Rachel Turba (rturba@ucla.edu) and Rachel Meyer (rsmeyer@ucla.edu).
#### Developed at UCLA for the University of California Conservation Consortium's CALeDNA Program

## Introduction
The toolkit is named for the iconic southern California island, Anacapa, that has significant cultural and biodiversity importance. The name derived from the Chumash word _Ennepah_ or _Anyapakh_ which translates to  "mirage island" (Bright, 2004; Gudde, 2010). Much like the name, using eDNA to monitor biodiversity seems like an illusion on the horizon, but like the real island, the __Anacapa__ toolkit can obtain true and quality results with full transparency of the caveats of eDNA. Here, we present __Anacapa__, an automated method to create custom reference databases and simultaneously analyze multiple metabarcoding reads produced by HiSeq and MiSeq Illumina sequence platforms, with a built-in exploration tool of the raw results output.

__Anacapa__ toolkit processes eDNA reads and assigns taxonomy using existing software or modifications to existing software. This modular toolkit is designed to analyze multiple samples and metabarcodes simultaneously from any Ilumina sequencing platform. __Anacapa__ accomplishes this in four steps: 1) building reference libraries using __CRUX__: Creating Reference libraries Using eXisting tools, 2) running quality control (QC) and assigning Amplicon Sequence Variants (ASV) using Dada2 (__Sequence QC and ASV Parsing__), 3) assigning taxonomy using bowtie2 and a bowtie2 specific Bayesian Least Common Ancestor (BLCA) (__Assignment__) and 4) Running exploratory data analysis and generating ecological diversity summary statistics for a set of samples (__ranacapa__). A significant advantage of the __Anacapa__ toolkit is that it does not require that paired reads overlap, or that both reads in a pair pass QC.  Taxonomy results are generated for all read types and the user can decide which read types they wish to retain for downstream analysis.

#### Step 1: CRUX: Creating Reference libraries Using eXisting tools
For full details on building reference libraries using CRUX, please refer to the following: https://github.com/limey-bean/CRUX_Creating-Reference-libraries-Using-eXisting-tools.

This first part of the toolkit generates reference libraries needed for taxonomic assignment using __CRUX__.  The output of __CRUX__ consists of two reference libraries, either unfiltered or filtered. Unfiltered libraries contain every dereplicated read found during the BLAST searches. The filtered library contains only reads with robust taxonomic assignments. Specifically we refer to robust taxonomic assignments as any reads that do not have the following in their taxonomic path: 'uncultured', 'environmental', 'sample', or 'NA;NA;NA;NA'. Prebuilt __CRUX__ reference libraries (12S - MiFish, 16S - EMP, 18S V4, 18S V8-9, 18S - EMP, PITS - Plant ITS2, CO1 and FITS - Fungal ITS [see Table 1] can be found at [link to dryad]. Each library contains unique metabarcode specific reads that correspond to NCBI accession version numbers. Libraries consist of fasta files, taxonomy files, and a bowtie2 index library.

<p align="center">
<img src="/figures-and-tables-for-the-Github/Table_1.png">
</p>

We acknowledge that users may wish to use their own custom sequences libraries to run __Anacapa__ or add additional custom sequences to a pre-made __CRUX__ reference library. Please refer to the __CRUX__ page (https://github.com/limey-bean/CRUX_Creating-Reference-libraries-Using-eXisting-tools) for instructions to create a library from a custom set of fasta formatted reads or add reads to a CRUX formatted reference library.

#### Step 2: Sequence QC and ASV Parsing using dada2
This next step of the toolkit aims to conduct standard sequence QC and then generate amplicon sequence variants (ASV) from Illumina data using **dada2** (Callahan et al. 2016). ASVs are a novel solution to identifying biologically informative unique sequences in metabarcoding samples that replaces the operational taxonomic unit (OTU) framework. Unlike OTUs which cluster sequences using an arbitrary sequence similarity (ex 97%), ASVs are unique sequence reads determined using Bayesian probabilities of known sequencing error. These unique sequences can be as little as 2 bp different, providing improved taxonomic resolution and an increase in observed diversity. Please see (Callahan et al. 2016, Amir et al. 2017) for further discussion.

An strong advantage of the __Anacapa__ toolkit is that is can simultaneously processes raw fastq reads for samples with single or multiple metabarcode targets generated on Illumina HiSeq and MiSeq machines. It is also not required that all samples contain reads for each metabarcode, thus allowing users to combine multiple projects or targets on the same sequencing run while only running the pipeline once.

__Anacapa__ takes raw Illumina fastq format reads and preprocesses them to assess file corruption (**md5sum**) and uncompresses (**gunzip**) and then renames the files for readability  readable. The QC portion of this script trims nextera and truseq adapters (**cutadapt**; Martin 2011), removes low quality reads **Fastx-toolkit**, and sorts reads by metabarcode primer sequence (**cutadapt**). Reads are trimmed using **cutadapt** (Martin 2011) to remove sequencing adapters from the 5' ends and sequencing adapters and primers from the 3' end of reads.  **Fastx-toolkit** (Gordon and Hannon, 2010) is then used to processed for quality control. Read are retained if they have a Q ≥ 35 and are at lease 100bp after adapter and 3' primer trimming. **cutadapt** is next used to sort reads by primer, and to trim additional basepairs from the end of read to increase quality going into **dada2**. Prior to running **dada2** a custom python script sorts reads into unpaired F, unpaired R and unmerged read files.  The files are passed separately into **dada2*** where they are denoised, dereplicated, merged (where possible), and  chimeric sequences removed from the data set.  


The input is raw Illumina metabarcode sequence data [\*.fastq.gz] reads and outputs are site frequency tables of ASVs and species count data for multiple samples and metabarcodes (ASV table). Successful implementation of this step requires: 1) raw illumina sequencing data and 2) a set of fasta formatted forward and reverse fasta format files that include the metabarcoding primers used to generate sequence data.

#### Step 3: Taxonomic Assignment using bowtie2 and BLCA
This next module of the pipeline assigns taxonomy to ASVs using **bowtie2** and a bowtie2 specific **Bayesian Least Common Ancestor** (**BLCA**) algorithm.

The **Anacapa** toolkit determines the best taxonomic hits for an ASV using **bowtie2** (Langmead and Salzberg, 2012). **Anacapa** considers paired merged, paired unmerged, and unpaired sequencing reads, and thus a fast and flexible read aligner, such as bowtie2, is required to handle all four read types. This script uses bowtie2, **CRUX** reference libraries, and **BLCA** to then assign taxonomy. All reads are first globally aligned against the **CRUX** database using **bowtie2**. Any reads that fail to align are then aligned locally. The best hits (the top 100 **bowtie2** returns) are then processed with **BLCA** script to assign taxonomy. The **bowtie2 BLCA** algorithm was adapted from https://github.com/qunfengdong/BLCA. **BLCA** uses pairwise sequence alignment to calculate sequence similarity between query sequences and reference library hits. Taxonomy is assigned based on the lowest common ancestor of multiple reference library hits for each query sequence. The reliability of each taxonomic assignment is then evaluated through bootstrap confidence scores [Gao et al. 2017].

Successful implementation of this script requires an ASV table (site frequency table) with the following columns (ASV number, Sequence, Samples). The output is multiple ASV tables for each metabarcode and both global and local alignments.


#### Step 4: ranacapa: Data exploration
This portion generates ecological diversity summary statistics for a set of samples.

The last step of the **Anacapa** Pipeline conducts exploratory data analysis to provide a first pass look at sequencing depth, taxonomic assignments, and generated data tables. This analysis is not meant for publication, but solely as a first stab at visualization of your data. This is helpful in identifying potential glaring errors or contamination, and identifying patterns worth investigating further through more robust analysis. We highly encourage data exploration before further analysis as different parameters within the **Anacapa** pipeline may produce differences in downstream results and these parameters will vary by project, stringency of taxonomic assignment, and users opinions. The exploratory_analysis.R script uses a variety of **R** packages, relying heavily on __phyloseq__, __vegan__, and __ggplot2__. See below for full list of R package dependencies and scripts. The output from reformat_summary_for_R.py is an ASV site frequency table with assigned taxonomy and is the input used for the exploratory_analysis.R script. In addition, the user supplies an input metadata table that only requires the first column be sample names. Users can include any type of metadata including categorical, continuous, and discete variables. The first step of the R script is to convert the input files into a **Phyloseq** class object. We then generate bar plots looking at total number of observed classes and relative abundance of each class. We then generate rarefaction curves, alpha diversity boxplots to observe total number of taxa and Shannon diversity, and alpha diversity statistics. In addition, we calculate jaccard and bray-curtis distance matrices and conduct NMDS ordination plots, network map, heat maps, and ward-linkage maps. Each of the above analyses are repeated with different grouping for each metadata column. In addition we conduct two betadiversity statistical tests, pairwise adonis and betadisp from the vegan package. Again each analyses is repeated across groupings of each metadata column.


## Required Programs and Dependencies
### Anacapa_db folder
* Four files:
 * **anacapa_QC_dada2.sh**
	 * _This script runs  QC and ASV generation to generate ASV tables_
 * **anacapa_bowtie2_blca.sh**
	 * _This script runs bowtie2 and BLCA to assign taxonomy_
 * **forward_primers.txt**
 * **reverse_primers.txt**
	 * _The two primer files are examples of how to format the primer forward and reverse input files.  It is_ **VERY IMPORTANT** _that you modify these files or make new files to reflect your data set!_
* Two folders
	* **adapters_and_PrimAdapt_rc/**
		* _The forward and reverse nextera adapters_
		* _The forward and reverse trueseq adapters_
	* **scripts/**
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
       * _These scripts have documentation in their headers. Please see each script to understand functionality._

### Dependencies

To run __Anacapa__, you need verify that the full path to each of the following programs is correctly indicated in the anacapa_config.sh file.  

1. __cutadapt__: http://cutadapt.readthedocs.io/en/stable/index.html

2. __fastxtoolkit__ (Version: 0.0.13) http://hannonlab.cshl.edu/fastx_toolkit/

3. __anaconda/python2-4.2__
	* make sure biopython is installed http://biopython.org/wiki/Packages
	* We recommend downloading using conda from conda

4. R 3.4.2
   * Cran Packages
	 * __ggplot2__
	 * __plyr__
	 * __dplyr__
	 * __seqRFLP__
	 * __reshape2__
	 * __tibble__
	 * __devtools__
	 * __Matrix__
	 * __mgcv__
	 * __readr__
	 * __stringr__
	 * __vegan__
	 * __plotly__
	 * __otparse__
	 * __ggrepel__
	 * __cluster__
 * __BIOCLITE__ Packages http://bioconductor.org/biocLite.R
	* __phyloseq__
	* __genefilter__
	* __impute__
	* __Biostrings__
 * __dada2__ https://github.com/benjjneb/dada2

5. __Bowtie2__: http://bowtie-bio.sourceforge.net/bowtie2/index.shtml
	* We recommend downloading using conda

6. muscle: https://www.drive5.com/muscle/downloads.htm
	* version muscle3.8.31
	* **__muscle must be installed within the anacapa_db folder__**


### CRUX Databases
Download taxonomy reference libraries from this google drive folder: https://drive.google.com/drive/folders/0BycoA83WF7aNOEFFV2Z6bC1GM1E?usp=sharing

Users can also make their own libraries using CRUX scripts. For example, Silva and greengeens libraries can easily be converted to __Anacapa__ compatible libraries. See https://github.com/limey-bean/CRUX_Creating-Reference-libraries-Using-eXisting-tools/tree/master/crux_release_V1_db/scripts for documentation.

Reference library folders must be transfered to the Anacapa_db folder.



## Running Anacapa

__Anacapa__ must be run in either local or default mode. Local mode is for personal computers and servers that are not Hoffman2 at UCLA. Default mode is for running on Hoffman2 at UCLA.



#### Preparing to Run Anacapa

Before running the __Anacapa__ toolkit you need to double check the anacapa_config.sh file and update the appropriate paths. For local mode set LOCALMODE=TRUE, CUTADAPT ="cutadapt",and  MUSCLE="muscle" ; replace all other values to "". Double check that all dependencies work in the terminal. This is the key for success.

#### How to run the QC / dada2 step:
```
sh ~/Anacapa_db/anacapa_QC_dada2.sh -i <input_dir> -o <out_dir> -d <database_directory> -u <hoffman_account_user_name> -f <fasta file of forward primers> -r <fasta file of reverse primers> -a <adapter type (nextera or truseq)>  -t <illumina run type HiSeq or MiSeq>
```

#### How to run the bowtie2 blca step:

```
# ~/Anacapa_db/anacapa_bowtie2_blca.sh -o <out_dir_for_anacapa_QC_run> -d <database_directory> -u <hoffman_account_user_name>
```
#### Hoffman Cluster UCLA
**_Hoffman users running the QC dada2 need to do the following before dada2 will run_**
```
qrsh
module load R/3.4.2
module load gcc/6.3.0
R
source("https://bioconductor.org/biocLite.R")
biocLite(suppressUpdates = FALSE)
```
When given this option "would you like to use a personal library" say "y"
This is the reason that it is not possible to install this in the R script.
```
biocLite("ShortRead", suppressUpdates = FALSE)
biocLite("devtools")
```
You may need to install MASS, mgcv, and rpart.
```
install.packages("MASS")
install.packages("mgcv")
install.packages("rpart")
```

Installing these dependencies may take a very long time so no worries...

Biopython must also be installed.

```
module load anconda
pip install biopython --user
```
"user" not your user name


This is a preliminary attempt a documentation...  To run the latest version of Anacapa using the 6 primer sets that are commonly used for CALeDNA projects, you will need to download the following:
* Download the Anacapa_db folder
* Download taxonomy reference libraries (updated 10-16-2017) from this google drive folder: https://drive.google.com/drive/folders/0BycoA83WF7aNOEFFV2Z6bC1GM1E?usp=sharing


### More to come, and it might be a bit buggy.
* if you choose to take this on...  Good Luck!



## References
Amir, A., McDonald, D., Navas-Molina, J.A., Kopylova, E., Morton, J.T., Xu, Z.Z., Kightley, E.P., Thompson, L.R., Hyde, E.R., Gonzalez, A. and Knight, R., 2017. Deblur rapidly resolves single-nucleotide community sequence patterns. MSystems, 2(2), pp.e00191-16.

Callahan, B. J., McMurdie, P. J., Rosen, M. J., Han, A. W., Johnson, A. J. A., & Holmes, S. P. (2016). DADA2: High-resolution sample inference from Illumina amplicon data. Nature Methods, 13(7), 581–583. doi:10.1038/nmeth.3869


Caporaso, J.G., Kuczynski, J., Stombaugh, J., Bittinger, K., Bushman, F.D., Costello, E.K., Fierer, N., Peña, A.G., Goodrich, J.K., Gordon, J.I. and Huttley, G.A., 2010. QIIME allows analysis of high-throughput community sequencing data. Nature methods, 7(5), pp.335-336.

Edgar, R.C., 2010. Search and clustering orders of magnitude faster than BLAST. Bioinformatics, 26(19), pp.2460-2461.

Gordon, A. and Hannon, G.J., 2010. Fastx-toolkit. FASTQ/A short-reads preprocessing tools (unpublished) http://hannonlab. cshl. edu/fastx_toolkit, 5.

Gudde, Erwin; William Bright (2004). California Place Names (Fourth ed.). University of California Press. p. 12. ISBN 0-520-24217-3.

Langmead, B. and Salzberg, S.L., 2012. Fast gapped-read alignment with Bowtie 2. Nature methods, 9(4), pp.357-359.

Martin, M., 2011. Cutadapt removes adapter sequences from high-throughput sequencing reads. EMBnet. journal, 17(1), pp.pp-10.

Team, R.C., 2000. R language definition. Vienna, Austria: R foundation for statistical computing.

Zhang, J., Kobert, K., Flouri, T. and Stamatakis, A., 2013. PEAR: a fast and accurate Illumina Paired-End reAd mergeR. Bioinformatics, 30(5), pp.614-620.
