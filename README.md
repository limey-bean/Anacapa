# Anacapa

### Anacapa_release_V1		09-07-2017
#### Written by Emily Curd (eecurd@g.ucla.edu), Jesse Gomer (jessegomer@gmail.com), Baochen Shi (biosbc@gmail.com), Zack Gold (zjgold@ucla.edu), Gaurav Kandlikar (gkandlikar@ucla.edu), and Rachel Meyer (rsmeyer@ucla.edu). 
#### Developed at UCLA for the University of California Conservation Consortium's CALeDNA Program

## Introduction
Anacapa Island's name is derived from the Chumash --Ennepah-- or --Anyapakh-- which means "mirage island" (Gudde and Bright, 2004). Much like Anacapa Island, using environmental DNA (eDNA) to uncover Biodiversity seems like an issusion on the horizon, however both are real phenomoena.  

Anacapa is an automated metabarcoding read processing pipeline.  It is designed to analyze multiple samples and metabarcodes simultaneously. The input is raw Illumina metabarcode sequence data and outputs are species count data for multiple samples and metabarcoes. Sucessful implementaion of Anacap requires: 1) raw illumina data, 2) primers for the metabrcodes of interest, 3) reference libraries for the metabrcodes of interest (for custom libraries check out CRUX), and 4) the dependencies indicated below. Anacapa takes raw Illumina fastq format reads and preprocesses them to assess file corruption (md5sum), uncompresses (gunzip), and renam files and reads within files.  Reads are next processed for quality control.  PEAR (Zhang et al., 2013) is used to merge quality trimmed reads (Q ≥ 30) of at least 125bp.  Sequencing adapters are removed and fastq files are converted to fasta files with cutadapt (Martin, 2011), and some reads are reverse complemented using Fastx-toolkit (Gordon and Hannon, 2010).  Metabarcoded reads generated using different primer sets are then sorted by primer set using split primers (https://github.com/jessegomer/Split_on_Primer). The resulting reads are assigned taxonomy using an open reference implementation of Uclust (Edgar, 2010) in Qiime (Caporaso et al., 2010).  





## References
Caporaso, J.G., Kuczynski, J., Stombaugh, J., Bittinger, K., Bushman, F.D., Costello, E.K., Fierer, N., Peña, A.G., Goodrich, J.K., Gordon, J.I. and Huttley, G.A., 2010. QIIME allows analysis of high-throughput community sequencing data. Nature methods, 7(5), pp.335-336.

Edgar, R.C., 2010. Search and clustering orders of magnitude faster than BLAST. Bioinformatics, 26(19), pp.2460-2461.

Gordon, A. and Hannon, G.J., 2010. Fastx-toolkit. FASTQ/A short-reads preprocessing tools (unpublished) http://hannonlab. cshl. edu/fastx_toolkit, 5.

Gudde, Erwin; William Bright (2004). California Place Names (Fourth ed.). University of California Press. p. 12. ISBN 0-520-24217-3.

Martin, M., 2011. Cutadapt removes adapter sequences from high-throughput sequencing reads. EMBnet. journal, 17(1), pp.pp-10.

Zhang, J., Kobert, K., Flouri, T. and Stamatakis, A., 2013. PEAR: a fast and accurate Illumina Paired-End reAd mergeR. Bioinformatics, 30(5), pp.614-620.
