#!/bin/bash

maindir=""
name=""
per=""

while getopts "m:n:p:" opt; do
    case $opt in
        m) maindir="$OPTARG" # path to desired Anacapa output
        ;;
        n) name="$OPTARG"  # path to Anacapa_db
        ;;
        p) per="$OPTARG"  # path to Anacapa_db
        ;;
    esac
done

# load modules
source /u/local/Modules/default/init/bash
module load ATS
module load qiime/1.8.0
module load bowtie2

#directory where everything is kept
#maindir=/u/project/rwayne/eecurd/eecurd/blca_test/  # <<<<<< pass in directory

#change directory into the scripts folder so hoffman can find python modules
cd ${maindir}scripts

#	get prefix and make it name
#	get test name 
#	get training name


#name=ncbi_16s_${per}                                   # <<<<<< pass in name
#made Edgar's training set an "all" (all reads possible) sequence set, because sometimes test reads were included with training reads.  I do a filter step to make sure that the test reads are removed from the training step.  I use Qiime to do this
testfa=${maindir}test_train/${name}_test.fasta


#####################  make directories
mkdir -p ${maindir}test_train/
traindir=${maindir}test_train/${name}
mkdir -p $traindir

mkdir -p ${maindir}namecounts
namecountsdir=${maindir}namecounts/${name}
mkdir -p $namecountsdir

mkdir -p ${maindir}bowtie2_libs/
bowtie2dir=${maindir}bowtie2_libs/${name}
mkdir -p $bowtie2dir

mkdir -p ${maindir}raw_output/
rawdir=${maindir}raw_output/${name}
mkdir -p $rawdir
raw=$rawdir

mkdir -p ${maindir}pred/
preddir=${maindir}pred/${name}
mkdir -p $preddir
pred=$preddir/${name}

mkdir -p ${maindir}stats/
statsdir=${maindir}stats/${name}
mkdir -p $statsdir


################# prepare files -> make taxonomy, files partition data, move things around, get ready for actual classification
#make taxonomy file for test data
grep -e ">" ${testfa} | awk 'sub(/^>/, "")'  >  ${maindir}test_train/${name}_orig_test_taxonomy.txt
#filer test data from full dataset to make a training dataset
filter_fasta.py -f ${maindir}test_train/${name}_all.fasta -o ${maindir}test_train/${name}/${name}_train.fasta -s ${maindir}test_train/${name}_orig_test_taxonomy.txt -n
#move file to appropriate test_train directory
mv ${maindir}test_train/${name}_* ${maindir}test_train/${name}
# make BLCA compatible taxonomy and fasta files for training vs test sets
python2 ${maindir}/scripts/convert_to_BLCA_format.py ${maindir}test_train/${name}/${name}_test.fasta ${maindir}test_train/${name}/${name}_blca_test.fasta ${maindir}test_train/${name}/${name}_blca_test_taxonomy.txt
python2 ${maindir}/scripts/convert_to_BLCA_format.py ${maindir}test_train/${name}/${name}_train.fasta ${maindir}test_train/${name}/${name}_blca_train.fasta ${maindir}test_train/${name}/${name}_blca_train_taxonomy.txt


################## make name count files.  These are critical for the stats!!!!
# get taxonomy from full dataset fasta files
grep -e ">" ${maindir}test_train/${name}/${name}_train.fasta | awk 'sub(/^>/, "")'  >  ${maindir}test_train/${name}/${name}_train_taxonomy.txt
#need this file in anacapa_blca_format REF_GI_343198724	Bacteria,Firmicutes,Bacilli,Bacillales,Listeriaceae,Listeria,Listeria_welshimeri
cat ${maindir}test_train/${name}/${name}_train_taxonomy.txt | sed 's/;tax=d:/'$'\t/g' | sed 's/[pcfogs]://g' | sed 's/;$//' | sed 's/gi_/REF_GI_/g' | sed 's/,/;/g' > ${maindir}test_train/${name}/${name}_anacapa_blca_train_taxonomy.txt

# get the frequency of occurrences of a taxonomic rank make a file 
python2 run_namecount.py ${maindir}test_train/${name}/${name}_train_taxonomy.txt ${namecountsdir} ${name}
# clean up the file by removing the damn space and also any lines without : and then remove temp files
cat ${namecountsdir}/${name}_namecount | sed 's/ //g' > ${namecountsdir}/${name}_namecount.txt
grep ':' ${namecountsdir}/${name}_namecount.txt > ${namecountsdir}/${name}_namecount.txt2; mv ${namecountsdir}/${name}_namecount.txt2 ${namecountsdir}/${name}_namecount.txt
rm ${namecountsdir}/${name}_namecount


################### run anacapa BLCA
# prepare for anacapa BLCA by making a bowtie2 database
bowtie2-build -f ${maindir}test_train/${name}/${name}_blca_train.fasta ${bowtie2dir}/${name}_bowtie2_index
# need biopython -> might be a waste of effort but Qiime alone did not seem to be working
module load anaconda
echo "global"
bowtie2 -x ${bowtie2dir}/${name}_bowtie2_index  -f -U ${maindir}test_train/${name}/${name}_blca_test.fasta -S ${maindir}test_train/${name}/${name}_test_end_to_end.sam --no-hd --no-sq --very-sensitive --end-to-end --no-unal -p 120 -k 50 --un ${maindir}test_train/${name}/${name}_test_end_to_end_reject.fasta
#unmerged pair reads local
echo "local"
bowtie2 -x ${bowtie2dir}/${name}_bowtie2_index  -f -U ${maindir}test_train/${name}/${name}_test_end_to_end_reject.fasta -S ${maindir}test_train/${name}/${name}_test_local.sam --no-hd --no-sq --very-sensitive --local --no-unal -p 120 -k 50 --un ${maindir}test_train/${name}/${name}_test_end_to_end_and_local_reject.fasta
### concat all of the sam files for blca
cat ${maindir}test_train/${name}/*.sam > ${maindir}test_train/${name}/${name}_bowtie2_all.sam
### run BLCA
python ${maindir}scripts/blca_from_bowtie.py -i ${maindir}test_train/${name}/${name}_bowtie2_all.sam -r ${maindir}test_train/${name}/${name}_anacapa_blca_train_taxonomy.txt -q ${maindir}test_train/${name}/${name}_blca_train.fasta -b .90 -l 0.85 -n 100 

################### Summary stats
# move a files around
python ${maindir}scripts/rename_anacapa_blca_to_blast_BLCA.py ${maindir}test_train/${name}/${name}_bowtie2_all.sam.blca.out	${raw}/raw_BLCA.out.${per}.tmp
cat ${raw}/raw_BLCA.out.${per}.tmp | sed 's/ //g' > ${raw}/raw_BLCA.out.${per}; rm ${raw}/raw_BLCA.out.${per}.tmp
# make a file that compares the actual taxonomy with the BLCA determined taxonomy
python2 ${maindir}scripts/blca2tab.py ${rawdir}/raw_BLCA.out.${per} ${maindir}test_train/${name}/${name}_test.fasta > $preddir/pred_BLCA.out.${per}
# make a stats file
for i in "s" "g" "f" "o" "c" "p"
do
	python2 ${maindir}scripts/taxbench.py $preddir/pred_BLCA.out.${per} anacapa_blca ${i} ${name} ${namecountsdir} >> $statsdir/stats_anacapa_BLCA.out.${per}
done