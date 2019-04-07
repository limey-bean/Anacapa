#!/bin/bash
# This script was designed to run a k-folds test on a specific dataset and then run the TAXXI framework (Edgar, R.C., 2018. Accuracy of taxonomy prediction for 16S rRNA and fungal ITS sequences. PeerJ, 6, p.e4652.) using the Anacapa Bowtie2-BLCA classifier.
# generate k-folds databases then use the /bin/bash Edgar_anacapa_BLCA_classifier_test_config.sh to test for TAXXI matrics
# example usage: sh K_folds_test.sh -f 2 -d ~/Anacapa_db/12S/12S_fasta_and_taxonomy -t ~/test12S -n 12S



database_dir=""  # path to the CRUX formated data base e.g. ~/Anacapa_db/12S/12S_fasta_and_taxonomy
kfold=""
name=""
ttdir=""

while getopts "d:f:n:t:" opt; do
    case $opt in
        d) database_dir="$OPTARG" # path to the CRUX formated data base e.g. ~/Anacapa_db/12S/12S_fasta_and_taxonomy
        ;;
        f) kfold="$OPTARG" # the number of folds
        ;;
        n) name="$OPTARG"  # name of the database to test e.g ${name}_test.fasta  e.g. for ncbi_16s_100_test.fasta the $name=ncbi_16s_100
        ;;
        t) ttdir="$OPTARG" # path to test and training directory
        ;;

    esac
done

####### randomly
#subsample your dataset
mkdir -p ${ttdir}/temp

cut -d "." -f 1 ${database_dir}/${name}_taxonomy.txt > ${ttdir}/${name}_accession.tmp.txt
sort -R ${ttdir}/${name}_accession.tmp.txt > ${ttdir}/${name}_accession_rand.txt
#split -${kfold} ${ttdir}/${name}_accession_rand.txt
rm ${ttdir}/${name}_accession.tmp.txt
total_lines=$(wc -l <${database_dir}/${name}_taxonomy.txt)
((lines_to_split = ( ( $total_lines / $kfold ) + 1 )))
#echo  $total_lines
#echo $lines_to_split
# Split the actual file, maintaining lines.
cd ${ttdir}/temp
split --lines=${lines_to_split} ${ttdir}/${name}_accession_rand.txt


i=1
for j in ${ttdir}/temp/*
do
  #echo ${j}
  grep -Fxvf ${j} ${ttdir}/${name}_accession_rand.txt > ${ttdir}/${name}_accession_all.txt
  grep -F -f ${j} ${database_dir}/${name}_taxonomy.txt > ${ttdir}/${name}_${i}_test_taxonomy.txt
  grep -F -f ${ttdir}/${name}_accession_all.txt ${database_dir}/${name}_taxonomy.txt > ${ttdir}/${name}_${i}_all_taxonomy.txt
  grep -F -f ${j} ${database_dir}/${name}_taxonomy.txt > ${ttdir}/${name}_${i}_test.fasta
  grep -F -f ${j} -A1 ${database_dir}/${name}_.fasta | sed 's/--//g' | sed -r '/^\s*$/d' > ${ttdir}/${name}_${i}_test.fasta
  grep -F -f ${ttdir}/${name}_accession_all.txt -A1 ${database_dir}/*.fasta | sed 's/--//g' | sed -r '/^\s*$/d' > ${ttdir}/${name}_${i}_all.fasta
i=$(expr $i + 1)
done

#rm ${ttdir}/${name}_accession_all.txt
#rm ${ttdir}/${name}_accession_rand.txt
rm -r ${ttdir}/temp/
