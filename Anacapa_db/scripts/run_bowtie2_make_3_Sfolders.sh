#! /bin/bash

### this script is run as follows
# sh ~/Anacapa_db/scripts/anacapa_release_V1.sh -o <working/out_dir> -d <database_directory> -n name  
OUT=""
DB=""
NAME=""

while getopts "o:d:n:" opt; do
    case $opt in
        o) OUT="$OPTARG" # path to desired Anacapa output
        ;;
        d) DB="$OPTARG"  # path to Anacapa_db
        ;;
        n) NAME="$OPTARG"  # need username for submitting sequencing job
        ;;
    esac
done
####################################script & software

##load module
source $DB/scripts/anacapa_vars.sh
source $DB/scripts/anacapa_config.sh

##load module
${MODULE_SOURCE} # use if you need to load modules from an HPC

${QIIME} #load qiime
${BOWTIE2} #load bowtie2

##### add the single and paired bowtie 2 files to different folders. Turn the following code into a for loop for the single bowtie2 reads



for pro in ${FILESTOPROCESS}
do
 echo "processing ${pro}"
 str1=${pro%/${NAME}*.fasta}
 folder=${str1#${OUT}/bowtie2_runs/${NAME}/}
 echo "processing ${folder}"
 echo "${pro}"
 

#############################################
# process the ${folder} reads first
#############################################

### Sort everything first with 99 and the base overhand amount
echo " "
echo "First run bowtie2 on ${NAME} reads using the ${NAME} 99 % reference library and allow the subject a maximum of ${BASEOVERHANG} bp overhang from the ends of the best matching reference(s)"
mkdir -p ${OUT}/bowtie2_runs/${NAME}/${folder}/Sort_${BASEOVERHANG}max_overhang  # folder that stores all of the isnitial split info...
echo "...bowtie2 running"
bowtie2 -x ${DB}/${NAME}/${NAME}_bowtie2_databases/${NAME}_bowtie2_.99/${NAME}_.99_bowtie2_index -f -U ${pro} --un ${OUT}/bowtie2_runs/${NAME}/${folder}/${NAME}_all.clean_${folder}_bowtie2_rejects.fasta -S ${OUT}/bowtie2_runs/${NAME}/${folder}/${NAME}_all.clean_${folder}_bowtie2.sam --no-hd --no-sq --very-sensitive-local  --no-unal -p 120  #run Bowtie2 on the entire dataset usng 99% ref lib
echo "...bowtie2 finished"
date
echo " "
mkdir -p ${OUT}/bowtie2_runs/${NAME}/${folder}/Sort_${BASEOVERHANG}max_overhang/bowtie2_at_.99       # get best hits for 99% at S25
python ${DB}/scripts/group_alignments_to_files_p_mod.py ${OUT}/bowtie2_runs/${NAME}/${folder}/Sort_${BASEOVERHANG}max_overhang/bowtie2_at_.99/ ${OUT}/bowtie2_runs/${NAME}/${folder}/${NAME}_all.clean_${folder}_bowtie2.sam ${BASEOVERHANG} .99 # sort reads into -> single best hits, multiple hits, to short.  Also sort my percent similarity to the 99% reference
echo "Summarizing single best hits for 99% with up to ${ove} basepair overhand "
mkdir -p ${OUT}/bowtie2_runs/${NAME}/${folder}/Sort_${BASEOVERHANG}max_overhang/single_best    # get best hits for 99% at S25
python ${DB}/scripts/summarize_bowtie2_hits_full_taxonomy.py ${DB}/${NAME}/${NAME}_final_database/${NAME}_taxonomy_.99.txt ${OUT}/bowtie2_runs/${NAME}/${folder}/Sort_${BASEOVERHANG}max_overhang/bowtie2_at_.99/bowtie2_good_hits.txt ${OUT}/bowtie2_runs/${NAME}/${folder}/Sort_${BASEOVERHANG}max_overhang/single_best/single_best${NAME}_.99_single_best.txt

###

### run bowtie2 at different overhang levels for reads with an overhand > base overhang for 99% ref
i="${BASEOVERHANG}"
for ove in ${OVERHANG}
do
 mkdir -p ${OUT}/bowtie2_runs/${NAME}/${folder}/plus${i}_unmapped_at_ends 
 cat ${OUT}/bowtie2_runs/${NAME}/${folder}/Sort_${i}max_overhang/bowtie2_at_.99/bowtie2_rejects_not_mapped_at_ends.*.fasta > ${OUT}/bowtie2_runs/${NAME}/${folder}/plus${i}_unmapped_at_ends/plus${i}_unmapped_at_ends.fasta
 rm ${OUT}/bowtie2_runs/${NAME}/${folder}/Sort_${i}max_overhang/bowtie2_at_.99/bowtie2_rejects_not_mapped_at_ends.*.fasta
 ###
 echo "...use Bowtie2 to sort for a maximum of ${ove} bp using the 99% reference database"
 bowtie2 -x ${DB}/${NAME}/${NAME}_bowtie2_databases/${NAME}_bowtie2_.99/${NAME}_.99_bowtie2_index -f -U ${OUT}/bowtie2_runs/${NAME}/${folder}/plus${i}_unmapped_at_ends/plus${i}_unmapped_at_ends.fasta --un ${OUT}/bowtie2_runs/${NAME}/${folder}/plus${i}_unmapped_at_ends/plus${i}_unmapped_at_ends_rejects.fasta -S ${OUT}/bowtie2_runs/${NAME}/${folder}/plus${i}_unmapped_at_ends/unmapped_at_ends.sam --no-hd --no-sq --very-sensitive-local  --no-unal -p 120  #run Bowtie2 on the entire dataset usng 99% ref lib
 echo "...bowtie2 finished"
 date
 echo " "
 mkdir -p ${OUT}/bowtie2_runs/${NAME}/${folder}/Sort_${ove}max_overhang/bowtie2_at_.99
 python ${DB}/scripts/group_alignments_to_files_p_mod.py ${OUT}/bowtie2_runs/${NAME}/${folder}/Sort_${ove}max_overhang/bowtie2_at_.99 ${OUT}/bowtie2_runs/${NAME}/${folder}/plus${i}_unmapped_at_ends/unmapped_at_ends.sam ${ove} .99 # sort reads into -> single best hits, multiple hits, to short.  Also sort my percent similarity to the 99% reference
 echo "Summarizing single best hits for 99% with up to ${ove} basepair overhang"
 mkdir -p ${OUT}/bowtie2_runs/${NAME}/${folder}/Sort_${ove}max_overhang/single_best     # get best hits for 99% at S25
 python ${DB}/scripts/summarize_bowtie2_hits_full_taxonomy.py ${DB}/${NAME}/${NAME}_final_database/${NAME}_taxonomy_.99.txt ${OUT}/bowtie2_runs/${NAME}/${folder}/Sort_${ove}max_overhang/bowtie2_at_.99/bowtie2_good_hits.txt ${OUT}/bowtie2_runs/${NAME}/${folder}/Sort_${ove}max_overhang/single_best/single_best${NAME}_99_single_best.txt
 date
 echo " "
 i="${ove}"
done
####

for hang in ${ALLOVER}
do 
 e="${BASESIM}"
 for perc in ${SIM}
 do
  if (( $(bc -l <<<"${perc} > ${TOOLOW}") )); then
  mkdir -p ${OUT}/bowtie2_runs/${NAME}/${folder}/Sort_${hang}max_overhang/bowtie2_at_${perc}
  mkdir -p ${OUT}/bowtie2_runs/${NAME}/${folder}/Sort_${hang}max_overhang/bowtie2_at_${perc}/fasta_to_process/
  p=${perc}
  #### acquire list of input files
  echo "--------------${hang}---------------------------${perc}---------------------"
  for str in `ls ${OUT}/bowtie2_runs/${NAME}/${folder}/Sort_${hang}max_overhang/bowtie2_at_${e}/bowtie2_rejects_multiple_hits*.fasta`
  do
   a=""
   b=""
   str1=${str%-*.fasta}
   j=${str1#${OUT}/bowtie2_runs/${NAME}/${folder}/Sort_${hang}max_overhang/bowtie2_at_${e}/bowtie2_rejects_multiple_hits}
   if (( $(bc -l <<<"${j} >= ${p}") )); then
    echo "${j}"
    a="${str}"
    echo "$a"
   fi
   ####
   if [ "x$a" != "x"  ] && [ "x$list_for_cat_del" != "x"  ];
   then
    list_for_cat_del="$a ${list_for_cat_del}"
   elif [ "x$a" != "x" ] && [ "x$list_for_cat_del" == "x"  ];
   then
    list_for_cat_del="$a"
   fi
   ###
   echo "a"
  done
  for st in `ls ${OUT}/bowtie2_runs/${NAME}/${folder}/Sort_${hang}max_overhang/bowtie2_at_${e}/bowtie2_rejects_low_percent_id*.fasta`
  do
   str2=${st%-*.fasta}
   k=${str2#${OUT}/bowtie2_runs/${NAME}/${folder}/Sort_${hang}max_overhang/bowtie2_at_${e}/bowtie2_rejects_low_percent_id}
   if (( $(bc -l <<<"${k} >= ${p}") ));
   then
    echo "${k}"
    b="${st}"
    echo "$b"
   fi
   ####
   if [ "x$b" != "x"  ] && [ "x$list_for_cat_del" != "x"  ];
   then
    list_for_cat_del="${b} ${list_for_cat_del}"
   elif [ "x$b" != "x" ] && [ "x$list_for_cat_del" == "x"  ];
   then
    list_for_cat_del="${b}"
   fi
   ###
   echo 
  done
  if [ "x${list_for_cat_del}" != "x"  ];
  then
   cat ${list_for_cat_del} > ${OUT}/bowtie2_runs/${NAME}/${folder}/Sort_${hang}max_overhang/bowtie2_at_${perc}/fasta_to_process/reads_from_${e}.fasta
   rm ${list_for_cat_del}
  fi
  list_for_cat_del=""
  p=""
  fi
 done
 e="${perc}"
 echo "${e}"
done  

##### run bowtie2 on the subfolders.....

move_fasta=""
for hang in ${ALLOVER}
do
 for perc in ${SIMBT}
 do
  if [ "x$move_fasta" != "x" ];
  then
    echo "${perc}"
    cat ${move_fasta} > ${OUT}/bowtie2_runs/${NAME}/${folder}/Sort_${hang}max_overhang/bowtie2_at_${perc}/fasta_to_process/from_${p}.fasta
    rm ${move_fasta}
  fi
  cat ${OUT}/bowtie2_runs/${NAME}/${folder}/Sort_${hang}max_overhang/bowtie2_at_${perc}/fasta_to_process/*.fasta > ${OUT}/bowtie2_runs/${NAME}/${folder}/Sort_${hang}max_overhang/bowtie2_at_${perc}/fasta_to_process/sam_input.fasta
  echo "Run bowtie2 on ${hang} with ${perc} ref library"
  bowtie2 -x ${DB}/${NAME}/${NAME}_bowtie2_databases/${NAME}_bowtie2_${perc}/${NAME}_${perc}_bowtie2_index -f -U ${OUT}/bowtie2_runs/${NAME}/${folder}/Sort_${hang}max_overhang/bowtie2_at_${perc}/fasta_to_process/sam_input.fasta --un ${OUT}/bowtie2_runs/${NAME}/${folder}/Sort_${hang}max_overhang/bowtie2_at_${perc}/fasta_to_process/sam_input_${hang}${perc}_rejects.fasta -S ${OUT}/bowtie2_runs/${NAME}/${folder}/Sort_${hang}max_overhang/bowtie2_at_${perc}/input_${hang}${perc}.sam --no-hd --no-sq --very-sensitive-local  --no-unal -p 120  #run Bowtie2 on the entire dataset usng 97% ref lib
  #sort sam to catergories
  python ${DB}/scripts/group_alignments_to_files_p_mod.py ${OUT}/bowtie2_runs/${NAME}/${folder}/Sort_${hang}max_overhang/bowtie2_at_${perc}/ ${OUT}/bowtie2_runs/${NAME}/${folder}/Sort_${hang}max_overhang/bowtie2_at_${perc}/input_${hang}${perc}.sam ${hang} ${perc} # sort reads into -> single best hits, multiple hits, to short.  Also sort my percent similarity to the 97% reference
  echo "Summarizing single best hits for ${perc} with up to ${hang} basepair overhang "
  python ${DB}/scripts/summarize_bowtie2_hits_full_taxonomy.py  ${DB}/${NAME}/${NAME}_final_database/${NAME}_taxonomy_${perc}.txt ${OUT}/bowtie2_runs/${NAME}/${folder}/Sort_${hang}max_overhang/bowtie2_at_${perc}/bowtie2_good_hits.txt ${OUT}/bowtie2_runs/${NAME}/${folder}/Sort_${hang}max_overhang/single_best/single_best${NAME}_${perc}_single_best.txt
  date
  echo " "
  move_fasta=""
  p="${perc}"
  for f in `ls ${OUT}/bowtie2_runs/${NAME}/${folder}/Sort_${hang}max_overhang/bowtie2_at_${perc}/*.fasta`
  do
   move_fasta="${f} ${move_fasta}"
  done
 #move rejected fasta down to the next percent lower file
 done
done

# now summarize this stuff!


mkdir -p ${OUT}/bowtie2_runs/${NAME}/${folder}/Sort_${hang}max_overhang/merged

#make biom tables
list=""
for hang in ${ALLOVER}
do
 for perc in ${FINSIM}
 do
  biom convert -i ${OUT}/bowtie2_runs/${NAME}/${folder}/Sort_${hang}max_overhang/single_best/single_best${NAME}_${perc}_single_best.txt -o ${OUT}/bowtie2_runs/${NAME}/${folder}/Sort_${hang}max_overhang/single_best/single_best${NAME}_${perc}_single_best.biom --table-type "otu table"
 done
 for f in `ls ${OUT}/bowtie2_runs/${NAME}/${folder}/Sort_${hang}max_overhang/single_best/*.biom`
  do
   list="${f} ${list}"
 done
 merge_otu_tables.py -i ${list} -o ${OUT}/bowtie2_runs/${NAME}/${folder}/Sort_${hang}max_overhang/merged/max_overhang_${hang}_sum_merged.biom
 biom convert -i ${OUT}/bowtie2_runs/${NAME}/${folder}/Sort_${hang}max_overhang/merged/max_overhang_${hang}_sum_merged.biom -o ${OUT}/bowtie2_runs/${NAME}/${folder}/Sort_${hang}max_overhang/merged/max_overhang_${hang}_sum_merged.txt -b --header-key="taxonomy" --output-metadata-id="Consensus Lineage"
 biom summarize_table -i ${OUT}/bowtie2_runs/${NAME}/${folder}/Sort_${hang}max_overhang/merged/max_overhang_${hang}_sum_merged.biom -o ${OUT}/bowtie2_runs/${NAME}/${folder}/Sort_${hang}max_overhang/merged/max_overhang_${hang}_sum_merged.biom.sum_table.txt
done

#make biom tables
for hang in ${ALLOVER}
do
list=""
mkdir -p ${OUT}/bowtie2_runs/${NAME}/${folder}/Sort_${hang}max_overhang/merged
 for perc in ${FINSIM}
 do
  echo "hi"
  biom convert -i ${OUT}/bowtie2_runs/${NAME}/${folder}/Sort_${hang}max_overhang/single_best/single_best${NAME}_${perc}_single_best.txt -o ${OUT}/bowtie2_runs/${NAME}/${folder}/Sort_${hang}max_overhang/single_best/single_best${NAME}_${perc}_single_best.biom --table-type="OTU table" --process-obs-metadata taxonomy
 done
 for f in `ls ${OUT}/bowtie2_runs/${NAME}/${folder}/Sort_${hang}max_overhang/single_best/*.biom`
 do
  if [ "x$f" != "x"  ] && [ "x$list" != "x"  ];
   then
    list="${f},${list}"
   elif [ "x$f" != "x" ] && [ "x$list" == "x"  ];
   then
    list="${f}"
   fi
 done
 echo "${list}"
 merge_otu_tables.py -i "${list}" -o ${OUT}/bowtie2_runs/${NAME}/${folder}/Sort_${hang}max_overhang/merged/max_overhang_${hang}_sum_merged.biom
 biom convert -i ${OUT}/bowtie2_runs/${NAME}/${folder}/Sort_${hang}max_overhang/merged/max_overhang_${hang}_sum_merged.biom -o ${OUT}/bowtie2_runs/${NAME}/${folder}/Sort_${hang}max_overhang/merged/max_overhang_${hang}_sum_merged.txt -b --header-key="taxonomy" --output-metadata-id="Consensus Lineage"
 biom summarize_table -i ${OUT}/bowtie2_runs/${NAME}/${folder}/Sort_${hang}max_overhang/merged/max_overhang_${hang}_sum_merged.biom -o ${OUT}/bowtie2_runs/${NAME}/${folder}/Sort_${hang}max_overhang/merged/max_overhang_${hang}_sum_merged.biom.sum_table.txt
done

# merge final set of biom tables

list=""
#make biom tables
mkdir -p ${OUT}/bowtie2_runs/${NAME}/${folder}/all_overhang_merged_results/
for hang in ${ALLOVER}
do
 for f in `ls ${OUT}/bowtie2_runs/${NAME}/${folder}/Sort_${hang}max_overhang/merged/max_overhang_${hang}_sum_merged.biom`
 do
   if [ "x$f" != "x"  ] && [ "x$list" != "x"  ];
   then
    list="${f},${list}"
   elif [ "x$f" != "x" ] && [ "x$list" == "x"  ];
   then
    list="${f}"
   fi
 done
echo "${list}"
done
merge_otu_tables.py -i "${list}" -o ${OUT}/bowtie2_runs/${NAME}/${folder}/all_overhang_merged_results/all_max_overhang_sum_merged.biom
biom convert -i ${OUT}/bowtie2_runs/${NAME}/${folder}/all_overhang_merged_results/all_max_overhang_sum_merged.biom -o ${OUT}/bowtie2_runs/${NAME}/${folder}/all_overhang_merged_results/all_max_overhang_sum_merged.txt -b --header-key="taxonomy" --output-metadata-id="Consensus Lineage"
biom summarize_table -i ${OUT}/bowtie2_runs/${NAME}/${folder}/all_overhang_merged_results/all_max_overhang_sum_merged.biom -o ${OUT}/bowtie2_runs/${NAME}/${folder}/all_overhang_merged_results/all_max_overhang_sum_merged.biom.sum_table.txt


done

#############################################
# process the paired reads
#############################################

###   no mixed!!!!!!!!!!

nextF=""
nextR=""
for hang in ${ALLOVER}
do 
 mkdir -p ${OUT}/bowtie2_runs/${NAME}/unassembled/Sort_${hang}max_overhang/
 mkdir -p ${OUT}/bowtie2_runs/${NAME}/unassembled/Sort_${hang}max_overhang/single_best 
 mkdir -p ${OUT}/bowtie2_runs/${NAME}/unassembled/Sort_${hang}max_overhang/merged
 for perc in ${FINSIM}
 do
  mkdir -p ${OUT}/bowtie2_runs/${NAME}/unassembled/Sort_${hang}max_overhang/bowtie2_at_${perc}
  mkdir -p ${OUT}/bowtie2_runs/${NAME}/unassembled/Sort_${hang}max_overhang/bowtie2_at_${perc}/fasta_to_process/
  if [ "x$nextF" != "x"  ];
  then
   mv $nextF ${OUT}/bowtie2_runs/${NAME}/unassembled/Sort_${hang}max_overhang/bowtie2_at_${perc}/fasta_to_process/${NAME}_pear_unassembled_F.fasta
   mv $nextR ${OUT}/bowtie2_runs/${NAME}/unassembled/Sort_${hang}max_overhang/bowtie2_at_${perc}/fasta_to_process/${NAME}_pear_unassembled_R.fasta
  fi
  echo "...bowtie2 running"
  bowtie2 -x ${DB}/${NAME}/${NAME}_bowtie2_databases/${NAME}_bowtie2_${perc}/${NAME}_${perc}_bowtie2_index -f -1 ${OUT}/bowtie2_runs/${NAME}/unassembled/Sort_${hang}max_overhang/bowtie2_at_${perc}/fasta_to_process/${NAME}_pear_unassembled_F.fasta  -2 ${OUT}/bowtie2_runs/${NAME}/unassembled/Sort_${hang}max_overhang/bowtie2_at_${perc}/fasta_to_process/${NAME}_pear_unassembled_R.fasta --un-conc ${OUT}/bowtie2_runs/${NAME}/unassembled/Sort_${hang}max_overhang/bowtie2_at_${perc}/fasta_to_process/${NAME}_all.clean_unassembled_bowtie2_rejects.fasta -S ${OUT}/bowtie2_runs/${NAME}/unassembled/Sort_${hang}max_overhang/bowtie2_at_${perc}/${NAME}_all.clean_unassembled_bowtie2.sam --no-hd --no-sq --very-sensitive-local --no-mixed --no-unal -p 120 #run Bowtie2 on the entire dataset usng 99% ref lib
  echo "...bowtie2 finished"
  date
  echo " "
  python ${DB}/scripts/group_alignments_to_files_p_mod.py ${OUT}/bowtie2_runs/${NAME}/unassembled/Sort_${hang}max_overhang/bowtie2_at_${perc}/ ${OUT}/bowtie2_runs/${NAME}/unassembled/Sort_${hang}max_overhang/bowtie2_at_${perc}/${NAME}_all.clean_unassembled_bowtie2.sam ${hang} ${perc} -p          # sort reads into -> single best hits, multiple hits, to short.  Also sort my percent similarity to the 99% reference
  echo "Summarizing single best hits for ${perc} with up to ${hang} basepair overhand "
   # get best hits for 99% at S25
  python ${DB}/scripts/summarize_bowtie2_hits_full_taxonomy.py ${DB}/${NAME}/${NAME}_final_database/${NAME}_taxonomy_${perc}.txt ${OUT}/bowtie2_runs/${NAME}/unassembled/Sort_${hang}max_overhang/bowtie2_at_${perc}/bowtie2_good_hits.txt ${OUT}/bowtie2_runs/${NAME}/unassembled/Sort_${hang}max_overhang/single_best/single_best${NAME}_${perc}_single_best.txt
  biom convert -i ${OUT}/bowtie2_runs/${NAME}/unassembled/Sort_${hang}max_overhang/single_best/single_best${NAME}_${perc}_single_best.txt -o ${OUT}/bowtie2_runs/${NAME}/unassembled/Sort_${hang}max_overhang/single_best/single_best${NAME}_${perc}_single_best.biom --table-type "otu table"
  cat ${OUT}/bowtie2_runs/${NAME}/unassembled/Sort_${hang}max_overhang/bowtie2_at_${perc}/*_forward.fasta > ${OUT}/bowtie2_runs/${NAME}/unassembled/Sort_${hang}max_overhang/bowtie2_at_${perc}/${NAME}_pear_unassembled_F.fasta
  rm ${OUT}/bowtie2_runs/${NAME}/unassembled/Sort_${hang}max_overhang/bowtie2_at_${perc}/*_forward.fasta
  nextF=${OUT}/bowtie2_runs/${NAME}/unassembled/Sort_${hang}max_overhang/bowtie2_at_${perc}/${NAME}_pear_unassembled_F.fasta
  cat ${OUT}/bowtie2_runs/${NAME}/unassembled/Sort_${hang}max_overhang/bowtie2_at_${perc}/*_reverse.fasta > ${OUT}/bowtie2_runs/${NAME}/unassembled/Sort_${hang}max_overhang/bowtie2_at_${perc}/${NAME}_pear_unassembled_R.fasta
  rm ${OUT}/bowtie2_runs/${NAME}/unassembled/Sort_${hang}max_overhang/bowtie2_at_${perc}/*_reverse.fasta
  nextR=${OUT}/bowtie2_runs/${NAME}/unassembled/Sort_${hang}max_overhang/bowtie2_at_${perc}/${NAME}_pear_unassembled_R.fasta
 done
done  



#merge biom tables for unasembled reads
for hang in ${ALLOVER}
do
list=""
mkdir -p ${OUT}/bowtie2_runs/${NAME}/unassembled/Sort_${hang}max_overhang/merged
 for f in `ls ${OUT}/bowtie2_runs/${NAME}/unassembled/Sort_${hang}max_overhang/single_best/*.biom`
 do
  if [ "x$f" != "x"  ] && [ "x$list" != "x"  ];
   then
    list="${f},${list}"
   elif [ "x$f" != "x" ] && [ "x$list" == "x"  ];
   then
    list="${f}"
   fi
 done
 echo "${list}"
 merge_otu_tables.py -i "${list}" -o ${OUT}/bowtie2_runs/${NAME}/unassembled/Sort_${hang}max_overhang/merged/max_overhang_${hang}_sum_merged.biom
 biom convert -i ${OUT}/bowtie2_runs/${NAME}/unassembled/Sort_${hang}max_overhang/merged/max_overhang_${hang}_sum_merged.biom -o ${OUT}/bowtie2_runs/${NAME}/unassembled/Sort_${hang}max_overhang/merged/max_overhang_${hang}_sum_merged.txt -b --header-key="taxonomy" --output-metadata-id="Consensus Lineage"
 biom summarize_table -i ${OUT}/bowtie2_runs/${NAME}/unassembled/Sort_${hang}max_overhang/merged/max_overhang_${hang}_sum_merged.biom -o ${OUT}/bowtie2_runs/${NAME}/unassembled/Sort_${hang}max_overhang/merged/max_overhang_${hang}_sum_merged.biom.sum_table.txt
done

# merge final set of assembled biom tables

list=""
#make biom tables
mkdir -p ${OUT}/bowtie2_runs/${NAME}/unassembled/all_overhang_merged_results/
for hang in ${ALLOVER}
do
 for f in `ls ${OUT}/bowtie2_runs/${NAME}/unassembled/Sort_${hang}max_overhang/merged/max_overhang_${hang}_sum_merged.biom`
 do
   if [ "x$f" != "x"  ] && [ "x$list" != "x"  ];
   then
    list="${f},${list}"
   elif [ "x$f" != "x" ] && [ "x$list" == "x"  ];
   then
    list="${f}"
   fi
 done
echo "${list}"
done
merge_otu_tables.py -i "${list}" -o ${OUT}/bowtie2_runs/${NAME}/unassembled/all_overhang_merged_results/all_max_overhang_sum_merged.biom
biom convert -i ${OUT}/bowtie2_runs/${NAME}/unassembled/all_overhang_merged_results/all_max_overhang_sum_merged.biom -o ${OUT}/bowtie2_runs/${NAME}/unassembled/all_overhang_merged_results/all_max_overhang_sum_merged.txt -b --header-key="taxonomy" --output-metadata-id="Consensus Lineage"
biom summarize_table -i ${OUT}/bowtie2_runs/${NAME}/unassembled/all_overhang_merged_results/all_max_overhang_sum_merged.biom -o ${OUT}/bowtie2_runs/${NAME}/unassembled/all_overhang_merged_results/all_max_overhang_sum_merged.biom.sum_table.txt
done


##### add the single and paired bowtie 2 files to different folders. Turn the following code into a for loop for the single bowtie2 reads

mkdir -p ${OUT}/bowtie2_runs/${NAME}/all_bowtie2_results_merged_all_overhang_all_percent_all_numbers_of_reads/
###merge all biom tables for all types of reads
type="unassembled discarded_R discarded_F assembled"
#make biom tables
for f in `ls ${OUT}/bowtie2_runs/${NAME}/*/all_overhang_merged_results/all_max_overhang_sum_merged.biom`
do
  if [ "x$f" != "x"  ] && [ "x$list" != "x"  ];
  then
   list="${f},${list}"
  elif [ "x$f" != "x" ] && [ "x$list" == "x"  ];
  then
   list="${f}"
  fi
echo "${list}"
done
merge_otu_tables.py -i "${list}" -o ${OUT}/bowtie2_runs/${NAME}/all_bowtie2_results_merged_all_overhang_all_percent_all_numbers_of_reads/all_bowtie2_results_merged_all_overhang_all_percent_all_numbers_of_reads.biom
biom convert -i ${OUT}/bowtie2_runs/${NAME}/all_bowtie2_results_merged_all_overhang_all_percent_all_numbers_of_reads/all_bowtie2_results_merged_all_overhang_all_percent_all_numbers_of_reads.biom -o ${OUT}/bowtie2_runs/${NAME}/all_bowtie2_results_merged_all_overhang_all_percent_all_numbers_of_reads/all_bowtie2_results_merged_all_overhang_all_percent_all_numbers_of_reads.txt -b --header-key="taxonomy" --output-metadata-id="Consensus Lineage"
biom summarize_table -i ${OUT}/bowtie2_runs/${NAME}/all_bowtie2_results_merged_all_overhang_all_percent_all_numbers_of_reads/all_bowtie2_results_merged_all_overhang_all_percent_all_numbers_of_reads.biom -o ${OUT}/bowtie2_runs/${NAME}/all_bowtie2_results_merged_all_overhang_all_percent_all_numbers_of_reads/all_bowtie2_results_merged_all_overhang_all_percent_all_numbers_of_reads.biom.sum_table.txt






































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































