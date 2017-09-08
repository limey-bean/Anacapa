#!/bin/bash

if [ $# -ne 4 ]; then
 echo "need 4 paras: <metabarcode_target> <out_dir> <ref_dir> <uclust_percent>";
 exit;

fi
####################################script & software
# need to fix the order at some point
script=./scripts
cutadapt=/u/local/apps/python/2.7.13/bin/cutadapt

##load module
source /u/local/Modules/default/init/bash

module load qiime


	echo "Pick_open_reference: mapping merged-reads to refs"
	echo "Up next $1"
	echo "assembled started"
	pick_open_reference_otus.py -i $2/primer_split/assembled/$1_F.fasta -o $2/Qiime_open_ref/$1/$1_assembled_open_ref_$4 -r $3/$1/$1_qiime.fasta -p $2/Qiime_open_ref/Params/Params$1_pick_open.txt --min_otu_size=1 -m uclust -s $4 --prefilter_percent_id=0.0 --suppress_align_and_tree  -a 
	echo "assembled completed"
	echo "unassembled started"
	pick_open_reference_otus.py -i $2/primer_split/unassembled_forward/$1_F.fasta,$2/primer_split/unassembled_reverse/$1_R.fasta -o $2/Qiime_open_ref/$1/$1_unassembled_open_ref_$4 -r $3/$1/$1_qiime.fasta -p $2/Qiime_open_ref/Params/Params$1_pick_open.txt --min_otu_size=1 -m uclust -s $4 --prefilter_percent_id=0.0 --suppress_align_and_tree -a
	date
	echo "unassembled completed"

	echo "discarded started"
	cat $2/primer_split/discarded/$1_F.fasta $2/primer_split/discarded/$1_R.fasta >> $2/primer_split/discarded/$1.fasta
	pick_open_reference_otus.py -i $2/primer_split/discarded/$1.fasta -o $2/Qiime_open_ref/$1/$1_discarded_open_ref_$4 -r $3/$1/$1_qiime.fasta -p $2/Qiime_open_ref/Params/Params$1_pick_open.txt --min_otu_size=1 -m uclust -s $4 --prefilter_percent_id=0.0 --suppress_align_and_tree -a
	date

	echo "discarded completed"

	echo "Time to summarize that we found"

	summarize_taxa.py  -i $2/Qiime_open_ref/$1/$1_assembled_open_ref_$4/otu_table_mc1_w_tax.biom -o $2/Qiime_open_ref/$1/$1_assembled_taxon_sum -L 3,4,5,6,7 --absolute_abundance
	summarize_taxa.py  -i $2/Qiime_open_ref/$1/$1_unassembled_open_ref_$4/otu_table_mc1_w_tax.biom -o $2/Qiime_open_ref/$1/$1_unassembled_taxon_sum -L 3,4,5,6,7  --absolute_abundance
	summarize_taxa.py  -i $2/Qiime_open_ref/$1/$1_discarded_open_ref_$4/otu_table_mc1_w_tax.biom -o $2/Qiime_open_ref/$1/$1_discarded_taxon_sum -L 3,4,5,6,7  --absolute_abundance


	cp $2/Qiime_open_ref/$1/$1_assembled_taxon_sum/otu_table_mc1_w_tax_L7.txt $2/taxon_summaries/$1/$1_assembled_taxon_sum.txt
	cp $2/Qiime_open_ref/$1/$1_unassembled_taxon_sum/otu_table_mc1_w_tax_L7.txt $2/taxon_summaries/$1/$1_unassembled_taxon_sum.txt
	cp $2/Qiime_open_ref/$1/$1_discarded_taxon_sum/otu_table_mc1_w_tax_L7.txt $2/taxon_summaries/$1/$1_discarded_taxon_sum.txt

	biom summarize_table -i $2/Qiime_open_ref/$1/$1_assembled_taxon_sum/otu_table_mc1_w_tax_L7.biom -o $2/taxon_summaries/$1/$1_assembled_taxon_biom_sum_table.txt
	biom summarize_table -i $2/Qiime_open_ref/$1/$1_unassembled_taxon_sum/otu_table_mc1_w_tax_L7.biom -o $2/taxon_summaries/$1/$1_unassembled_taxon_biom_sum_table.txt
	biom summarize_table -i $2/Qiime_open_ref/$1/$1_discarded_taxon_sum/otu_table_mc1_w_tax_L7.biom -o $2/taxon_summaries/$1/$1_discarded_taxon_biom_sum_table.txt


	mkdir $2/taxon_summaries/$1/merged
	merge_otu_tables.py -i $2/Qiime_open_ref/$1/$1_assembled_taxon_sum/otu_table_mc1_w_tax_L7.biom,$2/Qiime_open_ref/$1/$1_discarded_taxon_sum/otu_table_mc1_w_tax_L7.biom,$2/Qiime_open_ref/$1/$1_unassembled_taxon_sum/otu_table_mc1_w_tax_L7.biom  -o $2/taxon_summaries/$1/merged/$1_sum_merged.biom
	biom convert -i $2/taxon_summaries/$1/merged/$1_sum_merged.biom -o $2/taxon_summaries/$1/merged/$1_sum_merged.txt -b --header-key="taxonomy" --output-metadata-id="Consensus Lineage"

	biom summarize_table -i $2/taxon_summaries/$1/merged/$1_sum_merged.biom -o $2/taxon_summaries/$1/merged/$1_sum_merged.biom.sum_table.txt


	echo "Hope it looks good"

	

