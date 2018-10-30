#####################################################################
# Readme for the Edgar_anacapa_BLCA_classifier_test_config.sh script
#####################################################################

This script was designed to run the CVI analysis written by Edgar 2018 using the Anacapa Bowtie2-BLCA classifier.
  * Edgar, R.C., 2018. Accuracy of taxonomy prediction for 16S rRNA and fungal ITS sequences. PeerJ, 6, p.e4652.
  * https://drive5.com/taxxi/doc/index.html

The majority of scripts are directly copied from TAXXI
  * https://drive5.com/taxxi/doc/py_index.html

A few have been modified and or written explicitly for anacapa
for example:
  * run_namecount.py
  * rename_anacapa_blca_to_blast_BLCA.py
  * blca_from_bowtie.py
  * bcc_cuttoff.py
  * convert_to_BLCA_format.py
  * blca2tab_varaible_BCC.py

To modify this script to run on other platforms (not UCLA's Hoffman2 Cluster), change the ~/scripts/Edgar_anacapa_BLCA_classifier_test_config.sh file to match your computing environment.

You will also need a directory for the test and training datasets found at the TAXXI link: https://drive5.com/taxxi/doc/fasta_index.html
  * alternatively you can download the file in the github ~/Edgar_test_training_sets
  * The test set has the suffix *_test.fasta
  * The training set has the suffix *_all.fasta

example usage:
/bin/bash Edgar_anacapa_BLCA_classifier_test_config.sh -m ~/anacapa/bowtie2_blca_classifier_test -n ncbi_16s_100 -p 100 -t ~/anacapa/TAXXI_test_train_datasets
