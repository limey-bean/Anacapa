Run Examples for running the Anacapa QC and Anacapa Classifier 

1. Run the Anacapa QC example
This command runs the Anacapa QC pipeline with the included 12S_test_data. Edit the ANACAPA variable to point to your extracted anacapa folder.

$ export ANACAPA="/home/vagrant/anacapa"
$ /bin/bash $ANACAPA/Anacapa_db/anacapa_QC_dada2.sh -i $ANACAPA/Example_data/12S_example_anacapa_QC_dada2_and_BLCA_classifier/12S_test_data -o $ANACAPA/Example_data/12S_example_anacapa_QC_dada2_and_BLCA_classifier/12S_time_test -d $ANACAPA/Anacapa_db -f $ANACAPA/12S_test_data/forward.txt -r $ANACAPA/12S_test_data/reverse.txt -e $ANACAPA/Anacapa_db/metabarcode_loci_min_merge_length.txt -a nextera -t MiSeq -l
The expected results can be found in anacapa/Anacapa_test_data_expected_output_after_QC_dada2

Approximate time to run:

real	0m45.906s
user	0m43.568s
sys	0m1.396s


2. Run the Anacapa Classifier example
This command runs the Anacapa Classifier pipeline on the output of the QC pipeline. Edit the ANACAPA variable to point to your extracted anacapa folder.

$ export ANACAPA="/home/vagrant/anacapa"
$ /bin/bash $ANACAPA/Anacapa_db/anacapa_classifier.sh -o $ANACAPA/Example_data/12S_example_anacapa_QC_dada2_and_BLCA_classifier/12S_time_test -d $ANACAPA/Anacapa_db  -l
The expected results can be found in anacapa/Anacapa_test_data_expected_output_after_classifier

Approximate time to run:

real	0m19.467s
user	0m13.384s
sys	0m1.480s
If using slurm or qsub an example job file is available in the jobs/ folder of this repository.


For examples run in the containerized version of Anacapa please see https://github.com/datproject/anacapa-container