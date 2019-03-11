#!/bin/bash

 /bin/bash /home/vagrant/Anacapa-git/Anacapa_db/scripts/run_bowtie2_blca.sh -o /home/vagrant/12S_time_test -d /home/vagrant/Anacapa-git/Anacapa_db -m 12S -l -b 0.8 -p 0.8 -c /home/vagrant/Anacapa-git/Anacapa_db/scripts/BCC_default_cut_off.sh -n 100 -x 1 -f -2.5 -g -2 -k /home/vagrant/Anacapa-git/Anacapa_db/scripts/Hoffman2_HPC_header.sh

echo _END_ [run_bowtie2_blca.sh]