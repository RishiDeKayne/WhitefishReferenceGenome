#!/bin/bash
# Stefan Zoller, Genetic Diversity Centre, ETH Zurich, 2019

module load gcc/4.8.2 gdc snap/2006-07-28 blast/2.3.0 exonerate/2.4.0 genemark_es/4 hmmer/3.1 rmblast/2.6.0 repeatmasker/4.0.7 boost/1.55.0 sqlite3/3.11 gsl/1.16 bamtools/2.4.0
module load suitesparse/4.5.1 openblas/0.2.13_seq augustus/3.2.1 maker/2.31.9
module load perl/5.18.4  
source /cluster/apps/gdc/perl5/etc/bashrc

export OMP_NUM_THREADS=4

export AUGUSTUS_CONFIG_PATH="./local_augustus/config"

bsub -n 2 -W 120:00 -J "Maker.threaded" -R "rusage[mem=4000]" "maker -cpus 4 > maker.threaded.log 2> maker.threaded.err "


