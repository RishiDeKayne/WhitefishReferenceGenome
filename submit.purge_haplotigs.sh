#!/bin/bash

#BSUB -J "purge_haplotigs"
#BSUB -W 120:00 
#BSUB -n 16
#BSUB -R "rusage[mem=12000]"

module load gcc/4.8.2 gdc perl/5.18.4 samtools/1.8 bedtools/2.27.1 blast/2.3.0 python/2.7.11 r/3.2.2
source /cluster/apps/gdc/perl5/etc/bashrc

OMP_NUM_THREADS=16

export PATH="/cluster/project/gdc/shared/tools/purge_haplotigs/bin/:$PATH"
export PATH="/cluster/project/gdc/shared/tools/lastz1.0.4/bin:$PATH"    # 


#############################################################################################################
##
## Run each step indivdually. For example, uncomment step 1 (remove the # in front of "purge_haplotigs readhist ...") and submit it. 
## When done, comment step 1 (= add # in front), uncomment step 2, change options and submit 
## Step 2 needs options adjusted according to results from step 1, see also workflow file. 
##

### step 1, might need run-time > 1h -> -W 4:00
# purge_haplotigs  readhist  /cluster/scratch/rdekayne/PurgeHaplos_test/Minimap2/pb_minimapped_to_wtdbg2.run2.arrow.iter1.pilon.iter1.fasta.sorted.bam 

### step 2, runs very fast  -> change run-time  -W to 1:00
#  purge_haplotigs  contigcov  -i  /cluster/scratch/rdekayne/PurgeHaplos_test/Minimap2/wtdbg2_run2/pb_minimapped_to_wtdbg2.run2.arrow.iter1.pilon.iter1.fasta.sorted.bam.genecov   -l 10  -m 55  -h 120  

### step 3, might need a lot of time and memory, depending on genome size and mapped data
###  -> change -W to 72:00 and mem=45000  (maybe even higher)
# purge_haplotigs  purge -t 16  -g /cluster/project/gdc/shared/p298/Rishi/PrePurgeHaplotigGenomes/wtdbg2.run2.arrow.iter1.pilon.iter1.fasta -c /cluster/scratch/rdekayne/PurgeHaplos_test/Minimap2/wtdbg2_run2/coverage_stats.csv -b /cluster/scratch/rdekayne/PurgeHaplos_test/Minimap2/pb_minimapped_to_wtdbg2.run2.arrow.iter1.pilon.iter1.fasta.sorted.bam