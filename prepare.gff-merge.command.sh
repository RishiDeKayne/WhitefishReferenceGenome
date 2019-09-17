#!/bin/bash
# Stefan Zoller, Genetic Diversity Centre, ETH Zurich, 2019

ls -1 maker.*/*/*_datastore/*/*/*/*.gff > list.of.gffs.txt 

rm gff_merge.command.sh
echo '#!/bin/bash' > gff_merge.command.sh
echo  "module load gcc/4.8.2 gdc snap/2006-07-28 blast/2.3.0 exonerate/2.4.0 genemark_es/4 hmmer/3.1 rmblast/2.6.0 repeatmasker/4.0.7 boost/1.55.0 sqlite3/3.11 gsl/1.16 bamtools/2.4.0" >> gff_merge.command.sh
echo  "module load suitesparse/4.5.1 openblas/0.2.13_seq augustus/3.2.1 maker/2.31.9" >> gff_merge.command.sh
echo  "module load perl/5.18.4  " >> gff_merge.command.sh
echo  "source /cluster/apps/gdc/perl5/etc/bashrc" >> gff_merge.command.sh
echo  "" >> gff_merge.command.sh
echo  -n "gff3_merge -o allr3.gff " >> gff_merge.command.sh

for g in `cat list.of.gffs.txt `; do echo $g |tr "\n" " " ; done >> gff_merge.command.sh

chmod 755 gff_merge.command.sh


