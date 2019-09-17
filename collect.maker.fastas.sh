#!/bin/bash
# Stefan Zoller, Genetic Diversity Centre, ETH Zurich, 2019

module load gcc/4.8.2 gdc open_mpi/1.6.5 perl/5.18.4 snap/2006-07-28 blast/2.3.0 exonerate/2.4.0 genemark_es/4 hmmer/3.1 rmblast/2.2.28 repeatmasker/4.0.6 boost/1.55.0 sqlite3/3.11 gsl/1.16 bamtools/2.4.0 suitesparse/4.5.1 openblas/0.2.13_seq augustus/3.2.1 maker/2.31.9

for index in `ls -1 /cluster/scratch/rdekayne/Maker_wtdbg2_fullprot/threaded_r3/maker.*/*.maker.output/*index.log`
do 
	echo -n "."
	#ls -1 $index
	fasta_merge -d $index 
done
echo

rm -rf all.maker.proteins.fasta
for prot in *maker.proteins.fasta
do 
	cat $prot >> all.maker.proteins.fasta
done

rm -rf all.maker.transcripts.fasta
for trans in *maker.transcripts.fasta
do
	cat $trans >> all.maker.transcripts.fasta
done




