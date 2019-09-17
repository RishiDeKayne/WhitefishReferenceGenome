#!/bin/bash
# Stefan Zoller, Genetic Diversity Centre, ETH Zurich, 2019

for f in `cat scaffold.files.txt`
do 
	fbase=`basename $f .fasta`
	fbase=`basename $fbase .scaffold`
	i=$((i + 1))
	mkdir maker.$fbase
	cd maker.$fbase
	cp ../$f .
	cp ../maker_*.ctl .
	cp ../maker.run.cmd.sh . 
	export f=$f; perl -p -i -e 's/^genome=/genome=$ENV{f}/' maker_opts.ctl   
	./maker.run.cmd.sh 
	cd ..
done 

