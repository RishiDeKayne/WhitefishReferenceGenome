#! /usr/bin/perl -w

# Stefan Zoller, Genetic Diversity Centre, ETH Zurich, 2019

use File::Basename;

if (@ARGV < 2) {
        print "\nusage:   cov.per.window.pl  <bedtools.coverage.file>  <window-size> \n\n";
	print "examples: cov.per.window.pl file.bed 20  \n";
	print "          cov.per.windwo.pl file.bed.gz 20    # with a gzipped input file \n";
	print "          cov.per.window.pl file.bed.gz 20 | gzip > cov20.txt.gz    # this will create a gzipped output file\n\n"; 
	exit 1;
}

my $winsize = ($ARGV[1]);
my @tcov;

my($filename, $directories, $suffix) = fileparse($ARGV[0]);

if ( $filename =~ m/\.gz$/){
		open (FILE1, "gunzip -c $directories/$filename |"); 
} else {
        	open (FILE1, "$directories/$filename") or die "Cannot open $filename";
}

my $prevcontig="_na_";
my $start=0;
my $end=0;
my $contig="";
my $pos=0;
my $cov=0;
#my $type="tile";
my $type="sliding";

if( $type eq "tile") {
# this "tile" approach needs to be updated -> compare with "sliding"
   while (($line1 = <FILE1>)) {
	chomp $line1;
	#$lnum++;
	($contig, $pos, $cov)=split("\t", $line1);
	print "$contig, $pos, $cov\n";

	if($prevcontig eq "_na_") {
		$end=$pos;
		$prevcontig = $contig;
		#print "## $contig \n";
	}
	if($pos < int($winsize / 2)){
		$end=$pos;
		$sum=$sum+$cov;
		$nvals++;
		next;
	}
	if($prevcontig ne $contig) {
		printCov($prevcontig);
		$prevcontig = $contig;
		$end=$pos;
		$sum=0;
		$nvals=0;
		@tcov=();
		exit 1; 
	}
	$sum=$sum + $cov; 
	$nvals++;
	$end=$pos;
	#push (@tcov, $cov);
	if($nvals == $winsize){
		printCov($contig);
		$sum=0;
		$nvals=0;
		@tcov=();
	}
   }
}

if($type eq "sliding"){
   while (($line1 = <FILE1>)) {
	# should save the last "window-size" of values and stop printing when reaching end of contig.
	chomp $line1;
	($contig, $pos, $cov)=split("\t", $line1);
	#print "$contig, $pos, $cov\n";
	$nvals++;
	if($prevcontig eq "_na_") {
		$end=$pos;
		$prevcontig = $contig;
		$nvals=1;
		$sum=$cov;
		@tcov=();
		push (@tcov, $cov);
		next;
	}
	if($prevcontig ne $contig) {
		#printCov($prevcontig);
		$prevcontig = $contig;
		$end=$pos;
		$sum=$cov;
		$nvals=1;
		@tcov=();
		push (@tcov, $cov);
		#printCov($contig);
		next;
	}
	if($pos < int($winsize / 2)){
		$end=$pos;
		$sum=$sum+$cov;
		$nvals++;
		push (@tcov, $cov);
		#printCov($contig);
		next;
	}if($pos >= (int($winsize / 2))){
		$end=$pos;
		#print ">> $pos $winsize \n";
		$nvals--;
		printCov($contig);
		$rm = shift @tcov;
		$sum=$sum - $rm;
		$sum=$sum + $cov;
		push (@tcov, $cov);
		$end=$pos;
	}
   }
}	




sub printCov {
	my $cont = shift;
	my $mean = $sum / $nvals;
	my $mid = int ($end - ($winsize)/2);
	#print "   $cont\tpos: $pos\tmean_cov: $mean\tmid: $mid\tnvals: $nvals\twinsize: $winsize\n";
	print "$cont\t$pos\t$mean\n";
}


