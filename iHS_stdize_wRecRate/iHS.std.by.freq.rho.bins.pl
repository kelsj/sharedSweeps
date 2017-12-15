#!/usr/bin/perl
use strict; use warnings;

#standardize iHS scores by both rho and DAF bins

#input: pop, chrom#

#IMPT: autosomes have 4 rho bins, chrX has 3

if (scalar(@ARGV) != 2) {
	die "input: perl iHS.std.by.freq.rho.bins.pl POP chr#\n";
}
my $pop = $ARGV[0];
my $c = $ARGV[1];

#2D arrays mean & SD values
my @autoMean;
my @autoSD;
#get rho cutoffs used & store from stdization bin file
my @autoCutoffs;
open (my $autoStats, "<", "${pop}.new.std.bin.info.txt") or die "pop autosomes bins file\n";
#keep track of freq bin and last freq value
my $row = 0;
my $col = 0;
my $last = 0;
while (my $line = <$autoStats>) {
	chomp($line);
	my @stuff = split(/\s+/, $line);
	#skip header
	if ($line !~ m/^FREQ0/) {
		#first row
		if (($stuff[0]==0) && ($stuff[2]==0)) {
			$autoMean[$row][$col] = $stuff[5];
			$autoSD[$row][$col] = $stuff[6];
		}
		#new bin, advance row and reset columns & last
		elsif ($stuff[0] != $last) {
			$row++;
			$col = 0;
			$last = $stuff[0];
			$autoMean[$row][$col] = $stuff[5];
			$autoSD[$row][$col] = $stuff[6];			
		}
		#not new bin, advance column but not row
		else {
			$col++;
			$autoMean[$row][$col] = $stuff[5];
			$autoSD[$row][$col] = $stuff[6];
		}
		#if first line, store cutoffs for rho
		if ($stuff[0] == 0) {
			my @temp = ($stuff[2], $stuff[3]);
			push @autoCutoffs, [ @temp ];
		}
	}
}
close $autoStats;
#cutoffs format: D1 = pairs of thresh, D2=low,high
#mean/SD array format: D1=freq lower thresh, D2=rho bin


#for each SNP, calculate row# from freq0, find rho bin, then std with mean & var
open (my $out, ">", "${pop}chr${c}.new.std.ihs");
print $out "SNP\tPOS\tFREQ_1\tstd_iHS\tdensity max_gap\tngap\tL-EHH0\tL-EHH1\tR-EHH0\tR-EHH1\tL-EHH0_RHO\tR-EHH0_RHO\tL-EHH1_RHO\tR-EHH1_RHO\tWarnings\tLOC_RHO\n";
open (my $unstd, "<", "${pop}chr${c}.restd.w.rho.txt") or die "$pop chr$c new unstd\n";
while (my $line = <$unstd>) {
	chomp($line);
	my @stuff = split(/\s+/, $line);
	#skip header
	if ($line !~ m/^SNP/) {
		#calculate freq row pos
		my $rowNum = int($stuff[2]/0.02);
		#print "row $rowNum, freq $stuff[2], rho $stuff[3]\n";
		my $binMean;
		my $binSD;
		#get column for rho, then mean/sd
		my $colNum;
		ColLoop: for (my $t=0; $t<scalar(@autoCutoffs); $t++) {
			#if falls in this rho bin
			if (($stuff[3] >= $autoCutoffs[$t][0]) && ($stuff[3] < $autoCutoffs[$t][1])) {
				$colNum = $t;
				last ColLoop;
			}
		}
		$binMean = $autoMean[$rowNum][$colNum];
		$binSD = $autoSD[$rowNum][$colNum];
		#print "col $colNum, mean $binMean, sd $binSD\n";
		my $stdiHS = ($stuff[4]-$binMean)/$binSD;
		my $roundedStdiHS = sprintf("%.3f", $stdiHS);
		#print "$stuff[4] std: $roundedStdiHS\n";
		#print out in same order as original std ihs files, with local rho @ end to not change field positions
		#save rho & remove from array
		my $rho = $stuff[3];
		splice @stuff, 3, 1;
		#replace unstd with std iHS (NOW 3RD ARRAY ELEMENT)
		$stuff[3] = $roundedStdiHS;
		#print "@stuff\n";
		#print out
		for (my $i=0; $i<scalar(@stuff); $i++) {
			print $out "$stuff[$i]\t";
		}
		print $out "$rho\n";
	}
}
close $out;
close $unstd;
