#!/usr/bin/perl
use strict; use warnings;

#create new unstd iHS file with local rho estimate for stdization with DAF & rho

if (scalar(@ARGV) != 2) {
	die "input: perl makeInputFor.iHS.std.by.freq.rho.bins.pl POP chr#\n";
}
my $c = pop(@ARGV);
my $pop = pop(@ARGV);

#filter out SNPs with rho dist of 0 on both sides of one or other allele
#then create files for restandardization including SNPs with numeric iHS scores, and local rho estimates

open (my $out, ">", "${pop}chr${c}.restd.w.rho.txt");
print $out "SNP\tPOS\tFREQ_1\tLOC_RHO\tunstd_iHS\tdensity\tmax_gap\tngap\tL-EHH0\tL-EHH1\tR-EHH0\tR-EHH1\tL-EHH0_RHO\tR-EHH0_RHO\tL-EHH1_RHO\tR-EHH1_RHO\tWarnings\n";

#load local rho estimates into hash
my %rates;
open (my $rhoEst, "<", "${pop}chr${c}.local.rho.est.txt") or die "$pop $c local rho est\n";
while (my $line = <$rhoEst>) {
	chomp($line);
	my @stuff = split(/\s+/, $line);
	#skip header
	if ($line !~ m/^chr/) {
		$rates{$stuff[1]} = $stuff[4];
	}
}
close $rhoEst;

#go through ihs unstdized file
open (my $ihs, "<", "${pop}chr${c}.ihs") or die "$pop $c unstd ihs\n";
while (my $line = <$ihs>) {
	chomp($line);
	my @stuff = split(/\s+/, $line);
	#if iHS numeric
	if ($stuff[3] =~ m/^-?[0-9]+\.?[0-9]+$/) {
		#if SNP has local recomb rate
		if (exists $rates{$stuff[1]}) {
			#if both rho dists for each allele are not 0 (pos 11/12 13/14)
			if ((($stuff[11] != 0) || ($stuff[12] != 0)) && (($stuff[13] != 0) || ($stuff[14] != 0))) {
				#print out info
				print $out "$stuff[0]\t$stuff[1]\t$stuff[2]\t$rates{$stuff[1]}\t";
				for (my $i=3; $i<15; $i++) {
					print $out "$stuff[$i]\t";
				}
				print $out "$stuff[15]\n";
			}
		}
	}
}
close $ihs;
close $out;