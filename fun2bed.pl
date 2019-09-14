#!/usr/bin/perl -w

use strict;
use warnings;

###########################################################
## This program process the fuzznuc search result of PAM,
## generate the protospace positions and convert  to Bed 
## file of targeting region (protospacer).
##
##   Kabin Xie, 2018.8
###########################################################

die("Usage: divide.pl <in.fuzznuc> <out.bed> <sense==1/anti==0>\n") if (@ARGV != 3);

my ($chr, $line, $sense);
open(IN, $ARGV[0]) || die;
open(OUT, ">$ARGV[1]") || die;
$sense=$ARGV[2];

while($line=<IN>) {
	if ($line =~ m/^#/) {
		if($line =~ /Sequence:/) {	# retrieve chromosome ID (eq. id of each sequence)
		#$chr = substr($line, 12,10);      # specific for microbiome data
		my @splited_record = split(/\t|\s+/,$line);
		$chr = $splited_record[2];
		$chr =~ s/^\s+//;	#remove space (if existed!)
		$chr =~ s/\s+$//;
		}
		#print OUT $line;
	}
	
	else {
		#chomp $line;
		if ( ($line !~ /[a-z]/) || $line =~ /Start/) {next;} # skip empty line and header
		else{
			#my $t=substr($line, 0,8); # Get the start position of PAM
			my @record=split(/\t|\s+/, $line);
			my $t=$record[1];
	
			#print $record[1], $t, "\n";
			if(($sense and $t<23) or (not($sense) and $t<2)) {next;}  # make sure no minus coordinates
			else {
			 if($sense) {
			  print OUT $chr, "\t", $t-22, "\t", $t-2, "\t0\t0\t+", "\n";}  #sense 20bp PS sequence
			 else {
			 print OUT $chr, "\t", $t+2, "\t", $t+22, "\t0\t0\t-", "\n";} #antisense
			 }
			} 
		}
	}
	
close IN;
close OUT;
exit;

