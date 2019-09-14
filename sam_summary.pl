#!/usr/bin/perl -w

use strict;
use warnings;

##########################################################################
## This perl script will extract CIGAR and MD:Z from SAM file,         	##
## and then sort the OT hits to output file 1, and OT number           	##
## to out file 2.                                                      	##
## the off-target alignment meet one of following criteria             	##
##    (1) total unalignment < threshold                                	##
##    (2) total unalignment == threshold and no mismatch in seed region	##
##          							       	##
## seed region threshold x, means x bp next to PAM         		##
## 									##
## 									##
##   Kabin Xie, 2018.8.22						##
##########################################################################

die("Usage: sam_summary.pl seed_region unaln_threshold unaln_seed_threshold <in.sam> <out1_offtarget.tab> <out2_OT_no.tab> \n") if (@ARGV != 6);

my ($line, $seed_region);
my ($gRNA_id, $target_id, $cigar, $asi, $xm, $xo, $xg, $mdz);
my ($pre_gRNA_id, $min_unalign_no, $min_seed_unalign_no);
my ($unaln_threshold, $unaln_seed_threshold);
my ($ot_no);

open(IN, $ARGV[3]) 	   || die("Can't open the input sam file\n");
open(OUT1_OT, ">$ARGV[4]") || die("Can't open the output file 1\n");
open(OUT2_NO, ">$ARGV[5]") || die("Can't open the output file 1\n");

$seed_region          = $ARGV[0];
$unaln_threshold      = $ARGV[1];
$unaln_seed_threshold = $ARGV[2];

if($seed_region>20)           { die("err, the seed region is out of range \n");}
if($unaln_threshold>20)       { die("err, the unaln_threshold is out of range\n");}
if($unaln_seed_threshold >20) { die("err, the unaln_seed_threshold is out of range\n");}

$pre_gRNA_id 	     = "first";
$min_unalign_no      = 100;
$min_seed_unalign_no = 100;
$ot_no 		     = 0;
		
while($line=<IN>) {
	my @splited_record = split(/\t/, $line);
	my $unaligment;
	my $total_gap_no = 0;	

	$gRNA_id   = $splited_record[0];
	$target_id = $splited_record[2];
	

#Start process a new gRNA_id	
	if($pre_gRNA_id !~ $gRNA_id)  { 
		if($pre_gRNA_id =~ "first") { 	#print table header
		 	print OUT1_OT join("\t", ("gRNA_id", "Off_Target", "Total_unaln", "seed_unalign", "cigar", "nmd")), "\n";
		} else {			
			print OUT2_NO join("\t", ($pre_gRNA_id, $ot_no)), "\n";	
					
		}
		$pre_gRNA_id         = $gRNA_id;
		$min_unalign_no      = 100;
		$min_seed_unalign_no = 100;
		$ot_no               = 0;
	}

	$cigar = $splited_record[5];
	$xm    = $splited_record[13] =~ s/XM:i://r;  	# get the XM value
	$xo    = $splited_record[14] =~ s/XO:i://r;
	$xg    = $splited_record[15] =~ s/XG:i://r;
	$mdz   = $splited_record[17] =~ s/MD:Z://r;
	$asi   = $splited_record[11] =~ s/AS:i://r;
	
#Processing the CIGAR string
	my @letter= split(/[0-9]+/,$cigar);        #get the symbols
	my @num   = split(/[MDI]/, $cigar);        #get the numbers
	my $length= scalar@num;	
	splice(@letter,0,1); 
	
	if(scalar @letter != scalar@num) { 
		print $gRNA_id, " ", $target_id, " has error in cigar.\n";
		next;
	}
	
	# skip alignment with gaps next to PAM;
	if($letter[-1] =~ /[DI]/) { print $gRNA_id, " and ", $target_id, " ", $cigar, " have gaps at 3'end\n"; next; } 	
		
 	#Get the gaps in seed region
	my $seed_gap_no = 0;    	#number of gaps in seed region
	my $pos = 0;			#sequence coordinate number the N closest to PAM set to 1.
	#print $length;

	# start from the end of guide sequence, aka nucleotide closest to PAM
	for(my $i=$length-1;$i>=0;$i--) {	
		$pos=$pos+$num[$i];
		
		if($letter[$i]=~/D/ or $letter[$i]=~/I/) {	#gap region	
			$total_gap_no = $total_gap_no + $num[$i];
			if($pos<=$seed_region) {
				$seed_gap_no=$seed_gap_no+$num[$i]; #within seed region
			}	
 			else {
				if($pos-$num[$i]<=$seed_region) {	#cross the seed boundary					
					$seed_gap_no=$seed_gap_no+12-($pos-$num[$i]);
				}
				else {
					next;
				}
			} 
		}
		
	}

#processing the NM:D region
#get the mismatch (substitutions) in seed region
	my @mdz_snv   = split(/[0-9]+/,$mdz);
	my @mdz_match = split(/[A-Z]+|\^[A-Z]+/i, $mdz);
	splice(@mdz_snv,0,1);

	$pos   	       = 0;
	my $seed_mm_no = 0;
	#in case to skip hit has mis-match next to PAM, remove the # in next line
	#if($mdz_match[-1] == 0 ) {print $gRNA_id, " and ", $target_id, " ", $mdz, " have mismatches at 3' end\n"; next;  } 
	
	while (@mdz_match) {
		my $match_wid = pop @mdz_match; 
		my $snv = pop @mdz_snv;
		my $snv_wid = length($snv);
		
		$pos = $pos + $match_wid;
		#print $pos, "\t";
		if($snv) { 
			$pos=$pos+$snv_wid; 
			if($snv =~ /\^[A-Z]/i) { 
				$pos=$pos -1 ; # do not count symbol ^ in mdz
				next; 
			} else {
				if($pos<=$seed_region) { $seed_mm_no = $seed_mm_no+1;}
			}
		}
		
		
	}
	#print "\n";

	#if all record are needed, use the following output option. The min value could easily calculated using R
	#print join("\t", ($gRNA_id, $target_id, $unaligment, $seed_gap_no, $seed_mm_no, $cigar, $mdz)), "\n";
	
#Get the min_unalignment and min_seed_un value for each input guide
	$unaligment = $total_gap_no + $xm;
	print OUT1_OT join("\t", ($gRNA_id, $target_id, $unaligment, $seed_gap_no+$seed_mm_no, $cigar, $mdz)), "\n";
        #the target_id#No., while No. is the frequency of this sequence in RDP-16S rDNA 
	my @ot_no_this = split (/\#/, $target_id);
	if($unaligment < $unaln_threshold || 
	  ($unaligment ==$unaln_threshold && ($seed_gap_no + $seed_mm_no)< $unaln_seed_threshold)) {
		$ot_no = $ot_no+$ot_no_this[-1];
	}
	
#print $gRNA_id,"\t", $target_id, "\t", $unaligment, "\t",$asi,"\t", $seed_gap_no+$seed_mm_no, "\t", $cigar, "\t", $mdz, "\n";

}
	
print OUT2_NO join("\t", ($pre_gRNA_id, $ot_no)), "\n"; #print the last record	
close IN;
close OUT1_OT;
close OUT2_NO;
exit;

