#!/usr/bin/perl


use FileHandle;
use strict;
use Getopt::Long;


my $line;
my @array;
my $start_pos;
my $prev_gt;
my $counter;
my $start;
my $prev_chrom;
my $prev_pos;

my $args = scalar @ARGV;

my %options = (
    size	=> undef,
    infile 	=> undef,
    outfile	=> undef
    );
GetOptions(
    's=s'         => \$options{size},
    'i=s'	=> \$options{infile},
    'o=s'	=> \$options{outfile}
    );


my $in_file=$options{infile};
my $in_file_handle = new FileHandle;

if(defined($in_file)) {
	
	# check defined input file exists
	die("ERROR: Could not find input file ", $in_file, "\n") unless -e $in_file;
	
	$in_file_handle->open($in_file) or die("ERROR: Could not read from input file ", $in_file, "\n");
}

# no file specified - try to read data off command line
else {
	$in_file_handle = 'STDIN';
}



my $out_file=$options{infile}.".".$options{size}.".txt";   
#my $out_file="out.txt";
open(FILE, "s_4_1.recal.snps.all_conf_sites.gt.pass.txt");
my $out_file_handle = new FileHandle;
$out_file_handle->open(">$out_file") or die("ERROR: Could not write to output file ", $out_file, "\n");
while(<FILE>) {
	$line=$_;
	@array=split("\t",$line);
	my $gt=$array[2];
	my $geno = chomp($gt);
	#foreach $line ($lines){
	if ($prev_gt) {
		
		#if ($array[0] eq $prev_chrom && ($gt eq "0/0" || "1/1") && ($prev_gt eq "0/0" || "1/1"))  {
		unless ($gt eq "0/1" || $array[0] != $prev_chrom) {
			if ($counter eq 0 || !defined($counter)) {
				$start= $array[0].":".$array[1];
				$start_pos=$array[1];
			}
			$counter++;
			#f you want to print to STOUT, remove $out_file_handle from 2 following lines
			#print $out_file_handle $array[0]."\t".$array[1]."\t".$gt."\t".$counter."\n";
		}
		else {
			my $difference=$prev_pos-$start_pos;
			#unless ($counter <= $options{size}) {
			unless ($difference <= $options{size}) {
				print $out_file_handle "\n--END OF LROH--".$prev_gt."\tSTRETCH:\t".$counter."\tREGION:\t".$start."-".$prev_pos."\tLENGTH(bp):\t".$difference;
			}
			$counter=0;
			
		}
	}
	#else {
	#	$counter=0;
	#}
	#$prev_gt = $array[2];
	$prev_gt = $gt;
	$prev_chrom = $array[0];
	$prev_pos = $array[1];
	
}
close FILE;
