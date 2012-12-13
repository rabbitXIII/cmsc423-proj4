#!/usr/bin/perl -w

# Project Part 4
# Rohit Gopal
# CMSC 423 Fall 2012

use warnings;
use strict;
use Bio::SeqIO;
use Bio::Perl;

my $input_fasta = $ARGV[0];
my $seqio_obj;
my @seqs;
my $low = 999999999;
# consider intializing and doubling whenever reaching high utilization

# Validate the FASTA input file (use project 1)
$seqio_obj = Bio::SeqIO->new(-file => "$input_fasta", -format => "fasta");	

sub do_alignment {
	my ($id1, $id2) = @_;
	return 1 if $id1 < $id2;
}

while (my $seq = $seqio_obj->next_seq() ){
	$seqs[$seq->display_id] = $seq->seq;	
	$low = $seq->display_id if $seq->display_id < $low;
}
# build a hash of all k-mers in the set of sequences
# build clusters of hashes and test on each cluster

# run the local aligner on each pair of sequences that could overlap

for my $i ($low..$#seqs) {
	for my $j($i..$#seqs) { #only try to align each once
		if(do_alignment($i,$j)){
		#	print "align($seqs[$i],$seqs[$j]\n"; #Bio::Perl local alignment
		}
	}
}


# if S-W reports a proper dove-tail alignment, report the layout of the two sequences
# in an OVL formatted file
