#!/usr/bin/perl -w

# Project Part 4
# Rohit Gopal
# CMSC 423 Fall 2012

use warnings;
use strict;
use Bio::SeqIO;

my $input_fasta = $ARGV[0];

# consider intializing and doubling whenever reaching high utilization
my @fasta_sequences;

# taken from project 1
sub validate_and_store_fasta {
	my $filename = shift;
	my $handle;
	open($handle, $filename) or die "Error reading file: $filename.
					Please be sure that this is a valid file \n";
	my $index = 0;
	my $current_id = "";
	my $current_length = 0;
	my $new_id = "";
	while(<$handle>){
		chomp;
		if (/^>/) { #identifier line
			s/^>\s*//; #removes the > and whitespaces at start 
			s/^\s+//;  # and end of the identifier 
			s/ .+$//;
			$current_id = $_;
			$index++;
		} elsif($current_id ne ""){
			s/^\s+//; #removes whitespaces and such
			s/\s+$//;
			s/\s+//g;
			next if (/^$/); #only check basepairs
			my $bp = $_; 
			if ($bp =~ m/^[ACTG]/ and $current_id ne "") {
				$fasta_sequences[$current_id] .= $bp;
			} else {
				return 0;
			}
			$index++;
		} elsif($index == 0) {
			die "Error. Something is wrong with the file"
		}
	}
	1;
}

# Validate the FASTA input file (use project 1)

die "Invalid FASTA file: $input_fasta" if ! validate_and_store_fasta($input_fasta);

# build a hash of all k-mers in the set of sequences

# run the local aligner on each pair of sequences that could overlap

# if S-W reports a proper dove-tail alignment, report the layout of the two sequences
# in an OVL formatted file
