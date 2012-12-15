#!/usr/bin/perl -w

# Project Part 4
# Rohit Gopal
# CMSC 423 Fall 2012

use warnings;
use strict;
use Bio::SeqIO;
use Bio::Perl;
use List::Util qw[min max];

my $input_fasta = $ARGV[0];
my $input_hoxd = $ARGV[1];
my $output_file = $ARGV[2];
my $gap_open = $ARGV[3] || -2000;
my $gap_extend = $ARGV[4] || -200;
my $min_identity = $ARGV[5] || 98;
my $window = $ARGV[6] || 45;
my $max_hang = $ARGV[7] || 60;
my $seqio_obj;
my @seqs;
my $low = 999999999;
my %seq_hash;
my $output = "";
my %aligned;


my %hoxd_scores = (	A => {},
			T => {},
			C => {},
			G => {});

# consider intializing and doubling whenever reaching high utilization

# Validate the FASTA input file (use project 1)
$seqio_obj = Bio::SeqIO->new(-file => "$input_fasta", -format => "fasta");	

sub read_hoxd2 {
	my $hoxd_handle;
	open($hoxd_handle, $input_hoxd);
	my @hoxd_lines = <$hoxd_handle>;
	close($hoxd_handle);

	foreach(@hoxd_lines){ 
		if($_ =~ /([ACTGactg]),([ACTGactg])=(-?\d+)/){
			$hoxd_scores{lc($1)}{lc($2)} = $3; 	
			$hoxd_scores{lc($2)}{lc($1)} = $3;
			$hoxd_scores{uc($1)}{uc($2)} = $3; 	
			$hoxd_scores{uc($2)}{uc($1)} = $3;
		}
	}
}

sub do_alignment {
	my ($id1, $id2) = @_;
	return 1 if $id1 < $id2;
}

sub read_fasta {
	my $index = 1;
	while (my $seq = $seqio_obj->next_seq() ){
		$seqs[$index] = $seq->seq;	
		$index++;
	}
}

sub map_sequences {
	for my $i (1..$#seqs) {
		my $length = length($seqs[$i]);
		if ( $length < $window ) {
			next;
		}
		my $kmer = substr($seqs[$i],0, $window);
		if (exists($seq_hash{$kmer})) { 
			push (@{$seq_hash{$kmer}}, $i);
			
		} else {
			$seq_hash{$kmer} = [ $i ];
		}	
		for my $index (1..($length-$window)) {
			$kmer = substr($seqs[$i], $index, $window);
			if (exists($seq_hash{$kmer})) { 
				push (@{ $seq_hash{$kmer} }, $i);
			} else {
				$seq_hash{$kmer} = [ $i ];
			}	
		}
	}
}

sub align {
	my ( $rd1, $rd2 ) = @_;	
	my ( $seq1, $seq2 ) = ( $seqs[$rd1], $seqs[$rd2] );
	my  ( $ahg, $bhg ) = ( '\N' ) x 2;

	#perform the local alignment

		# we add one to make room for the initial row / col
	my $length1 = length($seq1);
	my $length2 = length($seq2);

		#doing the global alignment is easier with the strings being arrays

	my @seq1 = split(//, $seq1);
	my @seq2 = split(//, $seq2);


		#the traceback will have the following values
		#	1 	left		gap to the right
		#	2	middle		match / mismatch
		#	4	up		gap down
		#	3	left/middle	gap to right or match/mismatch
		#	5	left/up		gap right or gap down
		#	6	middle/up	match/mismatch or gap down
		#	7	all		all have same score
		# use the 1 2 4 method to keep track of traceback

	my $right_trace = 1;
	my $diag_trace = 2;
	my $down_trace = 4;

	my @matrix = undef;
	my @traceback = undef;
	my @matrix_x = undef;
	my @matrix_y = undef;

	#initialize the matrices
	# this could honestly just be done in a single matrix, but it's easier
	# to read
	for my $index (0..$length1) {
		for my $index2 (0..$length2) {
			$matrix[$index][$index2] = 0;
			$traceback[$index][$index2] = 0;
			$matrix_x[$index][$index2] = 0;
			$matrix_y[$index][$index2] = 0;
		}
	}	
		# O(n)
	for my $index (1..$length1) {
		$traceback[$index][0] =  $right_trace;
		$matrix_x[$index][0] = $gap_open + $index * $gap_extend;
	}

		# O(n)
	for my $index (1..$length2) {
		$traceback[0][$index] = $down_trace; 
		$matrix_y[0][$index] = $gap_open + $index * $gap_extend;
	}




	my @max = ({ i => 0,
			j => 0,
			score => 0 });

		#O(n^2) no way around this..? ;(
	for my $i (1..$length1) {
		for my $j (1..$length2) {
			$matrix[$i][$j] = $hoxd_scores{$seq1[$i-1]}{$seq2[$j-1]} + max( $matrix[$i-1][$j-1],
							 $matrix_x[$i-1][$j-1],
							 $matrix_y[$i-1][$j-1] );
			$matrix_x[$i][$j] = max( $gap_open + $gap_extend + $matrix[$i][$j-1],
						$gap_extend + $matrix_x[$i][$j-1],
						$gap_open + $gap_extend + $matrix_y[$i][$j-1] );
			$matrix_y[$i][$j] = max( $gap_open + $gap_extend + $matrix[$i-1][$j],
						$gap_extend + $matrix_y[$i-1][$j],
						$gap_open + $gap_extend + $matrix_x[$i-1][$j] );
			my $choice = max( $matrix[$i][$j], $matrix_x[$i][$j], $matrix_y[$i][$j]);
			if ($matrix[$i][$j] == $choice) {
				$traceback[$i][$j] = 2;
			} elsif ($matrix_x[$i][$j] == $choice) { 
				$traceback[$i][$j] = 4;
			} else {
				$traceback[$i][$j] = 1;	
			}
			# our max score can only end with a match, so we only check matrix
			if ( $max[0]{score} < $matrix[$i][$j] ) {
				@max = undef;
				@max = ( { i => $i,
					j => $j,
					score => $matrix[$i][$j]} );
			} elsif ($max[0]{score} == $matrix[$i][$j]){
				push(@max, { i => $i, j => $j, score => $matrix[$i][$j]});

			}
			
		}
	}
	#do the traceback for each valid sub section

	my ($i, $j, $score) = ($max[0]->{i}, $max[0]->{j}, $max[0]->{score});	
	my @pointers = undef; # tracks the indexes of the sequences
	my @result = undef;
	my $length_of_seq = 0;
	$pointers[1] = $i;
	$pointers[3] = $j;	
		#we need to print backwards into the array since we don't know how many gaps there are
	while ( ($matrix[$i][$j] != 0) && ($i != 0 || $j != 0)) {
		$result[1][$length_of_seq] = " ";
		if ($traceback[$i][$j] == 2) { #diagonal match / mismatch
			$result[1][$length_of_seq] = "|" if ($seq1[$i-1] eq $seq2[$j-1]);
			$result[0][$length_of_seq] = $seq1[--$i];
			$result[2][$length_of_seq] = $seq2[--$j];
		} elsif ($traceback[$i][$j] == 1) {  
			$result[0][$length_of_seq] = $seq1[--$i];
			$result[2][$length_of_seq] = "-";
		} else {
			$result[0][$length_of_seq] = "-";
			$result[2][$length_of_seq] = $seq2[--$j];
		}
		$length_of_seq++;
	}
	$pointers[0] = $i;
	$pointers[2] = $j;	

	#print out the score and alignment
	my $identities_count = 0;

	for my $i (1..$length_of_seq) {
		$identities_count++ if ($result[1][0-$i] eq "|");
	}
	my $longest = max($pointers[1] - $pointers[0], $pointers[3] - $pointers[2]);
	my $percent = $identities_count / ($longest) * 100;

	#not dove-tail alignment, so die
	if ( $pointers[0] != 0 && $pointers[2] != 0 ) {
		return 0;
	}
	
	if ($pointers[0] == 0 ){
		$ahg = $pointers[1] - length($seq1);
		$bhg = -$pointers[2];
	} else {
		$ahg = $pointers[0];
		$bhg = length($seq2) - $pointers[3];
		# not dove-tail alignment
	}

	if ($percent >= $min_identity && $ahg < $max_hang && $bhg < $max_hang) {
		# if within the constraints of alignment, return!
		# to prevent multiple copy entries
			# within % identity
			# ahg and bhg within supplied bounds
		# if S-W reports a proper dove-tail alignment, report the layout of the two sequences
		# in an OVL formatted file
			
		$aligned{$rd1}{$rd2} = 1;
		$aligned{$rd2}{$rd2} = 1;

		return { adj => 'N', rd1 => $rd1, rd2 => $rd2, scr => 0, ahg => $ahg, bhg => $bhg };
	}
	# else return bad
	return 0; 
}

sub local_align_sequences {
	foreach my $key ( keys %seq_hash ) {
		if($#{$seq_hash{$key}} < 2){ 
			next;
		}
		for my $i (0..( @{$seq_hash{$key}} -2 )){
			for my $j (($i+1)..( @{$seq_hash{$key}} -1 )) {
				if(! exists($aligned{$seq_hash{$key}[$i]}{$seq_hash{$key}[$j]})){
					if (my $res = align($seq_hash{$key}[$i], $seq_hash{$key}[$j])){
						$output .= "{OVL\nadj:$res->{adj}\nrds:$res->{rd1},$res->{rd2}\nscr:$res->{scr}\nahg:$res->{ahg}\nbhg:$res->{bhg}\n}\n";
					}
				}
			}
		}

	}
}

read_hoxd2();
read_fasta();


# build a hash of all k-mers in the set of sequences

map_sequences();

# run the local aligner on each pair of sequences that could overlap

local_align_sequences();

#print in OVL format 

open (MYFILE, ">>$output_file");
print MYFILE $output;
close (MYFILE);
