#!/usr/bin/perl -w

# Project Part 4
# Rohit Gopal
# CMSC 423 Fall 2012

use warnings;
use strict;
use Bio::SeqIO;
use Bio::Perl;
#use Bio::Tools::dpAlign;
#use Bio::AlignIO;

my $input_fasta = $ARGV[0];
my $input_hoxd = $ARGV[1];
my $seqio_obj;
my %seqs;
my %seq_hash;
my $window = 6;
my $output = "";
my @kees;
#my $matrix = Bio::Matrix::IO->new(-format => 'scoring', -file => '$input_hoxd');
#my $factory = new dpAlign(-matrix => $matrix,
				 #-alg => Bio::Tools::dpAlign::DPALIGN_LOCAL_MILLERMYERS);


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
	while (my $seq = $seqio_obj->next_seq() ){
		$seqs{$seq->display_id} = $seq->seq;	
	}
}

sub map_sequences {
	for my $i (0..$#kees) {
		my $length = length( $seqs{$kees[$i]} );
		if ( $length < $window ) {
			next;
		}
		my $kmer = substr($seqs{$kees[$i]},0, $window);
		if (exists($seq_hash{$kmer})) { 
			push (@{$seq_hash{$kmer}}, $kees[$i]);
			
		} else {
			$seq_hash{$kmer} = [ $kees[$i] ];
		}	
		for my $index (1..($length-$window)) {
			$kmer = substr($seqs{$kees[$i]}, $index, $window);
			if (exists($seq_hash{$kmer})) { 
				push (@{ $seq_hash{$kmer} }, $kees[$i]);
			} else {
				$seq_hash{$kmer} = [ $kees[$i] ];
			}	
		}
	}
}

sub align {
	my ( $rd1, $rd2 ) = @_;	
	my ( $seq1, $seq2 ) = ( $seqs{$rd1}, $seqs{$rd2} );
	my  ( $ahg, $bhg ) = ( '\N' ) x 2;

	# if within the constraints of alignment, return!
	return { adj => 'N', rd1 => $rd1, rd2 => $rd2, scr => 0, ahg => $ahg, bhg => $bhg };
	# else return bad
}

sub local_align_sequences {
	foreach my $key ( keys %seq_hash ) {
		if($#{$seq_hash{$key}} < 2){ 
			next;
		}
		for my $i (0..( @{$seq_hash{$key}} -2 )){
			for my $j (($i+1)..( @{$seq_hash{$key}} -1 )) {
				my $res = align($seq_hash{$key}[$i], $seq_hash{$key}[$j]);
				$output .= "{OVL\n adj:$res->{adj}\n rds:$res->{rd1},$res->{rd2}\n scr:$res->{scr}\n ahg:$res->{ahg}\n bhg:$res->{bhg}\n}\n";
			}
		}

	}
}

read_hoxd2();
read_fasta();

@kees = (keys %seqs);

# build a hash of all k-mers in the set of sequences

map_sequences();

# run the local aligner on each pair of sequences that could overlap

local_align_sequences();

# if S-W reports a proper dove-tail alignment, report the layout of the two sequences
# in an OVL formatted file

open (MYFILE, '>>output.txt');
print MYFILE $output;
close (MYFILE);
