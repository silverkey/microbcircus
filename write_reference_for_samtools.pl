#!/usr/bin/perl
use strict;
use warnings;
use Bio::SeqIO;

# In the pipeline you usually use pipes to apply samtools to the output of bwa at runtime.
# This might cause the rimotion of the header in the bam/sam output.
# With this script you create a table containg the two columns sequence_id and length
# required by samtools in case you run it on sam/bam files without a header

my $USAGE = "\n\tUSAGE: perl $0 [fasta with the genome] [name of the output file to generate]\n\n";
die $USAGE unless scalar(@ARGV) == 2;
die $USAGE unless -e $ARGV[0];
die "$USAGE\tThe names of the two parameters are equal!\n\n" if $ARGV[0] eq $ARGV[1];
die "$USAGE\tA file with the same name of the putput to generate already exists!\n\n" if -e $ARGV[1];

my $fasta = $ARGV[0];
my $out = $ARGV[1];
open(OUT,">$out");

my $seqio = Bio::SeqIO->new(-file => $fasta, -format => 'fasta');

while(my $seq = $seqio->next_seq) {
  print OUT $seq->id."\t".$seq->length."\n";
}

__END__
