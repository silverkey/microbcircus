#!/usr/bin/perl
use strict;
use warnings;
use Bio::Tools::GFF;
use Data::Dumper;

# This script can be run on the ensembl GFF3 annotation of the genome. It write in output
# a bed file with the locations of all the genes in order to be used with bedtools.

my $USAGE = "\n\tUSAGE: perl $0 [ensembl gff3 file] [name of the output file to generate]\n\n";
die $USAGE unless scalar(@ARGV) == 2;
die $USAGE unless -e $ARGV[0];
die "$USAGE\tThe names of the two parameters are equal!\n\n" if $ARGV[0] eq $ARGV[1];
die "$USAGE\tA file with the same name of the putput to generate already exists!\n\n" if -e $ARGV[1];
my $file = $ARGV[0];
my $out = $ARGV[1];
open(OUT,">$out");

my $gffio = Bio::Tools::GFF->new(-file => $file, -gff_version => 3);
while(my $feature = $gffio->next_feature()) {
  if($feature->primary_tag =~ /gene/) {
    my $chr = $feature->seq_id;
    my $start = $feature->start;
    my $end = $feature->end;
    my $strand = $feature->strand;
    my $id = ($feature->get_tag_values('gene_id'))[0];
    my $score = '.';
    my $name = $feature->has_tag('Name') ? ($feature->get_tag_values('Name'))[0] : 'NA';
    my $description = $feature->has_tag('description') ? ($feature->get_tag_values('description'))[0] : 'NA';
    print OUT join("\t",$chr,$start,$end,$id,$score,$strand,$name,$description)."\n";
  }
}

__END__
