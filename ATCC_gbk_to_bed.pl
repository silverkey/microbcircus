#!/usr/bin/perl
use strict;
use warnings;
use Bio::SeqIO;
use Data::Dumper;

# PARAMETERS TO SET UP
# -------------------------------------------------
# Genbank annotation from ATCC
my $gbk = 'Lactobacillus_acidophilus_ATCC_4356.gbk';
# -------------------------------------------------

# Output filename
my $bed = "$gbk.bed";
$bed =~ s/\.gbk\./\./;

open(OUT,">$bed");
my $seqio = Bio::SeqIO->new(-file => $gbk, -format => 'genbank');
while(my $seq = $seqio->next_seq) {
  for my $feature($seq->get_SeqFeatures) {
    if($feature->length < 1000000) {
      my $chr = $feature->seq_id;
      my $start = $feature->start;
      my $end = $feature->end;
      my $strand = $feature->strand;
      my $score = '.';
      my $id = $feature->has_tag('locus_tag') ? ($feature->get_tag_values('locus_tag'))[0] : 'NA';
      my $name = $feature->has_tag('gene') ? ($feature->get_tag_values('gene'))[0] : 'NA';
      my $description = $feature->has_tag('product') ? ($feature->get_tag_values('product'))[0] : 'NA';
      my $misc = $feature->has_tag('note') ? ($feature->get_tag_values('note'))[0] : 'NA';
      $description = $misc if $description eq 'NA' and $misc ne 'NA';
      print OUT join("\t",$chr,$start,$end,$id,$score,$strand,$name,$description)."\n";
    }
  }
}

__END__
