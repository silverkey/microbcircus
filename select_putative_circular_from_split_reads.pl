#!/usr/bin/perl
use strict;
use warnings;
use Bio::DB::Sam;
use Data::Dumper;

# ----------------------------------------------
# VARIABLES TO SET FOR FILENAMES AND PARAMETERS
# ----------------------------------------------
my $bwa = '/home/remo/bin/bwa/bwa';
my $genome = '/home/remo/data/ecoli/genome/lactobacillus/Lactobacillus_acidophilus_ncfm.ASM1198v1.dna.chromosome.Chromosome.fa';
my $fastq = 'SRR6329247.fastq.gz';
my $gene_bed = "lactobacillus_genes.bed";
my $gene_sorted_bed = "$gene_bed\_sorted.bed";

my $samtools = '/home/remo/bin/samtools-1.9/samtools';
my $reference = 'reference_genome_seqid_length.txt';
my $seqtk = '/home/remo/bin/seqtk-master/seqtk';
my $bedtools = '/home/remo/bin/bedtools.2.29.2';

my $run_commands = 0;
my $threads = 6;
my $experiment = 'SRR6329247';
my $min_leng_cutoff = 150;
my $max_leng_cutoff = 1500;
my $score_cutoff = 60;

# ----------------------------------------------
# OUTPUT FILENAMES
# ----------------------------------------------
my $bwa_out = "$experiment\_matches.sam";
my $samtools_out = "$experiment\_splitted.bam";
my $splitted_id = "$experiment\_splitted_id.txt";
my $splitted_bed = "$experiment\_splitted.bed";
my $splitted_sorted_bed = "$experiment\_splitted_sorted_bed";
my $splitted_circular = "$experiment\_reads_supporting_putative_circular.xls";
my $splitted_gene_overlap = "$experiment\_splitted_gene_overlap.txt";

# ----------------------------------------------
# LET'S THE CODE BEGIN!
# ----------------------------------------------

# Map on the genome with BWA and retain only all the matches, non-mapping reads are discarded piping to samtools
my $bwa_command = "$bwa mem -t $threads $genome $fastq \| $samtools view -F4 -o $bwa_out";
run_command($bwa_command) if $run_commands;
# Select only split reads using samtools
my $samtools_command = "$samtools view -b -t $reference -f 2048 -o $samtools_out $bwa_out";
run_command($samtools_command) if $run_commands;

# Load the .bam containing only the split reads
my $bam = Bio::DB::Sam->new(-bam => "$samtools_out",
                            -expand_flags  => 1,
                            -split_splices => 1);

# Create the iterator object to get every match from the .bam
my $matches = $bam->features(-type => 'match',
                             -iterator => 1);

my $readhref = {};
my $reshref = {};
my $overhref = {};

# Open the bed file to calculate the overlap with genes
open(BED,">$splitted_bed");

# Iterate for every match
while(my $match = $matches->next_seq) {
  my $read = $match->query->name;
  my $chr = $match->seq_id;
  my $start = $match->start;
  my $end = $match->end;
  my $strand = $match->strand == 1 ? '+' : '-';
  my $cigar = $match->cigar_str;
  my @cigar1 = split_cigar($cigar);
  next if $cigar1[0] eq 'NOTOK';

  my $score = $match->score;
  # The FIRST_MATE tag tells if the read is from the first fastq
  my $mate = $match->get_tag_values('FIRST_MATE') == 1 ? 1 : 2;
  # The SA tag contains information about the secondary match of the split read
  # in the form Chromosome,2002383,+,46M30S,60,0;
  my $splitstr = $match->get_tag_values('SA');
  my @split = split(',',$splitstr);
  my @cigar2 = split_cigar($split[3]);
  next if $cigar2[0] eq 'NOTOK';

  my $achr = $chr;
  my $astart = $start;
  my $aend = $end;
  my $astrand = $strand;
  my $ascore = $score;
  my $bchr = $split[0];
  my $bstart = $split[1];
  my $bend = $cigar2[1] eq 'M' ? $bstart+$cigar2[0] : $bstart+$cigar2[2];
  my $bstrand = $split[2];
  my $bscore = $split[4];

  # Select only those split reads that pass the stringent filtering
  # The 2 matches must be on the same reference sequence
  # Score must be [>=$score_cutoff] for both the primary and the secondary matches
  # The strand must be the same for both the matches
  # The cigars should be equal in terms of number and not equal in term of caracthers (e.g. 30M35H and 30S35M is OK while 30M35H and 30M35S is NOT OK)
  if($achr eq $bchr and $ascore>=$score_cutoff and $bscore>=$score_cutoff and $astrand eq $bstrand and 
     $cigar1[1] ne $cigar2[1] and $cigar1[3] ne $cigar2[3] and $cigar1[0]==$cigar2[0] and $cigar1[2]==$cigar2[2]) {

    # Now we need to select only the split that can suggest a backsplice site. The following workaround is done keeping in mind that
    # THE CIGAR IS REFERRED TO THE GENOME AND IS REVERSE COMPLEMENTED WHEN THE STRAND OF THE MATCH IS NEGATIVE

    my $abs_start;
    my $abs_end;

    # the strand is + and the M is the first so a is the first part of query and should come after b
    if($cigar1[1] eq 'M' and $strand eq '+') {
      if($astart > $bstart) {
        $abs_start = $bstart;
        $abs_end = $aend;
      }
    }
    # the strand is - and the M is the first so a is the second part of query and should come after b
    elsif($cigar1[1] eq 'M' and $strand eq '-') {
      if($astart > $bstart) {
        $abs_start = $bstart;
        $abs_end = $aend;
      }
    }
   # the strand is + and M is the second so a is the second part of query and should come before b
   elsif($cigar1[3] eq 'M' and $strand eq '+') {
      if($astart < $bstart) {
        $abs_start = $astart;
        $abs_end = $bend;
      }
    }
    # the strand is - and M is the second so a is the first part of query and should come before b
    elsif($cigar1[3] eq 'M' and $strand eq '-') {
      if($astart < $bstart) {
        $abs_start = $astart;
        $abs_end = $bend;
      }
    }

    # If abs_start is undefined means that the split does not suggest a backsplice but a simple splice
    next unless $abs_start;

    my $abs_length = $abs_end - $abs_start + 1;
    
    if($abs_length >= $min_leng_cutoff and $abs_length <= $max_leng_cutoff) {

      # Print all the info in output if the filters have been passed
      $readhref->{$read} ++;
      my $read_i = "$read\_".$readhref->{$read};
      $reshref->{$read}->{$read_i} = join("\t",$read,$read_i,$mate,$achr,$astart,$aend,$astrand,$cigar,@cigar1,$ascore,
                                                          $bchr,$bstart,$bend,$bstrand,$split[3],@cigar2,$bscore,
                                                          $abs_start,$abs_end,$abs_length);

      print BED join("\t",$chr,$abs_start,$abs_end,$read_i)."\n";
    }
  }
}

# Sort the bed files to use it in bedtools
my $bedtools_command_1 = "$bedtools sort -i $splitted_bed > $splitted_sorted_bed";
run_command($bedtools_command_1);
my $bedtools_command_2 = "$bedtools sort -i $gene_bed > $gene_sorted_bed";
run_command($bedtools_command_2);

# Calculate the overlap with genes
my $bedtools_command_3 = "$bedtools closest -d -t first -a $splitted_sorted_bed -b $gene_sorted_bed > $splitted_gene_overlap";
run_command($bedtools_command_3);
open(OVER,$splitted_gene_overlap);
while(my $row = <OVER>) {
  chomp($row);
  my @f = split("\t",$row);
  $overhref->{$f[3]} = join("\t",$f[7],$f[9],$f[10],$f[11],$f[12])
}


# Write the file containing the ids of the split reads
open(ID,">$splitted_id");
foreach my $key(keys %$readhref) {
  print ID "$key\n";
}

# Use seqtk to extract the reads corresponding to the split reads
my $selhref = {};
my @fastq = split(' ',$fastq);
for my $i(1..scalar(@fastq)) {
  my $fq = $fastq[$i-1];
  my $splitted_seq = "$fq\_splitted_seq.txt";
  my $seqtk_command = "$seqtk subseq -t $fq $splitted_id > $splitted_seq";
  run_command($seqtk_command);
  open(IN,"$splitted_seq");
  while(my $row = <IN>) {
    chomp($row);
    my @field = split("\t",$row);
    $selhref->{$field[0]}->{$i} = $field[2];
  }
}

# Write the output file
my $out = "$splitted_circular";
open(OUT,">$out");
# Prepare and print the header with the names of the rows to output
my $head = join("\t",'read','read_i','mate','achr','astart','aend','astrand','acigar','ac1',
                     'ac2','ac3','ac4','aclength','ascore','bchr','bstart','bend',
                     'bstrand','bcigar','bc1','bc2','bc3','bc4','bclength','bscore',
                     'abs_start','abs_end','abs_length','closest_id','closest_strand',
                     'closest_name','closest_description','closest_distance','read1','read2');
print OUT "$head\n";
foreach my $key (keys %$reshref) {
  foreach my $i (keys %{$reshref->{$key}}) {
    print OUT $reshref->{$key}->{$i}."\t";
    print OUT $overhref->{$i}."\t";
    print OUT $selhref->{$key}->{1};
    print OUT "\t".$selhref->{$key}->{2} if exists $selhref->{$key}->{2};
    print OUT "\tNA" if not exists $selhref->{$key}->{2};
    print OUT "\n";
  }
}

# Function to split the cigar and take only those with only 2 pieces (e.g. 30M35S is OK and 30M1D34S is NOT OK)
# It also calculates the length of the match by taking into account the digits associated to the M in the cigar
sub split_cigar {
  my $cigar = shift;
  if($cigar =~ /^(\d+)([MHS])(\d+)([MHS])$/) {
    my $length = $2 eq 'M' ? $1 : $3;
    return($1,$2,$3,$4,$length)
  }
  else {
    return('NOTOK')
  }
}

sub run_command {
  my $command = shift;
#  return unless $run_commands;
  print "\nLAUNCHING SYSTEM CALL:\n\t$command\n";
  system($command);
  die "ERROR using command:\n\t$command\:\n\t$!" unless $? == 0;
  print "$command DONE!\n";
}

__END__

#my @tags = $match->get_all_tags;
#foreach my $tag (@tags) {
#  print "$tag\n";
#}

#sel = subset(t, ascore==60 & bscore==60 & astrand==bstrand & ac2!=bc2 & ac4!=bc4 & ac1==bc1 & ac3==bc3)
#pos = subset(sel, astrand=='+' & aend-bstart>=150 & aend-bstart<=1500)
#neg = subset(sel, astrand=='-' & bend-astart>=150 & bend-astart<=1500)

/home/remo/bin/bwa/bwa mem -t 6 ../../genome/Escherichia_coli_str_k_12_substr_mg1655.ASM584v2.dna.chromosome.Chromosome.fa SRR441585_1.fastq.gz SRR441585_2.fastq.gz | /home/remo/bin/samtools-1.9/samtools view -F4 -o output
/home/remo/bin/samtools-1.9/samtools view -b -t reference -f 2048 -o splitted.bam output


 2109  /home/remo/bin/bedtools.2.29.2 sort -i SRR6329247_split_reads_id.bed > SRR6329247_split_reads_id_sorted.bed
 2110  /home/remo/bin/bedtools.2.29.2 sort -i lactobacillus_genes.bed > lactobacillus_genes_sorted.bed 
 2111  /home/remo/bin/bedtools.2.29.2 closest -d -t first -a SRR6329247_split_reads_id_sorted.bed -b lactobacillus_genes_sorted.bed
 2112  history 

