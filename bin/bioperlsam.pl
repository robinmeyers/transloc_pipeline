#!/usr/bin/perl

use strict;
use warnings;

use Bio::DB::Sam;

my $samfile = "/Volumes/AltLab/Translocation/RawData/Alt024-20130429/NewPipelineTest/results/AC004_Alt024/AC004_Alt024.bam";

my $sam = Bio::DB::Sam->new(-bam => "$samfile",
                            -fasta => "/Volumes/AltLab/Genomes/mm9/mm9.fa");


my $iterator = $sam->features(-iterator => 1);

while (my $feature = $iterator->next_seq) {
  my $feature2 = $iterator->next_seq if $feature->proper_pair;

  my $output = join("\t",)

}