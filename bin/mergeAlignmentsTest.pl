#!/usr/bin/perl

use strict;
use warnings;

use Bio::DB::Sam;

use Cwd qw(abs_path);
use FindBin;
use lib abs_path("$FindBin::Bin/../lib");

require "PerlSub.pl";
require "TranslocHelper.pl";

my $bam = "/Volumes/AltLab/Translocation/RawData/Alt024-20130429/NewPipelineTest/results/CC004_Alt024/CC004_Alt024.bam";
my $qname = "M01407:14:000000000-A3ERD:1:1101:3916:15223";
print $ENV{'GENOME_DB'}."/mm9/mm9.fa\n";

my $samobj = Bio::DB::Sam->new(-bam => $bam,
                                 -fasta => $ENV{'GENOME_DB'}."/mm9/mm9.fa",
                                 -expand_flags => 1);


my ($R1_aln, $R2_aln) = $samobj->get_features_by_name($qname);

merge_alignments($R1_aln,$R2_aln,$samobj);
