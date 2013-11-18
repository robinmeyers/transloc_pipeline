#!/usr/bin/env perl


use strict;
use warnings;
use Getopt::Long;
use Carp;
use IO::File;
use Text::CSV;
use threads;
use Interpolation 'arg:@->$' => \&argument;
use Time::HiRes qw(gettimeofday tv_interval);

use Cwd qw(abs_path);
use FindBin;
use lib abs_path("$FindBin::Bin/../lib");

my $GENOME_DB = $ENV{'GENOME_DB'};
defined $GENOME_DB or croak "Error: set environment variable GENOME_DB";

require "4CHelper.pl";
require "Perlsub.pl";
require "pslHelper.pl";

# Flush output after every write
select( (select(STDOUT), $| = 1 )[0] );
my $dir = "/Volumes/AltLab/4C/Results/Alt017-20130124/Sample/";
my $expt_id = "SJ001_Alt017";
my $expt_hash = {	psl 		=> "$dir/$expt_id/alignments/$expt_id.psl" ,
									refdb		=> "/Volumes/AltLab/Genomes/mm8/mm8_mask_chr12_113931642-113933108.2bit" ,
									raw			=> "/Volumes/AltLab/4C/RawData/Alt017-20130124/Final/sample/$expt_id.fasta" ,		
									redpsl	=> "$dir/$expt_id/alignments/red.psl" ,
									redppsl	=> "$dir/$expt_id/alignments/redp.psl" ,
									blupsl	=> "$dir/$expt_id/alignments/blu.psl" ,
									bluppsl	=> "$dir/$expt_id/alignments/blup.psl" ,
									exptdir	=> "$dir/$expt_id" };


make_tlxl($expt_id,$expt_hash);