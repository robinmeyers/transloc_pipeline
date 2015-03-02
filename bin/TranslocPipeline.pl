#!/usr/bin/env perl

use strict;
use warnings;
use Getopt::Long;
use Sys::Hostname;
use Carp;
use Switch;
use IO::Handle;
use IO::File;
use Text::CSV;
use File::Basename;
use File::Which;
use File::Copy;
use Bio::SeqIO;
use Bio::DB::Sam;
use List::Util qw(min max shuffle);
use List::MoreUtils qw(pairwise);
use Interpolation 'arg:@->$' => \&argument;
use Time::HiRes qw(gettimeofday tv_interval);
use Sys::Info;
use Sys::Info::Constants qw( :device_cpu );
use Data::Dumper;
use Cwd qw(abs_path);
use FindBin;
use lib abs_path("$FindBin::Bin/../lib");

require "PerlSub.pl";
require "TranslocHelper.pl";
require "TranslocFilters.pl";

# The pipeline searches for a genome build in a fasta file
# It expects this fasta file to live in the GENOME_DB
# e.g. GENOME_DB/mm9/mm9.fa
my $GENOME_DB = $ENV{'GENOME_DB'};
defined $GENOME_DB or croak "Error: set environment variable GENOME_DB";

# Bowtie2 searches for a reference index using this environment variable
# e.g. BOWTIE2_INDEXES/mm9.1.bt2
my $BOWTIE2_INDEXES = $ENV{'BOWTIE2_INDEXES'};
defined $BOWTIE2_INDEXES or croak "Error: set environment variable BOWTIE2_INDEXES";

# Flush output after every write
select( (select(STDOUT), $| = 1 )[0] );

##
## This program runs the main Translocation Pipeline
## It mainly consisists of alignments, junction detection,
## filtering, and various post-processing scripts
## 
## run with "--help" for usage information
##
## Robin Meyers

# Forward declarations
sub parse_command_line;
sub write_parameters_file;
sub align_to_breaksite;
sub align_to_adapter;
sub align_to_genome;
sub process_alignments;
sub mark_repeatseq_junctions;
sub mark_duplicate_junctions;
sub filter_junctions;
sub sort_junctions;
sub post_process_junctions;
sub clean_up;


my $commandline;
my $local_time;

# Global flags and arguments, 
# Set by command line arguments
my $read1;
my $read2;
my $workdir;
my $assembly;
my $brk_chr;
my $brk_start;
my $brk_end;
my $brk_strand;
my $mid_fa;
my $prim_fa;
my $break_fa;
my $breakcoord;
my $adapt_fa;
my $cut_fa;
my $threads = 1;
my $random_barcode;
my $repeatseq_bedfile;
my $simfile;


my $skip_alignment;
my $skip_process;
my $skip_dedup;
my $no_dedup;
my $no_clean;

our $debug_level = 0;

our $params = {};



# Bowtie parameters
$params->{match_award} = 2;
$params->{mismatch_pen} = "8,2";
$params->{n_base_pen} = 1;
$params->{read_gap_pen} = "8,4";
$params->{ref_gap_pen} = "8,4";
$params->{score_min} = "C,50";
$params->{D_effort} = 20;
$params->{R_effort} = 3;
$params->{seed_mismatch} = 0;
$params->{seed_length} = 20;
$params->{seed_interval} = "C,6";
$params->{breaksite_alignments} = 5;
$params->{genome_alignments} = 20;

$params->{force_bait} = 1;
# Bait alignment must be within this many bp if forcing
$params->{max_brkstart_dif} = 20;
# Maximum overlap of two adjacent alignments in OCS
# Calculated by intersection/union
$params->{max_overlap} = 0.5;
# OCS junction penalty
$params->{brk_pen} = 50;
# Max gap for concordant alignment (between R1 and R2)
# Need to check if it's the gap or the entire aln fragment
$params->{max_pe_gap} = 1500;
$params->{pe_pen} = 50;
$params->{max_dovetail} = 10;

# uncut filter params
$params->{max_bp_after_cutsite} = 10;

# misprimed filter params
$params->{min_bp_after_primer} = 10;

# freqcut filter params

# largegap filter params
$params->{max_largegap} = 30;

# mapqual filter params
$params->{mapq_ol_thresh} = 0.9;
$params->{mapq_score_thresh} = 10;
# $params->{mapq_mismatch_int} = 1.5;
# $params->{mapq_mismatch_coef} = 0.01;

# duplicate filter parameters
$params->{dedup_offset_dist} = 1;
$params->{dedup_break_dist} = 1;


# my $user_bowtie_opt = "";
# my $user_bowtie_adapter_opt = "";
# my $user_bowtie_breaksite_opt = "";

# Global variables
my @tlxl_header = tlxl_header();
my @tlx_header = tlx_header();
my @tlx_filter_header = tlx_filter_header();


my @dispatch_names = qw(unaligned baitonly uncut misprimed freqcut largegap mapqual breaksite sequential);

# This is a dispatch table
my %filter_dispatch;
@filter_dispatch{@dispatch_names} = map {eval '\&filter_'.$_ } @dispatch_names;

my @filters = (@dispatch_names,"repeatseq","duplicate");


#
# Start of Program
#

my $t0 = [gettimeofday];

parse_command_line;

my $bt2_break_opt = join(" ","--local --no-1mm-upfront",
                             "-D",$params->{D_effort},
                             "-R",$params->{R_effort},
                             "-N",$params->{seed_mismatch},
                             "-L",$params->{seed_length},
                             "-i",$params->{seed_interval},
                             "--ma",$params->{match_award},
                             "--mp",$params->{max_mismatch_pen}.",".$params->{min_mismatch_pen},
                             "--np",$params->{n_base_pen},
                             "--rdg",$params->{read_gap_open}.",".$params->{read_gap_ext},
                             "--rfg",$params->{ref_gap_open}.",".$params->{ref_gap_ext},
                             "--score-min",$params->{score_min},
                             "-k",$params->{breaksite_alignments},
                             "-p",$threads,
                             "--no-unal --reorder -t");

my $bt2_opt = join(" ","--local --no-1mm-upfront",
                             "-D",$params->{D_effort},
                             "-R",$params->{R_effort},
                             "-N",$params->{seed_mismatch},
                             "-L",$params->{seed_length},
                             "-i",$params->{seed_interval},
                             "--ma",$params->{match_award},
                             "--mp",$params->{max_mismatch_pen}.",".$params->{min_mismatch_pen},
                             "--np",$params->{n_base_pen},
                             "--rdg",$params->{read_gap_open}.",".$params->{read_gap_ext},
                             "--rfg",$params->{ref_gap_open}.",".$params->{ref_gap_ext},
                             "--score-min",$params->{score_min},
                             "-k",$params->{genome_alignments},
                             "-p",$threads,
                             "--reorder -t");

my $bt2_adapt_opt = join(" ","--local --no-1mm-upfront",
                             "-D",15,
                             "-R",2,
                             "-N",1,
                             "-L",10,
                             "-i","C,6",
                             "--ma",2,
                             "--mp","6,2",
                             "--np",1,
                             "--rdg","5,3",
                             "--rfg","5,3",
                             "--score-min","C,20",
                             "-p",$threads,
                             "--no-unal --reorder -t");

# my $default_bowtie_adapter_opt = "--very-sensitive-local -L 10 --ma 2 --mp 6,2 --np 1 --rdg 5,3 --rfg 5,3 --score-min C,20 --no-unal -p $bowtie_threads --reorder -t";
# my $default_bowtie_breaksite_opt = "--very-sensitive-local --ma 2 --mp 10,2 --np 2 --rdg 6,4 --rfg 6,4 --score-min C,50 -k 5 --no-unal -p $bowtie_threads  --reorder -t";
# my $default_bowtie_opt = "--very-sensitive-local --ma 2 --mp 10,2 --np 2 --rdg 6,4 --rfg 6,4 --score-min C,50 -k 20 -p $bowtie_threads --reorder -t";


# This may be too complicated -
# in the future perhaps just change it
# such that the user inputs full bowtie2 option string
# my $bt2_break_opt = manage_program_options($default_bowtie_breaksite_opt,$user_bowtie_breaksite_opt);
# my $bt2_adapt_opt = manage_program_options($default_bowtie_adapter_opt,$user_bowtie_adapter_opt);
# my $bt2_opt = manage_program_options($default_bowtie_opt,$user_bowtie_opt);


# croak "Error: cannot find match award in bowtie2 options" unless $bt2_opt =~ /-ma (\d+)/;
# my $match_award = $1;
# croak "Error: cannot find mismatch penalty in bowtie2 options" unless $bt2_opt =~ /-mp (\d+),\d+/;
# my $mismatch_penalty = $1;

# # I warned you
# carp "Warning: match award in bowtie2 does not equal OCS overlap penalty" unless $match_award eq $params->{overlap_mult};


my $expt = basename($workdir);
my $expt_stub = "$workdir/$expt";


my $assembly_fa = "$GENOME_DB/$assembly/$assembly.fa";
my $break_bt2idx;
if (defined $break_fa) {
  my @brk_file_path = parseFilename($break_fa);
  $break_bt2idx = join("",@brk_file_path[0..1]);
}
my @adapt_file_path = parseFilename($adapt_fa);
my $adapt_bt2idx = join("",@adapt_file_path[0..1]);

my $R1_brk_sam = "${expt_stub}_R1_brk.sam";
my $R2_brk_sam = "${expt_stub}_R2_brk.sam";
my $R1_sam = "${expt_stub}_R1.sam";
my $R2_sam = "${expt_stub}_R2.sam";
my $R1_adpt_sam = "${expt_stub}_R1_adpt.sam";
my $R2_adpt_sam = "${expt_stub}_R2_adpt.sam";

my $R1_brk_bam = "${expt_stub}_R1_brk.bam";
my $R2_brk_bam = "${expt_stub}_R2_brk.bam";
my $R1_bam = "${expt_stub}_R1.bam";
my $R2_bam = "${expt_stub}_R2.bam";
my $R1_adpt_bam = "${expt_stub}_R1_adpt.bam";
my $R2_adpt_bam = "${expt_stub}_R2_adpt.bam";

my ($R1_brk_samobj,$R2_brk_samobj,$R1_adpt_samobj,$R2_adpt_samobj,$R1_samobj,$R2_samobj);


my $tlxlfile = "${expt_stub}.tlxl";
my $tlxfile = "${expt_stub}.tlx";
my $mapqfile = "${expt_stub}_mapq.txt";

my ($tlxlfh,$tlxfh,$mapqfh);

my $statsfile = "${expt_stub}_stats.txt";
my $paramsfile = "${expt_stub}_params.txt";
my $dedup_output = "${expt_stub}_dedup.txt";


# Prepare breaksite hash
my $brksite = {chr => $brk_chr,
                start => $brk_start,
                end => $brk_end,
                strand => $brk_strand};

# Read in primer sequence
my $prim_io = Bio::SeqIO->new(-file => $prim_fa, -format => 'fasta');
$brksite->{primer} = $prim_io->next_seq();


# Add sequence to hash unless library has endogenous breaksite
if (defined $break_fa) {
  $brksite->{endogenous} = 0;
  my $brk_io = Bio::SeqIO->new(-file => $break_fa, -format => 'fasta');
  $brksite->{seq} = $brk_io->next_seq();
  $brksite->{len} = $brksite->{seq}->length;
  $brksite->{aln_name} = "Breaksite";
  $brksite->{aln_strand} = 1;
  $brksite->{primer_start} = index($brksite->{seq}->seq,$brksite->{primer}->seq)+1;
  $brksite->{breakcoord} = $breakcoord;
} else {
  $brksite->{endogenous} = 1;
  $brksite->{aln_name} = $brksite->{chr};
  $brksite->{aln_strand} = $brksite->{strand} eq "+" ? 1 : -1;
  $brksite->{primer_start} = $brksite->{strand} eq "+" ? $brksite->{start} : $brksite->{end} - 1;
  $brksite->{breakcoord} = $brksite->{strand} eq "+" ? $brksite->{end} - 1 : $brksite->{start};
}

# Calculate threshold for uncut/unjoin filter
$brksite->{uncut_threshold} = $brksite->{endogenous} ?
                                  ($brksite->{strand} eq "+" ?
                                    $brksite->{end} + $params->{max_bp_after_cutsite}:
                                    $brksite->{start} - $params->{max_bp_after_cutsite} - 1) :
                                  $brksite->{breakcoord} + $params->{max_bp_after_cutsite} + 1;

# Calculate thresholds for mispriming filters
croak "Error: could not find primer within breaksite sequence" unless $brksite->{primer_start} > 0;
$brksite->{misprimed_threshold} = $brksite->{primer_start} +
            $brksite->{aln_strand} * ($brksite->{primer}->length - 1 + $params->{min_bp_after_primer});

$params->{brksite} = $brksite;


# Read in adapter sequence
my $adpt_io = Bio::SeqIO->new(-file => $adapt_fa, -format => 'fasta');
my $adaptseq = $adpt_io->next_seq();

# Read in cutter sequence
my $cutseq;
if (defined $cut_fa) {
  my $cut_io = Bio::SeqIO->new(-file => $cut_fa, -format => 'fasta');
  $cutseq = $cut_io->next_seq();
}
$params->{cutter} = defined $cutseq ? $cutseq->seq : "";


write_parameters_file;



unless ($skip_alignment || $skip_process) {

  align_to_breaksite unless $brksite->{endogenous};
  align_to_genome unless -r $R1_bam && -r $R2_bam;
  align_to_adapter unless -r $R1_adpt_bam && -r $R2_adpt_bam;


} else {
  unless ($brksite->{endogenous}) {
    croak "Error: could not find breaksite alignment files when skipping align step" unless -r $R1_brk_bam && -r $R2_brk_bam;
  }
  croak "Error: could not find genome alignment files when skipping align step" unless -r $R1_bam && -r $R2_bam;
  croak "Error: could not find adapter alignment files when skipping align step" unless -r $R1_adpt_bam && -r $R2_adpt_bam;
}


unless ($skip_process) {
  $tlxlfh = IO::File->new(">$tlxlfile");
  $tlxfh = IO::File->new(">$tlxfile");
  $mapqfh = IO::File->new(">$mapqfile");

  $tlxlfh->print(join("\t", @tlxl_header)."\n");
  $tlxfh->print(join("\t", @tlx_header, @filters)."\n");

  my @mapq_header = qw(Qname Primary Read Rname Rstart Rend Strand Qstart Qend AS Cigar Overlap);
  push(@mapq_header,"Simulation") if defined $simfile;
  $mapqfh->print(join("\t", @mapq_header)."\n");
  $params->{mapqfh} = $mapqfh;
  $params->{mapq_header} = \@mapq_header;

  process_alignments;

  $tlxlfh->close;
  $tlxfh->close;

  mark_repeatseq_junctions;


  unless ($no_dedup) {

  #   mark_duplicate_junctions;

  }


} else {
  croak "Error: could not find tlx file when skipping processing step" unless -r $tlxfile;
}


filter_junctions;

# print stat results of filter

# sort_junctions;


# post_process_junctions;


my $t1 = tv_interval($t0);

printf("\nFinished all processes in %.2f seconds.\n", $t1);


# clean_up unless $no_clean;


#
# End of program
#


sub align_to_breaksite {
  print "\nRunning Bowtie2 alignment for $expt against breaksite sequence\n";

  System("bowtie2-build -q $break_fa $break_bt2idx");

  my $R1_brk_bt2_cmd = "bowtie2 $bt2_break_opt -x $break_bt2idx -U $read1 -S $R1_brk_sam";

  System($R1_brk_bt2_cmd);

  System("samtools view -bSh -o $R1_brk_bam $R1_brk_sam") unless sam_file_is_empty($R1_brk_sam);
  System("touch $R1_brk_bam",1);

  if (defined $read2) {
    my $R2_brk_bt2_cmd = "bowtie2 $bt2_break_opt -x $break_bt2idx -U $read2 -S $R2_brk_sam";

    System($R2_brk_bt2_cmd);

    System("samtools view -bSh -o $R2_brk_bam $R2_brk_sam") unless sam_file_is_empty($R2_brk_sam);
    System("touch $R2_brk_bam",1);
  }
}

sub align_to_adapter {
  print "\nRunning Bowtie2 alignment for $expt against adapter sequence\n";

  System("bowtie2-build -q $adapt_fa $adapt_bt2idx");

  my $R1_adpt_bt2_cmd = "bowtie2 $bt2_adapt_opt -x $adapt_bt2idx -U $read1 -S $R1_adpt_sam";

  System($R1_adpt_bt2_cmd);

  System("samtools view -bSh -o $R1_adpt_bam $R1_adpt_sam") unless sam_file_is_empty($R1_adpt_sam);
  System("touch $R1_adpt_bam",1);

  if (defined $read2) {
    my $R2_adpt_bt2_cmd = "bowtie2 $bt2_adapt_opt -x $adapt_bt2idx -U $read2 -S $R2_adpt_sam";

    System($R2_adpt_bt2_cmd);

    System("samtools view -bSh -o $R2_adpt_bam $R2_adpt_sam") unless sam_file_is_empty($R2_adpt_sam);
    System("touch $R2_adpt_bam",1);
  }
}

sub align_to_genome {

  print "\nRunning Bowtie2 alignment for $expt against $assembly genome \n";

  my $R1_bt2_cmd = "bowtie2 $bt2_opt -x $BOWTIE2_INDEXES/$assembly -U $read1 -S $R1_sam";

  System($R1_bt2_cmd);

  System("samtools view -bSh -o $R1_bam $R1_sam") unless sam_file_is_empty($R1_sam);
  System("touch $R1_bam",1);

  if (defined $read2) {
    my $R2_bt2_cmd = "bowtie2 $bt2_opt -x $BOWTIE2_INDEXES/$assembly -U $read2 -S $R2_sam";

    System($R2_bt2_cmd);

    System("samtools view -bSh -o $R2_bam $R2_sam") unless sam_file_is_empty($R2_sam);
    System("touch $R2_bam",1);
  }
  
}

sub process_alignments {

  debug_print("reading alignments",1,$expt);
  my ($R1_brk_iter,$R2_brk_iter);
  my ($R1_adpt_iter,$R2_adpt_iter);
  my ($R1_iter,$R2_iter);
  my ($simfh,$simcsv);

  my ($next_R1_brk_aln,$next_R2_brk_aln);
  my ($next_R1_adpt_aln,$next_R2_adpt_aln);
  my ($next_R1_aln,$next_R2_aln);
  my $next_sim_aln;

  unless ($brksite->{endogenous}) {
    $R1_brk_samobj = Bio::DB::Sam->new(-bam => $R1_brk_bam,
                                       -fasta => $break_fa,
                                       -expand_flags => 1);

    $R1_brk_iter = $R1_brk_samobj->get_seq_stream();
    $next_R1_brk_aln = wrap_alignment("R1",$R1_brk_iter->next_seq);
  }


  $R1_adpt_samobj = Bio::DB::Sam->new(-bam => $R1_adpt_bam,
                                      -fasta => $adapt_fa,
                                      -expand_flags => 1);
  $R1_adpt_iter = $R1_adpt_samobj->get_seq_stream();
  $next_R1_adpt_aln = wrap_alignment("R1",$R1_adpt_iter->next_seq);




  $R1_samobj = Bio::DB::Sam->new(-bam => $R1_bam,
                                 -fasta => $assembly_fa,
                                 -expand_flags => 1);

  $R1_iter = $R1_samobj->get_seq_stream();
  $next_R1_aln = wrap_alignment("R1",$R1_iter->next_seq);

  if (defined $simfile) {
    $simfh = IO::File->new("<$simfile");
    $simcsv = Text::CSV->new({sep_char => "\t"});
    $simcsv->column_names(qw(Qname Rname Junction Strand Length));
  }

  if (defined $read2) {

    unless ($brksite->{endogenous}) {
      $R2_brk_samobj = Bio::DB::Sam->new(-bam => $R2_brk_bam,
                                         -fasta => $break_fa,
                                         -expand_flags => 1);
    $R2_brk_iter = $R2_brk_samobj->get_seq_stream();
    $next_R2_brk_aln = wrap_alignment("R2",$R2_brk_iter->next_seq);

    }

    $R2_adpt_samobj = Bio::DB::Sam->new(-bam => $R2_adpt_bam,
                                        -fasta => $adapt_fa,
                                        -expand_flags => 1);
    $R2_adpt_iter = $R2_adpt_samobj->get_seq_stream();
    $next_R2_adpt_aln = wrap_alignment("R2",$R2_adpt_iter->next_seq);
    
    

    $R2_samobj = Bio::DB::Sam->new(-bam => $R2_bam,
                                   -fasta => $assembly_fa,
                                   -expand_flags => 1);

    
    $R2_iter = $R2_samobj->get_seq_stream();
    $next_R2_aln = wrap_alignment("R2",$R2_iter->next_seq);
  
  }

  if ($brksite->{endogenous}) {
    my $primer_check = $R1_samobj->seq($brksite->{chr},$brksite->{start},$brksite->{end}-1);
    $primer_check = reverseComplement($primer_check) if $brksite->{strand} eq "-";
    print $brksite->{primer}->seq ." $primer_check\n"; 
    croak "Error: primer sequence does not match reference genome"
        unless $brksite->{primer}->seq eq substr(uc($primer_check),0,$brksite->{primer}->length);
  }

  croak "Error: no R1 alignments" unless defined $next_R1_aln;


  while (defined $next_R1_aln) {


    # The strategy here is to read in all alignments from the genome, adapter,
    # and breaksite (for non-endogenous libraries) and pool them all into
    # an array for read 1, and read 2 if the library is paired-end

    my $qname = $next_R1_aln->{Qname};

    debug_print("initializing read - reading alignments",2,$qname);

    my %R1_alns_h;
    my %R2_alns_h;

    if (defined $simfile) {
      $next_sim_aln = $simcsv->getline_hr($simfh);
      $next_R1_aln->{Simulation} = alignment_matches_simulation($next_R1_aln,$next_sim_aln);
    } 

    $R1_alns_h{$next_R1_aln->{ID}} = $next_R1_aln;
    undef $next_R1_aln;

    while(my $aln = $R1_iter->next_seq) {
      $next_R1_aln = wrap_alignment("R1",$aln);
      last unless $next_R1_aln->{Qname} eq $qname;
      if (defined $simfile) {
        $next_R1_aln->{Simulation} = alignment_matches_simulation($next_R1_aln,$next_sim_aln);
      } 
      $R1_alns_h{$next_R1_aln->{ID}} = $next_R1_aln;
      undef $next_R1_aln;
    }

    if (defined $next_R2_aln && $next_R2_aln->{Qname} eq $qname) {

      if (defined $simfile) {
        $next_R2_aln->{Simulation} = alignment_matches_simulation($next_R2_aln,$next_sim_aln);
      } 

      $R2_alns_h{$next_R2_aln->{ID}} = $next_R2_aln;
      undef $next_R2_aln;
      while(my $aln = $R2_iter->next_seq) {
        $next_R2_aln = wrap_alignment("R2",$aln);
        last unless $next_R2_aln->{Qname} eq $qname;
        if (defined $simfile) {
          $next_R2_aln->{Simulation} = alignment_matches_simulation($next_R2_aln,$next_sim_aln);
        } 
        $R2_alns_h{$next_R2_aln->{ID}} = $next_R2_aln;
        undef $next_R2_aln;
      }
    }

    # Read in Breaksite alignments only if brksite is non-endogenous
    unless ($brksite->{endogenous}) {
      if (defined $next_R1_brk_aln && $next_R1_brk_aln->{Qname} eq $qname) {
        $R1_alns_h{$next_R1_brk_aln->{ID}} = $next_R1_brk_aln;
        undef $next_R1_brk_aln;
        while(my $aln = $R1_brk_iter->next_seq) {
          $next_R1_brk_aln = wrap_alignment("R1",$aln);
          last unless $next_R1_brk_aln->{Qname} eq $qname;
          $R1_alns_h{$next_R1_brk_aln->{ID}} = $next_R1_brk_aln;
          undef $next_R1_brk_aln;
        }
      }

      if (defined $next_R2_brk_aln && $next_R2_brk_aln->{Qname} eq $qname) {
        $R2_alns_h{$next_R2_brk_aln->{ID}} = $next_R2_brk_aln;
        undef $next_R2_brk_aln;
        while(my $aln = $R2_brk_iter->next_seq) {
          $next_R2_brk_aln = wrap_alignment("R2",$aln);
          last unless $next_R2_brk_aln->{Qname} eq $qname;
          $R2_alns_h{$next_R2_brk_aln->{ID}} = $next_R2_brk_aln;
          undef $next_R2_brk_aln;
        }
      }
    }


    # Read in Adapter alignments
    if (defined $next_R1_adpt_aln && $next_R1_adpt_aln->{Qname} eq $qname) {
      $R1_alns_h{$next_R1_adpt_aln->{ID}} = $next_R1_adpt_aln;
      undef $next_R1_adpt_aln;
      while(my $aln = $R1_adpt_iter->next_seq) {
        $next_R1_adpt_aln = wrap_alignment("R1",$aln);
        last unless $next_R1_adpt_aln->{Qname} eq $qname;
        $R1_alns_h{$next_R1_adpt_aln->{ID}} = $next_R1_adpt_aln;
        undef $next_R1_adpt_aln;
      }
    }

    if (defined $next_R2_adpt_aln && $next_R2_adpt_aln->{Qname} eq $qname) {
      $R2_alns_h{$next_R2_adpt_aln->{ID}} = $next_R2_adpt_aln;
      undef $next_R2_adpt_aln;
      while(my $aln = $R2_adpt_iter->next_seq) {
        $next_R2_adpt_aln = wrap_alignment("R2",$aln);
        last unless $next_R2_adpt_aln->{Qname} eq $qname;
        $R2_alns_h{$next_R2_adpt_aln->{ID}} = $next_R2_adpt_aln;
        undef $next_R2_adpt_aln;
      }
    }

    my $OCS = find_optimal_coverage_set(\%R1_alns_h,\%R2_alns_h);

    debug_print("processing OCS",2,$OCS->[0]->{R1}->{QnameShort});

    my $tlxls = create_tlxl_entries($OCS);

    my $tlxs = create_tlx_entries($tlxls, {genome => $R1_samobj,
                                           brk => $R1_brk_samobj,
                                           adpt => $R1_adpt_samobj} )  ;

    foreach my $tlx (@$tlxs) {
      my %filter_init;
      @filter_init{@filters} = (0) x @filters;
      $tlx->{filters} = {%filter_init};
    }


    # Bundle alignments and everything into a single object
    my $read_obj = {tlxs => $tlxs, tlxls => $tlxls, R1_alns => \%R1_alns_h, R2_alns => \%R2_alns_h};

    find_random_barcode($read_obj,$random_barcode);

    foreach my $filter (@dispatch_names) {
      $filter_dispatch{$filter}->($read_obj);
    }


    foreach my $tlxl (@{$read_obj->{tlxls}}) {
      write_entry($tlxlfh,$tlxl,\@tlxl_header);
    }
    foreach my $tlx (@{$read_obj->{tlxs}}) {  
      write_filter_entry($tlxfh,$tlx,\@tlx_header,\@filters);        
    }

  }

}

sub mark_repeatseq_junctions {

  (my $repeatseq_output = $tlxfile) =~ s/.tlx$/_repeatseq.tlx/;

  $repeatseq_bedfile = "$GENOME_DB/$assembly/annotation/repeatSeq.bed" unless defined $repeatseq_bedfile;

  return unless -r $repeatseq_bedfile;

  my $repeatseq_cmd = join(" ","$FindBin::Bin/TranslocRepeatSeq.pl",
                          $tlxfile,
                          $repeatseq_bedfile,
                          $repeatseq_output) ;

  System($repeatseq_cmd);

  rename $repeatseq_output, $tlxfile;

  return;

}

sub mark_duplicate_junctions {
  (my $duplicate_output = $tlxfile) =~ s/.tlx$/_duplicate.tlx.tmp/;

  my $duplicate_cmd = join(" ","$FindBin::Bin/../R/TranslocDedup.pl",
                          $tlxfile,
                          $duplicate_output,
                          "--cores $threads",
                          "--offset_dist",$params->{dedup_offset_dist},
                          "--break_dist",$params->{dedup_break_dist}) ;

  rename $duplicate_output, $tlxfile;

  return;

}

sub filter_junctions {
  (my $filter_output = $tlxfile) =~ s/.tlx$/_result.tlx/;

  my $filter_cmd = join(" ","$FindBin::Bin/TranslocFilter.pl",
                          $tlxfile,
                          $filter_output,
                          "--filters",
                          "\"f.unaligned=1",
                          "f.baitonly=1",
                          "f.uncut=1",
                          "f.misprimed=L".$params->{min_bp_after_primer},
                          "f.freqcut=1",
                          "f.largegap=G".$params->{max_largegap},
                          "f.mapqual=L".$params->{mapq_score_thresh},
                          "f.breaksite=1",
                          "f.sequential=1",
                          "f.repeatseq=1",
                          "f.duplicate=1\"") ;

  System($filter_cmd);

  return;
}

sub sort_junctions {
  my $tlxbak = "$tlxfile.bak";

  rename $tlxfile, $tlxbak;

  my @junctions;

  my $bak_tlxfh = IO::File->new("<$tlxbak");
  my $tlxfh = IO::File->new(">$tlxfile");

  my $csv = Text::CSV->new({sep_char => "\t"});
  my $header = $csv->getline($bak_tlxfh);
  $csv->column_names(@$header);
  $tlxfh->print(join("\t", @tlx_header)."\n");
  
  while (my $tlx = $csv->getline_hr($bak_tlxfh)) {
    push(@junctions, $tlx);
  }

  @junctions = sort {$a->{Rname} cmp $b->{Rname} || $a->{Junction} <=> $b->{Junction}} @junctions;

  foreach my $tlx (@junctions) {
    write_entry($tlxfh,$tlx,\@tlx_header);
  }

  unlink $tlxbak;

}

sub post_process_junctions {

  (my $html_reads = $tlxfile) =~ s/tlx$/html/;
  System("TranslocHTMLReads.pl $tlxfile $html_reads --primer ".$brksite->{primer}->seq." --adapter ".$adaptseq->seq);  

  if ($assembly =~ /^(hg\d+|mm\d+)$/) {

    (my $pdf_plot = $tlxfile) =~ s/tlx$/pdf/;

    System("$FindBin::Bin/../R/TranslocPlot.R $tlxfile $pdf_plot binsize=2000000 strand=2 assembly=$assembly " .
            "brkchr=$brk_chr brksite=" .($brk_strand eq "+" ? $brk_end : $brk_start) ." brkstrand=" . ($brk_strand eq "+" ? "1" : "-1") );


    (my $pdf_plot_breaksite = $pdf_plot) =~ s/\.pdf/_brksite.pdf/;

    System("$FindBin::Bin/../R/TranslocPlot.R $tlxfile $pdf_plot_breaksite strand=0 assembly=$assembly " .
          "brkchr=$brk_chr brksite=" .($brk_strand eq "+" ? $brk_end : $brk_start) ." brkstrand=" . ($brk_strand eq "+" ? "1" : "-1") .
          " chr=$brk_chr rmid=" .($brk_strand eq "+" ? $brk_end : $brk_start) ." rwindow=5000 binnum=100 plottype=linear" );
  }
}


sub write_parameters_file {
  my $paramsfh = IO::File->new(">$paramsfile");

  $paramsfh->print("$local_time\t$commandline\n\n");
  $paramsfh->print(Dumper($params)."\n");

  # $paramsfh->print(join("\n", map {join("\t",@$_)}
  #   (["Alignment Options"],
  #   ["Genome Bowtie2 Options", $bt2_opt],
  #   ["Breaksite Bowtie2 Options", $bt2_break_opt],
  #   ["Adapter Bowtie2 Options", $bt2_adapt_opt],
  #   [],
  #   ["Optimal Query Coverage Options"],
  #   ["Break Penalty",$Brk_pen_default],
  #   ["PE Gap Penalty",$PE_pen_default],
  #   ["Overlap Penalty",$OL_mult],
  #   ["Max Bait Offset",$maximum_brk_start_dif],
  #   ["Bait Offset Penalty",$Dif_mult],
  #   [],
  #   ["Filtering Options"],
  #   ["Max Bp After Cutsite", $max_bases_after_cutsite],
  #   ["Min Bp After Primer", $min_bases_after_primer],
  #   ["MapQuality Overlap Threshold", $mapq_ol_thresh],
  #   ["MapQuality Mismatch Threshold Intercept", $mapq_mismatch_thresh_int],
  #   ["MapQuality Mismatch Threshold Coefficient", $mapq_mismatch_thresh_coef],
  #   [],
  #   ["Dedup Options"],
  #   ["Max Prey Offset Distance",$dedup_offset_dist],
  #   ["Max Bait Junction Distance",$dedup_break_dist])));

}


sub clean_up {

  debug_print("Cleaning up extra files",1,$expt);
  unlink glob "${expt_stub}*.sam";

}

sub parse_command_line {
	my $help;

	usage() if (scalar @ARGV == 0);

  $commandline = join(" ",@ARGV);
  $local_time = localtime;

	my $result = GetOptions ( "read1=s" => \$read1,
                            "read2=s" => \$read2,
                            "workdir=s" => \$workdir,
                            "assembly=s" => \$assembly,
                            "chr=s" => \$brk_chr,
                            "start=i" => \$brk_start,
                            "end=i" => \$brk_end,
                            "strand=s" => \$brk_strand,
                            "mid=s" => \$mid_fa,
                            "primer=s" => \$prim_fa,
                            "breakseq=s" => \$break_fa,
                            "breaksite=i" => \$breakcoord,
                            "adapter=s" => \$adapt_fa,
                            "cutter=s" => \$cut_fa,
                            "threads=i" => \$threads,
                            "random-barcode=i" => \$random_barcode,
                            "skip-align" => \$skip_alignment,
                            "skip-process" => \$skip_process,
                            "skip-dedup" => \$skip_dedup,
                            "no-dedup" => \$no_dedup,
                            "no-clean" => \$no_clean,
                            "force-bait=i" => \$params->{force_bait},
                            "match-award=i" => \$params->{match_award},
                            "mismatch-pen=s" => \$params->{mismatch_pen},
                            "n-base-pen=i" => \$params->{n_base_pen},
                            "read-gap-pen=s" => \$params->{read_gap_pen},
                            "ref-gap-pen=s" => \$params->{ref_gap_pen},
                            "score-min=s" => \$params->{score_min},
                            "D-effort=i" => \$params->{D_effort},
                            "R-effort=i" => \$params->{R_effort},
                            "seed-mismatch=i" => \$params->{seed_mismatch},
                            "seed-length=i" => \$params->{seed_length},
                            "seed-interval=s" => \$params->{seed_interval},
                            "breaksite-alignments=i" => \$params->{breaksite_alignments},
                            "genome-alignments=i" => \$params->{genome_alignments},
                            "max-brkstart-dif=i" => \$params->{max_brkstart_dif},
                            "break-pen=i" => \$params->{brk_pen},
                            "max-pe-gap=i" => \$params->{max_pe_gap},
                            "pe-pen=i" => \$params->{pe_pen},
                            "max-dovetail=i" => \$params->{max_dovetail},
                            "max-uncut-bp=i" => \$params->{max_bp_after_cutsite},
                            "min-priming-bp=i" => \$params->{min_bp_after_primer},
                            "max-largegap=i" => \$params->{max_largegap},
                            "mapq-ol=f" => \$params->{mapq_ol_thresh},
                            "mapq-score=i" => \$params->{mapq_score_thresh},
                            "repeatseq-bed=s" => \$repeatseq_bedfile,
                            "dedup-offset-bp=i" => \$params->{dedup_offset_dist},
                            "dedup-bait-bp=i" => \$params->{dedup_break_dist},
                            "debug=i" => \$debug_level,
                            "simfile=s" => \$simfile,
				            				"help" => \$help
				            			) ;

	
	usage() if ($help);

  #Check options

  croak "Error: must specify --read1" unless defined $read1;
  croak "Error: cannot read read1 file" unless -r $read1;
  croak "Error: cannot read read2 file" if defined $read2 && ! -r $read2;
  croak "Error: working directory does not exist" unless (-d $workdir);

  croak "Error: priming-bp must be a positive integer" unless $params->{min_bp_after_primer} > 0;
  croak "Error: mapq-ol must be a fraction between 0 and 1" if $params->{mapq_ol_thresh} < 0 || $params->{mapq_ol_thresh} > 1;
  
  if ($params->{mismatch_pen} =~ /^(\d+),(\d+)$/) {
    $params->{max_mismatch_pen} = $1;
    $params->{min_mismatch_pen} = $2;
  } else {
    croak "Error: mismatch-pen not input in correct format 'max,min'";
  }

  if ($params->{read_gap_pen} =~ /^(\d+),(\d+)$/) {
    $params->{read_gap_open} = $1;
    $params->{read_gap_ext} = $2;
  } else {
    croak "Error: read-gap-pen not input in correct format 'open,ext'";
  }

  if ($params->{ref_gap_pen} =~ /^(\d+),(\d+)$/) {
    $params->{ref_gap_open} = $1;
    $params->{ref_gap_ext} = $2;
  } else {
    croak "Error: ref-gap-pen not input in correct format 'open,ext'";
  }
  
  


                            "max-mismatch-pen=i" => \$params->{max_mismatch_pen},
                            "min-mismatch-pen=i" => \$params->{min_mismatch_pen},
	exit unless $result;
}


sub usage()
{
print<<EOF;
TranslocPipeline.pl, by Robin Meyers, 2013


Usage: $0
        --read1 FILE --read2 FILE --workdir DIR --assembly (mm9|hg19|other)
        --chr chrN --start INT --end INT --strand (+|-) --MID FILE
        --primer FILE [--breaksite FILE] --adapter FILE --cutter FILE
        [--option VAL] [--flag] [--help]

Arguments set automatically by TranslocWrapper.pl:
$arg{"--read1","Input sequence file"}
$arg{"--read2","Input sequence file"}
$arg{"--workdir","Directory for results files"}
$arg{"--assembly","Genome assembly to align reads to"}
$arg{"--chr","Breaksite chromosome"}
$arg{"--start","Breaksite start"}
$arg{"--end","Breaksite end"}
$arg{"--strand","Breaksite orientation"}
$arg{"--mid","MID sequence"}
$arg{"--primer","Primer sequence"}
$arg{"--breakseq","Breaksite cassette sequence (non-endogenous only)"}
$arg{"--breaksite","Cutting coordinate on cassette (non-endogenous only)"}
$arg{"--adapter","Adapter sequence"}
$arg{"--cutter","Frequent cutter site sequence"}

All other arguments must be set manually with --pipeline-opt in TranslocWrapper.pl (defaults in parentheses):

$arg{"--threads","Number of threads to run bowtie/dedup on",$threads}
$arg{"--random-barcode","Set to integer N to find N-mer immediately before adapter"}
$arg{"--skip-align","Begin pipeline after alignment step - working directory must already have alignment files "}
$arg{"--skip-process","Begin pipeline after OCS and filtering steps - working directory must already have tlx files"}
$arg{"--skip-dedup","Begin pipeline after dedup step (post-processing only) "}
$arg{"--no-dedup","Do not run dedup filter"}
$arg{"--no-clean","Do not delete temp files at end of process"}
$arg{"--force-bait","Set to 0 to relax bait alignment restrictions",$params->{force_bait}}
$arg{"--simfile"," "}

Arguments sent to Bowtie2 alignment (see manual)
$arg{"--match-award","",$params->{match_award}}
$arg{"--mismatch-pen","",$params->{mismatch_pen}}
$arg{"--n-base-pen","",$params->{n_base_pen}}
$arg{"--read-gap-pen","",$params->{read_gap_pen}}
$arg{"--ref-gap-pen","",$params->{ref_gap_pen}}
$arg{"--score-min","",$params->{score_min}}
$arg{"--D-effort","",$params->{D_effort}}
$arg{"--R-effort","",$params->{R_effort}}
$arg{"--seed-mismatch","",$params->{seed_mismatch}}
$arg{"--seed-length","",$params->{seed_length}}
$arg{"--seed-interval","",$params->{seed_interval}}
$arg{"--breaksite-alignments","Number of breaksite alignments reported",$params->{breaksite_alignments}}
$arg{"--genome-alignments","Number of genome alignments reported",$params->{genome_alignments}}

OQC and Filtering parameters
$arg{"--max-brkstart-dif","Bait alignment forced to be within this distance of primer start site - not used if force-bait=0",$params->{max_brkstart_dif}}
$arg{"--max-overlap","Maximum OQC overlap between two segments",$params->{max_overlap}}
$arg{"--break-pen","Junction penalty incurred for every breakpoint in OQC",$params->{brk_pen}}
$arg{"--max-pe-gap","Maximum gap allowed between paired-end alignments to be considered concordant",$params->{max_pe_gap}}
$arg{"--pe-pen","Penalty incurred in OCS by gapped paired-end alignments (scaled to max-gap)",$params->{pe_pen}}
$arg{"--max-dovetail","Maximum dovetail allowed between concordant alignments",$params->{max_dovetail}}
$arg{"--max-uncut-bp","Maximum number of bases allowed to align after cutsite",$params->{max_bp_after_cutsite}}
$arg{"--min-priming-bp","Minimum number of bases in bait alignment after the primer",$params->{min_bp_after_primer}}
$arg{"--max-largegap","Maximum distance on query between bait and prey alignments",$params->{max_largegap}}
$arg{"--mapq-ol","Minimum overlapping fraction for mapq filter",$params->{mapq_ol_thresh}}
$arg{"--mapq-score","Minimum mapping quality score",$params->{mapq_score_thresh}}
$arg{"--repeatseq-bed","Location of repeatseq bedfile (automatically found for default genome assemblies)"}
$arg{"--dedup-offset-bp","Minimum offset distance between prey alignments in dedup filter",$params->{dedup_offset_dist}}
$arg{"--dedup-bait-bp","Minimum distance between bait alignments in dedup filter",$params->{dedup_break_dist}}

$arg{"--debug","Set to level 1-4 for increasing verbosity in log file",$debug_level}
$arg{"--help","This helpful help screen."}

EOF

exit 1;
}
