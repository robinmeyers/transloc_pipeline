#!/usr/bin/perl


use strict;
use warnings;
use Getopt::Long;
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

use List::Util qw(min max);
use Interpolation 'arg:@->$' => \&argument;
use Time::HiRes qw(gettimeofday tv_interval);
use Data::Dumper;

use threads;
use threads::shared;
use Thread::Queue;


use Cwd qw(abs_path);
use FindBin;
use lib abs_path("$FindBin::Bin/../lib");

require "PerlSub.pl";
require "TranslocHelper.pl";
require "TranslocFilters.pl";


# Flush output after every write
select( (select(STDOUT), $| = 1 )[0] );

##
## This program
## 
## run with "--help" for usage information
##
## Robin Meyers

# Forward declarations
sub parse_command_line;
sub align_to_breaksite;
sub align_to_adapter;
sub align_to_genome;
sub process_alignments;
sub find_optimal_coverage_set ($$);
sub process_optimal_coverage_set ($$$);
sub score_edge ($;$);
sub deduplicate_junctions;
sub post_process_junctions;



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
my $max_threads = 4;
my $max_queue_size = 5;
my $qblocksize = 500;

my $priming_bp = 12;

my $OL_mult = 2;
my $Dif_mult = 1;
my $Brk_pen_min = 10;
my $Brk_pen_power = 4;
my $Brk_pen_max = 60;
my $Brk_dist_max = 100000000;
my $Brk_pen_mult = ($Brk_pen_max-$Brk_pen_min)/(log10($Brk_dist_max)^$Brk_pen_power);

my $max_frag_len = 1500;
my $ol_thresh = 0.9;
my $score_thresh = 0.9;

my $user_bowtie_opt = "";
my $user_bowtie_adapter_opt = "";
my $user_bowtie_breaksite_opt = "";

my $use_current_tlx;

# Global variables
my @tlxl_header = tlxl_header();
my @tlx_header = tlx_header();
my @tlx_filter_header = tlx_filter_header();




my %stats :shared;
my $stats = \%stats;
share($stats);
$stats->{totalreads} = 0;
$stats->{aligned} = 0;
$stats->{junctions} = 0;
$stats->{junc_reads} = 0;
$stats->{mapqual} = 0;
$stats->{mapq_reads} = 0;
$stats->{priming} = 0;
$stats->{prim_reads} = 0;
$stats->{freqcut} = 0;
$stats->{freq_reads} = 0;
$stats->{breaksite} = 0;
$stats->{break_reads} = 0;
$stats->{splitjuncs} = 0;
$stats->{split_reads} = 0;
$stats->{dedup} = 0;
$stats->{final} = 0;



#
# Start of Program
#
parse_command_line;


# my $default_bowtie_adapter_opt = "--local -D 20 -R 3 -N 1 -L 6 -i C,4 --score-min C,30 -p $max_threads --no-unal --reorder -t";
# my $default_bowtie_breaksite_opt = "--local -D 20 -R 3 -N 1 -L 12 -i C,6 --score-min C,40 --mp 10,2 --rfg 10,10 --rdg 10,10 -p $max_threads -k 50 --reorder -t";
# my $default_bowtie_opt = "--local -D 20 -R 3 -N 0 -L 20 -i C,8 --score-min C,50 --mp 10,2 --rfg 10,2 --rdg 10,2 -p $max_threads -k 50 --no-unal --reorder -t";

my $default_bowtie_adapter_opt = "--local -D 20 -R 3 -N 1 -L 6 -i C,4 --score-min C,30 -p $max_threads --no-unal --reorder -t";
my $default_bowtie_breaksite_opt = "--local -D 20 -R 3 -N 1 -L 12 -i C,6 --score-min C,40 --mp 10,2 --rfg 10,2 --rdg 10,2 -p $max_threads -k 50 --reorder -t";
my $default_bowtie_opt = "--local -D 20 -R 3 -N 0 -L 20 -i C,8 --score-min C,50 --mp 10,2 --rfg 10,2 --rdg 10,2 -p $max_threads -k 50 --no-unal --reorder -t";


my $bt2_break_opt = manage_program_options($default_bowtie_breaksite_opt,$user_bowtie_breaksite_opt);
my $bt2_adapt_opt = manage_program_options($default_bowtie_adapter_opt,$user_bowtie_adapter_opt);
my $bt2_opt = manage_program_options($default_bowtie_opt,$user_bowtie_opt);

my $t0 = [gettimeofday];

#check_working_dir;

my $expt = basename($workdir);
my $expt_stub = "$workdir/$expt";


my $assembly_fa = $ENV{'GENOME_DB'}."/$assembly/$assembly.fa";
my $break_fa = "$workdir/misc/breaksite.fa";
my $prim_fa = "$workdir/misc/primer.fa";
my $adapt_fa = "$workdir/misc/adapter.fa";
my $cut_fa = "$workdir/misc/cutter.fa";
my $break_bt2idx = "$workdir/misc/breaksite";
my $adapt_bt2idx = "$workdir/misc/adapter";

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
my $filt_tlxfile = "${expt_stub}_filtered.tlx";
my $unjoin_tlxfile = "${expt_stub}_unjoined.tlx";
my ($tlxlfh,$tlxfh,$filt_tlxfh,$unjoin_tlxfh);

my $dedup_output = "${expt_stub}_dedup.txt";


my $brk_io = Bio::SeqIO->new(-file => $break_fa, -format => 'fasta');
my $breakseq = $brk_io->next_seq();
my $breaklen = length($breakseq->seq);

my $prim_io = Bio::SeqIO->new(-file => $prim_fa, -format => 'fasta');
my $primseq = $prim_io->next_seq();
my $primlen = length($primseq->seq);
my $priming_threshold = $primlen + $priming_bp;
my $primer_start = index($breakseq->seq,$primseq->seq)+1;
croak "Error: could not find primer within breaksite sequence" unless $primer_start > 0;

my $adpt_io = Bio::SeqIO->new(-file => $adapt_fa, -format => 'fasta');
my $adaptseq = $adpt_io->next_seq();
my $adaptlen = length($adaptseq->seq);

my $cut_io = Bio::SeqIO->new(-file => $cut_fa, -format => 'fasta');
my $cutseq = $cut_io->next_seq();

my $brk_hash = {chr => $brk_chr,
                start => $brk_start,
                end => $brk_end,
                strand => $brk_strand,
                len => $breaklen};



align_to_breaksite unless -r $R1_brk_bam && -r $R2_brk_bam;

align_to_adapter unless -r $R1_adpt_bam && -r $R2_adpt_bam;

align_to_genome unless -r $R1_bam && -r $R2_bam;


# my $aln_queue = Thread::Queue->new();
# my $ocs_queue = Thread::Queue->new();

# my $aln_reader_thr = threads->create('read_alignments');

# my @ocs_finder_thrs = ();
# foreach my $i (1..($max_threads-2)) {
#   my $thr = threads->create('find_optimal_coverage_set');
#   push(@ocs_finder_thrs,$thr);
# }

# # my $ocs_processer_thr = threads->create('process_OCS');

# $aln_reader_thr->join;
# foreach my $thr (@ocs_finder_thrs) {
#   $thr->join;
# }
# $ocs_processer_thr->join;

$tlxlfh = IO::File->new(">$tlxlfile");
$tlxfh = IO::File->new(">$tlxfile");
$filt_tlxfh = IO::File->new(">$filt_tlxfile");
$unjoin_tlxfh = IO::File->new(">$unjoin_tlxfile");

$tlxlfh->print(join("\t", @tlxl_header)."\n");
$tlxfh->print(join("\t", @tlx_header)."\n");
$filt_tlxfh->print(join("\t", @tlx_filter_header)."\n");
$unjoin_tlxfh->print(join("\t", @tlx_filter_header)."\n");

process_alignments;

$tlxlfh->close;
$tlxfh->close;
$filt_tlxfh->close;
$unjoin_tlxfh->close;

deduplicate_junctions;

post_process_junctions;

my $t1 = tv_interval($t0);

printf("\nFinished all processes in %.2f seconds.\n", $t1);

print("Some stats:\n".join("\t",$stats->{totalreads},
                                $stats->{aligned},
                                $stats->{junctions}." (".$stats->{junc_reads}.")",
                                $stats->{mapqual}." (".$stats->{mapq_reads}.")",
                                $stats->{priming}." (".$stats->{prim_reads}.")",
                                $stats->{freqcut}." (".$stats->{freq_reads}.")",
                                $stats->{breaksite}." (".$stats->{break_reads}.")",
                                $stats->{splitjuncs}." (".$stats->{split_reads}.")",
                                $stats->{dedup})."\n");



#
# End of program
#

sub align_to_breaksite {
  print "\nRunning Bowtie2 alignment for $expt against breaksite sequence\n";

  System("bowtie2-build -q $break_fa $break_bt2idx");

  my $R1_brk_bt2_cmd = "bowtie2 $bt2_break_opt -x $break_bt2idx -U $read1 -S $R1_brk_sam";

  System($R1_brk_bt2_cmd);

  System("samtools view -bS -o $R1_brk_bam $R1_brk_sam") unless sam_file_is_empty($R1_brk_sam);
  System("touch $R1_brk_bam",1);

  my $R2_brk_bt2_cmd = "bowtie2 $bt2_break_opt -x $break_bt2idx -U $read2 -S $R2_brk_sam";

  System($R2_brk_bt2_cmd);

  System("samtools view -bS -o $R2_brk_bam $R2_brk_sam") unless sam_file_is_empty($R2_brk_sam);
  System("touch $R2_brk_bam",1);

}

sub align_to_adapter {
  print "\nRunning Bowtie2 alignment for $expt against adapter sequence\n";

  System("bowtie2-build -q $adapt_fa $adapt_bt2idx");

  my $R1_adpt_bt2_cmd = "bowtie2 $bt2_adapt_opt -x $adapt_bt2idx -U $read1 -S $R1_adpt_sam";

  System($R1_adpt_bt2_cmd);

  System("samtools view -bS -o $R1_adpt_bam $R1_adpt_sam") unless sam_file_is_empty($R1_adpt_sam);
  System("touch $R1_adpt_bam",1);

  my $R2_adpt_bt2_cmd = "bowtie2 $bt2_adapt_opt -x $adapt_bt2idx -U $read2 -S $R2_adpt_sam";

  System($R2_adpt_bt2_cmd);

  System("samtools view -bS -o $R2_adpt_bam $R2_adpt_sam") unless sam_file_is_empty($R2_adpt_sam);
  System("touch $R2_adpt_bam",1);

}

sub align_to_genome {

  print "\nRunning Bowtie2 alignment for $expt against $assembly genome \n";

  my $R1_bt2_cmd = "bowtie2 $bt2_opt -x $assembly -U $read1 -S $R1_sam";

  System($R1_bt2_cmd);

  System("samtools view -bS -o $R1_bam $R1_sam") unless sam_file_is_empty($R1_sam);
  System("touch $R1_bam",1);

  my $R2_bt2_cmd = "bowtie2 $bt2_opt -x $assembly -U $read2 -S $R2_sam";

  System($R2_bt2_cmd);

  System("samtools view -bS -o $R2_bam $R2_sam") unless sam_file_is_empty($R2_sam);
  System("touch $R2_bam",1);

  
}

sub process_alignments {

  print "\nReading alignments\n";


  $R1_brk_samobj = Bio::DB::Sam->new(-bam => $R1_brk_bam,
                                     -fasta => $break_fa,
                                     -expand_flags => 1);

  $R2_brk_samobj = Bio::DB::Sam->new(-bam => $R2_brk_bam,
                                     -fasta => $break_fa,
                                     -expand_flags => 1);

  $R1_adpt_samobj = Bio::DB::Sam->new(-bam => $R1_adpt_bam,
                                     -fasta => $adapt_fa,
                                     -expand_flags => 1);

  $R2_adpt_samobj = Bio::DB::Sam->new(-bam => $R2_adpt_bam,
                                     -fasta => $adapt_fa,
                                     -expand_flags => 1);

  $R1_samobj = Bio::DB::Sam->new(-bam => $R1_bam,
                                 -fasta => $assembly_fa,
                                 -expand_flags => 1);

  $R2_samobj = Bio::DB::Sam->new(-bam => $R2_bam,
                                 -fasta => $assembly_fa,
                                 -expand_flags => 1);

  my $R1_brk_iter = $R1_brk_samobj->get_seq_stream();
  my $R2_brk_iter = $R2_brk_samobj->get_seq_stream();
  my $R1_adpt_iter = $R1_adpt_samobj->get_seq_stream();
  my $R2_adpt_iter = $R2_adpt_samobj->get_seq_stream();
  my $R1_iter = $R1_samobj->get_seq_stream();
  my $R2_iter = $R2_samobj->get_seq_stream();


  my ($next_R1_brk_aln,$next_R2_brk_aln);# :shared;
  my ($next_R1_adpt_aln,$next_R2_adpt_aln);# :shared;
  my ($next_R1_aln,$next_R2_aln);# :shared;

  $next_R1_brk_aln = wrap_alignment("R1",$R1_brk_iter->next_seq);
  $next_R2_brk_aln = wrap_alignment("R2",$R2_brk_iter->next_seq);
  $next_R1_adpt_aln = wrap_alignment("R1",$R1_adpt_iter->next_seq);
  $next_R2_adpt_aln = wrap_alignment("R2",$R2_adpt_iter->next_seq);
  $next_R1_aln = wrap_alignment("R1",$R1_iter->next_seq);
  $next_R2_aln = wrap_alignment("R2",$R2_iter->next_seq);

  croak "Error: no R1 alignments to breaksite" unless defined $next_R1_brk_aln;


  while (defined $next_R1_brk_aln) {

    # my @aln_block :shared;

    # for (my $i=0; $i < $qblocksize; $i++) {

    #   last unless defined $next_R1_brk_aln;

    my $qname = $next_R1_brk_aln->{Qname};

    my @R1_alns;# :shared;
    my @R2_alns;# :shared;

       
    push(@R1_alns,$next_R1_brk_aln);
    undef $next_R1_brk_aln;

    while(my $aln = $R1_brk_iter->next_seq) {
      $next_R1_brk_aln = wrap_alignment("R1",$aln);
      last unless $next_R1_brk_aln->{Qname} eq $qname;
      push(@R1_alns,$next_R1_brk_aln);
      undef $next_R1_brk_aln;
    }

    if (defined $next_R2_brk_aln && $next_R2_brk_aln->{Qname} eq $qname) {
      push(@R2_alns,$next_R2_brk_aln) unless $next_R2_brk_aln->{Unmapped};
      undef $next_R2_brk_aln;
      while(my $aln = $R2_brk_iter->next_seq) {
        $next_R2_brk_aln = wrap_alignment("R2",$aln);
        last unless $next_R2_brk_aln->{Qname} eq $qname;
        push(@R2_alns,$next_R2_brk_aln) unless $next_R2_brk_aln->{Unmapped};
        undef $next_R2_brk_aln;
      }
    }

    if (defined $next_R1_adpt_aln && $next_R1_adpt_aln->{Qname} eq $qname) {
      push(@R1_alns,$next_R1_adpt_aln) unless $next_R1_adpt_aln->{Unmapped};
      undef $next_R1_adpt_aln;
      while(my $aln = $R1_adpt_iter->next_seq) {
        $next_R1_adpt_aln = wrap_alignment("R1",$aln);
        last unless $next_R1_adpt_aln->{Qname} eq $qname;
        push(@R1_alns,$next_R1_adpt_aln) unless $next_R1_adpt_aln->{Unmapped};
        undef $next_R1_adpt_aln;
      }
    }

    if (defined $next_R2_adpt_aln && $next_R2_adpt_aln->{Qname} eq $qname) {
      push(@R2_alns,$next_R2_adpt_aln) unless $next_R2_adpt_aln->{Unmapped};
      undef $next_R2_adpt_aln;
      while(my $aln = $R2_adpt_iter->next_seq) {
        $next_R2_adpt_aln = wrap_alignment("R2",$aln);
        last unless $next_R2_adpt_aln->{Qname} eq $qname;
        push(@R2_alns,$next_R2_adpt_aln) unless $next_R2_adpt_aln->{Unmapped};
        undef $next_R2_adpt_aln;
      }
    }

    if (defined $next_R1_aln && $next_R1_aln->{Qname} eq $qname) {
      push(@R1_alns,$next_R1_aln) unless $next_R1_aln->{Unmapped};
      undef $next_R1_aln;
      while(my $aln = $R1_iter->next_seq) {
        $next_R1_aln = wrap_alignment("R1",$aln);
        last unless $next_R1_aln->{Qname} eq $qname;
        push(@R1_alns,$next_R1_aln) unless $next_R1_aln->{Unmapped};
        undef $next_R1_aln;
      }
    }

    if (defined $next_R2_aln && $next_R2_aln->{Qname} eq $qname) {
      push(@R2_alns,$next_R2_aln) unless $next_R2_aln->{Unmapped};
      undef $next_R2_aln;
      while(my $aln = $R2_iter->next_seq) {
        $next_R2_aln = wrap_alignment("R2",$aln);
        last unless $next_R2_aln->{Qname} eq $qname;
        push(@R2_alns,$next_R2_aln) unless $next_R2_aln->{Unmapped};
        undef $next_R2_aln;
      }
    }

    my $OCS = find_optimal_coverage_set(\@R1_alns,\@R2_alns);

    
    
    process_optimal_coverage_set($OCS,\@R1_alns,\@R2_alns);

    


  }

  #   while (1) {
  #     unless ($aln_queue->pending < $max_queue_size) {
  #       sleep 1;
  #       next;
  #     }
  #     my $t0 = [gettimeofday];
  #     $aln_queue->enqueue(\@aln_block);
  #     my $t1 = tv_interval($t0);
  #     print "queued alignments in $t1\n";

  #     last;
  #   }
  # }

  # $aln_queue->end;

}


sub find_optimal_coverage_set ($$) {

  print "finding OCS\n";
  # my $t0 = [gettimeofday];



  # while (my $aln_block = $aln_queue->dequeue) {

  #   my $t1 = tv_interval($t0);
  #   print "dequeued alignments in $t1\n";

  #   my @OCS_block :shared;

  #   foreach my $aln_item (@$aln_block) {

      # my ($R1_alns_ref,$R2_alns_ref) :shared;
      my $R1_alns_ref = shift;
      my $R2_alns_ref = shift;

      my @graph = ();
      my $OCS_ptr;


      my @unmapped_OCS = ( { R1 => { Qname => $R1_alns_ref->[0]->{Qname},
                                     Seq => $R1_alns_ref->[0]->{Seq},
                                     Qual => $R1_alns_ref->[0]->{Qual},
                                     Unmapped => 1 },
                             R2 => { Qname => $R1_alns_ref->[0]->{Qname},
                                     Seq => $R2_alns_ref->[0]->{Seq},
                                     Qual => $R2_alns_ref->[0]->{Qual},
                                     Unmapped => 1 } } );

      
      # if ($R1_alns_ref->[0]->{Unmapped}) {
      #   my @OCS_item;
      #   @OCS_item = (shared_clone(\@unmapped_OCS),$R1_alns_ref,$R2_alns_ref);
      #   push(@OCS_block,\@OCS_item);
      #   next;
      # }

      if ($R1_alns_ref->[0]->{Unmapped}) {
        return \@unmapped_OCS;
        next;
      }

      my @R1_alns = sort {$a->{Qstart} <=> $b->{Qstart}} @$R1_alns_ref;
      my @R2_alns = sort {$a->{Qstart} <=> $b->{Qstart}} @$R2_alns_ref;

      foreach my $R1_aln (@R1_alns) {

        # print $R1_aln->{Qname}." ".$R1_aln->{Unmapped}."\n";
        next if $R1_aln->{Unmapped};

        my $graphsize = scalar @graph;

        my $new_node = {R1 => $R1_aln};

        # my %new_node :shared;
        # $new_node{R1} = $R1_aln;
        # my $new_node = \%new_node;

        my $init_score = score_edge($new_node);
        $new_node->{score} = $init_score if defined $init_score;

        # print "not a possible initial node\n" unless defined $init_score;
        # print "initialized node to $init_score\n" if defined $init_score;

        my $nodenum = 1;
        foreach my $node (@graph) {
          # print "scoring edge against node $nodenum\n";
          my $edge_score = score_edge($node,$new_node);
          next unless defined $edge_score;
          # print "found edge score of $edge_score\n";

          if (! exists $new_node->{score} || $edge_score > $new_node->{score}) {
            $new_node->{score} = $edge_score;
            $new_node->{back_ptr} = $node;

            # print "setting back pointer to $nodenum\n";
          }
          $nodenum++;
        }

        if (defined $new_node->{score}) {
          push(@graph,$new_node) ;
          # print "pushed node into position ".scalar @graph." with score ".$new_node->{score}."\n";
          # print "setting OCS pointer to node\n" if ! defined $OCS_ptr || $new_node->{score} > $OCS_ptr->{score};
          $OCS_ptr = $new_node if ! defined $OCS_ptr || $new_node->{score} > $OCS_ptr->{score};
        }

        foreach my $R2_aln (@R2_alns) {

          next if $R2_aln->{Unmapped};

          # print "testing proper-pairedness with R2:\n";
          # print_aln($R2_aln_wrap);

          next unless pair_is_proper($R1_aln,$R2_aln,$max_frag_len);

          my $new_pe_node = {R1 => $R1_aln, R2 => $R2_aln};
          # my %new_pe_node :shared;
          # $new_pe_node{R1} = $R1_aln;
          # $new_pe_node{R2} = $R2_aln;
          # my $new_pe_node = \%new_pe_node;


          # print "pair is proper\n";
          my $init_score = score_edge($new_pe_node);
          $new_pe_node->{score} = $init_score if defined $init_score;
          # print "not a possible initial node\n" unless defined $init_score;
          # print "initialized node to $init_score\n" if defined $init_score;
          $nodenum = 1;
          if ($graphsize > 0) {
            foreach my $node (@graph[0..($graphsize-1)]) {
              # print "scoring edge against node $nodenum\n";

              my $edge_score = score_edge($node,$new_pe_node);
              next unless defined $edge_score;
              # print "found edge score of $edge_score\n";

              if (! defined $new_pe_node->{score} || $edge_score > $new_pe_node->{score}) {
                $new_pe_node->{score} = $edge_score;
                $new_pe_node->{back_ptr} = $node;

                # print "setting back pointer to $nodenum\n";

              }
              $nodenum++;
            }
          }
          if (defined $new_pe_node->{score}) {
            push(@graph,$new_pe_node);
            # print "pushed node into position ".scalar @graph." with score ".$new_pe_node->{score}."\n";
            # print "setting OCS pointer to node\n" if ! defined $OCS_ptr || $new_pe_node->{score} > $OCS_ptr->{score};
            $OCS_ptr = $new_pe_node if ! defined $OCS_ptr || $new_pe_node->{score} > $OCS_ptr->{score};
          }
        }
      }

      foreach my $R2_aln (@R2_alns) {

        # print "\nStarting test for R2:\n";
        # print_aln($R2_aln_wrap);
        next if $R2_aln->{Unmapped};

        my $new_node = {R2 => $R2_aln};
        # my %new_node :shared;
        # $new_node{R2} = $R2_aln;
        # my $new_node = \%new_node;

        my $nodenum = 1;
        foreach my $node (@graph) {
          # print "scoring edge against node $nodenum\n";
          my $edge_score = score_edge($node,$new_node);
          next unless defined $edge_score;
          # print "found edge score of $edge_score\n";

          if (! defined $new_node->{score} || $edge_score > $new_node->{score}) {
            $new_node->{score} = $edge_score;
            $new_node->{back_ptr} = $node;
            # print "setting back pointer to $nodenum\n";
          }
          $nodenum++;
        }

        if (defined $new_node->{score}) {
          push(@graph,$new_node) ;
          # print "pushed node into position ".scalar @graph." with score ".$new_node->{score}."\n";
          # print "setting OCS pointer to node\n" if ! defined $OCS_ptr || $new_node->{score} > $OCS_ptr->{score};
          $OCS_ptr = $new_node if ! defined $OCS_ptr || $new_node->{score} > $OCS_ptr->{score};
        }

      }

      unless (defined $OCS_ptr) {
        return \@unmapped_OCS;
        # my @OCS_item :shared;
        # @OCS_item = (shared_clone(\@unmapped_OCS),$R1_alns_ref,$R2_alns_ref);
        # push(@OCS_block,\@OCS_item);
        next;
      }

      my @OCS;# :shared;

      while (defined $OCS_ptr->{back_ptr}) {
        unshift(@OCS,$OCS_ptr);
        $OCS_ptr->{R1_Rgap} = find_genomic_distance($OCS_ptr->{back_ptr}->{R1},$OCS_ptr->{R1},$brk_hash)
          if defined $OCS_ptr->{R1} && defined $OCS_ptr->{back_ptr}->{R1};
        $OCS_ptr->{R2_Rgap} = find_genomic_distance($OCS_ptr->{back_ptr}->{R2},$OCS_ptr->{R2},$brk_hash)
          if defined $OCS_ptr->{R2} && defined $OCS_ptr->{back_ptr}->{R2};
        $OCS_ptr->{back_ptr}->{score} = $OCS_ptr->{score};
        $OCS_ptr = $OCS_ptr->{back_ptr};
      }
      unshift(@OCS,$OCS_ptr) if defined $OCS_ptr;


      return \@OCS;

      # my @OCS_item; # :shared;
      # @OCS_item = (\@OCS,$R1_alns_ref,$R2_alns_ref);

      # push(@OCS_block,\@OCS_item);      

  #   while (1) {

  #     unless ($ocs_queue->pending < $max_queue_size) {
  #       sleep 1;
  #       next;
  #     }
  #     $ocs_queue->enqueue(\@OCS_block);
  #     last;
  #   }

  #   $t0 = [gettimeofday]
  # }

  # $ocs_queue->end();

}




sub process_optimal_coverage_set ($$$) {

  # print "processing OCSs\n";

  # use Bio::DB::Fasta;

  # my $genome_obj = Bio::DB::Fasta->new($assembly_fa);
  # my $break_obj = Bio::DB::Fasta->new($break_fa);
  # my $adpt_obj = Bio::DB::Fasta->new($adapt_fa);

  

  # while (my $ocs_item = $ocs_queue->dequeue) {

  #   my ($OCS_ref,$R1_alns,$R2_alns) = @$ocs_item;
    
    my $OCS_ref = shift;
    my $R1_alns = shift;
    my $R2_alns = shift;


    # print "create tlxls\n";
    my $tlxls = create_tlxl_entries($OCS_ref);

    $stats->{totalreads}++;

    $stats->{aligned}++ unless $tlxls->[0]->{Unmapped};
     
    # print "create tlxs\n";
    create_tlx_entries($tlxls, { genome => $R1_samobj,
                                 brk => $R1_brk_samobj,
                                 adpt => $R1_adpt_samobj} )  ;

    # print "filter unjoined\n";
    my $junctions = filter_unjoined($tlxls);

    $stats->{junctions} += $junctions;
    $stats->{junc_reads}++ if $junctions > 0;

    # print "filter map quality\n";
    my $quality_maps = filter_mapping_quality($tlxls,$R1_alns,$R2_alns,
                                        $ol_thresh,$score_thresh,$max_frag_len);

    $stats->{mapqual} += $quality_maps;
    $stats->{mapq_reads}++ if $quality_maps > 0;

    # print "filter mispriming\n";
    my $correct_priming = filter_mispriming($tlxls,$priming_threshold);

    $stats->{priming} += $correct_priming;
    $stats->{prim_reads}++ if $correct_priming > 0;

    # print "filter frequent cutter\n";
    my $no_freq_cutter = filter_freq_cutter($tlxls,$cutseq->seq);

    $stats->{freqcut} += $no_freq_cutter;
    $stats->{freq_reads}++ if $no_freq_cutter > 0;

    # print "filter breaksite\n";
    my $outside_breaksite = filter_breaksite($tlxls);

    $stats->{breaksite} += $outside_breaksite;
    $stats->{break_reads}++ if $outside_breaksite > 0;


    # print "filter split juctions\n";
    my $primary_junction = filter_split_junctions($tlxls);

    $stats->{splitjuncs} += $primary_junction;
    $stats->{split_reads}++ if $primary_junction > 0;

    write_tlxls($tlxls);


  

}

sub write_tlxls ($) {
  my $tlxls = shift;
  print "writing stuff\n";
  foreach my $tlxl (@$tlxls) {
    write_entry($tlxlfh,$tlxl,\@tlxl_header);
    next unless defined $tlxl->{tlx};
    if (defined $tlxl->{tlx}->{Filter}) {
      if ($tlxl->{tlx}->{Filter} eq "Unjoined" || $tlxl->{tlx}->{Filter} eq "Unaligned") {
        write_entry($unjoin_tlxfh,$tlxl->{tlx},\@tlx_filter_header);
      } else {
        write_entry($filt_tlxfh,$tlxl->{tlx},\@tlx_filter_header);
      }
    } else {
      write_entry($tlxfh,$tlxl->{tlx},\@tlx_header);
    }
  }
}


sub score_edge ($;$) {
  my $node1 = shift;
  my $node2 = shift;

  my $score;

  if (defined $node2) {

    my $R1_Qgap = 0;
    my $R2_Qgap = 0;
    my $Rname1;
    my $Rname2;
    my $Strand1;
    my $Strand2;
    my $Junc1;
    my $Junc2;
    my $R1_Rdist;
    my $R2_Rdist;
    my $Len1 = 0;
    my $Len2 = 0;



    if (defined $node2->{R1}) {
      return undef unless defined $node1->{R1};
      return undef if $node1->{R1}->{Rname} eq "Adapter";
      return undef if $node2->{R1} == $node1->{R1};
      return undef unless $node2->{R1}->{Qend} > $node1->{R1}->{Qend};
      $Rname1 = $node1->{R1}->{Rname};
      $Strand1 = $node1->{R1}->{Strand};
      $Rname2 = $node2->{R1}->{Rname};
      $Strand2 = $node2->{R1}->{Strand};
      $R1_Qgap = $node2->{R1}->{Qstart} - $node1->{R1}->{Qend} - 1;
      $R1_Rdist = find_genomic_distance($node1->{R1},$node2->{R1},$brk_hash);
      $Len1 += $node1->{R1}->{Qend} - $node1->{R1}->{Qstart} + 1;
      $Len2 += $node2->{R1}->{Qend} - $node2->{R1}->{Qstart} + 1;
    } else {
      return undef unless defined $node1->{R2};
    }


    if (defined $node1->{R2}) {
      return undef unless defined $node2->{R2};
      return undef if $node1->{R2}->{Rname} eq "Adapter";
      return undef if $node2->{R2} == $node1->{R2};
      return undef unless $node2->{R2}->{Qend} > $node1->{R2}->{Qend};
      $Rname1 = $node1->{R2}->{Rname};
      $Strand1 = $node1->{R2}->{Strand};
      $Rname2 = $node2->{R2}->{Rname};
      $Strand2 = $node2->{R2}->{Strand};
      $R2_Qgap = $node2->{R2}->{Qstart} - $node1->{R2}->{Qend} - 1;
      $R2_Rdist = find_genomic_distance($node1->{R2},$node2->{R2},$brk_hash);
      $Len1 += $node1->{R2}->{Qend} - $node1->{R2}->{Qstart} + 1;
      $Len2 += $node2->{R2}->{Qend} - $node2->{R2}->{Qstart} + 1;
    } else {
      return undef unless defined $node2->{R1};
    }

    my $totalOverlap = -min($R1_Qgap,0) -min($R2_Qgap,0);
    return undef if $totalOverlap > 0.5 * $Len1 || $totalOverlap > 0.5 * $Len1;

    my $OL_correction;
    my $Qgap_pen;
    my $Rgap_pen;
    my $Brk_pen;

    if (defined $node2->{R1} && defined $node1->{R2}) {
      
      $Junc1 = $Strand1 == 1 ? max($node1->{R1}->{Rend},$node1->{R2}->{Rend}) : min($node1->{R1}->{Rstart},$node1->{R2}->{Rstart});
      $Junc2 = $Strand2 == 1 ? min($node2->{R1}->{Rstart},$node2->{R2}->{Rstart}) : max($node2->{R1}->{Rend},$node2->{R2}->{Rend});
      $OL_correction = $OL_mult * max(0,-$R1_Qgap) + $OL_mult * max(0,-$R2_Qgap);
      $Qgap_pen = $Dif_mult * abs($R1_Qgap - $R2_Qgap); #$Dif_mult * (max(abs($R1_Qgap),abs($R2_Qgap)) + abs($R1_Qgap - $R2_Qgap));
      $Rgap_pen = $Dif_mult * abs($node2->{R1}->{Rstart} - $node2->{R2}->{Rstart});
    } else {

      if (defined $node2->{R1}) {
        $Junc1 = $Strand1 == 1 ? $node1->{R1}->{Rend} : $node1->{R1}->{Rstart};
        $Junc2 = $Strand2 == 1 ? $node2->{R1}->{Rstart} : $node1->{R1}->{Rend};
        $OL_correction = $OL_mult * max(0,-$R1_Qgap);
        #$Qgap_pen = $Dif_mult * abs($R1_Qgap);
      } else {
        $Junc1 = $Strand1 == 1 ? $node1->{R2}->{Rend} : $node1->{R2}->{Rstart};
        $Junc2 = $Strand2 == 1 ? $node2->{R2}->{Rstart} : $node1->{R2}->{Rend};
        $OL_correction = $OL_mult * max(0,-$R2_Qgap);
        #$Qgap_pen = $Dif_mult * abs($R2_Qgap);
      }
      $Qgap_pen = 0;
      $Rgap_pen = 0;
    }

    

    if ($Rname2 eq "Adapter") {
      $Brk_pen = 0;
    # } elsif ($Rname1 eq "Breaksite" && $Rname2 eq "Breaksite" && $Strand1 == 1 && $Strand2 == 1 && $Junc1 < $Junc2) {
    #   $Brk_pen = 0;
    } else {

      my $R1_Brk_pen;
      my $R2_Brk_pen;

      if (defined $R1_Rdist) {
        $R1_Brk_pen = $R1_Rdist > 1 || $Strand1 ne $Strand2 ? $Brk_pen_min + $Brk_pen_mult * log10($R1_Rdist)^$Brk_pen_power : 0;
      } else {
        $R1_Brk_pen = $Brk_pen_max;
      }

      if (defined $R2_Rdist) {
        $R2_Brk_pen = $R2_Rdist > 1 || $Strand1 ne $Strand2 ? $Brk_pen_min + $Brk_pen_mult * log10($R2_Rdist)^$Brk_pen_power : 0;
      } else {
        $R2_Brk_pen = $Brk_pen_max;
      }

      $Brk_pen = min($R1_Brk_pen,$R2_Brk_pen,$Brk_pen_max);;
    }


    my $R1_AS = defined $node2->{R1} && $node2->{R1}->{AS};
    my $R2_AS = defined $node2->{R2} && $node2->{R2}->{AS};
    my $PEgap;
    my $PEgap_pen;

    if (defined $node2->{R1} && defined $node2->{R2}) {
      $PEgap = $node2->{R1}->{Strand} == 1 ? $node2->{R2}->{Rstart} - $node2->{R1}->{Rend} : $node2->{R1}->{Rstart} - $node2->{R1}->{Rend};
    }

    $PEgap_pen = defined $PEgap && $PEgap > 1 ? $Brk_pen_min + $Brk_pen_mult * log10($PEgap)^$Brk_pen_power : 0;

    $score = $node1->{score} + $R1_AS + $R2_AS - $PEgap_pen - $Brk_pen - $Qgap_pen - $OL_correction - $Rgap_pen;

    
  } else {

    return undef unless defined $node1->{R1};
    return undef unless $node1->{R1}->{Rname} eq "Breaksite";
    return undef unless $node1->{R1}->{Strand} == 1;

    my $brk_start_gap = $node1->{R1}->{Rstart} - $primer_start;

    my $R1_AS = $node1->{R1}->{AS};
    my $R2_AS = defined $node1->{R2} && $node1->{R2}->{AS};
    my $PEgap;
    my $PEgap_pen;


    if (defined $node1->{R1} && defined $node1->{R2}) {
      $PEgap = $node1->{R1}->{Strand} == 1 ? $node1->{R2}->{Rstart} - $node1->{R1}->{Rend} : $node1->{R1}->{Rstart} - $node1->{R1}->{Rend};
    }

    $PEgap_pen = defined $PEgap && $PEgap > 1 ? $Brk_pen_min + $Brk_pen_mult * log10($PEgap)^$Brk_pen_power : 0;

    $score = $R1_AS + $R2_AS - $PEgap_pen - $Dif_mult * $brk_start_gap;

  }

  return $score;

}

sub deduplicate_junctions {


  my $dedup_cmd = "$FindBin::Bin/../R/TranslocDedup.R $tlxfile $dedup_output cores=$max_threads";
  System($dedup_cmd);

  my $tlxbak = "$tlxfile.bak";
  rename $tlxfile, $tlxbak;

  my %hash;

  my $dedupfh = IO::File->new("<$dedup_output");

  while (my $read = $dedupfh->getline) {
    chomp($read);
    my @dup = split("\t",$read);
    $hash{$dup[0]} = 1;
  }

  $dedupfh->close;

  my $filt_tlxfh = IO::File->new(">>$filt_tlxfile");
  my $tlxfh = IO::File->new(">$tlxfile");
  my $bak_tlxfh = IO::File->new("<$tlxbak");

  my $csv = Text::CSV->new({sep_char => "\t"});
  my $header = $csv->getline($bak_tlxfh);
  $csv->column_names(@$header);
  $tlxfh->print(join("\t", @tlx_header)."\n");

  while (my $tlx = $csv->getline_hr($bak_tlxfh)) {
    if (exists $hash{$tlx->{Qname}}) {
      $tlx->{Filter} = "DeDup";
      write_entry($filt_tlxfh,$tlx,\@tlx_filter_header);
    } else {
      write_entry($tlxfh,$tlx,\@tlx_header);
      $stats->{dedup}++;

    }
  }

}

sub post_process_junctions {

  (my $html_reads = $tlxfile) =~ s/tlx$/html/;
  System("TranslocHTMLReads.pl $tlxfile $html_reads");  

}

sub parse_command_line {
	my $help;

	usage() if (scalar @ARGV == 0);


	my $result = GetOptions ( "read1=s" => \$read1,
                            "read2=s" => \$read2,
                            "assembly=s" => \$assembly,
                            "chr=s" => \$brk_chr,
                            "start=i" => \$brk_start,
                            "end=i" => \$brk_end,
                            "strand=s" => \$brk_strand,
                            "workdir=s" => \$workdir,
                            "threads=i" => \$max_threads,
                            "bt2opt=s" => \$user_bowtie_opt,
                            "bt2brkopt=s" => \$user_bowtie_breaksite_opt,
                            "usecurrtlx" => \$use_current_tlx,
				            				"help" => \$help
				            			) ;

	
	usage() if ($help);

  #Check options

  croak "Error: cannot read sequence files" unless (-r $read1 & -r $read2);
  croak "Error: working directory does not exist" unless (-d $workdir);
  croak "Error: pipeline requires at least 3 threads" if ($max_threads < 2);
  



	exit unless $result;
}


sub usage()
{
print<<EOF;
Title, by Robin Meyers, ddmonthyyyy

This program .


Usage: $0 arg1 arg2 arg3 ...
        [--option VAL] [--flag] [--help]

Arguments (defaults in parentheses):

$arg{"--read1","Input sequence file"}
$arg{"--read2","Input sequence file"}
$arg{"--workdir","Directory for results files - note: default Bowtie output goes to working directory"}
$arg{"--assembly","Genome assembly to align reads to"}
$arg{"--threads","Number of threads to run bowtie on","$max_threads"}
$arg{"--help","This helpful help screen."}


EOF

exit 1;
}
