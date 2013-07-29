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
use Bio::DB::Sam;
use List::Util qw(min max);
use Interpolation 'arg:@->$' => \&argument;
use Time::HiRes qw(gettimeofday tv_interval);

use Cwd qw(abs_path);
use FindBin;
use lib abs_path("$FindBin::Bin/../lib");

require "PerlSub.pl";
require "TranslocHelper.pl";


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
sub score_edge ($;$);



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
my $threads = 4;
my $priming_bp = 12;

my $OL_mult = 2;
my $Dif_mult = 1;
my $Brk_pen_mult = 5;
my $Brk_pen_max = 40;
my $max_frag_len = 1500;

my $user_bowtie_opt = "";
my $user_bowtie_adapter_opt = "";
my $user_bowtie_breaksite_opt = "";

# Global variables
my @tlxl_header = tlxl_header();
my @tlx_header = tlx_header();

my $stats = { totreads    => 0,
              aligned     => 0,
              alignment   => 0,
              orientation => 0,
              primaln     => 0,
              splitaln    => 0,
              mapqual     => 0,
              dedup       => 0,
              final       => 0 };


#
# Start of Program
#

parse_command_line;


my $default_bowtie_adapter_opt = "--local -D 20 -R 3 -N 1 -L 6 -i C,4 --score-min C,30 --mp 10 -p $threads --no-unal --reorder -t";
my $default_bowtie_breaksite_opt = "--local -D 20 -R 3 -N 1 -L 12 -i C,6 --score-min C,40 --mp 10 --rfg 12,12 --rdg 12,12 -p $threads -k 20 --reorder -t";
my $default_bowtie_opt = "--local -D 20 -R 3 -N 0 -L 20 -i C,8 --score-min C,50 --mp 10 --rfg 12,12 --rdg 12,12 -p $threads -k 10 --no-unal --reorder -t";

my $bt2_break_opt = manage_program_options($default_bowtie_breaksite_opt,$user_bowtie_breaksite_opt);
my $bt2_adapt_opt = manage_program_options($default_bowtie_adapter_opt,$user_bowtie_adapter_opt);
my $bt2_opt = manage_program_options($default_bowtie_opt,$user_bowtie_opt);

my $t0 = [gettimeofday];

#check_working_dir;

my $expt = basename($workdir);
my $expt_stub = "$workdir/$expt";
# my $tmpstitchdir = "$workdir/tmp_stitch";
# mkdir $tmpstitchdir;
# my $tmpbrkdir = "$workdir/tmp_brksite";
# mkdir $tmpbrkdir;
my $break_fa = "$workdir/misc/breaksite.fa";
my $prim_fa = "$workdir/misc/primer.fa";
my $adapt_fa = "$workdir/misc/adapter.fa";
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


# my $breaksite_sam = "${expt_stub}_breaksite.sam";
# my $breaksite_bam = "${expt_stub}_breaksite.bam";
# my $sam = "$expt_stub.sam";
# my $bam = "$expt_stub.bam";
# my $unaligned_tlxl = "${expt_stub}_unaligned.tlxl";
# my $unal_tlxlfh = IO::File->new(">$unaligned_tlxl");
# my $filtered_tlxl = "${expt_stub}_filtered.tlxl";
my $tlxlfile = "${expt_stub}.tlxl";
my $tlxlfh = IO::File->new(">$tlxlfile");

my $filt_tlxlfile = "${expt_stub}_filtered.tlxl";
my $filt_tlxlfh = IO::File->new(">$filt_tlxlfile");


# my $filtered_tlx = "${expt_stub}_filtered.tlx";
# my $tlx = "$expt_stub.tlx";
# my $tlxfh = IO::File->new(">$tlx");

# my $statsfile = "$expt_stub.stats";


my $brk_io = Bio::SeqIO->new(-file => $break_fa, -format => 'fasta');
my $breakseq = $brk_io->next_seq();
my $breaklen = length($breakseq->seq);

my $prim_io = Bio::SeqIO->new(-file => $prim_fa, -format => 'fasta');
my $primseq = $prim_io->next_seq();
my $primlen = length($primseq->seq);
my $priming_threshold = $primlen + $priming_bp;

my $adpt_io = Bio::SeqIO->new(-file => $adapt_fa, -format => 'fasta');
my $adaptseq = $adpt_io->next_seq();
my $adaptlen = length($adaptseq->seq);

my $brk_hash = {chr => $brk_chr,
                start => $brk_start,
                end => $brk_end,
                strand => $brk_strand,
                len => $breaklen};



align_to_breaksite unless -r $R1_brk_bam && -r $R2_brk_bam;

align_to_adapter unless -r $R1_adpt_bam && -r $R2_adpt_bam;

align_to_genome unless -r $R1_bam && -r $R2_bam;



# if ($read1 =~ s/\.gz//) {
#   System("gunzip -c $read1.gz > $read1") unless -r $read1;
# }
# if ($read2 =~ s/\.gz//) {
#   System("gunzip -c $read2.gz > $read2") unless -r $read2;
# }

process_alignments;

#filter_tlxl;


my $t1 = tv_interval($t0);

printf("\nFinished all processes in %.2f seconds.\n", $t1);

print("Some stats:\n".join("\t",$stats->{totreads},
                                $stats->{aligned},
                                $stats->{alignments},
                                $stats->{orientation},
                                $stats->{primaln},
                                $stats->{splitaln},
                                $stats->{mapqual},
                                $stats->{priming},
                                $stats->{final})."\n");


#
# End of program
#

sub align_to_breaksite {
  print "\nRunning Bowtie2 alignment for $expt against breaksite sequence\n";

  System("bowtie2-build -q $break_fa $break_bt2idx");

  my $R1_brk_bt2_cmd = "bowtie2 $bt2_break_opt -x $break_bt2idx -U $read1 -S $R1_brk_sam";

  System($R1_brk_bt2_cmd);

  System("samtools view -bS -o $R1_brk_bam $R1_brk_sam");

  my $R2_brk_bt2_cmd = "bowtie2 $bt2_break_opt -x $break_bt2idx -U $read2 -S $R2_brk_sam";

  System($R2_brk_bt2_cmd);

  System("samtools view -bS -o $R2_brk_bam $R2_brk_sam");

}

sub align_to_adapter {
  print "\nRunning Bowtie2 alignment for $expt against adapter sequence\n";

  System("bowtie2-build -q $adapt_fa $adapt_bt2idx");

  my $R1_adpt_bt2_cmd = "bowtie2 $bt2_adapt_opt -x $adapt_bt2idx -U $read1 -S $R1_adpt_sam";

  System($R1_adpt_bt2_cmd);

  System("samtools view -bS -o $R1_adpt_bam $R1_adpt_sam");

  my $R2_adpt_bt2_cmd = "bowtie2 $bt2_adapt_opt -x $adapt_bt2idx -U $read2 -S $R2_adpt_sam";

  System($R2_adpt_bt2_cmd);

  System("samtools view -bS -o $R2_adpt_bam $R2_adpt_sam");

}

sub align_to_genome {

  print "\nRunning Bowtie2 alignment for $expt against $assembly genome \n";

  my $R1_bt2_cmd = "bowtie2 $bt2_opt -x $assembly -U $read1 -S $R1_sam";

  System($R1_bt2_cmd);

  System("samtools view -bS -o $R1_bam $R1_sam");

  my $R2_bt2_cmd = "bowtie2 $bt2_opt -x $assembly -U $read2 -S $R2_sam";

  System($R2_bt2_cmd);

  System("samtools view -bS -o $R2_bam $R2_sam");

  
}

sub process_alignments {

  print "\nProcessing alignments and writing tlxl and tlx files\n";

  $tlxlfh->print(join("\t", @tlxl_header)."\n");
  $filt_tlxlfh->print(join("\t", @tlxl_header)."\n");


  my $R1_brk_samobj = Bio::DB::Sam->new(-bam => $R1_brk_bam,
                                     -fasta => $break_fa,
                                     -expand_flags => 1);

  my $R2_brk_samobj = Bio::DB::Sam->new(-bam => $R2_brk_bam,
                                     -fasta => $break_fa,
                                     -expand_flags => 1);

  my $R1_adpt_samobj = Bio::DB::Sam->new(-bam => $R1_adpt_bam,
                                     -fasta => $adapt_fa,
                                     -expand_flags => 1);

  my $R2_adpt_samobj = Bio::DB::Sam->new(-bam => $R2_adpt_bam,
                                     -fasta => $adapt_fa,
                                     -expand_flags => 1);

  my $R1_samobj = Bio::DB::Sam->new(-bam => $R1_bam,
                                 -fasta => $ENV{'GENOME_DB'}."/$assembly/$assembly.fa",
                                 -expand_flags => 1);

  my $R2_samobj = Bio::DB::Sam->new(-bam => $R2_bam,
                                 -fasta => $ENV{'GENOME_DB'}."/$assembly/$assembly.fa",
                                 -expand_flags => 1);




  my $R1_brk_iter = $R1_brk_samobj->get_seq_stream();
  my $R2_brk_iter = $R2_brk_samobj->get_seq_stream();
  my $R1_adpt_iter = $R1_adpt_samobj->get_seq_stream();
  my $R2_adpt_iter = $R2_adpt_samobj->get_seq_stream();
  my $R1_iter = $R1_samobj->get_seq_stream();
  my $R2_iter = $R2_samobj->get_seq_stream();

  my $next_R1_brk_aln = $R1_brk_iter->next_seq;
  my $next_R2_brk_aln = $R2_brk_iter->next_seq;
  my $next_R1_adpt_aln = $R1_adpt_iter->next_seq;
  my $next_R2_adpt_aln = $R2_adpt_iter->next_seq;
  my $next_R1_aln = $R1_iter->next_seq;
  my $next_R2_aln = $R2_iter->next_seq;

  croak "Error: no R1 alignments to breaksite" unless defined $next_R1_brk_aln;

  while (defined $next_R1_brk_aln) {

    my $qname = $next_R1_brk_aln->query->name;

    my @R1_alns = ();
    my @R2_alns = ();


    my $brk_unmapped = $next_R1_brk_aln->unmapped ? $next_R1_brk_aln : undef;

    push(@R1_alns,$next_R1_brk_aln);

    while($next_R1_brk_aln = $R1_brk_iter->next_seq) {
      last unless $next_R1_brk_aln->query->name eq $qname;
      push(@R1_alns,$next_R1_brk_aln);
    }

    if (defined $next_R2_brk_aln && $next_R2_brk_aln->query->name eq $qname) {
      push(@R2_alns,$next_R2_brk_aln);
      while($next_R2_brk_aln = $R2_brk_iter->next_seq) {
        last unless $next_R2_brk_aln->query->name eq $qname;
        push(@R2_alns,$next_R2_brk_aln);
      }
    }

    if (defined $next_R1_adpt_aln && $next_R1_adpt_aln->query->name eq $qname) {
      push(@R1_alns,$next_R1_adpt_aln);
      while($next_R1_adpt_aln = $R1_adpt_iter->next_seq) {
        last unless $next_R1_adpt_aln->query->name eq $qname;
        push(@R1_alns,$next_R1_adpt_aln);
      }
    }

    if (defined $next_R2_adpt_aln && $next_R2_adpt_aln->query->name eq $qname) {
      push(@R2_alns,$next_R2_adpt_aln);
      while($next_R2_adpt_aln = $R2_adpt_iter->next_seq) {
        last unless $next_R2_adpt_aln->query->name eq $qname;
        push(@R2_alns,$next_R2_adpt_aln);
      }
    }

    if (defined $next_R1_aln && $next_R1_aln->query->name eq $qname) {
      push(@R1_alns,$next_R1_aln);
      while($next_R1_aln = $R1_iter->next_seq) {
        last unless $next_R1_aln->query->name eq $qname;
        push(@R1_alns,$next_R1_aln);
      }
    }

    if (defined $next_R2_aln && $next_R2_aln->query->name eq $qname) {
      push(@R2_alns,$next_R2_aln);
      while($next_R2_aln = $R2_iter->next_seq) {
        last unless $next_R2_aln->query->name eq $qname;
        push(@R2_alns,$next_R2_aln);
      }
    }

    my $tlxls;

    unless (defined $brk_unmapped) {
      my ($OCS_ref,$R1_alns_ref,$R2_alns_ref) = find_optimal_coverage_set(\@R1_alns,\@R2_alns);

      $tlxls = create_tlxl_entries($OCS_ref);

      filter_unmapped($tlxls);


    } else {
      $tlxls = create_tlxl_entries([{R1 => {aln => $brk_unmapped}}]);
    }



    foreach my $tlxl (@$tlxls) {
      if (defined $tlxl->{Filter}) {
        write_tlxl_entry($filt_tlxlfh,$tlxl);
      } else {
        write_tlxl_entry($tlxlfh,$tlxl);
      }
    }

  }

}


sub find_optimal_coverage_set ($$) {

  my $R1_alns_ref = shift;
  my $R2_alns_ref = shift;

  my @graph = ();
  my $OCS_ptr;

  my @R1_aln_wraps = ();
  my @R2_aln_wraps = ();

  foreach my $R1_aln (@$R1_alns_ref) {
    next if $R1_aln->unmapped;
    push(@R1_aln_wraps,wrap_alignment($R1_aln,"R1"));
  }
  foreach my $R2_aln (@$R2_alns_ref) {
    next if $R2_aln->unmapped;
    push(@R2_aln_wraps,wrap_alignment($R2_aln,"R2"));
  }

  @R1_aln_wraps = sort {$a->{Qstart} <=> $b->{Qstart}} @R1_aln_wraps;
  @R2_aln_wraps = sort {$a->{Qstart} <=> $b->{Qstart}} @R2_aln_wraps;

  foreach my $R1_aln_wrap (@R1_aln_wraps) {

    print "\nStarting test for R1:\n";
    print_aln($R1_aln_wrap);

    my $graphsize = $#graph;

    my $new_node = {R1 => $R1_aln_wrap};
    
    my $init_score = score_edge($new_node);
    $new_node->{score} = $init_score if defined $init_score;

    print "not a possible initial node\n" unless defined $init_score;
    print "initialized node to $init_score\n" if defined $init_score;

    my $nodenum = 1;
    foreach my $node (@graph) {
      print "scoring edge against node $nodenum\n";
      my $edge_score = score_edge($node,$new_node);
      next unless defined $edge_score;
      print "found edge score of $edge_score\n";

      if (! exists $new_node->{score} || $edge_score > $new_node->{score}) {
        $new_node->{score} = $edge_score;
        $new_node->{back_ptr} = $node;
        print "setting back pointer to $nodenum\n";
      }
      $nodenum++;
    }

    if (defined $new_node->{score}) {
      push(@graph,$new_node) ;
      print "pushed node into position ".scalar @graph." with score ".$new_node->{score}."\n";
      print "setting OCS pointer to node\n" if ! defined $OCS_ptr || $new_node->{score} > $OCS_ptr->{score};
      $OCS_ptr = $new_node if ! defined $OCS_ptr || $new_node->{score} > $OCS_ptr->{score};
    }

    foreach my $R2_aln_wrap (@R2_aln_wraps) {

      print "testing proper-pairedness with R2:\n";
      print_aln($R2_aln_wrap);

      next unless pair_is_proper($R1_aln_wrap,$R2_aln_wrap,$max_frag_len);

      my $new_pe_node = {R1 => $R1_aln_wrap, R2 => $R2_aln_wrap};


      print "pair is proper\n";
      my $init_score = score_edge($new_pe_node);
      $new_pe_node->{score} = $init_score if defined $init_score;
      print "not a possible initial node\n" unless defined $init_score;
      print "initialized node to $init_score\n" if defined $init_score;
      $nodenum = 1;
      foreach my $node (@graph[0..$graphsize]) {
        print "scoring edge against node $nodenum\n";

        my $edge_score = score_edge($node,$new_pe_node);
        next unless defined $edge_score;
        print "found edge score of $edge_score\n";

        if (! defined $new_pe_node->{score} || $edge_score > $new_pe_node->{score}) {
          $new_pe_node->{score} = $edge_score;
          $new_pe_node->{back_ptr} = $node;
          print "setting back pointer to $nodenum\n";

        }
        $nodenum++;
      }
      if (defined $new_pe_node->{score}) {
        push(@graph,$new_pe_node);
        print "pushed node into position ".scalar @graph." with score ".$new_pe_node->{score}."\n";
        print "setting OCS pointer to node\n" if ! defined $OCS_ptr || $new_pe_node->{score} > $OCS_ptr->{score};
        $OCS_ptr = $new_pe_node if ! defined $OCS_ptr || $new_pe_node->{score} > $OCS_ptr->{score};
      }
    }
  }

  foreach my $R2_aln_wrap (@R2_aln_wraps) {

    print "\nStarting test for R2:\n";
    print_aln($R2_aln_wrap);


    my $new_node = {R2 => $R2_aln_wrap};

    my $nodenum = 1;
    foreach my $node (@graph) {
      print "scoring edge against node $nodenum\n";
      my $edge_score = score_edge($node,$new_node);
      next unless defined $edge_score;
      print "found edge score of $edge_score\n";

      if (! defined $new_node->{score} || $edge_score > $new_node->{score}) {
        $new_node->{score} = $edge_score;
        $new_node->{back_ptr} = $node;
        print "setting back pointer to $nodenum\n";
      }
      $nodenum++;
    }

    if (defined $new_node->{score}) {
      push(@graph,$new_node) ;
      print "pushed node into position ".scalar @graph." with score ".$new_node->{score}."\n";
      print "setting OCS pointer to node\n" if ! defined $OCS_ptr || $new_node->{score} > $OCS_ptr->{score};
      $OCS_ptr = $new_node if ! defined $OCS_ptr || $new_node->{score} > $OCS_ptr->{score};
    }

  }

  my @OCS = ();
  while (defined $OCS_ptr->{back_ptr}) {
    unshift(@OCS,$OCS_ptr);
    $OCS_ptr->{back_ptr}->{score} = $OCS_ptr->{score};
    $OCS_ptr = $OCS_ptr->{back_ptr};
  }
  unshift(@OCS,$OCS_ptr) if defined $OCS_ptr;
  return (\@OCS,\@R1_aln_wraps,\@R2_aln_wraps);

}





sub score_edge ($;$) {
  my $node1 = shift;
  my $node2 = shift;

  my $score;

  if (defined $node2) {

    my $R1_Qgap;
    my $R2_Qgap;
    my $Rname1;
    my $Rname2;
    my $Strand1;
    my $Strand2;
    my $Junc1;
    my $Junc2;
    my $R1_Rdist;
    my $R2_Rdist;



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
      $Junc1 = $Strand1 == 1 ? $node1->{R2}->{Rend} : $node1->{R2}->{Rstart};
      $Junc2 = $Strand1 == 1 ? $node1->{R2}->{Rstart} : $node1->{R2}->{Rend};
      $R2_Rdist = find_genomic_distance($node1->{R2},$node2->{R2},$brk_hash);
    } else {
      return undef unless defined $node2->{R1};
    }

    my $OL_correction;
    my $Qgap_pen;
    my $Rgap_pen;
    my $Brk_pen;

    if (defined $node2->{R1} && defined $node1->{R2}) {
      $Junc1 = $Strand1 == 1 ? max($node1->{R1}->{Rend},$node1->{R2}->{Rend}) : min($node1->{R1}->{Rstart},$node1->{R2}->{Rstart});
      $Junc2 = $Strand2 == 1 ? min($node2->{R1}->{Rstart},$node2->{R2}->{Rstart}) : max($node2->{R1}->{Rend},$node2->{R2}->{Rend});
      $OL_correction = $OL_mult * max(0,-$R1_Qgap) + $OL_mult * max(0,-$R2_Qgap);
      $Qgap_pen = $Dif_mult * (max(abs($R1_Qgap),abs($R2_Qgap)) + abs($R1_Qgap - $R2_Qgap));
      $Rgap_pen = $Dif_mult * abs($node2->{R1}->{Rstart} - $node2->{R2}->{Rstart});
    } else {

      if (defined $node2->{R1}) {
        $Junc1 = $Strand1 == 1 ? $node1->{R1}->{Rend} : $node1->{R1}->{Rstart};
        $Junc2 = $Strand2 == 1 ? $node2->{R1}->{Rstart} : $node1->{R1}->{Rend};
        $OL_correction = $OL_mult * max(0,-$R1_Qgap);
        $Qgap_pen = $Dif_mult * abs($R1_Qgap);
      } else {
        $Junc1 = $Strand1 == 1 ? $node1->{R2}->{Rend} : $node1->{R2}->{Rstart};
        $Junc2 = $Strand2 == 1 ? $node2->{R2}->{Rstart} : $node1->{R2}->{Rend};
        $OL_correction = $OL_mult * max(0,-$R2_Qgap);
        $Qgap_pen = $Dif_mult * abs($R2_Qgap);
      }
      $Rgap_pen = 0;
    }

    

    if ($Rname2 eq "Adapter") {
      $Brk_pen = 0;
    } else {
      my $R1_Brk_pen = defined $R1_Rdist ? $Brk_pen_mult * log10($R1_Rdist) : $Brk_pen_max;
      my $R2_Brk_pen = defined $R2_Rdist ? $Brk_pen_mult * log10($R2_Rdist) : $Brk_pen_max;
      $Brk_pen = min($R1_Brk_pen,$R2_Brk_pen);
    }

    if ($Rname1 eq "Breaksite" && $Rname2 eq "Breaksite" && $Strand1 == 1 && $Strand2 == 1 && $Junc1 < $Junc2) {
      $Brk_pen = 2;
      $Qgap_pen = 0;
    }

    my $R1_AS = defined $node2->{R1} && $node2->{R1}->{aln}->aux =~ /AS:i:(\d+)/ ? $1 : 0;
    my $R2_AS = defined $node2->{R2} && $node2->{R2}->{aln}->aux =~ /AS:i:(\d+)/ ? $1 : 0;
    my $PEgap;
    my $PEgap_pen;

    if (defined $node2->{R1} && defined $node2->{R2}) {
      $PEgap = $node2->{R1}->{Strand} == 1 ? $node2->{R2}->{Rstart} - $node2->{R1}->{Rend} : $node2->{R1}->{Rstart} - $node2->{R1}->{Rend};
    }

    $PEgap_pen = defined $PEgap && $PEgap > 0 ? $Brk_pen_mult * log10($PEgap) : 0;

    $score = $node1->{score} + $R1_AS + $R2_AS - $PEgap_pen - $Brk_pen - $Qgap_pen - $OL_correction - $Rgap_pen;

    
  } else {

    return undef unless defined $node1->{R1};
    return undef unless $node1->{R1}->{Rname} eq "Breaksite";
    return undef unless $node1->{R1}->{Strand} == 1;

    my $brk_start_gap = $node1->{R1}->{Rstart} - 1;

    my $R1_AS = $node1->{R1}->{aln}->aux =~ /AS:i:(\d+)/ ? $1 : 0;
    my $R2_AS = defined $node1->{R2} && $node1->{R2}->{aln}->aux =~ /AS:i:(\d+)/ ? $1 : 0;
    my $PEgap;
    my $PEgap_pen;


    if (defined $node1->{R1} && defined $node1->{R2}) {
      $PEgap = $node1->{R1}->{Strand} == 1 ? $node1->{R2}->{Rstart} - $node1->{R1}->{Rend} : $node1->{R1}->{Rstart} - $node1->{R1}->{Rend};
    }

    $PEgap_pen = defined $PEgap && $PEgap > 0 ? $Brk_pen_mult * log10($PEgap) : 0;

    $score = $R1_AS + $R2_AS - $PEgap_pen - $Dif_mult * $brk_start_gap;

  }

  return $score;

}


  # my $unal_tlxlfh = IO::File->new(">$unaligned_tlxl");
  # my $filt_tlxlfh = IO::File->new(">$filtered_tlxl");
  # my $tlxlfh = IO::File->new(">$tlxl");
  # my $filt_tlxfh = IO::File->new(">$filtered_tlx");
  # my $tlxfh = IO::File->new(">$tlx");

#   $unal_tlxlfh->print(join("\t", @tlxl_header)."\n");
#   $filt_tlxlfh->print(join("\t", @tlxl_header)."\n");
#   $tlxlfh->print(join("\t", @tlxl_header)."\n");
#   $filt_tlxfh->print(join("\t", @tlx_header)."\n");
#   $tlxfh->print(join("\t", @tlx_header)."\n");

#   my %alns;

#   my $iter = $samobj->get_seq_stream(-type=>'match');
#   while (my $aln = $iter->next_seq) {
#     my $qname = $aln->query->name;
#     if (defined $alns{$aln->query->name}) {
#       push(@{$alns{$qname}},$aln);
#     } else {
#       $alns{$qname} = [$aln];
#     }
#   }

#   while (my $brk_aln = $brk_iter->next_seq) {


#     $stats->{totreads}++;

#     my ($R1_brk_aln, $R2_brk_aln) = $brk_aln->get_SeqFeatures;
#     #my @q_alns = $samobj->get_features_by_name($R1_brk_aln->query->name);
#     my @q_alns = @{$alns{$R1_brk_aln->query->name}} if exists $alns{$R1_brk_aln->query->name};

#     if (scalar @q_alns > 0) {
#       $stats->{aligned}++;

#       my @alns = ();

#       my $i = 0;
#       my $j = 0;
#       while ($i < @q_alns) {

#         my $aln = $q_alns[$i++];
#         my ($R1_aln,$R2_aln);

#         if ($aln->get_tag_values('FIRST_MATE')) {
#           $R1_aln = $aln;
#           $R2_aln = $q_alns[$i++] unless ($R1_aln->munmapped);
#         } else {
#           $R2_aln = $aln;
#         }
#         $alns[$j++] = { R1_aln => $R1_aln, R2_aln => $R2_aln, mapq => undef, orientation => undef };
#       }
      
#       my @tlxl_objs = ();

#       foreach my $aln (@alns) {

#         $stats->{alignments}++;

#         my $tlxl_obj = create_tlxl_object($aln,$R1_brk_aln,$R2_brk_aln);

#         push(@tlxl_objs,$tlxl_obj);

#       }

#       my $primary_tlxl = find_optimal_coverage_set( \@tlxl_objs ) if @tlxl_objs > 1;

#       foreach my $tlxl_obj (@tlxl_objs) {
#         next if defined $tlxl_obj->{Filter};
#         $tlxl_obj->{tlx} = create_tlx_object($tlxl_obj,$samobj,$samobj_brk);
#       }

#       foreach my $tlxl_obj (@tlxl_objs) {
#         next if defined $tlxl_obj->{Filter};
#         filter_orientation($tlxl_obj);
#         $stats->{orientation}++ unless defined $tlxl_obj->{Filter};
#       }

#       foreach my $tlxl_obj (@tlxl_objs) {
#         next if defined $tlxl_obj->{Filter};
#         refine_breaksite_alignment($brksite,$tlxl_obj->{tlx},$tmpbrkdir);
#       }

#       # my @tlxl_objs_unfiltered = grep {! defined $_->{Filter}} @tlxl_objs;

#       # find_optimal_coverage_set( \@tlxl_objs_unfiltered ) if @tlxl_objs_unfiltered > 1;

#       foreach my $tlxl_obj (@tlxl_objs) {
#         unless (defined $tlxl_obj->{Filter}) {
#           $stats->{primaln}++;
#         } else {
#           $stats->{primaln}++ if $tlxl_obj->{Filter} eq "MapQuality";
#           $stats->{primaln}++ if $tlxl_obj->{Filter} eq "SplitAlignment";
#         }

#         unless (defined $tlxl_obj->{Filter}) {
#           $stats->{splitaln}++;
#         } else {
#           $stats->{splitaln}++ if $tlxl_obj->{Filter} eq "MapQuality";
#         }

#         $stats->{mapqual}++ unless defined $tlxl_obj->{Filter};
#       }


      

#       foreach my $tlxl_obj (@tlxl_objs) {
#         next if defined $tlxl_obj->{Filter};
#         filter_correct_priming($tlxl_obj,$priming_threshold);
#         $stats->{priming}++ unless defined $tlxl_obj->{Filter};
#       }

      


#       foreach my $tlxl_obj (@tlxl_objs) {
#         if (defined $tlxl_obj->{Filter}) {
#           write_tlxl_entry($filt_tlxlfh,$tlxl_obj);
#           write_tlx_entry($filt_tlxfh,$tlxl_obj->{tlx}) if defined $tlxl_obj->{tlx};
#         } else {
#           write_tlxl_entry($tlxlfh,$tlxl_obj);
#           write_tlx_entry($tlxfh,$tlxl_obj->{tlx});
#           $stats->{final}++;
#         }
#       }




#       #   unless ( defined $aln->{orientation} && $aln->{orientation} > 0 ) { 
#       #     create_tlxl_entry($tlxlfh,$aln,$R1_brk_aln,$R2_brk_aln,"BadOrientation");
#       #     next;
#       #   }

#       #   unless ( defined $aln->{mapq} && $aln->{mapq} > 0 ) {
#       #     create_tlxl_entry($tlxlfh,$aln,$R1_brk_aln,$R2_brk_aln,"MapQuality");
#       #     next;
#       #   }

#       #   create_tlxl_entry($tlxlfh,$aln,$R1_brk_aln,$R2_brk_aln,$aln->{orientation});

#       #   $stats->{mapqual}++;
#       #   create_tlx_entry($tlxfh,$aln,$R1_brk_aln,$R2_brk_aln,$samobj,$samobj_brk);

#       # }

#     } else {

#       my $tlxl_obj = create_tlxl_object(undef,$R1_brk_aln,$R2_brk_aln,"Unaligned");
#       write_tlxl_entry($unal_tlxlfh,$tlxl_obj);
#     }

#   }

#   $tlxlfh->close;
#   $tlxfh->close;

# }



# sub create_tlx_object ($$$) {
#   my $tlxl_obj = shift;
#   my $sam = shift;
#   my $sam_brk = shift;
  
#   my $R1_aln = $tlxl_obj->{R1_aln};
#   my $R2_aln = $tlxl_obj->{R2_aln};
#   my $R1_brk_aln = $tlxl_obj->{R1_brk_aln};
#   my $R2_brk_aln = $tlxl_obj->{R2_brk_aln};

#   my $R1_Qstart = $tlxl_obj->{R1_Qstart};
#   my $R1_Qend = $tlxl_obj->{R1_Qend};
#   my $R2_Qstart = $tlxl_obj->{R2_Qstart};
#   my $R2_Qend = $tlxl_obj->{R2_Qend};
#   my $R1_BrkQstart = $tlxl_obj->{R1_BrkQstart};
#   my $R1_BrkQend = $tlxl_obj->{R1_BrkQend};
#   my $R2_BrkQstart = $tlxl_obj->{R2_BrkQstart};
#   my $R2_BrkQend = $tlxl_obj->{R2_BrkQend};

#   my $R1_Rstart = $tlxl_obj->{R1_Rstart};
#   my $R1_Rend = $tlxl_obj->{R1_Rend};
#   my $R2_Rstart = $tlxl_obj->{R2_Rstart};
#   my $R2_Rend = $tlxl_obj->{R2_Rend};
#   my $R1_BrkRstart = $tlxl_obj->{R1_BrkRstart};
#   my $R1_BrkRend = $tlxl_obj->{R1_BrkRend};
#   my $R2_BrkRstart = $tlxl_obj->{R2_BrkRstart};
#   my $R2_BrkRend = $tlxl_obj->{R2_BrkRend};

#   my $R1_seq = $tlxl_obj->{R1_seq};
#   my $R2_seq = $tlxl_obj->{R2_seq};

#   my $tlx_obj = {};

#   $tlx_obj->{Qname} = $tlxl_obj->{Qname};
#   $tlx_obj->{Orientation} = $tlxl_obj->{Orientation};
#   $tlx_obj->{MapQ} = $tlxl_obj->{MapQ};
#   $tlx_obj->{Rname} = $tlxl_obj->{Rname};
#   $tlx_obj->{Strand} = $tlxl_obj->{Strand};


#   switch ($tlx_obj->{Orientation}) {
    
#     case 1 {
      
#       $tlx_obj->{Rstart} = $tlx_obj->{Strand} == 1 ? $R1_Rstart : $R2_Rstart;
#       $tlx_obj->{Rend} = $tlx_obj->{Strand} == 1 ? $R2_Rend : $R1_Rend;
#       $tlx_obj->{Junction} = $tlx_obj->{Strand} == 1 ? $tlx_obj->{Rstart} : $tlx_obj->{Rend};


#       my $left = substr($R1_seq,0,$R1_Qstart - 1);
#       my $mid = merge_alignments($R1_aln,$R2_aln,$sam);
#       $mid = $tlx_obj->{Strand} == 1 ? $mid : reverseComplement($mid);
#       my $right = substr($R2_seq,$R2_Qend);
#       $tlx_obj->{Seq} = $left . $mid . $right;
#       $tlx_obj->{Qlen} = length($tlx_obj->{Seq});


#       $tlx_obj->{BrkRstart} = $R1_BrkRstart;
#       $tlx_obj->{BrkQstart} = $R1_BrkQstart;
#       $tlx_obj->{BrkRend} = $R1_BrkRend;
#       $tlx_obj->{BrkQend} = $R1_BrkQend;

#       # if ($R2_BrkRend > $R1_BrkRend) {
#       #   $tlx_obj->{BrkRend} = $R2_BrkRend;

#       #   if ($tlx_obj->{Strand} == 1) {
#       #     $tlx_obj->{BrkQend} = $R1_Qstart + ( $R2_Rstart - $R1_Rstart ) - ( $R2_Qstart - $R2_BrkQend );
#       #   } else {
#       #     $tlx_obj->{BrkQend} = $R1_Qstart + ( $R1_Rend - $R2_Rend ) - ( $R2_Qstart - $R2_BrkQend );
#       #   }

#       # } else {
#       #   $tlx_obj->{BrkRend} = $R1_BrkRend;
#       #   $tlx_obj->{BrkQend} = $R1_BrkQend;
#       # }
   
#       $tlx_obj->{Qstart} = length($left) + 1;
#       $tlx_obj->{Qend} = length($left.$mid);



#     }

#     case 2 {
      
#       $tlx_obj->{Rstart} = $tlx_obj->{Strand} == 1 ? $R1_Rstart : $R2_Rstart;
#       $tlx_obj->{Rend} = $tlx_obj->{Strand} == 1 ? $R2_Rend : $R1_Rend;
#       $tlx_obj->{Junction} = $tlx_obj->{Strand} == 1 ? $tlx_obj->{Rstart} : $tlx_obj->{Rend};



#       my $left = substr($R1_seq,0,$R1_Qstart - 1);
#       my $mid = merge_alignments($R1_aln,$R2_aln,$sam);
#       $mid = $tlx_obj->{Strand} == 1 ? $mid : reverseComplement($mid);
#       my $right = substr(reverseComplement($R2_seq),$R2_Qend);
#       $tlx_obj->{Seq} = $left . $mid . $right;
#       $tlx_obj->{Qlen} = length($tlx_obj->{Seq});



#       $tlx_obj->{BrkRstart} = $R1_BrkRstart;
#       $tlx_obj->{BrkRend} = $R1_BrkRend;
#       $tlx_obj->{BrkQstart} = $R1_BrkQstart;
#       $tlx_obj->{BrkQend} = $R1_BrkQend;
#       $tlx_obj->{Qstart} = length($left) + 1;
#       $tlx_obj->{Qend} = length($left.$mid);

#     }

#     case 3 {
      
#       $tlx_obj->{Rstart} = $R2_Rstart;
#       $tlx_obj->{Rend} = $R2_Rend;
#       $tlx_obj->{Junction} = $tlx_obj->{Strand} == 1 ? $tlx_obj->{Rstart} : $tlx_obj->{Rend};


#       my ($ol_aln,$R1_frag,$R2_frag,$R1_OLstart,$R1_OLend,$R2_OLstart,$R2_OLend);

#       # unless ($R1_brk_aln->end > $R2_brk_aln->start) {
#       $ol_aln = sw_align_pairs($R1_brk_aln->query->seq,$R2_brk_aln->query->seq,$tmpstitchdir);
#       print "\n".$tlx_obj->{Qname}." - sw alignment\n";
#       $R1_frag = $ol_aln->get_seq_by_pos(1);
#       $R2_frag = $ol_aln->get_seq_by_pos(2);
#       $R1_OLstart = $R1_frag->location_from_column(1)->start;
#       $R1_OLend = $R1_frag->location_from_column($ol_aln->length())->start;
#       $R2_OLstart = $R2_frag->location_from_column(1)->start;
#       $R2_OLend = $R2_frag->location_from_column($ol_aln->length())->start;
#       # }
      
#       if ($ol_aln->length > $R1_BrkRend - $R2_BrkRstart && ($R1_OLend > $R1_BrkQend + 3 || $R2_OLstart < $R2_BrkQstart - 3)) {
#         print(join(" ","overlapped",$R1_OLstart,$R1_OLend,$R2_OLstart,$R2_OLend)."\n");
#         (my $cons = $ol_aln->consensus_string(100)) =~ s/\?/N/g;
#         print $cons."\n";

#         my $left = substr($R1_seq,0,$R1_OLstart - 1);
#         my $mid = $cons;
#         my $right = substr($R2_seq,$R2_OLend);
#         $tlx_obj->{Seq} = $left . $mid . $right;
#         $tlx_obj->{Qlen} = length($tlx_obj->{Seq});

#         $tlx_obj->{BrkRstart} = $R1_BrkRstart;
#         $tlx_obj->{BrkQstart} = $R1_BrkQstart;

#         $tlx_obj->{BrkRend} = $R2_BrkRend;
#         $tlx_obj->{BrkQend} = length($left.$mid) + $R2_BrkQend - $R2_OLend;


#         $tlx_obj->{Qstart} = length($left.$mid) + $R2_Qstart - $R2_OLend;
#         $tlx_obj->{Qend} = length($left.$mid) + $R2_Qend - $R2_OLend;

#       } else {
#         print(join(" ","merging",$R1_BrkQend,$R2_BrkQstart)."\n");
#         my $left = substr($R1_seq,0,$R1_BrkQstart - 1);
#         my $mid = merge_alignments($R1_brk_aln,$R2_brk_aln,$sam_brk);
#         my $right = substr($R2_seq,$R2_BrkQend);

#         $tlx_obj->{Seq} = $left . $mid . $right;
#         $tlx_obj->{Qlen} = length($tlx_obj->{Seq});


#         $tlx_obj->{BrkRstart} = $R1_BrkRstart;
#         $tlx_obj->{BrkQstart} = $R1_BrkQstart;

#         $tlx_obj->{BrkRend} = $R2_BrkRend;
#         $tlx_obj->{BrkQend} = length($left.$mid);
     
#         $tlx_obj->{Qstart} = $tlx_obj->{BrkQend} + $R2_Qstart - $R2_BrkQend;
#         $tlx_obj->{Qend} = $tlx_obj->{BrkQend} + $R2_Qend - $R2_BrkQend;

#       }


#     }

#     case 4 {
      
#       $tlx_obj->{Rstart} = $R1_Rstart;
#       $tlx_obj->{Rend} = $R1_Rend;
#       $tlx_obj->{Junction} = $tlx_obj->{Strand} == 1 ? $tlx_obj->{Rstart} : $tlx_obj->{Rend};

#       $tlx_obj->{Seq} = $R1_seq;
#       $tlx_obj->{BrkRstart} = $R1_BrkRstart;
#       $tlx_obj->{BrkQstart} = $R1_BrkQstart;
#       $tlx_obj->{BrkRend} = $R1_BrkRend;
#       $tlx_obj->{BrkQend} = $R1_BrkQend;
#       $tlx_obj->{Qstart} = $R1_Qstart;
#       $tlx_obj->{Qend} = $R1_Qend;

#       $tlx_obj->{Qlen} = length($tlx_obj->{Seq});


#     }

#   }

#   #refine_breaksite_alignment($brksite,$tlx_obj,$tmpbrkdir);

  
#   return $tlx_obj;

# }

# sub create_tlxl_object ($$$;$) {
#   my $aln = shift;
#   my $R1_brk_aln = shift;
#   my $R2_brk_aln = shift;
#   my $filter = shift;

#   my $R1_aln = $aln->{R1_aln};
#   my $R2_aln = $aln->{R2_aln};

#   my $tlxl_obj = {};

#   $tlxl_obj->{Qname} = $R1_brk_aln->query->name;
#   $tlxl_obj->{Filter} = $filter if defined $filter;
#   $tlxl_obj->{R1_aln} = $R1_aln;
#   $tlxl_obj->{R2_aln} = $R2_aln;
#   $tlxl_obj->{R1_brk_aln} = $R1_brk_aln;
#   $tlxl_obj->{R2_brk_aln} = $R2_brk_aln;


#   $tlxl_obj->{R1_seq} = $R1_brk_aln->query->seq->seq;
#   $tlxl_obj->{R2_seq} = $R2_brk_aln->reversed ? $R2_brk_aln->query->seq->seq : reverseComplement($R2_brk_aln->query->seq->seq);
#   $tlxl_obj->{R1_Qlen} = length($tlxl_obj->{R1_seq});
#   $tlxl_obj->{R2_Qlen} = length($tlxl_obj->{R2_seq});

#   unless ($R1_brk_aln->unmapped) {
#     $tlxl_obj->{R1_BrkQstart} = $R1_brk_aln->query->start;
#     $tlxl_obj->{R1_BrkQend} = $R1_brk_aln->query->end;
#     $tlxl_obj->{R1_BrkRstart} = $R1_brk_aln->start;
#     $tlxl_obj->{R1_BrkRend} = $R1_brk_aln->end;
#     $tlxl_obj->{R1_BrkCigar} = $R1_brk_aln->cigar_str;
#   }
#   unless ($R2_brk_aln->unmapped) {
#     $tlxl_obj->{R2_BrkQstart} = $R2_brk_aln->query->start;
#     $tlxl_obj->{R2_BrkQend} = $R2_brk_aln->query->end;
#     $tlxl_obj->{R2_BrkRstart} = $R2_brk_aln->start;
#     $tlxl_obj->{R2_BrkRend} = $R2_brk_aln->end;
#     $tlxl_obj->{R2_BrkCigar} = $R2_brk_aln->cigar_str;
#   }

#   if (defined $R1_aln) {
#     $tlxl_obj->{Rname} = $R1_aln->seq_id;
#     $tlxl_obj->{Strand} = $R1_aln->strand;
#     $tlxl_obj->{MapQ} = $R1_aln->qual;

#     $tlxl_obj->{R1_Qstart} = $R1_aln->reversed ? $R1_aln->l_qseq - $R1_aln->query->end + 1 : $R1_aln->query->start;
#     $tlxl_obj->{R1_Qend} = $R1_aln->reversed ? $R1_aln->l_qseq - $R1_aln->query->start + 1 : $R1_aln->query->end;
#     $tlxl_obj->{R1_Rstart} = $R1_aln->start;
#     $tlxl_obj->{R1_Rend} = $R1_aln->end;
#     $tlxl_obj->{R1_Cigar} = $R1_aln->cigar_str;
#   }

#   if (defined $R2_aln) {
#     if ($R2_aln->munmapped) {
#       $tlxl_obj->{Rname} = $R2_aln->seq_id;
#       $tlxl_obj->{Strand} = -1 * $R2_aln->strand;
#       $tlxl_obj->{MapQ} = $R2_aln->qual;
#     }

#     unless (defined $R1_aln && ! $R1_aln->proper_pair) {
#       $tlxl_obj->{R2_Qstart} = $R2_aln->reversed ? $R2_aln->query->start : $R2_aln->l_qseq - $R2_aln->query->end + 1;
#       $tlxl_obj->{R2_Qend} = $R2_aln->reversed ? $R2_aln->query->end : $R2_aln->l_qseq - $R2_aln->query->start + 1;
#       $tlxl_obj->{R2_Rstart} = $R2_aln->start;
#       $tlxl_obj->{R2_Rend} = $R2_aln->end;
#       $tlxl_obj->{R2_Cigar} = $R2_aln->cigar_str;
#     }
#   }


#   # Test for orientation #1
#   if ($R1_brk_aln->proper_pair && defined $R1_aln && $R1_aln->proper_pair) {

#     if ($tlxl_obj->{R1_BrkQstart} <= $tlxl_obj->{R1_Qstart} && $tlxl_obj->{R2_BrkQstart} <= $tlxl_obj->{R2_Qstart}) {
#       $tlxl_obj->{Orientation} = 1;
#     }

#   }

#   # Test for orientation #2
#   if ( ! defined $tlxl_obj->{Orientation} && ! $R1_brk_aln->unmapped && defined $R1_aln && $R1_aln->proper_pair) {

#     if ($tlxl_obj->{R1_BrkQstart} <= $tlxl_obj->{R1_Qstart}) {
#       $tlxl_obj->{Orientation} = 2;
#     }

#   }

#   # Test for orientation #3
#   if ( ! defined $tlxl_obj->{Orientation} && $R1_brk_aln->proper_pair && defined $R2_aln) {

#     if ($tlxl_obj->{R2_BrkQstart} <= $tlxl_obj->{R2_Qstart}) {
#       $tlxl_obj->{Orientation} = 3;
#     }

#   }

#   # Test for orientation #4
#   if ( ! defined $tlxl_obj->{Orientation} && ! $R1_brk_aln->unmapped && defined $R1_aln) {

#     unless (! $R2_brk_aln->unmapped) { # && $R2_brk_aln->start <= $R1_brk_aln->end) {

#       if ($tlxl_obj->{R1_BrkQstart} <= $tlxl_obj->{R1_Qstart}) {
#         $tlxl_obj->{Orientation} = 4;
#       }
#     }

#   }

#   $tlxl_obj->{Orientation} = 0 unless defined $tlxl_obj->{Orientation};


#   return $tlxl_obj;


# }

# sub write_tlx_entry ($$) {

#   my $fh = shift;
#   my $entry = shift;

#   my @tlx_header = tlx_header();

#   $fh->print(join("\t",map(check_undef($_,""),@{$entry}{@tlx_header}))."\n");

# }

sub write_tlxl_entry ($$) {

  my $fh = shift;
  my $entry = shift;

  my @tlxl_header = tlxl_header();

  $fh->print(join("\t",map(check_undef($_,""),@{$entry}{@tlxl_header}))."\n");

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
                            "threads=i" => \$threads,
                            "bt2opt=s" => \$user_bowtie_opt,
                            "bt2brkopt=s" => \$user_bowtie_breaksite_opt,
				            				"help" => \$help
				            			) ;
	
	usage() if ($help);

  #Check options

  croak "Error: cannot read sequence files" unless (-r $read1 & -r $read2);
  croak "Error: working directory does not exist" unless (-d $workdir);
  



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
$arg{"--threads","Number of threads to run bowtie on","$threads"}
$arg{"--help","This helpful help screen."}


EOF

exit 1;
}
