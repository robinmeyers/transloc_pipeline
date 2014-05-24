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
use List::Util qw(min max);
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


my $GENOME_DB = $ENV{'GENOME_DB'};
defined $GENOME_DB or croak "Error: set environment variable GENOME_DB";

my $BOWTIE2_INDEXES = $ENV{'BOWTIE2_INDEXES'};
defined $BOWTIE2_INDEXES or croak "Error: set environment variable BOWTIE2_INDEXES";





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
sub write_parameters_file;
sub align_to_breaksite;
sub align_to_adapter;
sub align_to_genome;
sub process_alignments;
sub find_optimal_coverage_set ($$);
sub process_optimal_coverage_set ($$$);
sub score_edge ($;$);
sub deduplicate_junctions;
sub sort_junctions;
sub post_process_junctions;
sub write_stats_file;
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
my $break_fa;
my $breakcoord;
my $prim_fa;
my $adapt_fa;
my $cut_fa;
my $bowtie_threads = 4;
my $dedup_threads = 4;

my $skip_alignment;
my $skip_process;
my $skip_dedup;
my $no_dedup;
my $no_clean;




# OL_mult must be set the same as --ma in bowtie2 options
my $OL_mult = 2;
my $maximum_brk_start_dif = 20;
my $Dif_mult = 2;
my $Brk_pen_default = 40;
my $max_frag_len = 1500;
my $PE_pen_default = 20;

my $min_bases_after_primer = 10;
my $max_bases_after_cutsite = 10;
my $mapq_ol_thresh = 0.90;
my $mapq_mismatch_thresh_int = 1.5;
my $mapq_mismatch_thresh_coef = 0.01;

my $dedup_offset_dist = 2;
my $dedup_break_dist = 2;





my $user_bowtie_opt = "";
my $user_bowtie_adapter_opt = "";
my $user_bowtie_breaksite_opt = "";

my $use_current_tlx;

# Global variables
my @tlxl_header = tlxl_header();
my @tlx_header = tlx_header();
my @tlx_filter_header = tlx_filter_header();


my $stats = {};
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
$stats->{sequentialjuncs} = 0;
$stats->{sequential_reads} = 0;
$stats->{dedup} = 0;
$stats->{final} = 0;


#
# Start of Program
#
parse_command_line;


if ($bowtie_threads == 0) {
  my $info = Sys::Info->new;
  my $cpu  = $info->device( 'CPU'  );
  $bowtie_threads = 2*$cpu->count;
}

# my $default_bowtie_adapter_opt = "--local -D 20 -R 3 -N 1 -L 10 -i C,6 --score-min C,20 -p $bowtie_threads --no-unal --reorder -t";
# my $default_bowtie_breaksite_opt = "--local -D 15 -R 2 -N 0 -L 20 -i C,8 --score-min C,50 --mp 10,2 --rfg 10,2 --rdg 10,2 -p $bowtie_threads -k 20 --no-unal --reorder -t";
# my $default_bowtie_opt = "--local -D 15 -R 2 -N 0 -L 20 -i C,8 --score-min C,50 --mp 10,2 --rfg 10,2 --rdg 10,2 -p $bowtie_threads -k 50 --reorder -t";

my $default_bowtie_adapter_opt = "--local -D 20 -R 3 -N 1 -L 10 -i C,6 --ma 2 --mp 6,2 --np 1 --rdg 5,3 --rfg 5,3 --score-min C,20 --no-unal -p $bowtie_threads --reorder -t";
my $default_bowtie_breaksite_opt = "--local -D 20 -R 3 -N 1 -L 20 -i C,8 --ma 2 --mp 10,6 --np 2 --rdg 6,4 --rfg 6,4 --score-min C,50 -k 5 --no-unal -p $bowtie_threads  --reorder -t";
my $default_bowtie_opt = "--local -D 20 -R 3 -N 1 -L 20 -i C,8 --ma 2 --mp 10,6 --np 2 --rdg 6,4 --rfg 6,4 --score-min C,50 -k 20 -p $bowtie_threads --reorder -t";

my $bt2_break_opt = manage_program_options($default_bowtie_breaksite_opt,$user_bowtie_breaksite_opt);
my $bt2_adapt_opt = manage_program_options($default_bowtie_adapter_opt,$user_bowtie_adapter_opt);
my $bt2_opt = manage_program_options($default_bowtie_opt,$user_bowtie_opt);


croak "Error: cannot find match award in bowtie2 options" unless $bt2_opt =~ /-ma (\d+)/;
my $match_award = $1;
croak "Error: cannot find mismatch penalty in bowtie2 options" unless $bt2_opt =~ /-mp (\d+),\d+/;
my $mismatch_penalty = $1;

carp "Warning: match award in bowtie2 does not equal OCS overlap penalty" unless $match_award eq $OL_mult;


my $t0 = [gettimeofday];

#check_working_dir;

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
my $filt_tlxfile = "${expt_stub}_filtered.tlx";
my $unjoin_tlxfile = "${expt_stub}_unjoined.tlx";
my ($tlxlfh,$tlxfh,$filt_tlxfh,$unjoin_tlxfh);

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
}

# Calculate threshold for uncut/unjoin filter
$brksite->{joining_threshold} = $brksite->{endogenous} ?
                                  ($brksite->{strand} eq "+" ?
                                    $brksite->{end} + $max_bases_after_cutsite :
                                    $brksite->{start} - $max_bases_after_cutsite ) :
                                  $brksite->{breakcoord} + $max_bases_after_cutsite;

# Calculate thresholds for mispriming filters
croak "Error: could not find primer within breaksite sequence" unless $brksite->{primer_start} > 0;
$brksite->{priming_threshold} = $brksite->{primer_start} + $brksite->{aln_strand} * ($brksite->{primer}->length - 1 + $min_bases_after_primer);

# Read in adapter sequence
my $adpt_io = Bio::SeqIO->new(-file => $adapt_fa, -format => 'fasta');
my $adaptseq = $adpt_io->next_seq();

# Read in cutter sequence
my $cutseq;
if (defined $cut_fa) {
  my $cut_io = Bio::SeqIO->new(-file => $cut_fa, -format => 'fasta');
  $cutseq = $cut_io->next_seq();
}


write_parameters_file;



unless ($skip_alignment || $skip_process || $skip_dedup) {

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


unless ($skip_process || $skip_dedup) {
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

} else {
  croak "Error: could not find tlx file when skipping processing step" unless -r $tlxfile;
}

unless ($skip_dedup || $no_dedup) {

  deduplicate_junctions;

} else {
  $stats->{dedup} = $stats->{sequentialjuncs};
  croak "Error: could not find tlx file when skipping dedup step" unless -r $tlxfile;
}

sort_junctions;


write_stats_file unless ($skip_process || $skip_dedup);


post_process_junctions;


my $t1 = tv_interval($t0);

printf("\nFinished all processes in %.2f seconds.\n", $t1);

print("\nStats\n".join("\n","Total Reads: ".$stats->{totalreads},
                          "Aligned: ".$stats->{aligned},
                          "Junctions: ".$stats->{junctions}." (".$stats->{junc_reads}.")",
                          "Priming: ".$stats->{priming}." (".$stats->{prim_reads}.")",
                          "FrequentCutter: ".$stats->{freqcut}." (".$stats->{freq_reads}.")",
                          "MapQuality: ".$stats->{mapqual}." (".$stats->{mapq_reads}.")",
                          "Breaksite: ".$stats->{breaksite}." (".$stats->{break_reads}.")",
                          "SequentialJuncs: ".$stats->{sequentialjuncs}." (".$stats->{sequential_reads}.")",
                          "Dedup: ".$stats->{dedup})."\n");

clean_up unless $no_clean;


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

  print "\nReading alignments\n";
  my ($R1_brk_iter,$R2_brk_iter);
  my ($R1_adpt_iter,$R2_adpt_iter);
  my ($R1_iter,$R2_iter);
  my ($next_R1_brk_aln,$next_R2_brk_aln);
  my ($next_R1_adpt_aln,$next_R2_adpt_aln);
  my ($next_R1_aln,$next_R2_aln);

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

    my @R1_alns = ();
    my @R2_alns = ();

       
    push(@R1_alns,$next_R1_aln);
    undef $next_R1_aln;

    while(my $aln = $R1_iter->next_seq) {
      $next_R1_aln = wrap_alignment("R1",$aln);
      last unless $next_R1_aln->{Qname} eq $qname;
      push(@R1_alns,$next_R1_aln);
      undef $next_R1_aln;
    }

    if (defined $next_R2_aln && $next_R2_aln->{Qname} eq $qname) {
      push(@R2_alns,$next_R2_aln);
      undef $next_R2_aln;
      while(my $aln = $R2_iter->next_seq) {
        $next_R2_aln = wrap_alignment("R2",$aln);
        last unless $next_R2_aln->{Qname} eq $qname;
        push(@R2_alns,$next_R2_aln);
        undef $next_R2_aln;
      }
    }

    # Read in Breaksite alignments only if brksite is non-endogenous
    unless ($brksite->{endogenous}) {
      if (defined $next_R1_brk_aln && $next_R1_brk_aln->{Qname} eq $qname) {
        push(@R1_alns,$next_R1_brk_aln);
        undef $next_R1_brk_aln;
        while(my $aln = $R1_brk_iter->next_seq) {
          $next_R1_brk_aln = wrap_alignment("R1",$aln);
          last unless $next_R1_brk_aln->{Qname} eq $qname;
          push(@R1_alns,$next_R1_brk_aln);
          undef $next_R1_brk_aln;
        }
      }

      if (defined $next_R2_brk_aln && $next_R2_brk_aln->{Qname} eq $qname) {
        push(@R2_alns,$next_R2_brk_aln);
        undef $next_R2_brk_aln;
        while(my $aln = $R2_brk_iter->next_seq) {
          $next_R2_brk_aln = wrap_alignment("R2",$aln);
          last unless $next_R2_brk_aln->{Qname} eq $qname;
          push(@R2_alns,$next_R2_brk_aln);
          undef $next_R2_brk_aln;
        }
      }
    }


    # Read in Adapter alignments
    if (defined $next_R1_adpt_aln && $next_R1_adpt_aln->{Qname} eq $qname) {
      push(@R1_alns,$next_R1_adpt_aln);
      undef $next_R1_adpt_aln;
      while(my $aln = $R1_adpt_iter->next_seq) {
        $next_R1_adpt_aln = wrap_alignment("R1",$aln);
        last unless $next_R1_adpt_aln->{Qname} eq $qname;
        push(@R1_alns,$next_R1_adpt_aln);
        undef $next_R1_adpt_aln;
      }
    }

    if (defined $next_R2_adpt_aln && $next_R2_adpt_aln->{Qname} eq $qname) {
      push(@R2_alns,$next_R2_adpt_aln);
      undef $next_R2_adpt_aln;
      while(my $aln = $R2_adpt_iter->next_seq) {
        $next_R2_adpt_aln = wrap_alignment("R2",$aln);
        last unless $next_R2_adpt_aln->{Qname} eq $qname;
        push(@R2_alns,$next_R2_adpt_aln);
        undef $next_R2_adpt_aln;
      }
    }

    

    # print "\nbefore ocs ". Dumper(\@R2_alns) if @R2_alns < 2;


    my $OCS = find_optimal_coverage_set(\@R1_alns,\@R2_alns);

    
    
    process_optimal_coverage_set($OCS,\@R1_alns,\@R2_alns);

    


  }


}


sub find_optimal_coverage_set ($$) {

  # print "finding OCS\n";
  # my $t0 = [gettimeofday];



  my $R1_alns_ref = shift;
  my $R2_alns_ref = shift;

  my @graph = ();
  my $OCS_ptr;

  my @R1_alns = sort {$a->{Qstart} <=> $b->{Qstart} || $a->{Rstart} <=> $b->{Rstart}} @$R1_alns_ref;
  my @R2_alns = sort {$a->{Qstart} <=> $b->{Qstart} || $a->{Rstart} <=> $b->{Rstart}} @$R2_alns_ref;

  # print "\nafter sort ".Dumper(\@R2_alns) if @R2_alns < 2;

  foreach my $R1_aln (@R1_alns) {

    next if $R1_aln->{Unmapped};

    my $graphsize = scalar @graph;

    my $new_node = {R1 => $R1_aln};



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
      # print $R1_aln->{Qname}."\n" if ! defined $R2_aln->{Unmapped};
      # print Dumper($R2_aln) if $R2_aln->{Unmapped};

      next unless defined $R2_aln && ! $R2_aln->{Unmapped};

      # print "testing proper-pairedness with R2:\n";
      # print_aln($R2_aln_wrap);

      next unless pair_is_proper($R1_aln,$R2_aln,$max_frag_len);

      my $new_pe_node = {R1 => $R1_aln, R2 => $R2_aln};



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
    next unless defined $R2_aln && ! $R2_aln->{Unmapped};

    my $new_node = {R2 => $R2_aln};


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
    my @unmapped_OCS = ( { R1 => { Qname => $R1_alns_ref->[0]->{Qname},
                                   Seq => $R1_alns_ref->[0]->{Seq},
                                   Qual => $R1_alns_ref->[0]->{Qual},
                                   Unmapped => 1 },
                           R2 => { Qname => $R2_alns_ref->[0]->{Qname},
                                   Seq => $R2_alns_ref->[0]->{Seq},
                                   Qual => $R2_alns_ref->[0]->{Qual},
                                   Unmapped => 1 } } );
    return \@unmapped_OCS;

    next;
  }

  my @OCS;

  while (defined $OCS_ptr->{back_ptr}) {
    unshift(@OCS,$OCS_ptr);
    $OCS_ptr->{R1_Rgap} = find_genomic_distance($OCS_ptr->{back_ptr}->{R1},$OCS_ptr->{R1},$brksite)
      if defined $OCS_ptr->{R1} && defined $OCS_ptr->{back_ptr}->{R1};
    $OCS_ptr->{R2_Rgap} = find_genomic_distance($OCS_ptr->{back_ptr}->{R2},$OCS_ptr->{R2},$brksite)
      if defined $OCS_ptr->{R2} && defined $OCS_ptr->{back_ptr}->{R2};
    $OCS_ptr->{back_ptr}->{score} = $OCS_ptr->{score};
    $OCS_ptr = $OCS_ptr->{back_ptr};
  }
  unshift(@OCS,$OCS_ptr) if defined $OCS_ptr;


  return \@OCS;


}




sub process_optimal_coverage_set ($$$) {

  
    
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
  my $junctions = filter_unjoined($tlxls,$brksite);

  $stats->{junctions} += $junctions;
  $stats->{junc_reads}++ if $junctions > 0;

  # print "filter mispriming\n";
  my $correct_priming = filter_mispriming($tlxls,$brksite);

  $stats->{priming} += $correct_priming;
  $stats->{prim_reads}++ if $correct_priming > 0;

  # print "filter frequent cutter\n";
  my $no_freq_cutter = filter_freq_cutter($tlxls,$cutseq);

  $stats->{freqcut} += $no_freq_cutter;
  $stats->{freq_reads}++ if $no_freq_cutter > 0;

  # print "filter map quality\n";
  my $quality_maps = filter_mapping_quality($tlxls,$R1_alns,$R2_alns,
                                      $mapq_ol_thresh,$mapq_mismatch_thresh_int,$mapq_mismatch_thresh_coef,
                                      $max_frag_len,$match_award,$mismatch_penalty);

  $stats->{mapqual} += $quality_maps;
  $stats->{mapq_reads}++ if $quality_maps > 0;


  # print "filter breaksite\n";
  my $outside_breaksite = filter_breaksite($tlxls);

  $stats->{breaksite} += $outside_breaksite;
  $stats->{break_reads}++ if $outside_breaksite > 0;


  # print "filter sequential juctions\n";
  my $primary_junction = filter_sequential_junctions($tlxls);

  $stats->{sequentialjuncs} += $primary_junction;
  $stats->{sequential_reads}++ if $primary_junction > 0;


  write_tlxls($tlxls);


}

sub write_tlxls ($) {
  my $tlxls = shift;
  # print "writing stuff\n";
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
    my $OL_correction;

    if ( defined $node1->{R2} ) {
      return undef if defined $node2->{R1};
      return undef unless defined $node2->{R2};
      return undef if $node1->{R2}->{Rname} eq "Adapter";
      return undef if $node2->{R2} == $node1->{R2};
      return undef unless $node2->{R2}->{Qstart} > $node1->{R2}->{Qstart};
      return undef unless $node2->{R2}->{Qend} > $node1->{R2}->{Qend};
      $Rname1 = $node1->{R2}->{Rname};
      $Strand1 = $node1->{R2}->{Strand};
      $Rname2 = $node2->{R2}->{Rname};
      $Strand2 = $node2->{R2}->{Strand};
      # if ($Rname1 eq "Breaksite" && $Rname2 eq "Breaksite" && $Strand1 == 1 && $Strand2 == 1) {
      #   return undef unless $node2->{R2}->{Rstart} > $node1->{R2}->{Rstart};
      #   return undef unless $node2->{R2}->{Rend} > $node1->{R2}->{Rend};
      # }
      $R2_Qgap = $node2->{R2}->{Qstart} - $node1->{R2}->{Qend} - 1;
      $R2_Rdist = find_genomic_distance($node1->{R2},$node2->{R2},$brksite);
      $Len1 += $node1->{R2}->{Qend} - $node1->{R2}->{Qstart} + 1;
      $Len2 += $node2->{R2}->{Qend} - $node2->{R2}->{Qstart} + 1;
      $Junc1 = $Strand1 == 1 ? $node1->{R2}->{Rend} : $node1->{R2}->{Rstart};
      $Junc2 = $Strand2 == 1 ? $node2->{R2}->{Rstart} : $node1->{R2}->{Rend};
      $OL_correction = $OL_mult * max(0,-$R2_Qgap);
    } else {
      return undef unless defined $node2->{R1};
      return undef if $node1->{R1}->{Rname} eq "Adapter";
      return undef if $node2->{R1} == $node1->{R1};
      return undef unless $node2->{R1}->{Qstart} > $node1->{R1}->{Qstart};
      return undef unless $node2->{R1}->{Qend} > $node1->{R1}->{Qend};
      $Rname1 = $node1->{R1}->{Rname};
      $Strand1 = $node1->{R1}->{Strand};
      $Rname2 = $node2->{R1}->{Rname};
      $Strand2 = $node2->{R1}->{Strand};
      # if ($Rname1 eq "Breaksite" && $Rname2 eq "Breaksite" && $Strand1 == 1 && $Strand2 == 1) {
      #   return undef unless $node2->{R1}->{Rstart} > $node1->{R1}->{Rstart};
      #   return undef unless $node2->{R1}->{Rend} > $node1->{R1}->{Rend};
      # }
      $R1_Qgap = $node2->{R1}->{Qstart} - $node1->{R1}->{Qend} - 1;
      $R1_Rdist = find_genomic_distance($node1->{R1},$node2->{R1},$brksite);
      $Len1 += $node1->{R1}->{Qend} - $node1->{R1}->{Qstart} + 1;
      $Len2 += $node2->{R1}->{Qend} - $node2->{R1}->{Qstart} + 1;
      $Junc1 = $Strand1 == 1 ? $node1->{R1}->{Rend} : $node1->{R1}->{Rstart};
      $Junc2 = $Strand2 == 1 ? $node2->{R1}->{Rstart} : $node1->{R1}->{Rend};
      $OL_correction = $OL_mult * max(0,-$R1_Qgap);

    }

    # my $totalOverlap = -min($R1_Qgap,0) -min($R2_Qgap,0);
    # # return undef if $totalOverlap > $ol_thresh * $Len1 || $totalOverlap > $ol_thresh * $Len2;
    # return undef if $totalOverlap > $Len1 - 20  || $totalOverlap > $Len2 - 20;



    

    my $Brk_pen;


    if ($Rname2 eq "Adapter") {
      
      $Brk_pen = 0;

    } elsif (defined $R1_Rdist) {

      $Brk_pen = $R1_Rdist > 1 || $R1_Rdist < 0 ? $Brk_pen_default : 0;

    } elsif (defined $R2_Rdist) {
      $Brk_pen = $R2_Rdist > 1 || $R2_Rdist < 0 ? $Brk_pen_default : 0;

    } else {

      $Brk_pen = $Brk_pen_default;

    }



    my $R1_AS = defined $node2->{R1} ? $node2->{R1}->{AS} : 0;
    my $R2_AS = defined $node2->{R2} ? $node2->{R2}->{AS} : 0;
    my $PEgap;
    my $PEgap_pen;

    if (defined $node2->{R1} && defined $node2->{R2}) {
      $PEgap = $node2->{R1}->{Strand} == 1 ? $node2->{R2}->{Rstart} - $node2->{R1}->{Rend} : $node2->{R1}->{Rstart} - $node2->{R1}->{Rend};
    }

    # $PEgap_pen = defined $PEgap && $PEgap > 1 ? $Brk_pen_min + $Brk_pen_mult * log10($PEgap)**$Brk_pen_power : 0;
    $PEgap_pen = defined $PEgap && $PEgap > 1 ? $PE_pen_default : 0;
    # print $node1->{R1}->{Qname}." - $PEgap - $PEgap_pen\n" if defined $PEgap;

    $score = $node1->{score} + $R1_AS + $R2_AS - $PEgap_pen - $Brk_pen - $OL_correction;

    # my $qname = defined $node1->{R1} ? $node1->{R1}->{Qname} : $node1->{R2}->{Qname};
    # print "$qname PE: $score\n" if $Rname2 eq "Adapter" && defined $node2->{R1} && defined $node2->{R2};
    # print "$qname SE1: $score\n" if $Rname2 eq "Adapter" && (defined $node2->{R1} && !defined $node2->{R2});
    # print "$qname SE2: $score\n" if $Rname2 eq "Adapter" && (!defined $node2->{R1} && defined $node2->{R2});

    
  } else {

    return undef unless defined $node1->{R1};
    return undef unless $node1->{R1}->{Rname} eq $brksite->{aln_name};
    return undef unless $node1->{R1}->{Strand} == $brksite->{aln_strand};

    my $brk_start_dif = $brksite->{aln_strand} == 1 ?
                          abs($node1->{R1}->{Rstart} - $brksite->{primer_start}) :
                          abs($node1->{R1}->{Rend} - $brksite->{primer_start}) ;

    return undef unless $brk_start_dif < $maximum_brk_start_dif;


    my $R1_AS = $node1->{R1}->{AS};
    my $R2_AS = defined $node1->{R2} ? $node1->{R2}->{AS} : 0;
    my $PEgap;
    my $PEgap_pen;


    if (defined $node1->{R1} && defined $node1->{R2}) {
      $PEgap = $node1->{R1}->{Strand} == 1 ? $node1->{R2}->{Rstart} - $node1->{R1}->{Rend} : $node1->{R1}->{Rstart} - $node1->{R1}->{Rend};
    }

    # $PEgap_pen = defined $PEgap && $PEgap > 1 ? $Brk_pen_min + $Brk_pen_mult * log10($PEgap)**$Brk_pen_power : 0;
    $PEgap_pen = defined $PEgap && $PEgap > 1 ? $PE_pen_default : 0;
    # print $node1->{R1}->{Qname}." - $PEgap - $PEgap_pen\n" if defined $PEgap;

    $score = $R1_AS + $R2_AS - $PEgap_pen - $Dif_mult * $brk_start_dif;

  }

  return $score;

}

sub deduplicate_junctions {



  my $dedup_cmd = join(" ","$FindBin::Bin/../R/TranslocDedup.R",
                          $tlxfile,
                          $dedup_output,
                          "cores=$dedup_threads",
                          "offset_dist=$dedup_offset_dist",
                          "break_dist=$dedup_break_dist") ;

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

  unlink $tlxbak;

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

    System("TranslocPlot.R $tlxfile $pdf_plot binsize=2000000 strand=2 assembly=$assembly " .
            "brkchr=$brk_chr brksite=" .($brk_strand eq "+" ? $brk_end : $brk_start) ." brkstrand=" . ($brk_strand eq "+" ? "1" : "-1") );


    (my $pdf_plot_breaksite = $pdf_plot) =~ s/\.pdf/_brksite.pdf/;

    System("TranslocPlot.R $tlxfile $pdf_plot_breaksite strand=0 assembly=$assembly " .
          "brkchr=$brk_chr brksite=" .($brk_strand eq "+" ? $brk_end : $brk_start) ." brkstrand=" . ($brk_strand eq "+" ? "1" : "-1") .
          " chr=$brk_chr rmid=" .($brk_strand eq "+" ? $brk_end : $brk_start) ." rwindow=5000 binnum=100 plottype=linear" );
  }
}

sub write_parameters_file {
  my $paramsfh = IO::File->new(">$paramsfile");

  $paramsfh->print("$local_time\t$commandline\n\n");

  $paramsfh->print(join("\n", map {join("\t",@$_)}
    (["Alignment Options"],
    ["Genome Bowtie2 Options", $bt2_opt],
    ["Breaksite Bowtie2 Options", $bt2_break_opt],
    ["Adapter Bowtie2 Options", $bt2_adapt_opt],
    [],
    ["Optimal Query Coverage Options"],
    ["Break Penalty",$Brk_pen_default],
    ["PE Gap Penalty",$PE_pen_default],
    ["Overlap Penalty",$OL_mult],
    ["Max Bait Offset",$maximum_brk_start_dif],
    ["Bait Offset Penalty",$Dif_mult],
    [],
    ["Filtering Options"],
    ["Max Bp After Cutsite", $max_bases_after_cutsite],
    ["Min Bp After Primer", $min_bases_after_primer],
    ["MapQuality Overlap Threshold", $mapq_ol_thresh],
    ["MapQuality Mismatch Threshold Intercept", $mapq_mismatch_thresh_int],
    ["MapQuality Mismatch Threshold Coefficient", $mapq_mismatch_thresh_coef],
    [],
    ["Dedup Options"],
    ["Max Prey Offset Distance",$dedup_offset_dist],
    ["Max Bait Junction Distance",$dedup_break_dist])));

}

sub write_stats_file {

  my $statsfh = IO::File->new(">$statsfile");

  $statsfh->print(join("\t","TotalReads",
                            "Aligned",
                            "Junctions",
                            "Priming",
                            "FrequentCutter",
                            "MappingQuality",                            
                            "Breaksite",
                            "SequentialJuncs",
                            "DeDup")."\n");

  $statsfh->print(join("\t",$stats->{totalreads},
                            $stats->{aligned},
                            $stats->{junctions}." (".$stats->{junc_reads}.")",
                            $stats->{priming}." (".$stats->{prim_reads}.")",
                            $stats->{freqcut}." (".$stats->{freq_reads}.")",
                            $stats->{mapqual}." (".$stats->{mapq_reads}.")",                            
                            $stats->{breaksite}." (".$stats->{break_reads}.")",
                            $stats->{sequentialjuncs}." (".$stats->{sequential_reads}.")",
                            $stats->{dedup})."\n");

  $statsfh->close;

}

sub clean_up {

  print "\nCleaning up\n";

  unlink glob "${expt_stub}*.sam";

}

sub parse_command_line {
	my $help;

	usage() if (scalar @ARGV == 0);

  $commandline = join(" ",@ARGV);
  $local_time = localtime;

	my $result = GetOptions ( "read1=s" => \$read1,
                            "read2=s" => \$read2,
                            "assembly=s" => \$assembly,
                            "chr=s" => \$brk_chr,
                            "start=i" => \$brk_start,
                            "end=i" => \$brk_end,
                            "strand=s" => \$brk_strand,
                            "workdir=s" => \$workdir,
                            "mid=s" => \$mid_fa,
                            "primer=s" => \$prim_fa,
                            "breakseq=s" => \$break_fa,
                            "breaksite=i" => \$breakcoord,
                            "adapter=s" => \$adapt_fa,
                            "cutter=s" => \$cut_fa,
                            "threads-bt=i" => \$bowtie_threads,
                            "threads-dedup=i" => \$dedup_threads,
                            "skip-align" => \$skip_alignment,
                            "skip-process" => \$skip_process,
                            "skip-dedup" => \$skip_dedup,
                            "no-dedup" => \$no_dedup,
                            "mapq-ol=f" => \$mapq_ol_thresh,
                            "mapq-mm-int=f" => \$mapq_mismatch_thresh_int,
                            "mapq-mm-coef=f" => \$mapq_mismatch_thresh_coef,
                            "priming-bp=i" => \$min_bases_after_primer,
                            # "bowtie-opt=s" => \$user_bowtie_opt,
				            				"help" => \$help
				            			) ;

	
	usage() if ($help);

  #Check options

  croak "Error: must specify --read1" unless defined $read1;
  croak "Error: cannot read read1 file" unless -r $read1;
  croak "Error: cannot read read2 file" if defined $read2 && ! -r $read2;
  croak "Error: working directory does not exist" unless (-d $workdir);

  croak "Error: priming-bp must be a positive integer" unless $min_bases_after_primer > 0;
  croak "Error: mapq-ol must be a fraction between 0 and 1" if $mapq_ol_thresh < 0 || $mapq_ol_thresh > 1;
  # croak "Error: mapq-score must be a fraction between 0 and 1" if $mapq_score_thresh < 0 || $mapq_score_thresh > 1;
  
	exit unless $result;
}


sub usage()
{
print<<EOF;
TranslocPipeline.pl, by Robin Meyers, 2013


Usage: $0
        --read1 FILE --read2 FILE --workdir DIR --assembly (mm9|hg19)
        --chr chrN --start INT --end INT --strand (+|-) --MID FILE
        --primer FILE [--breaksite FILE] --adapter FILE --cutter FILE
        [--option VAL] [--flag] [--help]

Arguments set automatically by TranslocWrapper.pl:
$arg{"--read1","Input sequence file"}
$arg{"--read2","Input sequence file"}
$arg{"--workdir","Directory for results files - note: default Bowtie output goes to working directory"}
$arg{"--assembly","Genome assembly to align reads to"}
$arg{"--chr"," "}
$arg{"--start"," "}
$arg{"--end"," "}
$arg{"--strand"," "}
$arg{"--mid"," "}
$arg{"--primer"," "}
$arg{"--breakseq"," "}
$arg{"--breaksite"," "}
$arg{"--adapter"," "}
$arg{"--cutter"," "}


Arguments set manually with --pipeline-opt in TranslocWrapper.pl (defaults in parentheses):
$arg{"--threads-bt","Number of threads to run bowtie on",$bowtie_threads}
$arg{"--threads-dedup","Number of threads to run dedup on",$dedup_threads}
$arg{"--skip-align"," "}
$arg{"--skip-process"," "}
$arg{"--skip-dedup"," "}
$arg{"--no-dedup"," "}
$arg{"--mapq-ol","",$mapq_ol_thresh}
$arg{"--priming-bp","",$min_bases_after_primer}


$arg{"--help","This helpful help screen."}


EOF



exit 1;
}
