#!/usr/bin/perl


use strict;
use warnings;
use Getopt::Long;
use Carp;
use IO::Handle;
use IO::File;
use Text::CSV;
use File::Basename;
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
sub align_to_genome;
sub make_tlxl;


# Global flags and arguments, 
# Set by command line arguments
my $read1;
my $read2;
my $workdir;
my $assembly;
my $threads = 4;

my $user_bowtie_opt = "";
my $user_bowtie_breaksite_opt = "";

# Global variables 



#
# Start of Program
#

parse_command_line;


my $default_bowtie_breaksite_opt = "--very-sensitive-local --no-discordant -k 3 --score-min C,50 --mp 12,4 --rdg 10,6 --rfg 10,6 --gbar 1 -p $threads -t";
my $default_bowtie_opt = "--very-sensitive-local --score-min G,20,6 -p $threads -k 5 --no-unal --gbar 1 --dovetail -t";

my $bt2_breaksite_opt = manage_program_options($default_bowtie_breaksite_opt,$user_bowtie_breaksite_opt);
my $bt2_opt = manage_program_options($default_bowtie_opt,$user_bowtie_opt);

my $t0 = [gettimeofday];

#check_working_dir;

my $expt = basename($workdir);
my $breaksite_fa = "$workdir/misc/breaksite.fa";
my $bt2_breaksite_idx = "$workdir/misc/breaksite";
my $breaksite_R1_sam = "$workdir/${expt}_R1_breaksite.sam";
my $breaksite_R1_bam = "$workdir/${expt}_R1_breaksite.bam";
my $breaksite_R2_sam = "$workdir/${expt}_R2_breaksite.sam";
my $breaksite_R2_bam = "$workdir/${expt}_R2_breaksite.bam";
my $sam = "$workdir/$expt.sam";
my $bam = "$workdir/$expt.bam";
my $tlxl = "$workdir/$expt.tlxl";


align_to_breaksite;

align_to_genome; # unless (-r $bam);

if ($read1 =~ s/\.gz//) {
  System("gunzip -c $read1.gz > $read1");
}
if ($read2 =~ s/\.gz//) {
  System("gunzip -c $read2.gz > $read2");
}

make_tlxl;


my $t1 = tv_interval($t0);

printf("\nFinished all processes in %.2f seconds.\n", $t1);


#
# End of program
#

sub align_to_breaksite {
  print "\nRunning Bowtie2 alignment for $expt against breaksite locus\n";

  System("bowtie2-build -q $breaksite_fa $bt2_breaksite_idx");

  my $breaksite_R1_bt2_cmd = "bowtie2 $bt2_breaksite_opt --norc -x $bt2_breaksite_idx -U $read1 -S $breaksite_R1_sam";

  System($breaksite_R1_bt2_cmd);

  System("samtools view -bS -o $breaksite_R1_bam $breaksite_R1_sam");

  my $breaksite_R2_bt2_cmd = "bowtie2 $bt2_breaksite_opt --nofw -x $bt2_breaksite_idx -U $read2 -S $breaksite_R2_sam";

  System($breaksite_R2_bt2_cmd);

  System("samtools view -bS -o $breaksite_R2_bam $breaksite_R2_sam");
}

sub align_to_genome {

  print "\nRunning Bowtie2 alignment for $expt against genome $assembly\n";

  my $bt2_cmd = "bowtie2 $bt2_opt -x $assembly -1 $read1 -2 $read2 -S $sam";

  System($bt2_cmd);

  System("samtools view -bS -o $bam $sam");


}

sub make_tlxl {

  # my %R1_idx = %{make_raw_fastq_idx($read1)};
  # my %R2_idx = %{make_raw_fastq_idx($read2)};

  my $samobj_R1_brk = Bio::DB::Sam->new(-bam => $breaksite_R1_bam,
                                     -fasta => $breaksite_fa,
                                     -expand_flags => 1);

  my $samobj_R2_brk = Bio::DB::Sam->new(-bam => $breaksite_R2_bam,
                                     -fasta => $breaksite_fa,
                                     -expand_flags => 1);

  my $samobj = Bio::DB::Sam->new(-bam => $bam,
                                 -fasta => $ENV{'BOWTIE2_INDEXES'}."$assembly.fa",
                                 -expand_flags => 1);

  my $R1_brk_iter = $samobj_R1_brk->get_seq_stream;
  my $R2_brk_iter = $samobj_R2_brk->get_seq_stream;



  my $samiterator = $samobj->features(-iterator => 1);

  my $tlxlfh = IO::File->new(">$tlxl");
  $tlxlfh->print(join("\t",qw(Qname Rname Rstart Rend Strand Pair MapQ Qstart Qend BrkQstarts BrkQends BrkRstarts BrkRends R1 R2))."\n");

  my $R1_brk_aln = $R1_brk_iter->next_seq;
  my $R2_brk_aln = $R2_brk_iter->next_seq;


  my @R1_aln;
  my @R2_aln;

  do {

    my $qname = $R1_brk_aln->query->name;
    my $next_R1_brk_aln;
    while (1) {
      $next_R1_brk_aln = $R1_brk_iter->next_seq;
      if (defined $next_R1_brk_aln && $next_R1_brk_aln->name->query eq $qname) {
        #compare qstarts
      } else {
        last;
      }
    }
    my $next_R2_brk_aln
    while (1)









    $R1_brk_aln = $next_R1_brk_aln;
    $R2_brk_aln = $next_R2_brk_aln;

  } while (defined $R1_brk_aln);

  while (my $aln = $samiterator->next_seq) {
    my $qname = $aln->query->name;
    my $aln2 = $samiterator->next_seq if $aln->proper_pair;


    my $bqstart = 0;
    my $bqend = 0;
    my $brstart = 0;
    my $brend = 0;
    my @brstarts = ();
    my @brends = ();
    my @bqstarts = ();
    my @bqends = ();
    



    # my $R1fh = IO::File->new("<$read1");
    # my $R2fh = IO::File->new("<$read2");
    # my $R1idxfh = IO::File->new("<$read1.idx");
    # my $R2idxfh = IO::File->new("<$read2.idx");

    # my $R1seq = line_with_index($R1fh,$R1idxfh,$R1_idx{$qname});
    # my $R2seq = line_with_index($R2fh,$R2idxfh,$R2_idx{$qname});
    # chomp($R1seq);
    # chomp($R2seq);


    my $rname = $aln->seq_id;
    my $strand = $aln->strand;
    my $mapq = $aln->qual;
    my $pair;
    my $rstart;
    my $rend;
    my $qstart;
    my $qend;
    my $R1seq = "";
    my $R2seq = "";
    

    if ($aln->proper_pair) {

      $pair = 0;
      $rstart = min($aln->start,$aln2->start);
      $rend = max($aln->end,$aln2->end);
      $qstart = $aln->reversed ? $aln->l_qseq - $aln->query->end + 1 : $aln->query->start;
      $qend = $aln->reversed ? $aln->l_qseq - $aln->query->start + 1 : $aln->query->end;
      $R1seq = $aln->reversed ? reverseComplement($aln->query->seq->seq) : $aln->query->seq->seq;
      $R2seq = $aln2->reversed ? $aln2->query->seq->seq : reverseComplement($aln2->query->seq->seq);

      my @brk = $samobj_R1_brk->get_feature_by_name($qname);

      my $i = 0;
      while ($i < scalar @brk) {
        #if ($brk[$i]->get_tag_values('FIRST_MATE')) {
          push(@brstarts,$brk[$i]->start);
          push(@brends,$brk[$i]->end);
          push(@bqstarts,$brk[$i]->query->start);
          push(@bqends,$brk[$i]->query->end);

          if ( $bqstart==0 || $bqstart > $brk[$i]->query->start) {
            $brstart  = $brk[$i]->start;
            $brend = $brk[$i]->end;
            $bqstart = $brk[$i]->query->start;
            $bqend = $brk[$i]->query->end;
          }
        #}
        $i++;
      }


    } elsif ($aln->get_tag_values('FIRST_MATE')) {
      
      $pair = 1;
      $rstart = $aln->start;
      $rend = $aln->end;
      $qstart = $aln->reversed ? $aln->l_qseq - $aln->query->end + 1 : $aln->query->start;
      $qend = $aln->reversed ? $aln->l_qseq - $aln->query->start + 1 : $aln->query->end;
      $R1seq = $aln->reversed ? reverseComplement($aln->query->seq->seq) : $aln->query->seq->seq;

      my @brk = $samobj_R1_brk->get_feature_by_name($qname);


      my $i = 0;
      while ($i < scalar @brk) {
        #if ($brk[$i]->get_tag_values('FIRST_MATE')) {
          push(@brstarts,$brk[$i]->start);
          push(@brends,$brk[$i]->end);
          push(@bqstarts,$brk[$i]->query->start);
          push(@bqends,$brk[$i]->query->end);
          if ( $bqstart==0 || $bqstart > $brk[$i]->query->start) {
            $brstart  = $brk[$i]->start;
            $brend = $brk[$i]->end;
            $bqstart = $brk[$i]->query->start;
            $bqend = $brk[$i]->query->end;
          }
        #}
        $i++;
      }

    } else {

      $pair = 2;
      $rstart = $aln->start;
      $rend = $aln->end;
      $qstart = $aln->reversed ? $aln->query->start : $aln->l_qseq - $aln->query->end + 1;
      $qend = $aln->reversed ? $aln->query->end : $aln->l_qseq - $aln->query->start + 1;
      $R2seq = $aln->reversed ? $aln->query->seq->seq : reverseComplement($aln->query->seq->seq);

      my @brk = $samobj_R2_brk->get_feature_by_name($qname);

      my $i = 0;
      while ($i < scalar @brk) {
        #if ($brk[$i]->get_tag_values('SECOND_MATE')) {
          push(@brstarts,$brk[$i]->start);
          push(@brends,$brk[$i]->end);
          push(@bqstarts,$brk[$i]->query->start);
          push(@bqends,$brk[$i]->query->end);
          if ( $bqstart==0 || $bqstart > $brk[$i]->query->start) {
            $brstart  = $brk[$i]->start;
            $brend = $brk[$i]->end;
            $bqstart = $brk[$i]->query->start;
            $bqend = $brk[$i]->query->end;
          }
        #}
        $i++;
      }

    }


    $tlxlfh->print(join("\t", $qname,
                              $rname,
                              $rstart,
                              $rend,
                              $strand,
                              $pair,
                              $mapq,
                              $qstart,
                              $qend,
                              $bqstart,
                              $bqend,
                              $brstart,
                              $brend,
                              # join(",",@bqstarts),
                              # join(",",@bqends),
                              # join(",",@brstarts),
                              # join(",",@brends),
                              $R1seq,
                              $R2seq)."\n");
    
  }


}

sub parse_command_line {
	my $help;

	usage() if (scalar @ARGV == 0);

	my $result = GetOptions ( "read1=s" => \$read1,
                            "read2=s" => \$read2,
                            "assembly=s" => \$assembly,
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
