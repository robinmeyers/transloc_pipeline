#!/usr/bin/perl


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

require "TranslocHelper.pl";
require "PerlSub.pl";

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
sub read_in_meta_file;
sub check_existance_of_files;
sub process_experiment ($);


# Global flags and arguments, 
# Set by command line arguments
my $meta_file;
my $seqdir;
my $outdir;
my $pipeline_threads = 2;
my $expt_threads = 4;
my $use_current_tlx;

my $bsub;
my $user_bsub_opt = "";
my $default_bsub_opt = "-n 4 -q short -W 12:00";
my $user_bowtie_opt = "";
my $user_bowtie_breaksite_opt = "";

# Global variabless
my %meta;
my %stats;


#
# Start of Program
#

parse_command_line;

my $t0 = [gettimeofday];

read_in_meta_file;

check_existance_of_files;

prepare_reference_genomes (\%meta);

if (defined $bsub) {
  foreach my $expt_id (sort keys %meta) {
    
    process_experiment($expt_id);
    
  }
} else {

  my @threads = ();

  foreach my $expt_id (sort keys %meta) {

      while (1) {

      # joins any threads if possible
          foreach my $thr (@threads) {
              $thr->join() if $thr->is_joinable();
          }

          my @running = threads->list(threads::running);
          
          # if there are open threads, create a new one, push it onto list, and exit while loop
          if (scalar @running < $pipeline_threads) {
              my $thr = threads->create( sub {
                          
                          process_experiment($expt_id);
                          
                      });
              push(@threads,$thr);
              sleep(1);
              last;
          }
          sleep(1);
      } 
  }

  # waits for all threads to finish
  while( scalar threads->list(threads::all) > 0) {
      for my $thr (@threads) {
          $thr->join() if $thr->is_joinable;
      }
      sleep(1);
  }
}

my $t1 = tv_interval($t0);

printf("\nFinished all processes in %.2f seconds.\n", $t1);


#
# End of program
#

sub process_experiment ($) {

	my $expt_id = shift;

  my $expt_hash = $meta{$expt_id};

  my $t0_expt = [gettimeofday];
  print "\nStarting $expt_id\n";
  unless (-d $expt_hash->{exptdir}) {
    mkdir $expt_hash->{exptdir} or croak "Error: cannot create experiment directory";
  }

	prepare_working_directory($expt_id);

  my $assembly = $expt_hash->{mask} =~ /\S/ ? $expt_hash->{mask_assembly} : $expt_hash->{assembly};


  my $tl_cmd = join(" ","TranslocPipeline.pl --workdir",$expt_hash->{exptdir},
                    "--assembly",$assembly,"--chr",$expt_hash->{chr},"--start",$expt_hash->{start},"--end",$expt_hash->{end},"--strand",$expt_hash->{strand},"--threads $expt_threads --bt2opt \"$user_bowtie_opt\" --bt2brkopt \"$user_bowtie_breaksite_opt\"",
                    "--read1",$expt_hash->{R1});
  $tl_cmd .= " --read2 " . $expt_hash->{R2} if defined $expt_hash->{R2};


  $tl_cmd .= " --usecurrtlx" if defined $use_current_tlx;

  my $bsubopt = manage_program_options($default_bsub_opt,$user_bsub_opt);

  my $log = $expt_hash->{exptdir} . "/$expt_id.log";

  $tl_cmd = join(" ","bsub",$bsubopt,"-J $expt_id -o $log",$tl_cmd) if defined $bsub;

  $tl_cmd .= " > $log 2>&1" unless defined $bsub;

  System($tl_cmd);

  my $t1 = tv_interval($t0_expt);
  printf("\nFinished %s in %.2f seconds.\n", $expt_id,$t1);

}

sub read_in_meta_file {
	System("perl -pi -e 's/\\r/\\n/g' $meta_file",1);

	print "\nReading in meta file...\n";

	my $metafh = IO::File->new("<$meta_file");
	my $csv = Text::CSV->new({sep_char => "\t"});
	my $header = $csv->getline($metafh);
	$csv->column_names( map { lc } @$header );

	while (my $expt = $csv->getline_hr($metafh)) {

		my $expt_id = $expt->{experiment};
		$meta{$expt_id} = $expt;
		$meta{$expt_id}->{exptdir} = "$outdir/$expt_id";

	}
	print join("\t",qw(Experiment Researcher Chr Start End Strand))."\n";
	foreach my $expt (sort keys %meta) {
    print join("\t",$expt,$meta{$expt}->{researcher},
                          $meta{$expt}->{chr},
                          $meta{$expt}->{start},
                          $meta{$expt}->{end},
                          $meta{$expt}->{strand})."\n";
  }
}

sub prepare_working_directory ($) {


  my $expt_id = shift;
  my $expt_hash = $meta{$expt_id};

  my $miscdir = $expt_hash->{exptdir} . "/misc";

  unless (-d $miscdir) {
    mkdir $miscdir or croak "Error: could not create misc directory for $expt_id";
  } 
  $expt_hash->{breakfa} = "$miscdir/breaksite.fa";
  $expt_hash->{primfa} = "$miscdir/primer.fa";
  $expt_hash->{adaptfa} = "$miscdir/adapter.fa";
  $expt_hash->{midfa} = "$miscdir/mid.fa";
  $expt_hash->{cutfa} = "$miscdir/cutter.fa";

  my $brkfh = IO::File->new(">".$expt_hash->{breakfa}) or croak "Error: could not write to breaksite fasta file";
  $brkfh->print(">Breaksite\n");
  $brkfh->print(uc($expt_hash->{breaksite})."\n");
  $brkfh->close;

  my $fprimfh = IO::File->new(">".$expt_hash->{primfa}) or croak "Error: could not write to forward primer fasta file";
  $fprimfh->print(">Primer\n");
  $fprimfh->print(uc($expt_hash->{primer})."\n");
  $fprimfh->close;

  my $rprimfh = IO::File->new(">".$expt_hash->{adaptfa}) or croak "Error: could not write to reverse primer fasta file";
  $rprimfh->print(">Adapter\n");
  $rprimfh->print(uc($expt_hash->{adapter})."\n");
  $rprimfh->close;

  my $midfh = IO::File->new(">".$expt_hash->{midfa}) or croak "Error: could not write to mid fasta file";
  $midfh->print(">MID\n");
  $midfh->print(uc($expt_hash->{mid})."\n");
  $midfh->close;

  my $cutfh = IO::File->new(">".$expt_hash->{cutfa}) or croak "Error: could not write to frequent cutter fasta file";
  $cutfh->print(">Cutter\n");
  $cutfh->print(uc($expt_hash->{cutter})."\n");
  $cutfh->close;


}

sub check_existance_of_files {
	print "\nSearching for files\n";
	foreach my $expt_id (sort keys %meta) {
		my $base = $seqdir."/".$expt_id;

    if (-r $base."_R1.fq.gz") {
      $meta{$expt_id}->{R1} = $base."_R1.fq.gz";
      $meta{$expt_id}->{R2} = $base."_R2.fq.gz" if -r $base."_R2.fq.gz";
    } elsif  (-r $base.".fq.gz") {
      $meta{$expt_id}->{R1} = $base.".fq.gz";
    } else {
      croak "Error: could not locate read 1 file for $expt_id in $seqdir";
    }

  }
}

sub parse_command_line {
	my $help;

	usage() if (scalar @ARGV == 0);

	my $result = GetOptions ( 
														"othreads=i" => \$pipeline_threads,
                            "ithreads=i" => \$expt_threads,
                            "usecurrtlx" => \$use_current_tlx,
														"bowtie2opt=s" => \$user_bowtie_opt,
														"help" => \$help

				            			);
	
	usage() if ($help);

	croak "Error: not enough input arguments" if (scalar @ARGV < 3);

	$meta_file = shift(@ARGV);
	$seqdir = shift(@ARGV);
	$outdir = shift(@ARGV);

  #Check options

  croak "Error: cannot find $meta_file" unless (-r $meta_file);
  croak "Error: input directory $seqdir does not exist" unless (-d $seqdir);
  unless (-d $outdir) {
  	mkdir $outdir or croak "Error: output directory $outdir does not exist and cannot be created";
  }



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

$arg{"metafile","File containing meta data for one experiment per row - follow correct format"}
$arg{"seqdir","Directory containing all input sequence files"}
$arg{"outdir","Directory for results files"}
$arg{"--threads","Number of threads to run bowtie on","$pipeline_threads"}
$arg{"--bowtie2opt"," "}
$arg{"--help","This helpful help screen."}


EOF

exit 1;
}
