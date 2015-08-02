#!/usr/bin/env perl

use strict;
use warnings;
use Getopt::Long;
use Carp;
use File::Spec;
use IO::File;
use Text::CSV;
use Bio::DB::Fasta;
use threads;
use Interpolation 'arg:@->$' => \&argument;
use Time::HiRes qw(gettimeofday tv_interval);
use IPC::System::Simple qw(system capture);

use Cwd qw(abs_path);
use FindBin;
use lib abs_path("$FindBin::Bin/../lib");

require "TranslocSub.pl";

my $GENOME_DB = $ENV{'GENOME_DB'};
defined $GENOME_DB or croak "Error: set environment variable GENOME_DB";


# Flush output after every write
select( (select(STDOUT), $| = 1 )[0] );

##
## This program reads in a metadata file and for each library
## starts a separate perl thread to execute the full
## TranslocPipeline script.
## 
## run with "--help" for usage information
##
## Robin Meyers

# Forward declarations
sub parse_command_line;
sub read_in_meta_file;
sub check_existance_of_files;
sub process_experiment ($);
sub usage ($);


# Global flags and arguments, 
# Set by command line arguments
my $meta_file;
my $seqdir;
my $outdir;
my $which;
my $pipeline_threads = 2;
my $pipeline_opt;
my $print_only;
my $simulate;
my $bsub;
my $bsub_opt = "";

our $debug_level = 0;

# Global variabless
my %meta;

#
# Start of Program
#

parse_command_line;

my $t0 = [gettimeofday];

# Parse through 'which' variable to determine which libraries to run
# @which array will contain index of each library to run
my @which = ();
if (defined $which) {
  my @tmpwhich = split(",",$which);
  foreach (@tmpwhich) {
    if (/(\d+)-(\d+)/) {
      # add range of indices
      push(@which,$1..$2);
    } elsif (/(\d+)/) {
      # add single index
      push(@which,$1);
    } else {
      croak "Error: invalid 'which' statement";
    }
  }
}

read_in_meta_file;

unless (defined $print_only) {

  check_existance_of_files;

  if (defined $bsub) {
    # No need to start threads in bsub option
    # Just submit each job with short sleep between

    foreach my $expt_id (sort keys %meta) {
      
      process_experiment($expt_id);
      sleep(1);
      
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

                debug_print("creating new thread",1);
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
}
my $t1 = tv_interval($t0);

printf("\nFinished all processes in %.2f seconds.\n", $t1);


#
# End of program
#

# Main subroutine to launch job
sub process_experiment ($) {

	my $expt_id = shift;

  my $expt_hash = $meta{$expt_id};

  my $t0_expt = [gettimeofday];
  print "\nStarting $expt_id\n";
  unless (-d $expt_hash->{exptdir}) {
    mkdir $expt_hash->{exptdir} or croak "Error: cannot create experiment directory";
  }

	prepare_working_directory($expt_id);

  my $tl_cmd = join(" ","TranslocPipeline.pl --workdir",$expt_hash->{exptdir},
                    "--read1",$expt_hash->{R1});

  $tl_cmd = join(" ", $tl_cmd, "--read2", $expt_hash->{R2}) if defined $expt_hash->{R2};
                
  $tl_cmd = join(" ", $tl_cmd,
                    "--assembly", $expt_hash->{assembly},
                    "--chr", $expt_hash->{chr},
                    "--start",$expt_hash->{start},
                    "--end",$expt_hash->{end},
                    "--strand",$expt_hash->{strand},
                    "--primer",$expt_hash->{primfa},
                    "--adapter",$expt_hash->{adaptfa});

  $tl_cmd = join(" ", $tl_cmd, "--mid", $expt_hash->{midfa}) if -r $expt_hash->{midfa};
  $tl_cmd = join(" ", $tl_cmd, "--breakseq", $expt_hash->{breakfa}) if -r $expt_hash->{breakfa};
  $tl_cmd = join(" ", $tl_cmd, "--breaksite", $expt_hash->{breaksite}) if $expt_hash->{breaksite};
  $tl_cmd = join(" ", $tl_cmd, "--cutter", $expt_hash->{cutfa}) if -r $expt_hash->{cutfa};

  $tl_cmd = join(" ", $tl_cmd, "--simfile", $expt_hash->{simfile}) if defined $expt_hash->{simfile};

  $tl_cmd = join(" ", $tl_cmd, $pipeline_opt) if defined $pipeline_opt;

  $tl_cmd = join(" ", $tl_cmd, "--debug", $debug_level) if $debug_level > 0;



  my $log = $expt_hash->{exptdir} . "/$expt_id.log";

  if (defined $bsub) {
    $tl_cmd = join(" ","bsub",$bsub_opt,"-J",$expt_hash->{library},"-o $log -N",$tl_cmd);
  } else {
    $tl_cmd .= " >> $log 2>&1";
  }

  System($tl_cmd);

  my $t1 = tv_interval($t0_expt);
  printf("\nFinished %s in %.2f seconds.\n", $expt_id,$t1);

}

sub read_in_meta_file {
	System("perl -pi -e 's/\\r/\\n/g' $meta_file",1);

	debug_print("reading in meta file",1);

  print join("\t",qw(. Library Sequencing Researcher Genome Chr Start End Strand))."\n";


	my $metafh = IO::File->new("<$meta_file");
	my $csv = Text::CSV->new({sep_char => "\t"});
	my $header = $csv->getline($metafh);
	$csv->column_names( map { lc } @$header );

  my $i = 0;
	while (my $expt = $csv->getline_hr($metafh)) {
    next unless $expt->{library} =~ /\S/;
    $i++;
    
    print join("\t",$i,$expt->{library},
                      $expt->{sequencing},
                      $expt->{researcher},
                      $expt->{assembly},
                      $expt->{chr},
                      $expt->{start},
                      $expt->{end},
                      $expt->{strand})."\n";


    unless ($print_only) {
      if (@which > 0) {

        # Only proceed with checking metadata and 
        # adding to metahash if library is in @which array
        next unless grep { $_ == $i } @which;
      }

      check_validity_of_metadata($expt);

  		my $expt_id = $expt->{library} . "_" . $expt->{sequencing};
      croak "Error: Experiment ID $expt_id is already taken" if exists $meta{$expt_id};
  		$meta{$expt_id} = $expt;
  		$meta{$expt_id}->{exptdir} = "$outdir/$expt_id";
    }

	}
}

sub check_validity_of_metadata ($) {
  my $expt = shift;

  debug_print("checking validity of metadata",1,$expt->{library});

  # Convert all sequence fields to uppercase only
  $expt->{breakseq} = uc($expt->{breakseq});
  $expt->{mid} = uc($expt->{mid});
  $expt->{primer} = uc($expt->{primer});
  $expt->{adapter} = uc($expt->{adapter});
  $expt->{cutter} = uc($expt->{cutter});


  my $assembly = $expt->{assembly};
  my $assembly_fa = "$GENOME_DB/$assembly/$assembly.fa";

  croak "Metadata error: could not find genome assembly $assembly_fa" unless -r $assembly_fa;

  my $assembly_obj = Bio::DB::Fasta->new($assembly_fa);
  my @chrlist = $assembly_obj->get_all_ids;

  croak "Metadata error: chr must be valid" unless grep { $_ eq $expt->{chr} } @chrlist; 
  croak "Metadata error: end must not be less than start" if $expt->{end} < $expt->{start};
  croak "Metadata error: strand must be one of + or -" unless grep { $_ eq $expt->{strand} } qw(+ -);

  if ($expt->{breakseq} ne "") {
    croak "Metadata error: breaksite sequence contains non AGCT characters" unless $expt->{breakseq} =~ /^[AGCTagct]+$/;
    croak "Metadata error: breaksite must be defined if breakseq is" unless $expt->{breaksite};
    croak "Metadata e for rror: breaksite cannot be greater than the length of breakseq" if $expt->{breaksite} > length($expt->{breakseq});
    croak "Metadata error: primer sequence not found in breakseq"
      unless index($expt->{breakseq},$expt->{primer}) > -1;
  } else {
    my $genome_seq = $expt->{strand} eq "+" ?
                      $assembly_obj->seq($expt->{chr},$expt->{start},$expt->{start}+length($expt->{primer})-1) :
                      reverseComplement($assembly_obj->seq($expt->{chr},$expt->{end}-length($expt->{primer}),$expt->{end}-1));
    croak "Metadata error: primer sequence does not match reference genome" 
      unless $expt->{primer} eq uc($genome_seq);
  }

  
  croak "Metadata error: primer sequence contains non AGCT characters" unless $expt->{primer} =~ /^[AGCT]+$/;
  croak "Metadata error: adapter sequence contains non AGCT characters" unless $expt->{adapter} =~ /^[AGCT]+$/;
  croak "Metadata error: cutter sequence contains non AGCT characters" unless $expt->{cutter} =~ /^[AGCT]*$/;
  croak "Metadata error: MID sequence contains non AGCT characters" unless $expt->{mid} =~ /^[AGCT]*$/;

}

sub prepare_working_directory ($) {

  my $expt_id = shift;
  my $expt_hash = $meta{$expt_id};

  debug_print("creating working directory",1,$expt_id);


  my $seqdir = File::Spec->catfile($expt_hash->{exptdir},
                                    $expt_id . "_sequences");

  unless (-d $seqdir) {
    mkdir $seqdir or croak "Error: could not create sequenes directory for $expt_id";
  } 
  $expt_hash->{breakfa} = "$seqdir/breaksite.fa";
  $expt_hash->{primfa} = "$seqdir/primer.fa";
  $expt_hash->{adaptfa} = "$seqdir/adapter.fa";
  $expt_hash->{midfa} = "$seqdir/mid.fa";
  $expt_hash->{cutfa} = "$seqdir/cutter.fa";

  if ($expt_hash->{breakseq} =~ /\S/) {
    my $brkfh = IO::File->new(">".$expt_hash->{breakfa}) or croak "Error: could not write to breaksite fasta file";
    $brkfh->print(">Breaksite\n");
    $brkfh->print($expt_hash->{breakseq}."\n");
    $brkfh->close;
  }

  my $primfh = IO::File->new(">".$expt_hash->{primfa}) or croak "Error: could not write to primer fasta file";
  $primfh->print(">Primer\n");
  $primfh->print($expt_hash->{primer}."\n");
  $primfh->close;

  my $adptfh = IO::File->new(">".$expt_hash->{adaptfa}) or croak "Error: could not write to adapter fasta file";
  $adptfh->print(">Adapter\n");
  $adptfh->print($expt_hash->{adapter}."\n");
  $adptfh->close;

  if ($expt_hash->{mid} =~ /\S/) {
    my $midfh = IO::File->new(">".$expt_hash->{midfa}) or croak "Error: could not write to mid fasta file";
    $midfh->print(">MID\n");
    $midfh->print($expt_hash->{mid}."\n");
    $midfh->close;
  }

  if ($expt_hash->{cutter} =~ /\S/) {
    my $cutfh = IO::File->new(">".$expt_hash->{cutfa}) or croak "Error: could not write to frequent cutter fasta file";
    $cutfh->print(">Cutter\n");
    $cutfh->print($expt_hash->{cutter}."\n");
    $cutfh->close;
  }


}

sub check_existance_of_files {

	debug_print("searching for input files",1);

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

    if (defined $simulate) {
      if (-r $base.".txt") {
        $meta{$expt_id}->{simfile} = $base.".txt";  
      } else {
        croak "Error could not locate simulated file $base.txt";
      }
    }

  }
}

sub parse_command_line {

  debug_print("parsing command line",1);

	my $help;

	usage(0) if (scalar @ARGV == 0);

	my $result = GetOptions ( "which=s" => \$which,
                            "bsub" => \$bsub,
                            "bsub-opt=s" => \$bsub_opt,
														"threads=i" => \$pipeline_threads,
                            "pipeline-opt=s" => \$pipeline_opt,
                            "print" => \$print_only,
                            "debug=i" => \$debug_level,
                            "simulate" => \$simulate,
														"help" => \$help

				            			);
	
	usage(1) if ($help);



  #Check options


  if ($print_only) {
    if (scalar @ARGV < 1) {
      carp "Error: print option requires metadata argument";
      usage(0);
    }
  } else {
     if (scalar @ARGV < 3) {
      carp "Error: not enough input arguments";
      usage(0);
    }
  }

  $meta_file = shift(@ARGV);

  croak "Error: cannot find $meta_file" unless (-r $meta_file);

  unless ($print_only) {

    $seqdir = shift(@ARGV);
    $outdir = shift(@ARGV);

    croak "Error: input directory $seqdir does not exist" unless (-d $seqdir);
    unless (-d $outdir) {
    	mkdir $outdir or croak "Error: output directory $outdir does not exist and cannot be created";
    }
  }



	exit unless $result;
}


sub usage ($)
{

my $full_help = shift;

print<<EOF;
TranslocWrapper, by Robin Meyers, 2013

Usage: $0 metafile seqdir outdir
        [--option VAL] [--flag] [--help]

Print: $0 metafile --print

Arguments (defaults in parentheses):

$arg{"metafile","File containing meta data for one experiment per row - follow correct format"}
$arg{"seqdir","Directory containing all input sequence files"}
$arg{"outdir","Directory for results output"}
$arg{"--which","Only run specific jobs, numbered by order in metafile"}
$arg{"--bsub","Submit as LSF jobs"}
$arg{"--bsub-opt","Specify bsub options different from default",$bsub_opt}
$arg{"--threads","Number of libraries to run at once",$pipeline_threads}
$arg{"--pipeline-opt","Specify pipeline options - see below"}
$arg{"--print","Do not execute jobs, only print libraries found in metafile"}
$arg{"--debug","Print debug info with verbosity at specified level, 1-4"}
$arg{"--simulate","Indicate that these are simulations"}
$arg{"--help","This helpful help screen."}

--------------------------------------------

EOF

system("$FindBin::Bin/TranslocPipeline.pl --help") if $full_help;

exit 1;
}
