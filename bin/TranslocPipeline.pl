#!/usr/bin/perl


use strict;
use warnings;
use Getopt::Long;
use Carp;
use IO::File;
use Text::CSV;
use Interpolation 'arg:@->$' => \&argument;
use Time::HiRes qw(gettimeofday tv_interval);

my $TOOLS = $ENV{'TOOLS'};
my $ANNOT = $ENV{'ANNOT'};

require "$TOOLS/Perlsub.pl";


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


# Global flags and arguments, 
# Set by command line arguments
my $meta_file;
my $indir;
my $outdir;
my $assembly;
my $trim3 = 16;
my $mismatch = 1;
my $report = 1;
my $max_threads = 8;
my $skipmap;
my $skipsam;
my $skipbgbw;

# Global variables 
my %expt_hash;


#
# Start of Program
#

parse_command_line;

my $t0 = [gettimeofday];

read_in_meta_file;

check_existance_of_files;

#my @threads = ();
#
#foreach my $expt (sort keys %expts) {
#
#    while (1) {
#
#    # joins any threads if possible
#        foreach my $thr (@threads) {
#            $thr->join() if $thr->is_joinable();
#        }
#
#        my @running = threads->list(threads::running);
#        
#        # if there are open threads, create a new one, push it onto list, and exit while loop
#        if (scalar @running < $max_threads) {
#            my $thr = threads->create( sub {
#                        my $t0_expt = [gettimeofday];
#                        print "\nStarting $expt\n";
#                        process_experiment($expt);
#                        my $t1 = tv_interval($t0_expt);
#                        printf("\nFinished %s with %d reads in %.2f seconds.\n", $expt, $expts{$expt}->{rmdups},$t1);
#                    });
#            push(@threads,$thr);
#            sleep(1);
#            last;
#        }
#        sleep(5);
#    } 
#}
#
## waits for all threads to finish
#while( scalar threads->list(threads::all) > 0) {
#    for my $thr (@threads) {
#        $thr->join() if $thr->is_joinable;
#    }
#    sleep(5);
#}

my $t1 = tv_interval($t0);

printf("\nFinished all processes in %.2f seconds.\n", $t1);


#
# End of program
#

sub process_experiment ($) {
	my $expt = shift;
		
}



sub create_sequence_files;

sub parse_command_line {
	my $help;

	usage() if (scalar @ARGV == 0);

	my $result = GetOptions ( "trim3=i" => \$trim3,
														"mismatch=i" => \$mismatch,
														"report=i" => \$report,
														"threads=i" => \$threads,
														"skipmap" => \$skipmap,
														"skipsam" => \$skipsam,
														"skipbgbw" => \$skipbgbw,
														"help" => \$help
														) ;
				
	usage() if ($help);

	croak "Error: not enough input arguments" if (scalar @ARGV < 3);

	$meta_file = shift(@ARGV);
	$indir = shift(@ARGV);
	$outdir = shift(@ARGV);
				
	#Check options

	croak "Error: cannot read meta file $meta_file" unless (-r $meta_file);
	unless (-d $outdir) {
		System("mkdir $outdir") or croak "Error: output directory $outdir does not exist and cannot be created";
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

$arg{"fastq","Input sequence file"}
$arg{"outdir","Directory for results files - note: default Bowtie output goes to working directory"}
$arg{"assembly","Genome assembly to align reads to"}
$arg{"--trim3","Number of basepairs to trim off of the 3' end of each read",$trim3}
$arg{"--mismatch","Number of mismatches allowed in the alignment",$mismatch}
$arg{"--report","Bowtie -m options: Number of alignments to report (reads with more equivalent alignments are suppressed)",$report}
$arg{"--threads","Number of threads to run bowtie on","$threads"}
$arg{"--skipmap","Skip reporting output in default Bowtie format"}
$arg{"--skipsam","Skip reporting output in SAM format and all downstream steps (BAM,bedgraph,bigwig)"}
$arg{"--skipbgbw","Skip converting BAM format to bedgraph and bigwig formats"}
$arg{"--help","This helpful help screen."}


EOF

exit 1;
}
