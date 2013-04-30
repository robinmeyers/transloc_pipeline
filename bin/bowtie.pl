#!/usr/bin/perl


use strict;
use warnings;
use Getopt::Long;
use Carp;
use Interpolation 'arg:@->$' => \&argument;
use Time::HiRes qw(gettimeofday tv_interval);

use Cwd qw(abs_path);
use FindBin;
use lib abs_path("$FindBin::Bin/../lib");

my $ANNOT = $ENV{'ANNOT'};
my $BOWTIE_IDX = $ENV{'BOWTIE_IDX'};

require "Perlsub.pl";


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


# Global flags and arguments, 
# Set by command line arguments
my $fastq;
my $outdir;
my $assembly;
my $trim3 = 16;
my $mismatch = 1;
my $report = 1;
my $threads = 8;
my $skipmap;
my $skipsam;
my $skipbgbw;

# Global variables 


#
# Start of Program
#

parse_command_line;

my $t0 = [gettimeofday];

print "\nRunning Bowtie alignment for $fastq against $assembly assembly\n";

my $index = "$BOWTIE_IDX/$assembly";
my $chrsize = "$ANNOT/chrSize.$assembly.txt";

my ($path,$name,$ext) = parseFilename($fastq);

my $sam = "$outdir/${name}.sam";
my $bam = "$outdir/${name}.bam";
my $sort = "$outdir/${name}_sort";
my $bg_p = "$outdir/${name}_p.bg";
my $bg_m = "$outdir/${name}_m.bg";
my $bw_p = "$outdir/${name}_p.bw";
my $bw_m = "$outdir/${name}_m.bw";

unless (defined $skipmap) {
	my $cmd = "bowtie $index -q $fastq --refout -v $mismatch -m $report -3 $trim3 --best --strata -t -a -p $threads";
	System($cmd) or croak "Error: cannot complete bowtie command to default output format";
}


unless (defined $skipsam) {
	my $cmd = "bowtie $index -q $fastq -S $sam -v $mismatch -m $report -3 $trim3 --best --strata -t -a -p $threads";
	System($cmd) or croak "Error: cannot complete bowtie command to SAM format output";

	$cmd = "samtools view -S -b -o $bam $sam";
	System($cmd) or croak "Error: cannot complete conversion to BAM format";

	$cmd = "samtools sort $bam $sort";
	System($cmd) or croak "Error: cannot complete sorting BAM file";

	unless (defined $skipbgbw) {
		$cmd = "time genomeCoverageBed -ibam $sort.bam -g $chrsize -strand + -bg > $bg_p";
		System($cmd) or croak "Error: cannot complete conversion to bedgraph format for plus strand";

		$cmd = "time genomeCoverageBed -ibam $sort.bam -g $chrsize -strand - -bg > $bg_m";
		System($cmd) or croak "Error: cannot complete conversion to bedgraph format for minus strand";

		$cmd = "time bedGraphToBigWig $bg_p $chrsize $bw_p";
		System($cmd) or croak "Error: cannot complete conversion to bigwig format for plus strand";

		$cmd = "time bedGraphToBigWig $bg_m $chrsize $bw_m";
		System($cmd) or croak "Error: cannot complete conversion to bigwig format for minus strand";
	}
}

my $t1 = tv_interval($t0);

printf("\nFinished all processes in %.2f seconds.\n", $t1);


#
# End of program
#



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

	$fastq = shift(@ARGV);
	$outdir = shift(@ARGV);
	$assembly = shift(@ARGV);

  #Check options

  croak "Error: cannot read input file" unless (-r $fastq);
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
