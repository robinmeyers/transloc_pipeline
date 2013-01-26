#!/usr/bin/perl


use strict;
use warnings;
use Getopt::Long;
use Carp;
use Interpolation 'arg:@->$' => \&argument;
use Time::HiRes qw(gettimeofday tv_interval);

my $TOOLS = $ENV{'TOOLS'};
my $ANNOT = $ENV{'ANNOT'};
my $BLAT_DB = $ENV{'BLAT_DB'};

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


# Global flags and arguments, 
# Set by command line arguments
my $assembly;
my $query;
my $output;
my $blatopt = "";
my $outType;
my $noHead;
my $mask;
my $tileSize;
my $stepSize;
my $oneOff;
my $minMatch;
my $minScore;
my $minIdentity;
my $extendThroughN; 

# Global variables 


#
# Start of Program
#

parse_command_line;

my $t0 = [gettimeofday];

System("blat $BLAT_DB/$assembly.2bit $query $output $blatopt") or croak "Error: blat alignment failed";

my $t1 = tv_interval($t0);

printf("\nFinished alignment to $assembly genome in %.2f seconds.\n", $t1);


#
# End of program
#



sub parse_command_line {
	my $help;

	usage() if (scalar @ARGV == 0);

	my $result = GetOptions ( "blatopt=s" => \$blatopt,
														"outType=s" => \$outType,
														"noHead" => \$noHead,
														"mask=s" => \$mask,
														"tileSize=i" => \$tileSize,
														"stepSize=i" => \$stepSize,
														"oneOff" => \$oneOff,
														"minMatch=i" => \$minMatch,
														"minScore=i" => \$minScore,
														"minIdentity=i" => \$minIdentity,
														"extendThroughN" => \$extendThroughN,
				            				"help" => \$help
				            			) ;
	
	usage() if ($help);

	if ($outType) {
		$blatopt .= " -out=$outType";
	}
	if ($mask) {
		$blatopt .= " -mask=$mask";
	}
	if ($tileSize) {
		$blatopt .= " -tileSize=$tileSize";
	}
	if ($stepSize) {
		$blatopt .= " -stepSize=$stepSize";
	}
	if ($oneOff) {
		$blatopt .= " -oneOff=$oneOff";
	}
	if ($minMatch) {
		$blatopt .= " -minMatch=$minMatch";
	}
	if ($minScore) {
		$blatopt .= " -minScore=$minScore";
	}
	if ($minIdentity) {
		$blatopt .= " -minIdentity=$minIdentity";
	}
	if ($extendThroughN) {
		$blatopt .= " -extendThroughN";
	}

	croak "Error: not enough input arguments" if (scalar @ARGV < 3);

	my $database = shift(@ARGV);
	my $query = shift(@ARGV);
	my $output = shift(@ARGV);

  #Check options

  croak if 0;



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

$arg{"database","Reference fasta file"}
$arg{"query","Query fasta file"}
$arg{"output","Output file"}
$arg{"--blatopt","Directly input blat options or use shortcut options below", $blatopt}
$arg{"--outType","Output file format","psl"}
$arg{"--noHead","Supress .psl header"}
$arg{"--mask","Mask out repeats - see blat options for full description","none"}
$arg{"--tileSize","Size of match that triggers alignment","11"}
$arg{"--stepSize","Spacing between tiles","tileSize"}
$arg{"--oneOff","Allow one mismatch in tile"}
$arg{"--minMatch","Number of tile matches","2"}
$arg{"--minScore","Minimum score - matches minus mismatches minus gap penalties","30"}
$arg{"--minIdentity","Minimum sequence identity (in percent)","90"}
$arg{"--extendThroughN","Allow extension of alignment through large blocks of N"}
$arg{"--help","This helpful help screen."}


EOF

exit 1;
}
