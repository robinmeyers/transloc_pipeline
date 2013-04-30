#!/usr/bin/perl


use strict;
use warnings;
use Getopt::Long;
use Carp;
use IO::File;
use Text::CSV;
use Interpolation 'arg:@->$' => \&argument;
use Time::HiRes qw(gettimeofday tv_interval);

use Cwd qw(abs_path);
use FindBin;
use lib abs_path("$FindBin::Bin/../lib");

my $ANNOT = $ENV{'ANNOT'};
my $BOWTIE_IDX = $ENV{'BOWTIE_IDX'};
my $GENOME_DB = $ENV{'GENOME_DB'};

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
my $tlx;
my $outdir;
my $assembly;
my $bgbw;

# Global variables 


#
# Start of Program
#

parse_command_line;

my $t0 = [gettimeofday];

my ($path,$name,$ext) = parseFilename($tlx);

$path = defined $outdir ? $outdir : $path;

my $bed = "$path/$name.bed";

my $tlxio = IO::File->new("<$tlx");
my $bedio = IO::File->new(">$bed");
my $csv = Text::CSV->new({sep_char => "\t"});
my $header = $csv->getline($tlxio);
$csv->column_names(@$header);

$bedio->print("track name=$name");

while (my $tl = $csv->getline_hr($tlxio)) {

	$bedio->print(join("\t", 	$tl->{Tname} ,
														$tl->{Tstart} ,
														$tl->{Tend} ,
														$tl->{Qname} ,
														0 ,
														$tl->{Strand} )."\n");

}

$bedio->close;
$tlxio->close;

if ($bgbw) {

	System("sort -k1,1 -k2,2n $bed > $bed.tmp");
	System("mv $bed.tmp $bed");

	my $bg_p = "$path/${name}_p.bg";
	my $bg_m = "$path/${name}_m.bg";
	my $bg = "$path/${name}.bg";
	my $bw_p = "$path/${name}_p.bw";
	my $bw_m = "$path/${name}_m.bw";
	my $bw = "$path/${name}.bw";

	my $chrsize = "$GENOME_DB/$assembly/annotation/ChromInfo.txt";


	my $cmd = "time genomeCoverageBed -i $bed -g $chrsize -strand + -bg > $bg_p";
	System($cmd) or croak "Error: cannot complete conversion to bedgraph format for plus strand";

	$cmd = "time genomeCoverageBed -i $bed -g $chrsize -strand - -bg > $bg_m";
	System($cmd) or croak "Error: cannot complete conversion to bedgraph format for minus strand";

	$cmd = "time genomeCoverageBed -i $bed -g $chrsize -bg > $bg";
	System($cmd) or croak "Error: cannot complete conversion to bedgraph format for combined strands";

	$cmd = "time bedGraphToBigWig $bg_p $chrsize $bw_p";
	System($cmd) or croak "Error: cannot complete conversion to bigwig format for plus strand";

	$cmd = "time bedGraphToBigWig $bg_m $chrsize $bw_m";
	System($cmd) or croak "Error: cannot complete conversion to bigwig format for minus strand";

	$cmd = "time bedGraphToBigWig $bg $chrsize $bw";
	System($cmd) or croak "Error: cannot complete conversion to bigwig format for combined strands";



}


my $t1 = tv_interval($t0);

printf("\nFinished all processes in %.2f seconds.\n", $t1);


#
# End of program
#



sub parse_command_line {
	my $help;

	usage() if (scalar @ARGV == 0);

	my $result = GetOptions ( "outdir=s" => \$outdir,
														"assembly=s" => \$assembly,
                            "bgbw" => \$bgbw,
				            				"help" => \$help
				            			) ;
	
	usage() if ($help);

	croak "Error: not enough input arguments" if (scalar @ARGV < 1);

	$tlx = shift(@ARGV);
	
  #Check options

  croak "Error: cannot read input file" unless (-r $tlx);
	if (defined $outdir) {
	  unless (-d $outdir) {
	  	System("mkdir $outdir") or croak "Error: output directory $outdir does not exist and cannot be created";
	  }
	}
	croak "Error: must define assembly if using --bgbw option" if defined $bgbw && !defined $assembly;


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

$arg{"tlx","Input tlx file"}
$arg{"--outdir","Directory for results files"}
$arg{"--assembly","Genome assembly to align reads to"}
$arg{"--bgbw","Convert bed format to bedgraph and bigwig formats"}
$arg{"--help","This helpful help screen."}


EOF

exit 1;
}
