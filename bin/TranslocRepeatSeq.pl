#!/usr/bin/env perl


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


require "TranslocHelper.pl";
require "PerlSub.pl";

# Flush output after every write
select( (select(STDOUT), $| = 1 )[0] );

##
## This program r
## run with "--help" for usage information
##
## Robin Meyers

# Forward declarations
sub parse_command_line;

# Global flags and arguments, 
# Set by command line arguments
my $tlxfile;
my $bedfile;
my $output;


# Global variabless
my %filtered_reads;
my $repeatseq_output;
#
# Start of Program
#

parse_command_line;

# Run Rscript to print list of junctions to filter
($repeatseq_output = $output) =~ s/\.tlx$/.txt/;

my $repeatseq_cmd = join(" ","$FindBin::Bin/../R/TranslocRepeatSeq.R",
                        $tlxfile,
                        $bedfile,
                        $repeatseq_output);

System($repeatseq_cmd);

# Read in filtered junctions
open JUNC, "<", $repeatseq_output;
while (<JUNC>) {
  chomp;
  my @read = split("\t");
  $filtered_reads{$read[0]}->{$read[1]} = 1;
}
close JUNC;


# Read in tlxfile line by line, writing to output

my $infh = IO::File->new("<$tlxfile");
my $outfh = IO::File->new(">$output");
my $csv = Text::CSV->new({sep_char => "\t"});
my $header = $csv->getline($infh);
$csv->column_names(@$header);

$outfh->print(join("\t",@$header)."\n");

my $junc_id;
my $qname;

while (my $tlx = $csv->getline_hr($infh)) {

  if (! defined $qname || $tlx->{Qname} ne $qname) {
    $junc_id = 1;
    $qname = $tlx->{Qname}
  } else {
    $junc_id++;
  }

  if (defined $filtered_reads{$tlx->{Qname}}->{$junc_id}) {
    $tlx->{repeatseq} = 1;
    print "found one\n";
  }
  $outfh->print(join("\t",@{$tlx}{@$header})."\n");
}


#
# End of program
#


sub parse_command_line {

  my $help;

  usage() if (scalar @ARGV == 0);

  my $result = GetOptions ( 
                            "help" => \$help

                          );
  
  usage() if ($help);



  #Check options

  usage() if (scalar @ARGV < 3);

  $tlxfile = shift(@ARGV);
  $bedfile = shift(@ARGV);
  $output = shift(@ARGV);

  croak "Error: cannot read tlxfile" unless -r $tlxfile;
  croak "Error: cannot read bedfile" unless -r $bedfile;
  croak "Error: tlxfile must have .tlx extension" unless $tlxfile =~ /\.tlx$/;


  exit unless $result;
}


sub usage()
{
print<<EOF;
TranslocRepeatSeq, by Robin Meyers, 2013

Usage: $0 tlxfile bedfile output
        [--option VAL] [--flag] [--help]


Arguments (defaults in parentheses):

$arg{"tlxfile","File containing meta data for one experiment per row - follow correct format"}
$arg{"bedfile","File containing meta data for one experiment per row - follow correct format"}
$arg{"output","File containing meta data for one experiment per row - follow correct format"}

$arg{"--help","This helpful help screen."}

--------------------------------------------

EOF

exit 1;
}
