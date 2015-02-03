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
my $output;
my $filters;


# Global variabless
my %filtered_reads;
my $filter_output;
#
# Start of Program
#

parse_command_line;

# Run Rscript to print list of junctions to filter
($filter_output = $output) =~ s/\.tlx$/.txt/;

my $filter_cmd = join(" ","$FindBin::Bin/../R/TranslocFilter.R",
                        $tlxfile,
                        $filter_output,
                        $filters);

System($filter_cmd);

# Read in filtered junctions
open JUNC, "<", $filter_output;
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


while (my $tlx = $csv->getline_hr($infh)) {
  my $qname = $tlx->{Qname};
  my $junc_id = $tlx->{JuncID};
  
  if (defined $filtered_reads{$qname}->{$junc_id}) {
    $outfh->print(join("\t",@{$tlx}{@$header})."\n");
  }
}


#
# End of program
#


sub parse_command_line {

  my $help;

  usage() if (scalar @ARGV == 0);

  my $result = GetOptions ( "filters=s" => \$filters,
                            "help" => \$help

                          );
  
  usage() if ($help);



  #Check options

  usage() if (scalar @ARGV < 2);

  $tlxfile = shift(@ARGV);
  $output = shift(@ARGV);

  croak "Error: cannot read tlxfile" unless -r $tlxfile;
  croak "Error: tlxfile must have .tlx extension" unless $tlxfile =~ /\.tlx$/;


  exit unless $result;
}


sub usage()
{
print<<EOF;
TranslocFilter, by Robin Meyers, 2013

Usage: $0 tlxfile output
        [--option VAL] [--flag] [--help]


Arguments (defaults in parentheses):

$arg{"tlxfile"," "}
$arg{"output"," "}
$arg{"--filters"," "}

$arg{"--help","This helpful help screen."}

--------------------------------------------

EOF

system("$FindBin::Bin/../R/TranslocFilter.R");

exit 1;
}
