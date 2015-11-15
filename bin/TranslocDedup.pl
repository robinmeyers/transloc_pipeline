#!/usr/bin/env perl


use strict;
use warnings;
use Getopt::Long;
use Carp;
use IO::File;
use Text::CSV;
use Interpolation 'arg:@->$' => \&argument;
use Time::HiRes qw(gettimeofday tv_interval);
use IPC::System::Simple qw(system capture);

use Cwd qw(abs_path);
use FindBin;
use lib abs_path("$FindBin::Bin/../lib");


require "TranslocSub.pl";


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
my $cores = 1;
my $offset_dist = 0;
my $break_dist = 0;

# Global variabless
my %deduped_reads;
my $dedup_output;

#
# Start of Program
#

parse_command_line;

# Run Rscript to print list of junctions to filter
($dedup_output = $output) =~ s/\.tlx$/.txt/;


my $dedup_cmd = join(" ","$FindBin::Bin/../R/TranslocDedup.R",
                        $tlxfile,
                        $dedup_output,
                        "--offset.dist", $offset_dist,
                        "--break.dist", $break_dist,
                        "--cores", $cores);

System($dedup_cmd);



# Read in filtered junctions
open JUNC, "<", $dedup_output;
while (<JUNC>) {
  chomp;
  my @read = split("\t");
  $deduped_reads{$read[0]} = 1;
}
close JUNC;

# Read in tlxfile line by line, writing to output
my $infh = IO::File->new("<$tlxfile");
my $outfh = IO::File->new(">$output");
my $csv = Text::CSV->new({sep_char => "\t"});


my $header = $csv->getline($infh);
$outfh->print(join("\t",@$header)."\n");
$csv->column_names(@$header);


while (my $tlx = $csv->getline_hr($infh)) {
  if (defined $deduped_reads{$tlx->{Qname}}) {
    $tlx->{duplicate} = 1;
  }
  $outfh->print(join("\t",@{$tlx}{@$header})."\n");
}


#
# End of program
#


sub parse_command_line {

  debug_print("parsing command line",1);

  my $help;

  usage() if (scalar @ARGV == 0);

  my $result = GetOptions ( 
                            "offset_dist=i" => \$offset_dist,
                            "break_dist=i" => \$break_dist,
                            "cores=i" => \$cores,
                            "help" => \$help

                          );
  
  usage() if ($help);



  #Check options

  usage() if (scalar @ARGV < 2);

  $tlxfile = shift(@ARGV);
  $output = shift(@ARGV);

  croak "Error: cannot read tlxfile" unless -r $tlxfile;
  
  exit unless $result;
}


sub usage()
{
print<<EOF;
TranslocRepeatSeq, by Robin Meyers, 2013

Usage: $0 tlxfile output
        [--option VAL] [--flag] [--help]


Arguments (defaults in parentheses):

$arg{"tlxfile","File containing meta data for one experiment per row - follow correct format"}
$arg{"output","File containing meta data for one experiment per row - follow correct format"}

$arg{"--help","This helpful help screen."}

--------------------------------------------

EOF

exit 1;
}
