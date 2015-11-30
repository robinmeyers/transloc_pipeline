#!/usr/bin/env perl


use strict;
use warnings;
use Getopt::Long;
use Carp;
use Bio::SeqIO;
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
## This program
## 
## run with "--help" for usage information
##
## Robin Meyers

# Forward declarations
sub parse_command_line;
sub read_in_unmatched;
sub sort_and_print_result;


# Global flags and arguments, 
# Set by command line arguments
my $fastq;
my $barcodes_to_print = 25;
my $barcode_length = 10;
my $reads = 10000;

# Global variables 
my $barcodes = {};

#
# Start of Program
#

parse_command_line;

my $t0 = [gettimeofday];

read_in_unmatched;

sort_and_print_result;

my $t1 = tv_interval($t0);

printf("\nFinished all processes in %.2f seconds.\n", $t1);


#
# End of program
#

sub read_in_unmatched {

  my $gunzippath = capture("which gunzip");
  chomp $gunzippath;
  my $file = $fastq =~ /\.gz$/ ? "gunzippath -c $fastq |" : $fastq;

  my $fh = Bio::SeqIO->new(-file => $file,
                           -format => 'fastq');

  my $i = 0;

  while (my $seq = $fh->next_seq) {

    $i++;
    last if ($reads > 0 && $i > $reads);


    my $barcode = $seq->subseq(1,$barcode_length);

    if (exists $barcodes->{$barcode}) {
      $barcodes->{$barcode}++;
    } else {
      $barcodes->{$barcode} = 1;
    }

  }


  $fh->close;
}

sub sort_and_print_result {

  my @sorted_barcodes = sort { $barcodes->{$b} <=> $barcodes->{$a} } keys %$barcodes;


  my $i = 0;
  foreach my $barcode (@sorted_barcodes) {
    $i++;
    last if $i > $barcodes_to_print;
    print $barcode . "\t" . $barcodes->{$barcode} . "\n";
  }
}

sub parse_command_line {
  my $help;

  usage() if (scalar @ARGV == 0);

  my $result = GetOptions ( "print=i" => \$barcodes_to_print,
                            "length=i" => \$barcode_length,
                            "reads=i" => \$reads,
                            "help" => \$help
                          ) ;
  
  usage() if ($help);

  croak "Error: not enough input arguments" if (scalar @ARGV < 1);

  $fastq = shift(@ARGV);
  
  #Check options

  croak "Error: cannot read fastq file" unless (-r $fastq);
  

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

$arg{"fastq","Input fastq file"}
$arg{"--print","Print top n frequent barcodes",$barcodes_to_print}
$arg{"--length","Length of barcodes in bp to search for",$barcode_length}
$arg{"--reads","Number of reads to check - default 0 reads whole file",$reads}
$arg{"--help","This helpful help screen."}


EOF

exit 1;
}