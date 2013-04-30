#!/usr/bin/perl


use strict;
use warnings;
use Getopt::Long;
use Carp;
use Interpolation 'arg:@->$' => \&argument;
use Time::HiRes qw(gettimeofday tv_interval);
use List::Util qw(first);
use File::Slurp;
use Math::Round;


use Cwd qw(abs_path);
use FindBin;
use lib abs_path("$FindBin::Bin/../lib");

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
sub calculate_information_content;
sub search_database;

# Global flags and arguments, 
# Set by command line arguments
my $query;
my $database;
my $output;
my $score_thresh = 0;

# Global variables 
my $query_length;
my @sequence;
my $query_num = 0;
my $max_score = 0;
my $pen_thresh = 0;

#
# Start of Program
#

parse_command_line;

my $t0 = [gettimeofday];

print "\nCalculating information content from sequences in $query\n";

calculate_information_content;

$pen_thresh = $max_score - $score_thresh;

search_database if ($score_thresh > 0);



my $t1 = tv_interval($t0);

printf("\nFinished all processes in %.2f seconds.\n", $t1);


#
# End of program
#

sub calculate_information_content {
	my $qfh = IO::File->new("<$query") or croak "Error: could not open $query for reading";

	while ( my ($seqid,$seq) = read_fasta($qfh) ) {
		unless (defined $query_length) {
			$query_length = length($seq);
			# initialize @sequence hashes
			foreach my $i (0..($query_length-1)) {
				$sequence[$i] = {A => 0, C => 0, G => 0, T => 0, N => 0};
			}
		} else {
			croak "Error: all query sequences must be the same length" unless $query_length == length($seq);
		}
		foreach my $i (0..($query_length-1)) {
			my $base = uc(substr($seq,$i,1));
			$base = "N" unless first { $base eq $_ } ("A","C","G","T");
			$sequence[$i]->{$base} ++;
		}
		$query_num++;
	}

	print(join("\t","Pos","A","C","G","T","N")."\n");

	foreach my $i (0..($query_length-1)) {
		my $uncertainty = 0;
		my $max_freq = 0;
		foreach my $b (keys %{$sequence[$i]}) {
			my $freq = $sequence[$i]->{$b}/$query_num;
			$uncertainty += -$freq * log2($freq) if $freq > 0;
			$max_freq = $freq if $freq > $max_freq;
		}
		$sequence[$i]->{info} = 2 - $uncertainty;
		$sequence[$i]->{max} = $sequence[$i]->{info} * $max_freq;
		$max_score += $sequence[$i]->{max};
		print(join("\t",$i+1,$sequence[$i]->{A},$sequence[$i]->{C},$sequence[$i]->{G},$sequence[$i]->{T},$sequence[$i]->{N},$sequence[$i]->{info},$sequence[$i]->{max})."\n");
	}

	print "\nMaximum score is $max_score over $query_length bases\n";

	$qfh->close;

}


sub search_database {

	my $dbfh = IO::File->new("<$database") or croak "Error: could not open $database for reading";
	my $bedfh = IO::File->new(">$output") or croak "Error: could not open $output for writing";
	my $total_hits = 0;
	while ( my ($seqid,$seq) = read_fasta($dbfh) ) {
		my $chr_hits = 0;

		$seq = uc($seq);
		$seq =~ s/[^ACGTN]/N/g;
		my $db_len = length($seq)-$query_length+1;
		my $subseq = substr($seq,0,$query_length,"");
		print "Searching $seqid for matches...   0%";
		my $percentdone = 0;
		foreach my $i (1..$db_len) {

			if ($i/$db_len*100 >= $percentdone + 1) {
        $percentdone++;
        printf("\b\b\b\b%3d%%",$percentdone);
      }


			print "$subseq\n" if ($query_length != length($subseq));
			my $p_score = score_sequence($subseq);
			my $m_score = score_sequence(reverseComplement($subseq));
			$bedfh->print(join("\t",$seqid,$i,$i+$query_length-1,$p_score,"hit".(++$chr_hits+$total_hits),"+",$subseq)."\n") if ($p_score >= $score_thresh);
			$bedfh->print(join("\t",$seqid,$i,$i+$query_length-1,$m_score,"hit".(++$chr_hits+$total_hits),"-",reverseComplement($subseq))."\n") if ($m_score >= $score_thresh);
			$subseq = substr($subseq,1).substr($seq,0,1,"");
		}
		print("\b\b\b\b100%%\n");

		print "Found $chr_hits in $seqid\n";
		$total_hits += $chr_hits;
	}

	$dbfh->close;
	$bedfh->close;
}

sub score_sequence ($) {
	my $seq = shift;
	my @seq = split("",$seq);
	my $penalty = 0;
	foreach my $i (0..($query_length-1)) {
		my $base_score = $sequence[$i]->{$seq[$i]}/$query_num * $sequence[$i]->{info};
		$penalty += $sequence[$i]->{max} - $base_score;
		return(0) if ($penalty > $pen_thresh);
	}
	return($max_score-$penalty);
}


sub parse_command_line {
	my $help;

	usage() if (scalar @ARGV == 0);

	my $result = GetOptions ( "score=f" => \$score_thresh,
				            				"help" => \$help
				            			) ;
	
	usage() if ($help);

	croak "Error: not enough input arguments" if (scalar @ARGV < 3);

	$query = shift(@ARGV);
	$database = shift(@ARGV);
	$output = shift(@ARGV);

  #Check options

  croak "Error: cannot read input file" unless (-r $query);



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

$arg{"query","Fasta formatted file with sequences to search"}
$arg{"database","Either a genome assembly or fasta file with sequences to be searched"}
$arg{"output","File path to write bed file to"}
$arg{"--score","Minimum score of a sequence match to be written to bedfile - leave blank or set to 0 to only calculate information content",$score_thresh}
$arg{"--help","This helpful help screen."}


EOF

exit 1;
}
