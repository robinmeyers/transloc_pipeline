#!/usr/bin/perl


use strict;
use warnings;
use Getopt::Long;
use Carp;
use Interpolation 'arg:@->$' => \&argument;
use Time::HiRes qw(gettimeofday tv_interval);
use threads;

use Cwd qw(abs_path);
use FindBin;
use lib abs_path("$FindBin::Bin/../lib");


require "Perlsub.pl";

my $ANNOT = $ENV{'ANNOT'};


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
sub bin_bowtie_tags;
sub process_chromosome ($);


# Global flags and arguments, 
# Set by command line arguments
my $indir;
my $outdir;
my $nonmapdir;
my $assembly = "mm9";
my $binsize = 500;
my $max_threads = 4;

# Global variables 
my $refgene;
my $chrsize;

#
# Start of Program
#

parse_command_line;

$refgene = "$ANNOT/refGene.$assembly.txt";
$chrsize = "$ANNOT/chrSize.$assembly.txt";

my %chrsize = readChromsizeFile($chrsize);

my @tagfiles = listFilesInDir($indir);
my @binfiles = ();

bin_bowtie_tags;



#
# End of program
#

sub bin_bowtie_tags {

	my $t0 = [gettimeofday];
	my @threads = ();

	foreach my $tagfile (sort @tagfiles) {
		next unless ($tagfile =~ /\.map$/);

		
		
		while (1) {

			# joins any threads if possible
			foreach my $thr (@threads) {
				$thr->join() if $thr->is_joinable();
			}

			my @running = threads->list(threads::running);
			
			# if there are open threads, create a new one, push it onto list, and exit while loop
			if (scalar @running < $max_threads) {
				my $thr = threads->create( sub {
							my $t0_chr = [gettimeofday];

							bin_bowtie_tags_from_file($tagfile);
							
							my $t1 = tv_interval($t0_chr);
							printf("\nFinished chr%s in %.2f seconds.\n", $chrname, $t1);
						});
				push(@threads,$thr);
				sleep(1);
				last;
			}
			sleep(5);
		}
	
	}

	# waits for all threads to finish
	while( scalar threads->list(threads::all) > 0) {
		for my $thr (@threads) {
			$thr->join() if $thr->is_joinable;
		}
		sleep(5);
	}

	my $total = "$outdir/wholeGenome.${binsize}kb-bins.txt";
	System("cat @binfiles | perl -p -e 'if (\$.>1){s/^Tname.*\\n//}' > $total");

	my $elapsed = tv_interval($t0);

	printf("\nFinished binning all reads in %.2f seconds.\n", $elapsed);

}

sub process_chromosome ($) {
	my $tagfile = shift;
	my $chrname = substr($tagfile,-6,2) + 1;
	if ($assembly =~ /mm/) {
		next if ($chrname == 22);
		$chrname =~ s/20/X/;
		$chrname =~ s/21/Y/;
	} elsif ($assembly =~ /hg/) {
		next if ($chrname == 25);
		$chrname =~ s/23/X/;
		$chrname =~ s/24/Y/;
	}
	my $chrlen = $chrsize{"chr$chrname"};
	my $smalltagfile = "$indir/chr$chrname.txt";
	unless (-r $smalltagfile) {
		print "\nConverting $tagfile to smaller $smalltagfile\n";
		open TAG, "<", $tagfile or croak "Error: cannot open $tagfile for reading";
		open STAG, ">", $smalltagfile or croak "Error: cannot open $smalltagfile for writing";
		while (<TAG>) {
			next unless /\S/;
			chomp;
			my @tag = split("\t");
			print STAG join("\t",$tag[2],$tag[1],$tag[3]);
			print STAG "\n";
		}
	}
	
	my $output = "$outdir/chr${chrname}-bins.txt";
	push(@binfiles,$output);
	
	my $nonmapfile = "$nonmapdir/chr$chrname-nonmappable.bed";


	print "\nChromsome $chrname\nLength: $chrlen\n";
	my $Rcmd = "Rscript $FindBin::Bin/../R/BinBowtieTagsToGenome.R $smalltagfile $output $nonmapfile $chrlen binsize=$binsize";
	System($Rcmd) or croak "Error: failed on chr$chrname during binBowtieTags.R call";
}


sub parse_command_line {
	my $help;

	usage() if (scalar @ARGV==0);

	my $result = GetOptions ( "assembly=s" => \$assembly,
                            "binsize=i" => \$binsize,
                            "threads=i" => \$max_threads,
				            "help" => \$help
				            ) ;
	
	usage() if ($help);

    $indir = shift(@ARGV);
    $outdir = shift(@ARGV);
    $nonmapdir = shift(@ARGV);

    #Check options

    croak "Error: input directory $indir does not exist" unless (-d $indir);
    unless (-d $outdir) {
        croak "Error: output directory $outdir does not exist and cannot be created" unless System("mkdir $outdir");
    }
    croak "Error: nonmappable directory $nonmapdir does not exist" unless (-d $nonmapdir);

    croak "Error: unrecognizable assembly $assembly" unless $assembly =~ /^(mm8|mm9|hg18|hg19)$/;


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

$arg{"indir","File path containing input file for each chromosome"}
$arg{"outdir","File path to put output files"}
$arg{"nonmapdir","File path containing bedfiles of all nonmappable basepairs in genome"}
$arg{"--assembly","Genome assembly of organism",$assembly}
$arg{"--binsize","Divide genome into this size bins (in kilobases)",$binsize}
$arg{"--threads","Number of threads to run program on",$max_threads}
$arg{"--help","This helpful help screen."}


EOF

exit 1;
}
