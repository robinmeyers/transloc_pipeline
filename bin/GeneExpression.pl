#!/usr/bin/perl


use strict;
use warnings;
use Getopt::Long;
use Carp;
use Interpolation 'arg:@->$' => \&argument;
use Time::HiRes qw(gettimeofday tv_interval);
use threads;
use threads::shared;

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
sub calculate_gene_activity;


# Global flags and arguments, 
# Set by command line arguments
my $indir;
my $outdir;
my $nonmapdir;
my $assembly = "mm9";
my $tagsize = 35;
my $max_threads = 2;
my $nodes = 8;
my $lambda = 0.1;
my $alpha = 2;
my $skipprofile;
my $upstream = 5;
my $downstream = 10;
my $binsize = 50;

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
my @summaryfiles = ();
my @profiles = ();

calculate_gene_activity;





#
# End of program
#

sub calculate_gene_activity {

	my $t0 = [gettimeofday];
	my @threads = ();

	foreach my $tagfile (sort @tagfiles) {
		next unless ($tagfile =~ /\.map$/);
		
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
		
		my $summary = "$outdir/chr${chrname}-genes.txt";
		my $profile;
		push(@summaryfiles,$summary);

		unless (defined $skipprofile) {
			$profile = "$outdir/chr${chrname}-profile.txt";
			push(@profiles,$profile);
		} else {
			$profile = "NA";
		}

		my $nonmapfile = "$nonmapdir/chr$chrname-nonmappable.bed";
		

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
							print "\nChromsome $chrname\nLength: $chrlen\n";
							my $Rcmd = "$FindBin::Bin/../R/GeneExpression.R $smalltagfile $summary $refgene $nonmapfile"
								. " tagsize=$tagsize lambda=$lambda alpha=$alpha nodes=$nodes";
							$Rcmd = $Rcmd . " profile_file=$profile upstream=$upstream downstream=$downstream binsize=$binsize"
								unless (defined $skipprofile);
							System($Rcmd) or croak "Error: failed on chr$chrname during GeneExpression.R call";
							my $elapsed = tv_interval($t0_chr);
							printf("\nFinished chr%s in %.2f seconds.\n", $chrname, $elapsed);
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

	my $total = "$outdir/totalGeneExpression.txt";
	System("cat @summaryfiles | perl -p -e 'if (\$.>1){s/^Gene.*\\n//}' > $total");

	unless (defined $skipprofile) {
		my $total_profile = "$outdir/totalGeneProfile.txt";
		System("cat @profiles | perl -p -e 'if (\$.>1){s/^Gene.*\\n//}' > $total_profile");
	}

	my $elapsed = tv_interval($t0);

	printf("\nFinished calculating gene expression in %.2f seconds.\n", $elapsed);

}


sub parse_command_line {
	my $help;

	usage() if (scalar @ARGV==0);

	my $result = GetOptions ( "assembly=s" => \$assembly,
                            "threads=i" => \$max_threads,
                            "nodes=i" => \$nodes,
                            "lambda=f" => \$lambda,
                            "alpha=f" => \$alpha,
                            "tagsize=i" => \$tagsize,
                            "skipprofile" => \$skipprofile,
                            "upstream=i" => \$upstream,
                            "downstream=i" => \$downstream,
                            "binsize=i" => \$binsize,
                            "help" => \$help
				            ) ;
	
	usage() if ($help);
	usage() if (scalar @ARGV < 3);

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
$arg{"--tagsize","Length of tags in bp",$tagsize}
$arg{"--lambda","Background tag density in tags/kb (on both strands)",$lambda}
$arg{"--alpha","Over-dispersion parameter from negative binomial background distribution model",$alpha}
$arg{"--threads","Number of perl threads to run on - number of chromosomes analyzed at once",$max_threads}
$arg{"--nodes","Computes activity on clusters with this many nodes in the R program - number of genes analyzed at once", $nodes}
$arg{"--skipprofile","Skip computation of high-res gene expression profile - Subsequent options are for this module"}
$arg{"--upstream","Kilobases upstream from TSS to start binning",$upstream}
$arg{"--downstream","Kilobases downstream from TSS to stop binning",$downstream}
$arg{"--binsize","Bin size in bp for high-res profile",$binsize}
$arg{"--help","This helpful help screen."}


EOF

exit 1;
}
