#!/usr/bin/perl


use strict;
use warnings;
use Getopt::Long;
use Carp;
use IO::File;
use Text::CSV;
use threads;
use Interpolation 'arg:@->$' => \&argument;
use Time::HiRes qw(gettimeofday tv_interval);

my $TOOLS = $ENV{'TOOLS'};
my $ANNOT = $ENV{'ANNOT'};

require "$TOOLS/Perlsub.pl";
#require "$TOOLS/PSLhelpers.pl";

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
sub read_in_meta_file;
sub check_existance_of_files;
sub process_experiment;


# Global flags and arguments, 
# Set by command line arguments
my $meta_file;
my $indir;
my $outdir;
my $max_threads = 4;

# Global variables
my %expt_hash;
my %stats; 
my $genome2bit;

#
# Start of Program
#

parse_command_line;

my $t0 = [gettimeofday];

read_in_meta_file;

check_existance_of_files;

my @threads = ();

foreach my $expt_id (sort keys %expt_hash) {

    while (1) {

    # joins any threads if possible
        foreach my $thr (@threads) {
            $thr->join() if $thr->is_joinable();
        }

        my @running = threads->list(threads::running);
        
        # if there are open threads, create a new one, push it onto list, and exit while loop
        if (scalar @running < $max_threads) {
            my $thr = threads->create( sub {
                        my $t0_expt = [gettimeofday];
                        print "\nStarting $expt_id\n";
                        process_experiment($expt_id);
                        my $t1 = tv_interval($t0_expt);
                        printf("\nFinished %s with %d reads in %.2f seconds.\n", $expt_id, $stats{$expt_id}->{final},$t1);
                    });
            push(@threads,$thr);
            sleep(1);
            last;
        }
        sleep(1);
    } 
}

# waits for all threads to finish
while( scalar threads->list(threads::all) > 0) {
    for my $thr (@threads) {
        $thr->join() if $thr->is_joinable;
    }
    sleep(1);
}

my $t1 = tv_interval($t0);

printf("\nFinished all processes in %.2f seconds.\n", $t1);


#
# End of program
#

sub process_experiment ($) {
#	create_sequence_files;
#
#	align_to_sequence_files;
#
#	align_to_genome;
#
# combine_to_tlx; 
#
# filter_reads;
}

sub read_in_meta_file {
	System("perl -pi -e 's/\\r/\\n/g' $meta_file");

	print "\nReading in meta file...\n";

	my $meta = IO::File->new("<$meta_file");
	my $csv = Text::CSV->new({sep_char => "\t"});
	my $header_ref = $csv->getline($meta);
	my @header = @$header_ref;
	$csv->column_names(@header);

	while (my $row_ref = $csv->getline_hr($meta)) {

		$expt_hash{$row_ref->{experiment}."_".$row_ref->{seqrun}} = $row_ref;
	}
	#print join("\t",@header)."\n";
	#foreach my $expt (sort keys %expt_hash) {
	#	my $chr = $expt_hash{$expt}->{Chr};
	#	my $brksite = $expt_hash{$expt}->{Brksite};
	#	my $strand = $expt_hash{$expt}->{Strand};
	#	my %hash = %{$expt_hash{$expt}};
	#	print join("\t", @hash{@header} )."\n";
	#}

}

sub check_existance_of_files {
	print "\nSearching for files...\n";
	foreach my $expt_id (sort keys %expt_hash) {
		my $file = $indir."/".$expt_id;
		my @exts = qw(.fa .fasta .fq .fastq);
		foreach my $ext (@exts) {
			if (-r $file.$ext) {
				if ($ext =~ /q/) {
					(my $next = $ext) =~ s/q/a/;
					print "Converting $file to fasta format\n";
					System("fastq_to_fasta -Q33 -n -i $file$ext -o $file$next") or croak "Error: could not execute fastq_to_fastq";
					$expt_hash{$expt_id}->{file} = $file.$next;
				} else {
					$expt_hash{$expt_id}->{file} = $file.$ext;
				}
				last;
			}
		}
		croak "Error: Could not locate reads file $file in $indir" unless (defined $expt_hash{$expt_id}->{file});
	}
	print "Done.\n";
}

sub create_sequence_files {

}

sub parse_command_line {
	my $help;

	usage() if (scalar @ARGV == 0);

	my $result = GetOptions ( 
														"threads" => \$max_threads ,
														"help" => \$help

				            			);
	
	usage() if ($help);

	croak "Error: not enough input arguments" if (scalar @ARGV < 3);

	$meta_file = shift(@ARGV);
	$indir = shift(@ARGV);
	$outdir = shift(@ARGV);

  #Check options

  croak "Error: cannot find $meta_file" unless (-r $meta_file);
  croak "Error: input directory $indir does not exist" unless (-d $indir);
  unless (-d $outdir) {
  	System("mkdir -p $outdir") or croak "Error: output directory $outdir does not exist and cannot be created";
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

$arg{"metafile","File containing meta data for one experiment per row - follow correct format"}
$arg{"indir","Directory containing all input sequence files"}
$arg{"outdir","Directory for results files"}
$arg{"--threads","Number of threads to run bowtie on","$max_threads"}
$arg{"--help","This helpful help screen."}


EOF

exit 1;
}
