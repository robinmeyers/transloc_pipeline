#!/usr/bin/env perl


use strict;
use warnings;
use Getopt::Long;
use Carp;
use IO::Handle;
use IO::File;
use Text::CSV;
use Interpolation 'arg:@->$' => \&argument;
use POSIX qw(ceil floor);
use threads;
use threads::shared;
use Time::HiRes qw(gettimeofday tv_interval);


use Cwd qw(abs_path);
use FindBin;
use lib abs_path("$FindBin::Bin/../lib");

require "PerlSub.pl";
require "PipelineHelper.pl";



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
sub create_barcode_file;
sub check_existance_of_files;
sub update_stats_from_multx;
sub update_stats_from_trim;
sub process_experiment ($);
sub write_stats_file;
sub clean_up;

# Global flags and arguments, 
# Set by command line arguments
my $read1;
my $read2;
my $indir;
my $outdir;
my $meta_file;
my $max_threads = 4;
my $bc_len = 10;
my $bc_mismatch = 2;
my $adapt_max_dif = 20;
my $minlen = 30;
my $join;
my $join_max_dif = 10;
my $join_min_ol = 30;
my $skipclean;
my $forward_adapter = "AGATCGGAAGAGCGGTTCAG";
my $reverse_adapter = "AGATCGGAAGAGCGTCGTGT";
my $adapter_fa = "$FindBin::Bin/../ref/IlluminaAdapters-PE.fa";


# Global variables 
my %stats :shared;
my %meta;
my $unmatched = 0;
my $totalreads = 0;
my @bc_output;
my $bcfile;



#
# Start of Program
#

my $t0 = [gettimeofday];

parse_command_line;

read_in_meta_file;

unless (defined $indir) {
	create_barcode_file;
	mkdir "$outdir/multx";
	my $multx_cmd = "fastq-multx -m $bc_mismatch -x -b -B $bcfile $read1 $read2 -o $outdir/multx/%_R1.fq.gz $outdir/multx/%_R2.fq.gz";
	@bc_output = Capture($multx_cmd) or croak "Error: fastq-multx failed";
	print "\n@bc_output\n";
} else {
	check_existance_of_files;
}

update_stats_from_multx;


mkdir "$outdir/trim";
mkdir "$outdir/join";
mkdir "$outdir/logs";

my @threads = ();

foreach my $expt (sort keys %meta) {

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
						print "\nStarting $expt\n";
						process_experiment($expt);
						my $t1_expt = tv_interval($t0_expt);
						printf("\nFinished %s with %d reads in %.2f seconds.\n", $expt, $meta{$expt}->{multx},$t1_expt);
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

unless (defined $indir) {}


update_stats_from_trim;


write_stats_file;

clean_up unless $skipclean;

my $t1 = tv_interval($t0);
printf("\nFinished pre-preprocessing all experiments in %.2f seconds.\n", $t1);


#
# End of program
#


sub read_in_meta_file {
	print "\nReading in metadata file\n";

	System("perl -pi -e 's/\\r/\\n/g' $meta_file",1) or croak "Error: problem converting line feeds to unix format";

	my $metafh = IO::File->new("<$meta_file") or croak "Error: could not read meta file $meta_file";
	my $csv = Text::CSV->new({sep_char => "\t"});
	my $header = $csv->getline($metafh);
	$csv->column_names( map { lc } @$header );


	while (my $expt = $csv->getline_hr($metafh)) {

		next unless $expt->{library} =~ /\S/;

		my $name = $expt->{library} . "_" . $expt->{sequencing};

		$meta{$name} = $expt;
		
		my %temp_stats_hash :shared;
		$temp_stats_hash{totreads} = 0;
		
		$stats{$name} = \%temp_stats_hash;
	}

	$metafh->close;
}

sub create_barcode_file {
	$bcfile = "$outdir/barcodes.txt";
	print "\nWriting barcodes to $bcfile \n";
	my $bcfh = IO::File->new(">$bcfile");

	foreach my $expt (sort keys %meta) {
		$meta{$expt}->{R1} = "$outdir/multx/${expt}_R1.fq.gz";
		$meta{$expt}->{R2} = "$outdir/multx/${expt}_R2.fq.gz";
		my $barcode = substr($meta{$expt}->{mid}.$meta{$expt}->{primer},0,$bc_len);
		$bcfh->print(join("\t",$expt,$barcode)."\n");
		print "$expt $barcode\n";
	}

	$bcfh->close;
}

sub check_existance_of_files {
	print "\nSearching for files\n";
	foreach my $expt (sort keys %meta) {
		$meta{$expt}->{R1} = "$indir/${expt}_R1.fq.gz";
		$meta{$expt}->{R2} = "$indir/${expt}_R2.fq.gz";
		croak "Error: Could not locate $expt read 1 file" unless ($meta{$expt}->{R1});
		croak "Error: Could not locate $expt read 2 file " unless ($meta{$expt}->{R2});

	}
}

sub update_stats_from_multx {
	my $flag = 0;
	unless (defined $indir) {
		foreach (@bc_output) {
			if (/^Id/) {
				$flag = 1;
				next;
			}
			next unless $flag;
			my @row = split("\t");
			if ($row[0] =~ /unmatched/) {
				$unmatched = $row[1];
			} elsif ($row[0] =~ /total/) {
				$totalreads = $row[1];
			} else {
				$meta{$row[0]}->{multx} = $row[1];
			}
		}
	} else {
		foreach my $expt (sort keys %meta) {
			my @fqstats = Capture("fastq-stats ".$meta{$expt}->{R1},1);
			foreach (@fqstats) {
				chomp;
				if ( /reads\s+(\d+)/ ) {
					$meta{$expt}->{multx} = $1;
					last;
				}
			}
		}
	}
}

sub update_stats_from_trim {
	foreach my $expt (sort keys %meta) {
		my @fqstats = Capture("fastq-stats $outdir/trim/${expt}_R1.fq.gz",1);
		foreach (@fqstats) {
			chomp;
			if ( /reads\s+(\d+)/ ) {
				$meta{$expt}->{trim} = $1;
				last;
			}
		}
	}
}

sub process_experiment ($) {
	my $expt = shift;

	my $logfile = "$outdir/logs/$expt.log";

	System("echo \"Pre-processing $expt\" > $logfile",1);


	$meta{$expt}->{R1trim} = "$outdir/trim/${expt}_R1.fq.gz";
	$meta{$expt}->{R2trim} = "$outdir/trim/${expt}_R2.fq.gz";

	System("echo \"Running SeqPrep on $expt\" >> $logfile",1);

	my $t0_trim = [gettimeofday];
						

	my $trim_cmd = join(" ","SeqPrep -f",$meta{$expt}->{R1},"-r",$meta{$expt}->{R2},
													"-1",$meta{$expt}->{R1trim},"-2",$meta{$expt}->{R2trim},
													"-L $minlen -A $forward_adapter -B $reverse_adapter >> $logfile 2>&1");

	System("echo \"$trim_cmd\" >> $logfile",1);

	System($trim_cmd,1) or croak "Error: failed running SeqPrep on $expt";

	my $t1_trim = tv_interval($t0_trim);

	System("echo \"Finished fastq-trim on $expt in $t1_trim seconds\" >> $logfile",1);

	if ($join) {
		$meta{$expt}->{R1join} = "$outdir/join/${expt}_R1.fq.gz";
		$meta{$expt}->{R2join} = "$outdir/join/${expt}_R2.fq.gz";
		$meta{$expt}->{join} = "$outdir/join/${expt}_join.fq.gz";

		System("echo \"Running fastq-join on $expt\" >> $logfile",1);

		my $t0_join = [gettimeofday];

		my $join_cmd = join(" ","fastq-join -p $join_max_dif -m $join_min_ol",$meta{$expt}->{R1trim},$meta{$expt}->{R2trim},
			"-o",$meta{$expt}->{R1join},"-o",$meta{$expt}->{R2join},"-o",$meta{$expt}->{join},">> $logfile");

		System("echo \"$join_cmd\" >> $logfile",1);

		System($join_cmd,1) or croak "Error: failed running fastq-join on $expt";

		my $t1_join = tv_interval($t0_join);

		System("echo \"Finished fastq-join on $expt in $t1_join seconds\" >> $logfile",1);

	}

}


sub write_stats_file {
	print "\nWriting stats file\n";
	my $statsfh = IO::File->new(">$outdir/preprocess_stats.txt");
	my @header = qw(Expt Multx Trim);
	$statsfh->print(join("\t",@header)."\n");
	foreach my $expt (sort keys %meta) {
		my @row = ($expt,
			$meta{$expt}->{multx},
			$meta{$expt}->{trim});
		$statsfh->print(join("\t",@row)."\n");
	}
	unless (defined $indir) {
		$statsfh->print("Unmatched (PhiX)\t$unmatched\n");
		$statsfh->print("Total\t$totalreads\n");
	}
}

sub clean_up {

	print "\nCleaning up\n";

	System("rm $outdir/multx/*");

	unless ($join) {
		System("mv $outdir/trim/* $outdir/" );
	} else {
		System("mv $outdir/join/* $outdir/");
		System("rm $outdir/trim/*")
	}

}

sub parse_command_line {
  my $help;

  usage() if (scalar @ARGV==0);

  my $result = GetOptions (
  			"read1=s" => \$read1,
  			"read2=s" => \$read2,
  			"indir=s" => \$indir,
				"threads=i" => \$max_threads,
				"bclen=i" => \$bc_len,
				"bcmismatch=i" => \$bc_mismatch,
				"adapt_max_dif=i" => \$adapt_max_dif,
				"join" => \$join,
				"join_max_dif=i" => \$join_max_dif,
				"join_min_ol=i" => \$join_min_ol,
				"minlen=i" => \$minlen,
				"skipclean" => \$skipclean,
				"help" => \$help
		  ) ;
  
  usage() if ($help);
	
	$meta_file = shift(@ARGV);
	$outdir = shift(@ARGV);

	#Check options
	if (defined $indir) {
		croak "Error: cannot find input directory $indir" unless (-d $indir);
		croak "Error: cannot define both input directory and non-de-multiplexed reads" if (defined $read1 || defined $read2);		
	} else {
		croak "Error: cannot find read 1 $read1 does not exist" unless (-r $read1);
		croak "Error: cannot find read 1 $read2 does not exist" unless (-r $read2);
	}
	unless (-d $outdir) {
		mkdir $outdir or croak "Error: output directory $outdir neither exists nor could be created";
	}
	croak "Error: could not find or read meta file $meta_file" unless (-r $meta_file);

  exit unless $result;
}


sub usage()
{
print<<EOF;
Title, by Robin Meyers, ddmonthyyyy

This program .


Usage: $0 <meta_file> <outdir> (--read1 FIL --read2 FIL | --indir DIR) [--opts N] [--help]

Arguments (defaults in parentheses):


$arg{"meta_file","List of files with primer sequences"}
$arg{"outdir","Output directory"}
$arg{"--read1","Illumina read 1 file for de-multiplexing"}
$arg{"--read2","Illumina read 2 file for de-multiplexing"}
$arg{"--indir","Input directory - if specified assumes reads are already de-multiplexed"}
$arg{"--threads","Number of threads to run program on - experiments processed simultaneously",$max_threads}
$arg{"--bclen","Number of bases to use from primer file for de-multiplexing",$bc_len}
$arg{"--bcmismatch","Number of mismatches allowed in de-multiplexing",$bc_mismatch}
$arg{"--adapter","Fasta file of adapter sequences"}
$arg{"--adapt_max_dif","Maximum percent difference for match with adapter",$adapt_max_dif}
$arg{"--join","Stitch reads together using fastq-join"}
$arg{"--join_max_dif","Maximum percent difference between reads for quality-trimmed stitching",$join_max_dif}
$arg{"--join_min_ol","Minimum basepair overlap between two reads for quality-trimmed stitching",$join_min_ol}
$arg{"--minlength","Minimum length of quality trimmed reads",$minlen}
$arg{"--help","This helpful help screen."}


EOF

exit 1;
}
