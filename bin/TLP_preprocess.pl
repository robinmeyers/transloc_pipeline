#!/usr/bin/env perl


use strict;
use warnings;
use Getopt::Long;
use Carp;
use Interpolation 'arg:@->$' => \&argument;
use POSIX qw(ceil floor);
use IPC::System::Simple qw(capture);
use threads;
use threads::shared;
use Time::HiRes qw(gettimeofday tv_interval);

my $TOOLS = $ENV{'TOOLS'};
my $ANNOT = $ENV{'ANNOT'};
my $TLPDIR = $ENV{'TLPDIR'};

require "$TOOLS/Perlsub.pl";


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
sub read_in_primer_file;
sub check_existance_of_files;
sub update_stats_from_barcode;
sub process_experiment ($);
sub trim_sequence_to_primer ($);
sub stitch_primer_results ($);
sub trim_noprimer_results ($);
sub revert_qualtrimmed ($);
sub concat_qualtrimmed ($);
sub write_primer_stats_file;
sub write_full_stats_file;

# Global flags and arguments, 
# Set by command line arguments
my $read1;
my $read2;
my $indir;
my $outdir;
my $primer_file;
my $max_threads = 4;
my $dup_threads = 4;
my $bc_len = 10;
my $bc_mismatch = 2;
my $prim_mismatch = 0.2;
my $partial = 1;
my $dont_take_blue_rc;
my $read_length = 250;
my $binsize = 5;
my $prim_max_dif = 10;
my $prim_min_ol = 25;
my $qual_max_dif = 10;
my $qual_min_ol = 10;
my $quality = 10;
my $window = 10;
my $R1minlen = 30;
my $R2minlen = 6;
my $revert;
my $concat;
my $keepdups;

# Global variables 
my %expts :shared;
my $unmatched = 0;
my $totalreads = 0;
my @bc_output;



#
# Start of Program
#

parse_command_line;

System("perl -pi -e 's/\\r/\\n/g' $primer_file") or croak "Error: problem converting line feeds to unix format";
System("cat $primer_file");
unless (defined $indir) {
	my $bc_cmd = "$TLPDIR/TLP_barcode_splitter.pl --mate1 $read1 --mate2 $read2 --bcfile $primer_file --bol --bclen $bc_len --mismatches $bc_mismatch  --partial $partial --prefix $outdir --suffix .fastq";
	print "\n$bc_cmd\n";
	@bc_output = capture($bc_cmd) or croak "Error: TLP_barcode_splitter.pl failed";
	print "\n@bc_output\n";
}

read_in_primer_file;

check_existance_of_files;

update_stats_from_barcode;

my $PrimerDir = "$outdir/Primer/";
my $PrimerStitchDir = "$outdir/PrimerStitch/";
my $QualDir = "$outdir/Quality/";
my $QualStitchDir = "$outdir/QualityStitch/";
my $FinalDir = "$outdir/Final/";

System("mkdir $PrimerDir") or croak "Error: cannot make directory $PrimerDir";
System("mkdir $PrimerStitchDir") or croak "Error: cannot make directory $PrimerStitchDir";
System("mkdir $QualDir") or croak "Error: cannot make directory $QualDir";
System("mkdir $QualStitchDir") or croak "Error: cannot make directory $QualStitchDir";
System("mkdir $FinalDir") or croak "Error: cannot make directory $FinalDir";

my @threads = ();

foreach my $expt (sort keys %expts) {

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
						my $t1 = tv_interval($t0_expt);
						printf("\nFinished %s with %d reads in %.2f seconds.\n", $expt, $expts{$expt}->{rmdups},$t1);
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


write_primer_stats_file;

write_full_stats_file;

#
# End of program
#


sub read_in_primer_file {
	print "\nReading in primer file...";

	open PRIM, "<", $primer_file or croak "Error: Could not open primer file $primer_file for reading";
	while (<PRIM>) {
		next unless /\S/;
		chomp;
		my @row = split("\t");
		unless ($dont_take_blue_rc) {
			$row[2] = reverseComplement($row[2]);
		}



		my @temp_R1_array :shared;
		@temp_R1_array = (0) x (ceil($read_length/$binsize)+1);
		my @temp_R2_array :shared;
		@temp_R2_array = (0) x (ceil($read_length/$binsize)+1);
		
		my %temp_expt_hash :shared;
		$temp_expt_hash{R1} = $row[0]."_R1.fastq";
		$temp_expt_hash{R2} = $row[0]."_R2.fastq";
		$temp_expt_hash{red} = uc($row[1]);
		$temp_expt_hash{blue} = uc($row[2]);
		$temp_expt_hash{R1len} = \@temp_R1_array;
		$temp_expt_hash{R2len} = \@temp_R2_array;
		$temp_expt_hash{totreads} = 0;
		$temp_expt_hash{blueprim} = 0;
		$temp_expt_hash{redprim} = 0;
		$temp_expt_hash{primstitch} = 0;
		$temp_expt_hash{primstitchlen} = 0;
		$temp_expt_hash{primstitchstd} = 0;
		$temp_expt_hash{noprimstitch} = 0;
		$temp_expt_hash{qualstitch} = 0;
		$temp_expt_hash{qualstitchlen} = 0;
		$temp_expt_hash{qualstitchstd} = 0;
		$temp_expt_hash{noqualstitch} = 0;
		$temp_expt_hash{rmdups} = 0;
		$temp_expt_hash{reads_with_primer_file} = "";

		#%temp_expt_hash = {
		#    R1 => $row[0]."_R1.fastq",
		#    R2 => $row[0]."_R2.fastq",
		#    red => uc($row[1]),
		#    blue => uc($row[2]),
		#    R1len => \@temp_R1_array,
		#    R2len => \@temp_R2_array,
		#    totreads => 0,
		#    blueprim => 0,
		#    redprim => 0,
		#    primstitch => 0,
		#    primstitchlen => 0,
		#    primstitchstd => 0,
		#    noprimstitch => 0,
		#    qualstitch => 0,
		#    qualstitchlen => 0,
		#    qualstitchstd => 0,
		#    noqualstitch => 0 };

		my $expt = $row[0];
		$expts{$expt} = \%temp_expt_hash;
	}
	close PRIM;
	print "Done.\n";
	foreach my $expt (sort keys %expts) {
		printf ("%s\t%s\t%s\t%s\n", $expts{$expt}->{R1},$expts{$expt}->{R2},$expts{$expt}->{red},$expts{$expt}->{blue});
	}

}


sub check_existance_of_files {
	print "\nSearching for files...";
	foreach my $expt (sort keys %expts) {
		my $R1_file = defined $indir ? "$indir/".$expts{$expt}->{R1} : "$outdir/".$expts{$expt}->{R1};
		my $R2_file = defined $indir ? "$indir/".$expts{$expt}->{R2} : "$outdir/".$expts{$expt}->{R2};
		croak "Error: Could not locate read 1 file $R1_file" unless (-r $R1_file);
		croak "Error: Could not locate read 2 file $R2_file" unless (-r $R2_file);
	}
	print "Done.\n";
}

sub update_stats_from_barcode {
	my $flag = 0;
	unless (defined $indir) {
		foreach (@bc_output) {
			if (/^Barcode/) {
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
				$expts{$row[0]}->{totreads} = $row[1];
			}
		}
	} else {
		foreach my $expt (sort keys %expts) {
			my $linect = capture("wc -l $indir/".$expts{$expt}->{R1});
			chomp($linect);
			$linect =~ s/\s*(\d+).*/$1/ ;
			$expts{$expt}->{totreads} = $linect/4 ;
		}
	}
}

sub process_experiment ($) {
	my $expt = shift;

	trim_sequence_to_primer($expt);
	stitch_primer_results($expt);
	trim_noprimer_results($expt);
	if (defined $revert) {
		revert_qualtrimmed ($expt);
	}
	if (defined $concat) {
		concat_qualtrimmed ($expt);
	}
	print "\n$expt: Concatenating outputs of different streams\n";
	my $out1 = $PrimerStitchDir . $expt . ".join.fastq";
	my $out2 = $QualStitchDir . $expt . ".join.fastq";
	my $out3 = (defined $concat) ? ($QualStitchDir . $expt . ".fastq") : ($QualStitchDir . $expt . ".un1.fastq");
	System("cat $out1 $out2 $out3 > $FinalDir/$expt.fastq") or croak "Error: could not concatenate output streams";
	my $cd_hit_cmd = "$TOOLS/cd-hit-dup-para.pl $FinalDir/$expt.fastq $FinalDir/$expt.fastq --tmpdir $FinalDir/tmp_$expt --div 8 --threads $dup_threads > $FinalDir/${expt}_cdhit.log";
	unless (defined $keepdups) {
		System($cd_hit_cmd) or croak "Error: cd-hit-dup-para.pl failed on $expt";
	}
	$expts{$expt}->{rmdups} = (split(' ',capture("wc -l $FinalDir/$expt.fastq")))[0]/4;
	
}


sub trim_sequence_to_primer ($) {
	my $expt = shift;
	
	my $R1_file = $expts{$expt}->{R1};
	my $R2_file = $expts{$expt}->{R2};
	my $red = $expts{$expt}->{red};
	my $blue = $expts{$expt}->{blue};

	my $red_rc = reverseComplement($red);
	my $blue_rc = reverseComplement($blue);

	$expts{$expt}->{reads_with_primer_file} = "$PrimerDir/${expt}_reads_with_3primer.txt";
	my $rwpfh = IO::File->new(">".$expts{$expt}->{reads_with_primer_file}) or croak "Error: could not write to".$expts{$expt}->{reads_with_primer_file};

	my $R1_in = defined $indir ? IO::File->new("<$indir/$R1_file") : IO::File->new("<$outdir/$R1_file");
	croak "Error: could not open $R1_file for reading" unless $R1_in->opened;
	my $R1_prim_out = IO::File->new(">$PrimerDir/$R1_file");
	croak "Error: could not open $PrimerDir/$R1_file for writing" unless $R1_prim_out->opened;
	my $R2_in = defined $indir ? IO::File->new("<$indir/$R2_file") : IO::File->new("<$outdir/$R2_file");
	croak "Error: could not open $R2_file for reading" unless $R2_in->opened;
	my $R2_prim_out = IO::File->new(">$PrimerDir/$R2_file");
	croak "Error: could not open $PrimerDir/$R2_file for writing" unless $R2_prim_out->opened;


	print "\n$expt: Trimming to 3' primers\n";

	my ($red_mismatch, $blue_mismatch);
	if ($prim_mismatch > 0 && $prim_mismatch < 1) {
		$red_mismatch = floor($prim_mismatch * length($red_rc));
		$blue_mismatch = floor($prim_mismatch * length($blue_rc));
	} else {
		$red_mismatch = $prim_mismatch;
		$blue_mismatch = $prim_mismatch;
	}

	while (my @R1_seq = read_fastq($R1_in)) {
		
		my $R1_pos = 0;
		while($R1_pos <= length($R1_seq[1]) - length($blue_rc)) {
			last if (mismatch_count($blue_rc,substr($R1_seq[1],$R1_pos,length($blue_rc))) <= ceil(length($blue_rc)/2) &&
								levenshtein($blue_rc,substr($R1_seq[1],$R1_pos,length($blue_rc))) <= $blue_mismatch);
			$R1_pos++;
		}

		#while (1) {
		#    $R1_pos = index($R1_seq[1],$blue_rc,$R1_pos);
		#    last if ($R1_pos < 0);
		#    push(@R1_positions,$R1_pos);
		#    $R1_pos++;
		#}

		my @R2_seq = read_fastq($R2_in);
		my $R2_pos = 0;
		while($R2_pos <= length($R2_seq[1]) - length($red_rc)) {
			last if (mismatch_count($red_rc,substr($R2_seq[1],$R2_pos,length($red_rc))) <= ceil(length($red_rc)/2) && 
								levenshtein($red_rc,substr($R2_seq[1],$R2_pos,length($red_rc))) <= $red_mismatch);
			$R2_pos++;
		}


		#while (1) {
		#    $R2_pos = index($R2_seq[1],$red_rc,$R2_pos);
		#    last if ($R2_pos < 0);
		#    push(@R2_positions,$R2_pos);
		#    $R2_pos++;
		#}

		if ($R1_pos <= length($R1_seq[1]) - length($blue_rc)) {   # && @R2_positions > 0
			$expts{$expt}->{blueprim}++; 
			$R1_seq[1] = substr($R1_seq[1],0,$R1_pos+length($blue_rc));
			$R1_seq[3] = substr($R1_seq[3],0,$R1_pos+length($blue_rc));
			$rwpfh->print($R1_seq[0]."\n");
		}
		if ($R2_pos <= length($R2_seq[1]) - length($red_rc)) {
			$expts{$expt}->{redprim}++;
			$R2_seq[1] = substr($R2_seq[1],0,$R2_pos+length($red_rc));
			$R2_seq[3] = substr($R2_seq[3],0,$R2_pos+length($red_rc));
		}
		
		write_fastq($R1_prim_out,@R1_seq);
		write_fastq($R2_prim_out,@R2_seq);
		$expts{$expt}->{R1len}->[floor(length($R1_seq[1])/$binsize)]++;
		$expts{$expt}->{R2len}->[floor(length($R2_seq[1])/$binsize)]++;

	}

$R1_in->close;
$R1_prim_out->close;
$R2_in->close;
$R2_prim_out->close;
$rwpfh->close;

print "\n$expt primer trim stats:\nR1: ".$expts{$expt}->{blueprim}."\nR2: ".$expts{$expt}->{redprim}."\n";

}

sub stitch_primer_results ($) {
	my $expt = shift;
	
	my $R1_file = $PrimerDir . $expts{$expt}->{R1};
	my $R2_file = $PrimerDir . $expts{$expt}->{R2};
	my $output = $PrimerStitchDir . $expt . ".%.fastq";

	print "\n$expt: Stitching primer trimmed reads\n";
	unless (-z $R1_file || -z $R2_file) {
		my $stitch_cmd = "fastq-join -p $prim_max_dif -m $prim_min_ol $R1_file $R2_file -o $output";
		print "$stitch_cmd\n";
		my @stitch_output = capture($stitch_cmd);
		print "\n$expt primer stitch stats:\n@stitch_output\n";
		
		foreach (@stitch_output) {
			if (/Total joined: (\d+)/) {
				$expts{$expt}->{primstitch} = $1;
				$expts{$expt}->{noprimstitch} = $expts{$expt}->{totreads} - $1;
			} elsif (/Average join len: ([\d\.]+)/) {
				$expts{$expt}->{primstitchlen} = $1;
			} elsif (/Stdev join len: ([\d\.]+)/) {
				$expts{$expt}->{primstitchstd} = $1;
			}
		}
	} else {
		(my $un1 = $output) =~ s/%/un1/;
		(my $un2 = $output) =~ s/%/un2/;
		(my $join = $output) =~ s/%/join/;
		my $cmd = "touch $join $un1 $un2";
		System($cmd);
		$expts{$expt}->{primstitch} = 0;
		$expts{$expt}->{noprimstitch} = $expts{$expt}->{totreads};
		$expts{$expt}->{primstitchlen} = 0;
		$expts{$expt}->{primstitchstd} = 0;
	}



}


sub trim_noprimer_results ($) {
	my $expt = shift;

	my $un1 = $PrimerStitchDir . $expt . ".un1.fastq";
	my $un2 = $PrimerStitchDir . $expt . ".un2.fastq";

	print "\n$expt: Trimming by quality scores\n";
	System("$TLPDIR/TLP_quality_trimmer.pl $un1 $un2 $QualDir --quality $quality --window $window --r1minlength $R1minlen --r2minlength $R2minlen");

	$un1 =~ s/$PrimerStitchDir/$QualDir/;
	$un2 =~ s/$PrimerStitchDir/$QualDir/;

	my $R1_file = $QualDir . $expts{$expt}->{R1};
	my $R2_file = $QualDir . $expts{$expt}->{R2};

	System("mv $un1 $R1_file") or croak "Error: there was a problem renaming files";
	System("mv $un2 $R2_file") or croak "Error: there was a problem renaming files";

	my $output = $QualStitchDir . $expt . ".%.fastq";


	print "\n$expt: Stitching quality trimmed reads\n";
	unless (-z $R1_file || -z $R2_file) {
		my $stitch_cmd = "fastq-join -p $qual_max_dif -m $qual_min_ol $R1_file $R2_file -o $output";
		print "$stitch_cmd\n";
		my @stitch_output = capture($stitch_cmd);
		print "\n$expt quality stitch stats:\n@stitch_output\n";
		
		foreach (@stitch_output) {
			if (/Total joined: (\d+)/) {
				$expts{$expt}->{qualstitch} = $1;
				$expts{$expt}->{noqualstitch} = $expts{$expt}->{noprimstitch} - $1;
			} elsif (/Average join len: ([\d\.]+)/) {
				$expts{$expt}->{qualstitchlen} = $1;
			} elsif (/Stdev join len: ([\d\.]+)/) {
				$expts{$expt}->{qualstitchstd} = $1;
			}
		}
	} else {
		(my $un1 = $output) =~ s/%/un1/;
		(my $un2 = $output) =~ s/%/un2/;
		(my $join = $output) =~ s/%/join/;
		my $cmd = "touch $join $un1 $un2";
		System($cmd);
		$expts{$expt}->{qualstitch} = 0;
		$expts{$expt}->{noqualstitch} = $expts{$expt}->{noprimstitch};
		$expts{$expt}->{qualstitchlen} = 0;
		$expts{$expt}->{qualstitchstd} = 0;
	}

}

sub concat_qualtrimmed ($) {
	my $expt = shift;
	print "\n$expt: Concatentating unstitched reads\n";

	my $infile1 = "$QualStitchDir/$expt.un1.fastq";
	my $infile2 = "$QualStitchDir/$expt.un2.fastq";

	my $concat_len = 2 * $read_length + $concat;

	my $cmd = "$TLPDIR/TLP_fastx_joiner.pl $infile1 $infile2 --length $concat_len --noconcat ".$expts{$expt}->{reads_with_primer_file};
	System($cmd) or croak "Error: could not concatenate unstitched reads";
}

sub revert_qualtrimmed ($) {
	my $expt = shift;
	
	print "\n$expt: Reverting unstitched reads to length after primer trim step\n";

	my $reffile = "$QualStitchDir/$expt.join.fastq";
	my $out1file = "$QualStitchDir/$expt.un1.fastq";
	my $out2file = "$QualStitchDir/$expt.un2.fastq";
	my $in1file = "$PrimerStitchDir/$expt.un1.fastq";
	my $in2file = "$PrimerStitchDir/$expt.un2.fastq";

	my $ref = IO::File->new("<$reffile");
	croak "Error: could not open $reffile for reading" unless $ref->opened;
	my $out1 = IO::File->new(">$out1file");
	croak "Error: could not open $out1file for writing" unless $out1->opened;
	my $out2 = IO::File->new(">$out2file");
	croak "Error: could not open $out2file for writing" unless $out2->opened;
	my $in1 = IO::File->new("<$in1file");
	croak "Error: could not open $in1file for reading" unless $in1->opened;
	my $in2 = IO::File->new("<$in2file");
	croak "Error: could not open $in2file for reading" unless $in2->opened;
	
	my %qualstitch_hash;
	while (my @joinseq = read_fastq($ref)) {
		$qualstitch_hash{$joinseq[0]} = 1;
	}

	while (my @seq = read_fastq($in1)) {
		write_fastq($out1,@seq) unless (defined $qualstitch_hash{$seq[0]});
	}

	while (my @seq = read_fastq($in2)) {
		write_fastq($out2,@seq) unless (defined $qualstitch_hash{$seq[0]});
	}

	$ref->close;
	$out1->close;
	$out2->close;
	$in1->close;
	$in2->close;

}

	
sub write_primer_stats_file {
	print "\nWriting primer stats file...";
	open STATS, ">", "$outdir/primer-trim-stats.txt" or croak "Error: Can not open primer stats file for writing";
	my @header = ("Bin");
	foreach my $expt (sort keys %expts) {
		push(@header,$expts{$expt}->{R1},$expts{$expt}->{R2});
	}
	print STATS join("\t",@header)."\n";
	foreach my $bin (0 .. ($read_length/$binsize) ) {
		my @row = ($bin*$binsize);
		foreach my $expt (sort keys %expts) {
			push(@row,$expts{$expt}->{R1len}->[$bin],$expts{$expt}->{R2len}->[$bin]);
		}
		print STATS join("\t",@row)."\n";
	}
	print "Done\n";
}

sub write_full_stats_file {
	print "\nWriting full stats file...";
	open STATS, ">", "$outdir/full-stats.txt" or croak "Error: Can not open total stats file for writing";
	my @header = qw(Expt De-multiplex R1_3Primer R2_3Primer Stitch StitchLen StitchStd NoStitch QualStitch QualStitchLen QualStitchStd NoQualStitch RmDups);
	print STATS join("\t",@header)."\n";
	foreach my $expt (sort keys %expts) {
		my @row = ($expt,
			$expts{$expt}->{totreads},
			$expts{$expt}->{blueprim},
			$expts{$expt}->{redprim},
			$expts{$expt}->{primstitch},
			$expts{$expt}->{primstitchlen},
			$expts{$expt}->{primstitchstd},
			$expts{$expt}->{noprimstitch},
			$expts{$expt}->{qualstitch},
			$expts{$expt}->{qualstitchlen},
			$expts{$expt}->{qualstitchstd},
			$expts{$expt}->{noqualstitch},
			$expts{$expt}->{rmdups} );
		print STATS join("\t",@row)."\n";
	}
	unless (defined $indir) {
		print STATS "Unmatched (PhiX)\t$unmatched\n";
		print STATS "Total\t$totalreads\n";
	}
	print "Done\n";
}
	

sub parse_command_line {
  my $help;

  usage() if (scalar @ARGV==0);

  my $result = GetOptions (
  			"read1=s" => \$read1,
  			"read2=s" => \$read2,
  			"indir=s" => \$indir,
				"threads=i" => \$max_threads,
				"dupthreads=i" => \$dup_threads, 
				"bclen=i" => \$bc_len,
				"bcmismatch=i" => \$bc_mismatch,
				"bcpartial=i" => \$partial,
				"primmismatch=f" => \$prim_mismatch,
				"dont_take_blue_rc" => \$dont_take_blue_rc,
				"readlength=i" => \$read_length,
				"binsize=i" => \$binsize,  
				"max_dif_prim=i" => \$prim_max_dif,
				"min_ol_prim=i" => \$prim_min_ol,
				"max_dif_qual=i" => \$qual_max_dif,
				"min_ol_qual=i" => \$qual_min_ol,
				"quality=i" => \$quality,
				"window=i" => \$window,
				"r1minlength=i" => \$R1minlen,
				"r2minlength=i" => \$R2minlen,
				"revert" => \$revert,
				"concat=i" => \$concat,
				"keepdups" => \$keepdups,
				"help" => \$help
		  ) ;
  
  usage() if ($help);
	
	$primer_file = shift(@ARGV);
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
		System("mkdir $outdir") or croak "Error: output directory $outdir neither exists nor could be created";
	}
	$outdir .= "/" unless ($outdir =~ /\/$/);
	croak "Error: could not find or read primer file $primer_file" unless (-r $primer_file);

  exit unless $result;
}


sub usage()
{
print<<EOF;
Title, by Robin Meyers, ddmonthyyyy

This program .


Usage: $0 primer_file outdir (--read1 FILE --read2 FILE | --indir DIR)[--opts N] [--help]

Arguments (defaults in parentheses):


$arg{"primer_file","List of files with primer sequences"}
$arg{"outdir","Output directory"}
$arg{"--read1","Illumina read 1 file for de-multiplexing"}
$arg{"--read2","Illumina read 2 file for de-multiplexing"}
$arg{"--indir","Input directory - if specified assumes reads are already de-multiplexed"}
$arg{"--threads","Number of threads to run program on - experiments processed simultaneously",$max_threads}
$arg{"--dupthreads","Number of threads to run each experiment through cd-hit-dup-para",$dup_threads}
$arg{"--bclen","Number of bases to use from primer file for de-multiplexing",$bc_len}
$arg{"--bcmismatch","Number of mismatches allowed in de-multiplexing",$bc_mismatch}
$arg{"--partial","Allow partial overlap of barcode by N bases - see TLP_barcode_splitter.pl", $partial}
$arg{"--primmismatch","When searching for 3' primer: maximum number of mismatches if integer; maximum percent mismatches if 0<x<1",$prim_mismatch}
$arg{"--max_dif_prim","Maximum percent difference between reads for primer-trimmed stitching",$prim_max_dif}
$arg{"--min_ol_prim","Minimum basepair overlap between two reads for primer-trimmed stitching",$prim_min_ol}
$arg{"--max_dif_qual","Maximum percent difference between reads for quality-trimmed stitching",$qual_max_dif}
$arg{"--min_ol_qual","Minimum basepair overlap between two reads for quality-trimmed stitching",$qual_min_ol}
$arg{"--quality","Minimum quality score threshold",$quality}
$arg{"--window","Size of window (in bp) to slide across read for trimming by quality",$window}
$arg{"--r1minlength","Minimum length of quality trimmed reads",$R1minlen}
$arg{"--r2minlength","Will not quality trim read 2 sequences below this length",$R2minlen}
$arg{"--dont_take_blue_rc","Use this flag when the blue primer is in the barcode file as if it were read on read2"}
$arg{"--readlength","Length of Illumina reads",$read_length}
$arg{"--binsize","Bin size to group trimmed read lengths in stats file",$binsize}
$arg{"--revert","Revert qual-trimmed unstitched reads to original length"}
$arg{"--concat","Join qual-trimmed unstitched reads, padded with Ns to length two times raw reads + this"}
$arg{"--keepdups","Don't run cd-hit-dup step at end of pre-processing pipeline"}
$arg{"--help","This helpful help screen."}


EOF

exit 1;
}
