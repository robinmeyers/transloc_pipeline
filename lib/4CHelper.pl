sub create_sequence_files ($$$) {

	my $outdir = shift;
	my $expt_id = shift;
	my $expt_hash = shift;

	my $seqdir = "$outdir/$expt_id/sequences";

	System("mkdir -p $seqdir") or croak "Error: could not create sequences directory for $expt_id";

	my $redpf = "$seqdir/redp.fa";
	my $redcf = "$seqdir/redc.fa";
	my $blupf = "$seqdir/blup.fa";
	my $blucf = "$seqdir/bluc.fa";

	my $redpio = IO::File->new(">$redpf") or croak "Error: could not open $redpf for writing";
	print $redpio "> $expt_id Red Primer\n";
	print $redpio uc($expt_hash->{redp});
	$redpio->close;

	my $redcio = IO::File->new(">$redcf") or croak "Error: could not open $redcf for writing";
	print $redcio "> $expt_id Red Continued (MID + Primer + Cont)\n";
	print $redcio uc( $expt_hash->{mid} . $expt_hash->{redp} . $expt_hash->{redc} );
	$redcio->close;

	my $blupio = IO::File->new(">$blupf") or croak "Error: could not open $blupf for writing";
	print $blupio "> $expt_id Blu Primer\n";
	print $blupio uc($expt_hash->{blup});
	$blupio->close;

	my $blucio = IO::File->new(">$blucf") or croak "Error: could not open $blucf for writing";
	print $blucio "> $expt_id Blu Continued (Cont + Primer)\n";
	print $blucio uc( $expt_hash->{bluc} . $expt_hash->{blup} );
	$blucio->close;


}
















1;
