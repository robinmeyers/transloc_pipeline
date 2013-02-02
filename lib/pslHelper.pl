use IO::File;
use List::MoreUtils qw(pairwise);



sub compressHeader ($$) {

	my $psl = shift;
	my $rmbak = shift;
	my $bak = "$psl.bak";

	system("mv $psl $bak");

	my $pslx = 1 if ($psl =~ /pslx$/);

	my $pslio = IO::File->new(">$psl");
	my $bakio = IO::File->new("<$bak");

	my $nrow = 0;
	my (@header1,@header2);

	while ( defined($_ = $bakio->getline) ) {
		chomp;
		next unless /\S/;
		next if /^psLayout/;
		next if /^-------/;
		$pslio->print($_."\n") if (/mis-match/); #header has already been compressed
		if (/^match/) {
			@header1 = split("\t");
			if ($pslx) {
				push(@header1,"qSeq","tSeq");
			}
			foreach(@header1) {s/\s//g};
			next;
		}
		if (/^\s+match/) {
			@header2 = split("\t");
			foreach(@header2) {s/\s//g};
			$pslio->print( join("\t", pairwise {$a.$b} @header1,@header2)."\n" );
			next;
		}
		$pslio->print($_."\n");
	}
	system("rm $bak") if $rmbak;

}

1;
