# Source for mostly smaller subroutines

use Switch::Plain;

### Program execution

sub debug_print ($$;$) {

  my $msg = shift;
  my $code = shift;
  my $id = shift;

  return if ! defined $main::debug_level || $main::debug_level < $code;

  my $new_msg = $code . "|";
  $new_msg = $new_msg . $id . "|" if defined $id;
  $new_msg = $new_msg . $msg . "\n";

  print $new_msg;

}

sub argument {
	my $var = shift;
	my $description = shift;
	my $default = shift;

	return sprintf("  \%\-16s - %s",
		(defined $default ? sprintf("%s (%s)",$var,$default) : sprintf("%s",$var)),
		$description);
}

sub parseFilename ($) {
	my $fullname = shift;
	my ($name, $path, $ext) = fileparse($fullname, qr/\.\w{2,5}$/);
	return ($path, $name, $ext);
}

sub System ($;$) {
  my $cmd = shift;
  my $quiet = shift;
  $quiet = 0 unless defined $quiet;
  print "$cmd\n" unless $quiet;
  my $status = system($cmd);
  return !$status;
}

sub Capture ($;$) {
  my $cmd = shift;
  my $quiet = shift;
  $quiet = 0 unless defined $quiet;
  print "$cmd\n" unless $quiet;
  my @output = capture($cmd);
  return @output;
}

sub check_undef ($;$) {
  my $val = shift;
  my $replace = shift;
  $replace = defined $replace ? $replace : "";
  $val = defined $val ? $val : $replace;
  return $val;
}

### File formatting

sub tlx_header {
  my @tlx_header = ( qw(Qname JuncID Rname Junction Strand Rstart Rend),
                      qw(B_Rname B_Rstart B_Rend B_Strand B_Qstart B_Qend),
                      qw(Qstart Qend Qlen B_Cigar Cigar Seq J_Seq Barcode) );
  return (@tlx_header);
}

sub tlx_filter_header {
  my @tlx_header = tlx_header();
  splice(@tlx_header,1,0,("Filter"));
  return(@tlx_header);
}

sub tlxl_header {
  my @tlxl_header = ( qw(Qname OCS_score Rname Strand),
                      qw(R1_Rstart R1_Rend R1_Qstart R1_Qend R1_Qlen),
                      qw(R2_Rstart R2_Rend R2_Qstart R2_Qend R2_Qlen),
                      qw(R1_Cigar R2_Cigar) );
  return (@tlxl_header);
}

sub write_filter_entry ($$$$) {
  my $fh = shift;
  my $entry = shift;
  my $header = shift;
  my $filters = shift;

  my @part1 = map {check_undef($_,"")} @{$entry}{@$header};
  my @part2 = map {check_undef($_,"")} @{$entry->{filters}}{@$filters};

  $fh->print(join("\t",@part1,@part2)."\n");

}

sub write_entry ($$$) {

  my $fh = shift;
  my $entry = shift;
  my $header = shift;


  $fh->print(join("\t",map(check_undef($_,""),@{$entry}{@$header}))."\n");
  
}


### Math

sub mean {
	return sum(@_)/@_;
}


sub log2 ($) {
	my $n = shift;
	return log($n)/log(2);
}

sub log10 ($) {
	my $n = shift;
	return log($n)/log(10);
}

### Sequence manipulation


sub reverseComplement ($) {
	my $seq = shift;
	(my $rc = reverse($seq)) =~ tr/ACGTacgtNn/TGCAtgcaNn/;
	return $rc;
}


sub nearby_quals ($$$) {
  my $qual_ref = shift;
  my $pos = shift;
  my $mar = shift;

  return mean(@$qual_ref[max(0,$pos-$mar)..min($pos+$mar,$#{$qual_ref})]);
}

### Alignment subroutines

sub sam_file_is_empty ($) {
  my $file = shift;
  my $fh = IO::File->new("<$file");

  while (my $line = $fh->getline) {
    next unless $line =~ /\S/;
    next if $line =~ /^@/;
    $fh->close;
    return(0);
  }
  $fh->close;
  return(1);

}


sub stringify_alignment ($) {
  my $aln = shift;

  return(join(",",$aln->{Rname},
                  $aln->{Rstart},
                  $aln->{Rend},
                  $aln->{Qstart},
                  $aln->{Qend},
                  $aln->{Strand}));
}


sub aln_reference_length ($) {
  my $aln = shift;
  return $aln->{Rend} - $aln->{Rstart} + 1;
}

sub aln_reference_overlap ($$) {
  my $aln1 = shift;
  my $aln2 = shift;

  return 0 unless $aln1->{Rname} eq $aln2->{Rname};

  my $overlap_length = min($aln1->{Rend},$aln2->{Rend}) - max($aln1->{Rstart},$aln2->{Rstart}) + 1;

  if ($overlap_length > 0) {
    return $overlap_length;
  } else {
    return 0;
  }
}

sub aln_query_overlap ($$) {
  my $aln1 = shift;
  my $aln2 = shift;

  my $overlap_length = min($aln1->{Qend},$aln2->{Qend}) - max($aln1->{Qstart},$aln2->{Qstart}) + 1;

  if ($overlap_length > 0) {
    return $overlap_length;
  } else {
    return 0;
  }
}

sub aln_query_length ($) {
  my $aln = shift;
  return $aln->{Qend} - $aln->{Qstart} + 1;
}


sub calculate_fraction_overlap ($$;$$) {
  my $aln1_R1 = shift;
  my $aln2_R1 = shift;
  my $aln1_R2 = shift;
  my $aln2_R2 = shift;

  my $aln1_length = 0;
  my $overlap_length = 0;

  $aln1_length += aln_query_length($aln1_R1);
  my $ol = aln_query_overlap($aln1_R1,$aln2_R1);
  $overlap_length += $ol if $ol > 0;

  if (defined $aln1_R2) {

    $aln1_length += aln_query_length($aln1_R2);
    my $ol = aln_query_overlap($aln1_R2,$aln2_R2);
    $overlap_length += $ol if $ol > 0;
  }

  return($overlap_length/$aln1_length);

}


sub pair_is_proper ($$) {
  my $R1 = shift;
  my $R2 = shift;
  # my $max_frag_len = shift;
  # my $max_dovetail = 10;

  return 0 unless $R1->{Rname} eq $R2->{Rname};
  return 0 unless $R1->{Strand} == $R2->{Strand};


  if ($R1->{Strand} == 1) {
    return 0 unless $R1->{Rstart} < $R2->{Rend};
    return 0 if $R2->{Rstart} - $R1->{Rend} - 1 > $main::params->{max_pe_gap};
    return 0 if $R1->{Rstart} > $R2->{Rstart} + $main::params->{max_dovetail};
    return 0 if $R1->{Rend} > $R2->{Rend} + $main::params->{max_dovetail};
  } else {
    return 0 unless $R2->{Rstart} < $R1->{Rend};
    return 0 if $R1->{Rstart} - $R2->{Rend} - 1 > $main::params->{max_pe_gap};
    return 0 if $R2->{Rstart} > $R1->{Rstart} + $main::params->{max_dovetail};
    return 0 if $R2->{Rend} > $R1->{Rend} + $main::params->{max_dovetail};
  }

  return 1;
}


### TLX Junction subroutines


sub is_a_junction ($) {
  my $tlx = shift;
  if (! defined $tlx->{Rname} || $tlx->{Rname} eq "" || $tlx->{Rname} eq "Adapter") {
    return(0);
  } else {
    return(1);
  }
}

### CIGAR subroutines

sub cigar_array_to_string ($) {
  my $array_ref = shift;
  my $string = join("",map {join("", reverse @$_)} @$array_ref);
}

sub cigar_string_to_array ($) {
  my $string = shift;
  my @array;
  while ($string =~ /([0-9]+)([SMXIDN])/g) {
    push(@array,[$2,$1]);
  }
  return \@array;
}


sub soft_clip_cigar ($) {
  my $cigar_ref = shift;
  my @cigar = @$cigar_ref;
  if ($cigar[0]->[0] eq "S") {
    shift(@cigar);
  }
  if ($cigar[$#cigar]->[0] eq "S") {
    pop(@cigar);
  }
  return \@cigar;
}

sub expand_cigar_array ($) {
  my $compact_cigar = shift;
  my @expand_cigar;
  foreach my $i (@$compact_cigar) {
    # next if $i->[0] eq "S";
    push(@expand_cigar, ($i->[0]) x $i->[1]);
  }
  return \@expand_cigar;
}

sub compact_cigar_array ($) {
  my $expand_cigar = shift;
  my @compact_cigar = ();

  my $prev_code = $expand_cigar->[0];
  my $count = 1;
  return [[$prev_code,$count]] if $#$expand_cigar < 1;

  foreach my $i (1..$#$expand_cigar) {
    if ($expand_cigar->[$i] eq $prev_code) {
      $count++;
    } else {
      push(@compact_cigar,[$prev_code,$count]);
      $prev_code = $expand_cigar->[$i];
      $count = 1;
    }
  }

  push(@compact_cigar,[$prev_code,$count]);

  return \@compact_cigar;
}


sub remap_cigar ($$$) {

  my $cigar_array_ref = shift;
  my $qseq = shift;
  my $rseq = shift;

  my @Rseq = split("",uc($rseq));
  my @Qseq = split("",uc($qseq));

  my $Qpos = 1;
  my $Rpos = 1;
  # print Dumper($cigar_array_ref);
  my @old_cigar = @{ expand_cigar_array($cigar_array_ref) };
  my @new_cigar;

  debug_print("remapping a cigar string\nold: @old_cigar\nref: $rseq\nqry: $qseq\n",3);

  while ($Qpos <= @Qseq) {
    # print "$Qpos\n";
    my $c = shift(@old_cigar);
    sswitch ($c) {
      case 'S': {
        push(@new_cigar,"S");
        $Qpos++;
      }
      case 'M':
      case 'X': {
        if ($Qseq[$Qpos-1] eq $Rseq[$Rpos-1]) {
          push(@new_cigar,"M");
        } else {
          push(@new_cigar,"X")
        }
        $Qpos++;
        $Rpos++;
      }
      case 'D': {
        push(@new_cigar,"D");
        $Rpos++;
      }
      case 'I': {
        push(@new_cigar,"I");
        $Qpos++;
      }
    }
  }
  return(compact_cigar_array(\@new_cigar));
}



### Filtering subroutines

sub filter_entire_read ($$;$) {
  my $tlxs = shift;
  my $filter = shift;
  my $value = shift;
  $value = 1 unless defined $value;

  filter_remainder_of_read($tlxs,$filter,0,$value);

}

sub filter_remainder_of_read ($$$;$) {
  my $tlxs = shift;
  my $filter = shift;
  my $i = shift;
  my $value = shift;
  $value = 1 unless defined $value;

  foreach my $tlx (@$tlxs[$i..$#$tlxs]) {
    $tlx->{filters}->{$filter} = $value;
  }
}

1;