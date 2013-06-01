use strict;
use warnings;
use Switch;
use Bio::Factory::EMBOSS;
use List::Util qw(min max);

sub manage_program_options ($$) {

  my $defaultoptstr = shift;
  my $useroptstr = shift;

  return $defaultoptstr unless $useroptstr =~ /\S/;
  return $useroptstr unless $defaultoptstr =~ /\S/;

  my @defaultopt = split(/\s+/,$defaultoptstr);
  my @useropt = split(/\s+/,$useroptstr);

  my %opt_hash;
  my $optkey;
  my $optval;

  my @return_opt = ();

  my $i = 0;

  while ($i < @defaultopt) {
    if ( defined $defaultopt[$i+1] && $defaultopt[$i+1] !~ /^-/ ) {
      $opt_hash{$defaultopt[$i]} = $defaultopt[$i+1];
      $i += 2;
    } else {
      $opt_hash{$defaultopt[$i]} = "";
      $i += 1;
    }
  }

  $i = 0;

  while ($i < scalar @useropt) {
    if ( defined $useropt[$i+1] && $useropt[$i+1] !~ /^-/ ) {
      if ($useropt[$i+1] eq "OFF") {
        delete $opt_hash{$useropt[$i]};
      } else {
        $opt_hash{$useropt[$i]} = $useropt[$i+1];
      }
      $i += 2;
    } else {
      $i += 1;
    }
  }
  
  foreach $optkey (sort keys %opt_hash) {
    $optval = $opt_hash{$optkey};
    push(@return_opt,$optkey,$optval);
  }

  return(join(" ",@return_opt));


}



sub find_orientations ($$$) {
  my $alns = shift;
  my $R1_brk_aln = shift;
  my $R2_brk_aln = shift;

  foreach my $aln (@$alns) {

    my $R1_aln = $aln->{R1_aln};
    my $R2_aln = $aln->{R2_aln};

    my ($R1_Qstart,$R1_Qend,$R2_Qstart,$R2_Qend);

    if (defined $R1_aln) {
      $R1_Qstart = $R1_aln->reversed ? $R1_aln->l_qseq - $R1_aln->query->end + 1 : $R1_aln->query->start;
      $R1_Qend = $R1_aln->reversed ? $R1_aln->l_qseq - $R1_aln->query->start + 1 : $R1_aln->query->end;
    }

    if (defined $R2_aln) {
      $R2_Qstart = $R2_aln->reversed ? $R2_aln->query->start : $R2_aln->l_qseq - $R2_aln->query->end + 1;
      $R2_Qend = $R2_aln->reversed ? $R2_aln->query->end : $R2_aln->l_qseq - $R2_aln->query->start + 1;
    }

    # Test for orientation #1
    if ($R1_brk_aln->proper_pair && defined $R1_aln && $R1_aln->proper_pair) {

      if ($R1_brk_aln->query->end < $R1_Qend && $R2_brk_aln->query->end < $R2_Qend) {
        $aln->{orientation} = 1;
        next;
      }

    }

    # Test for orientation #2
    if ( ! $R1_brk_aln->unmapped && defined $R1_aln && $R1_aln->proper_pair) {

      if ($R1_brk_aln->query->end < $R1_Qend) {
        $aln->{orientation} = 2;
        next;
      }

    }

    # Test for orientation #3
    if ( $R1_brk_aln->proper_pair && defined $R2_aln) {

      if ($R2_brk_aln->query->end < $R2_Qend) {
        $aln->{orientation} = 3;
        next;
      }

    }

    # Test for orientation #4
    if ( ! $R1_brk_aln->unmapped && defined $R1_aln) {

      next if (! $R2_brk_aln->unmapped && $R2_brk_aln->start <= $R1_brk_aln->end);

      if ($R1_brk_aln->query->end < $R1_Qend) {
        $aln->{orientation} = 4;
        next;
      }

    }

    $aln->{orientation} = 0;
  }

}

sub filter_mapquals ($$$) {
  my $alns = shift;
  my $R1_brk_aln = shift;
  my $R2_brk_aln = shift;

  my $threshold = 5;

  my $mapq;

  foreach my $aln (@$alns) {

    next unless defined $aln->{orientation} && $aln->{orientation} > 0;

    my $R1_aln = $aln->{R1_aln};
    my $R2_aln = $aln->{R2_aln};

    if (defined $mapq) {
      $aln->{mapq} = 0;
      next;
    }

    if (defined $R1_aln && $R1_aln->qual > $threshold) {
      $mapq = $R1_aln->qual;
      $aln->{mapq} = $R1_aln->qual;
    } elsif (defined $R2_aln && $R2_aln->qual > $threshold) {
      $mapq = $R2_aln->qual;
      $aln->{mapq} = $R2_aln->qual;
    } else {
      $aln->{mapq} = 0;
    }

  }

}

sub merge_alignments ($$$) {
  my $aln1 = shift;
  my $aln2 = shift;
  my $ref = shift;

  if ($aln1->reversed) {
    ($aln1,$aln2) = ($aln2,$aln1);
  }

  my $rname = $aln1->seq_id;
  my $rstart = $aln1->start;
  my $rend = $aln2->end;

  my @refseq = split("",$ref->seq($rname,$rstart,$rend));

  my $cigar1ref = $aln1->cigar_array;
  my @cigar1 = ();
  foreach my $i (@$cigar1ref) {
    next if $i->[0] eq "S";
    push(@cigar1,split("",($i->[0] x $i->[1])));
  }

  my $cigar2ref = $aln2->cigar_array;
  my @cigar2 = ();
  foreach my $i (@$cigar2ref) {
    next if $i->[0] eq "S";
    push(@cigar2,split("",($i->[0] x $i->[1])));
  }


  my @qseq1 = split("",substr($aln1->qseq,$aln1->query->start-1,$aln1->query->length));
  my @qseq2 = split("",substr($aln2->qseq,$aln2->query->start-1,$aln2->query->length));

  my @qual1 = @{$aln1->qscore}[$aln1->query->start-1..$aln1->query->end-1];
  my @qual2 = @{$aln2->qscore}[$aln2->query->start-1..$aln2->query->end-1];
  
  print $aln1->start ." ".$aln1->end."\n";
  print join("",@cigar1)."\n";
  print join("",@qseq1)."\n";
  print join(" ",@qual1)."\n\n";
  print $aln2->start ." ".$aln2->end."\n";
  print join("",@cigar2)."\n";
  print join("",@qseq2)."\n";
  print join(" ",@qual2)."\n\n";


  print "$rstart $rend\n";
  print join("",@refseq)."\n";


  my $pos1 = 0;
  my $cpos1 = 0;
  my $pos2 = 0;
  my $cpos2 = 0;
  my $rpos = $rstart;

  my @consensus = ();
  my @quality = ();

  while ( $rpos <= $aln1->end && $rpos < $aln2->start  && $rpos <= $rend) {

    my $c1 = $cigar1[$cpos1++];

    switch ($c1) {
      case "D" { $rpos++; }
      case "I" { push(@consensus,$qseq1[$pos1++]); }
      case "M" { push(@consensus,$qseq1[$pos1++]); $rpos++; }
      else { croak "Error: unrecognized character in CIGAR string $c1"; }
    }

  }

  while ( $rpos > $aln1->end && $rpos < $aln2->start  && $rpos <= $rend) {
    push(@consensus,lc($refseq[$rpos++ - $rstart]));
  }

  while ($rpos <= $aln1->end && $rpos >= $aln2->start && $rpos <= $rend) {

    my $c1 = $cigar1[$cpos1++];
    my $c2 = $cigar2[$cpos2++];

    print join("\t",$rpos,$c1,$c2,$pos1,$pos2)."\n";

    switch ($c1) {
      case "D" {
        switch ($c2) {
          case "D" { $rpos++; }
          case "I" {
            if (mean(@qual1[max(0,$pos1-2)..min($pos1+2,$#qual1)]) >=
                mean(@qual2[max(0,$pos2-2)..min($pos2+2,$#qual2)])) {
              $pos2++;              
              $cpos1--;
            } else {
              push(@consensus,$qseq2[$pos2++]);
              $cpos1--;
            }
           }
          case "M" { 
            if (mean(@qual1[max(0,$pos1-2)..min($pos1+2,$#qual1)]) >=
                mean(@qual2[max(0,$pos2-2)..min($pos2+2,$#qual2)])) {
              $rpos++;
            } else {
              push(@consensus,$qseq2[$pos2++]);
              $rpos++;
            }
          }
          else { croak "Error: unrecognized character '$c2' in CIGAR string"; }
        }
      }
      case "I" {
        switch ($c2) {
          case "D" { $pos1++; $cpos2--;
            if (mean(@qual1[max(0,$pos1-2)..min($pos1+2,$#qual1)]) >=
                mean(@qual2[max(0,$pos2-2)..min($pos2+2,$#qual2)])) {
              push(@consensus,$qseq1[$pos1++]);
              $cpos2--;
            } else {
              $pos1++;
              $cpos2--;
            }
          }
          case "I" {
            if (mean(@qual1[max(0,$pos1-2)..min($pos1+2,$#qual1)]) >=
                mean(@qual2[max(0,$pos2-2)..min($pos2+2,$#qual2)])) {
              push(@consensus,$qseq1[$pos1++]);
              $pos2++;
            } else {
              push(@consensus,$qseq2[$pos2++]);
              $pos1++;
            }
          }
          case "M" {
            if (mean(@qual1[max(0,$pos1-2)..min($pos1+2,$#qual1)]) >=
                mean(@qual2[max(0,$pos2-2)..min($pos2+2,$#qual2)])) {
              push(@consensus,$qseq1[$pos1++]);
              $cpos2--;
            } else {
              $pos1++;
              $cpos2--;
            }
          }
          else { croak "Error: unrecognized character '$c2' in CIGAR string"; }
        }
      }
      case "M" {
        switch ($c2) {
          case "D" {
            if (mean(@qual1[max(0,$pos1-2)..min($pos1+2,$#qual1)]) >=
                mean(@qual2[max(0,$pos2-2)..min($pos2+2,$#qual2)])) {
              push(@consensus,$qseq1[$pos1++]);
              $rpos++;
            } else {
              $rpos++;
            }
          }
          case "I" {
            if (mean(@qual1[max(0,$pos1-2)..min($pos1+2,$#qual1)]) >=
                mean(@qual2[max(0,$pos2-2)..min($pos2+2,$#qual2)])) {
              $pos2++;
              $cpos1--;
            } else {
              push(@consensus,$qseq2[$pos2++]);
              $cpos1--;
            }
          }
          case "M" { 
            if (mean(@qual1[max(0,$pos1-2)..min($pos1+2,$#qual1)]) >=
                mean(@qual2[max(0,$pos2-2)..min($pos2+2,$#qual2)])) {
              push(@consensus,$qseq1[$pos1++]);
              $pos2++;
              $rpos++;
            } else {
              push(@consensus,$qseq2[$pos2++]);
              $pos1++;
              $rpos++;
            }
          }
          else { croak "Error: unrecognized character '$c2' in CIGAR string" ; }
        }
      }
      else { croak "Error: unrecognized character '$c1' in CIGAR string"; }
    }

  }

  while ( $rpos >= $aln2->start && $rpos <= $rend ) {
    
    my $c2 = $cigar2[$cpos2++];

    switch ($c2) {
      case "D" { $rpos++; }
      case "I" { push(@consensus,$qseq2[$pos2++]); }
      case "M" { push(@consensus,$qseq2[$pos2++]); $rpos++; }
      else { croak "Error: unrecognized character in CIGAR string $c2"; }
    }
  } 

  print "\n".join("",@consensus)."\n\n\n";

  return(join("",@consensus));

}

sub sw_align_pairs ($$$) {
  my $R1_aln = shift;
  my $R2_aln = shift;
  my $tmpdir = shift;

  my $f = Bio::Factory::EMBOSS -> new();
  my $water = $f->program('water');

  my $R1_seq = $R1_aln->query->seq;
  my $R2_seq = $R2_aln->query->seq;

  (my $R1_id = $R1_seq->id) =~ s/\W/_/g;

  my $tmpout = join("/",$tmpdir,$R1_id).".water";

  $water->run({ -asequence => $R1_seq,
                -bsequence => $R2_seq,
                -gapopen   => '4.0',
                -gapextend => '2.0',
                -datafile  => 'EDNASIMPLE4',
                -outfile   => $tmpout });

  return "";



}


sub tlx_header {
  my @tlx_header = ( qw(Qname Rname Junction Strand Rstart Rend),
                      qw(BrkRstart BrkRend BrkQstart BrkQend Qstart Qend),
                      qw(Qlen MapQ Orientation Seq) );
  return (@tlx_header);
}


sub tlxl_header {
  my @tlxl_header = ( qw(Qname Filter MapQ Rname Strand),
                      qw(R1_BrkRstart R1_BrkRend R1_Rstart R1_Rend),
                      qw(R1_BrkQstart R1_BrkQend R1_Qstart R1_Qend),
                      qw(R2_BrkRstart R2_BrkRend R2_Rstart R2_Rend),
                      qw(R2_BrkQstart R2_BrkQend R2_Qstart R2_Qend),
                      qw(R1_BrkCigar R1_Cigar R2_BrkCigar R2_Cigar R1_seq R2_seq) );
  return (@tlxl_header);
}




sub make_raw_fastq_idx ($) {
  my $fq = shift;

  

  my $fqfh = IO::File->new("<$fq");
  my $idx = "$fq.idx";
  my $idxfh = IO::File->new("+>$idx");

  build_index($fqfh,$idxfh);

  seek($fqfh,0,0);


  my %fqidx;
  my $line_no = 1;

  while ( my @seq = read_fastq($fqfh) ) {
    (my $qname) = $seq[0] =~ /^\s*(\S+)\s*/;
    croak "Error: multiple sequences with same ID in $fq" if defined $fqidx{$qname};
    $fqidx{$qname} = $line_no + 1;
    $line_no += 4;
  }

  return(\%fqidx);
}


1;
