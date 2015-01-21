use strict;
use warnings;



sub filter_unjoined ($$) {
  my $tlxls = shift;
  my $brksite = shift;

  my $filter;

  my $junctions = 0;

  return 0 unless defined $tlxls->[0]->{R1_ID};

  if (@$tlxls < 2) {
    $filter = "Unjoined";
  } else {

    foreach my $tlxl (@$tlxls[1..$#{$tlxls}]) {
      last if $tlxl->{Rname} eq "Adapter";
      next if defined $tlxl->{R1_Rgap} && $tlxl->{R1_Rgap} >= 0 && $tlxl->{R1_Rgap} < 3;
      next if defined $tlxl->{R2_Rgap} && $tlxl->{R2_Rgap} >= 0 && $tlxl->{R2_Rgap} < 3;
      $junctions++;
    }

    $filter = "Unjoined" unless $junctions;

    if ($brksite->{aln_strand} == 1) {
      $filter = "Unjoined" if $tlxls->[1]->{tlx}->{B_Rend} > $brksite->{joining_threshold};
    } else {
      $filter = "Unjoined" if $tlxls->[1]->{tlx}->{B_Rstart} < $brksite->{joining_threshold};
    }
  }

  $junctions = 0;

  foreach my $tlxl (@$tlxls) {
    if (defined $tlxl->{tlx} && ! defined $tlxl->{tlx}->{Filter}) {
      if (defined $filter) {
        $tlxl->{tlx}->{Filter} = $filter;
      } else {
        $junctions++;
      }
    }
  }

  return $junctions;

}


sub filter_mapping_quality ($$$$$$$$$$){

  my $tlxls_ref = shift;
  my $R1_alns_ref = shift;
  my $R2_alns_ref = shift;
  my $ol_thresh = shift;
  my $mismatch_thresh_int = shift;
  my $mismatch_thresh_coef = shift;
  my $max_frag_length = shift;
  my $match_award = shift;
  my $mismatch_penalty = shift;
  my $mapqfh = shift;

  my @tlxls = @$tlxls_ref;
  my @R1_alns = @$R1_alns_ref;
  my @R2_alns = @$R2_alns_ref;


  my $filter;

  my $quality_maps = 0;

  return 0 unless defined $tlxls[0]->{R1_ID};


  TLXL: foreach my $tlxl (@tlxls) {
    if (defined $filter) {
      $tlxl->{tlx}->{Filter} = $filter if defined $tlxl->{tlx} && ! defined $tlxl->{tlx}->{Filter};
      next TLXL;
    }
    last TLXL if $tlxl->{Rname} eq "Adapter";
    my @R1_OL;
    my @R2_OL;
    if (defined $tlxl->{R1_ID}) {
      my $tlxl_R1_length = $tlxl->{R1_Qend} - $tlxl->{R1_Qstart} + 1;
      my $tlxl_R1_score = $tlxl->{R1_AS};
      my $score_difference_thresh = ($match_award + $mismatch_penalty) * 
                                    ($mismatch_thresh_int + $mismatch_thresh_coef * $tlxl_R1_length);
      # print(join("\t",$tlxl_R1_length,$tlxl_R1_score,$score_difference_thresh)."\n");
      foreach my $R1_aln (@R1_alns) {
        next unless defined $R1_aln->{ID} && $R1_aln->{ID} ne $tlxl->{R1_ID};
        next if $tlxl->{Rname} eq "Breaksite" && $R1_aln->{Rname} eq "Breaksite";
        my $overlap = min($tlxl->{R1_Qend},$R1_aln->{Qend}) - max($tlxl->{R1_Qstart},$R1_aln->{Qstart}) + 1;
        my $length = $R1_aln->{Qend} - $R1_aln->{Qstart} + 1;
        my $score = $R1_aln->{AS};

        if ($overlap >= $ol_thresh * ($tlxl->{R1_Qend} - $tlxl->{R1_Qstart} + 1) 
              && $score >= $tlxl_R1_score - $score_difference_thresh) {
          push (@R1_OL,$R1_aln);
        }
      }
    }

    if (defined $tlxl->{R2_ID}) {
      my $tlxl_R2_length = $tlxl->{R2_Qend} - $tlxl->{R2_Qstart} + 1;
      my $tlxl_R2_score = $tlxl->{R2_AS};
      my $score_difference_thresh = ($match_award + $mismatch_penalty) * 
                                    ($mismatch_thresh_int + $mismatch_thresh_coef * $tlxl_R2_length);
      foreach my $R2_aln (@R2_alns) {
        next unless defined $R2_aln->{ID} && $R2_aln->{ID} ne $tlxl->{R2_ID};
        next if $tlxl->{Rname} eq "Breaksite" && $R2_aln->{Rname} eq "Breaksite";
        my $overlap = min($tlxl->{R2_Qend},$R2_aln->{Qend}) - max($tlxl->{R2_Qstart},$R2_aln->{Qstart}) + 1;
        my $length = $R2_aln->{Qend} - $R2_aln->{Qstart} + 1;
        my $score = $R2_aln->{AS};


        if ($overlap >= $ol_thresh * ($tlxl->{R2_Qend} - $tlxl->{R2_Qstart} + 1) 
              && $score >= $tlxl_R2_score - $score_difference_thresh) {
          push (@R2_OL,$R2_aln);
        }
      }
    }

    if (defined $tlxl->{R1_ID} && $tlxl->{R2_ID}) {
      ALN_PAIR: foreach my $R1_aln (@R1_OL) {
        foreach my $R2_aln (@R2_OL) {
          if (pair_is_proper($R1_aln,$R2_aln,$max_frag_length)) {
            $filter = "MappingQuality";
            $mapqfh->print(join("\t",$tlxl->{Qname},
                                     $tlxl->{Rname},
                                     $tlxl->{R1_Rstart},
                                     $tlxl->{R1_Rend},
                                     $tlxl->{Strand},
                                     $tlxl->{R1_Qstart},
                                     $tlxl->{R1_Qend},
                                     $tlxl->{R1_AS},
                                     $tlxl->{R1_Cigar},
                                     $tlxl->{Rname},
                                     $tlxl->{R2_Rstart},
                                     $tlxl->{R2_Rend},
                                     $tlxl->{Strand},
                                     $tlxl->{R2_Qstart},
                                     $tlxl->{R2_Qend},
                                     $tlxl->{R2_AS},
                                     $tlxl->{R2_Cigar})."\n");
            $mapqfh->print(join("\t",$R1_aln->{Qname},
                                     $R1_aln->{Rname},
                                     $R1_aln->{Rstart},
                                     $R1_aln->{Rend},
                                     $R1_aln->{Strand},
                                     $R1_aln->{Qstart},
                                     $R1_aln->{Qend},
                                     $R1_aln->{AS},
                                     $R1_aln->{Cigar},
                                     $R2_aln->{Rname},
                                     $R2_aln->{Rstart},
                                     $R2_aln->{Rend},
                                     $R2_aln->{Strand},
                                     $R2_aln->{Qstart},
                                     $R2_aln->{Qend},
                                     $R2_aln->{AS},
                                     $R2_aln->{Cigar})."\n");
            last ALN_PAIR;
          }
        }
      }
    } else {
      if (scalar @R1_OL > 0 || scalar @R2_OL > 0) {
        $filter = "MappingQuality";
        if (scalar @R1_OL > 0) {
          $mapqfh->print(join("\t",$tlxl->{Qname},
                                       $tlxl->{Rname},
                                       $tlxl->{R1_Rstart},
                                       $tlxl->{R1_Rend},
                                       $tlxl->{Strand},
                                       $tlxl->{R1_Qstart},
                                       $tlxl->{R1_Qend},
                                       $tlxl->{R1_AS},
                                       $tlxl->{R1_Cigar})."\n");
          foreach my $aln (@R1_OL) {
            $mapqfh->print(join("\t", $aln->{Qname},
                                      $aln->{Rname},
                                      $aln->{Rstart},
                                      $aln->{Rend},
                                      $aln->{Strand},
                                      $aln->{Qstart},
                                      $aln->{Qend},
                                      $aln->{AS},
                                      $aln->{Cigar})."\n");
          }
        } else {
          $mapqfh->print(join("\t",$tlxl->{Qname},
                                       "",
                                       "",
                                       "",
                                       "",
                                       "",
                                       "",
                                       "",
                                       "",
                                       $tlxl->{Rname},
                                       $tlxl->{R2_Rstart},
                                       $tlxl->{R2_Rend},
                                       $tlxl->{Strand},
                                       $tlxl->{R2_Qstart},
                                       $tlxl->{R2_Qend},
                                       $tlxl->{R2_AS},
                                       $tlxl->{R2_Cigar})."\n");
          foreach my $aln (@R2_OL) {
            $mapqfh->print(join("\t", $aln->{Qname},
                                      "",
                                      "",
                                      "",
                                      "",
                                      "",
                                      "",
                                      "",
                                      "",
                                      $aln->{Rname},
                                      $aln->{Rstart},
                                      $aln->{Rend},
                                      $aln->{Strand},
                                      $aln->{Qstart},
                                      $aln->{Qend},
                                      $aln->{AS},
                                      $aln->{Cigar})."\n");
          }
        }
      }
    }

    if (defined $tlxl->{tlx} && ! defined $tlxl->{tlx}->{Filter}) {
      if (defined $filter) {
        $tlxl->{tlx}->{Filter} = $filter;
      } else {
        $quality_maps++;
      }
    }
  }

  return $quality_maps;

}


sub filter_mispriming ($$) {
  my $tlxls = shift;
  my $brksite = shift;

  my $filter;
  my $priming = 0;

  return 0 unless defined $tlxls->[0]->{R1_ID};

  $filter = "Mispriming" if $brksite->{aln_strand} == 1 && $tlxls->[0]->{R1_Rend} < $brksite->{priming_threshold};
  $filter = "Mispriming" if $brksite->{aln_strand} == -1 && $tlxls->[0]->{R1_Rstart} > $brksite->{priming_threshold};

  foreach my $tlxl (@$tlxls) {
    if (defined $tlxl->{tlx} && ! defined $tlxl->{tlx}->{Filter}) {
      if (defined $filter) {
        $tlxl->{tlx}->{Filter} = $filter;
      } else {
        $priming++;
      }
    }
  }

  return $priming;

}

sub filter_freq_cutter ($$) {
  my $tlxls = shift;
  my $cutter = shift;

  my $filter;

  my $no_cutter = 0;

  return 0 unless defined $tlxls->[0]->{R1_ID};

  foreach my $tlxl (@$tlxls) {
    if (defined $tlxl->{tlx} && ! defined $tlxl->{tlx}->{Filter}) {
      if (defined $filter) {
        $tlxl->{tlx}->{Filter} = $filter;
        next;
      }


      if (defined $cutter && $cutter->seq =~ /\S/) {
        if (uc($tlxl->{tlx}->{J_Seq}) =~ $cutter->seq || substr($tlxl->{tlx}->{Seq},0,$tlxl->{tlx}->{Qstart}+4) =~ $cutter->seq) {
          $filter = "FreqCutter";
          $tlxl->{tlx}->{Filter} = $filter;
          next;
        }
      }

      $no_cutter++;
    }
  }

  return $no_cutter;
}

sub filter_breaksite ($) {
  my $tlxls = shift;

  my $filter;

  my $outside_breaksite = 0;

  return 0 unless defined $tlxls->[0]->{R1_ID};
  return 0 unless defined $tlxls->[1];

  foreach my $i (0..$#{$tlxls}) {

    my $tlxl = $tlxls->[$i];

    if (defined $tlxl->{tlx} && ! defined $tlxl->{tlx}->{Filter}) {
      if (defined $filter) {
        $tlxl->{tlx}->{Filter} = $filter;
        next;
      }



      if (($i == 1 && $tlxl->{tlx}->{Rname} eq "Breaksite") ||
          ($i == 1 && defined $tlxl->{R1_Rgap} && $tlxl->{R1_Rgap} >=0 && $tlxl->{R1_Rgap} < 10) ||
          ($i == 1 && defined $tlxl->{R2_Rgap} && $tlxl->{R2_Rgap} >=0 && $tlxl->{R2_Rgap} < 10)) {
        $filter = "Breaksite";
        $tlxl->{tlx}->{Filter} = $filter;
        next;
      }

      $outside_breaksite++;
    }
  }

  return $outside_breaksite;
}

sub filter_sequential_junctions ($) {
  my $tlxls = shift;

  my $filter;

  my $primary_junction = 0;

  return 0 unless defined $tlxls->[0]->{R1_ID};

  foreach my $tlxl (@$tlxls) {

    if (defined $tlxl->{tlx} && ! defined $tlxl->{tlx}->{Filter}) {
      if (defined $filter) {
        $tlxl->{tlx}->{Filter} = $filter;
        next;
      }

      $primary_junction++;
      $filter = "SequentialJunction";


    }
  }

  return $primary_junction;



}


1;