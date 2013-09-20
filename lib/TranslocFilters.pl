use strict;
use warnings;



sub filter_unjoined ($) {
  my $tlxls = shift;

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


sub filter_mapping_quality ($$$$$$){

  my $tlxls_ref = shift;
  my $R1_alns_ref = shift;
  my $R2_alns_ref = shift;
  my $ol_thresh = shift;
  my $score_thresh = shift;
  my $max_frag_length = shift;

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
      my $R1_AS = $tlxl->{R1_AS};
      foreach my $R1_aln (@R1_alns) {
        next unless defined $R1_aln->{ID} && $R1_aln->{ID} ne $tlxl->{R1_ID};
        next if $tlxl->{Rname} eq "Breaksite" && $R1_aln->{Rname} eq "Breaksite";
        my $overlap = min($tlxl->{R1_Qend},$R1_aln->{Qend}) - max($tlxl->{R1_Qstart},$R1_aln->{Qstart}) + 1;
        my $score = $R1_aln->{AS};

        if ($overlap > $ol_thresh * ($tlxl->{R1_Qend} - $tlxl->{R1_Qstart} + 1) 
              && $score > $score_thresh * $R1_AS) {
          push (@R1_OL,$R1_aln);
        }
      }
    }

    if (defined $tlxl->{R2_ID}) {
      my $R2_AS = $tlxl->{R2_AS};
      foreach my $R2_aln (@R2_alns) {
        next unless defined $R2_aln->{ID} && $R2_aln->{ID} ne $tlxl->{R2_ID};
        next if $tlxl->{Rname} eq "Breaksite" && $R2_aln->{Rname} eq "Breaksite";
        my $overlap = min($tlxl->{R2_Qend},$R2_aln->{Qend}) - max($tlxl->{R2_Qstart},$R2_aln->{Qstart}) + 1;
        my $score = $R2_aln->{AS};
        if ($overlap > $ol_thresh * ($tlxl->{R2_Qend} - $tlxl->{R2_Qstart} + 1) 
              && $score > $score_thresh * $R2_AS) {
          push (@R2_OL,$R2_aln);
        }
      }
    }

    if (defined $tlxl->{R1_Rname} && $tlxl->{R2_Rname}) {
      ALN_PAIR: foreach my $R1_aln (@R1_OL) {
        foreach my $R2_aln (@R2_OL) {
          if (pair_is_proper($R1_aln,$R2_aln,$max_frag_length)) {
            $filter = "MappingQuality";
            last ALN_PAIR;
          }
        }
      }
    } else {
      if (scalar @R1_OL > 0 || scalar @R2_OL > 0) {
        $filter = "MappingQuality";
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
  my $priming_threshold = shift;

  my $filter;
  my $priming = 0;

  return 0 unless defined $tlxls->[0]->{R1_ID};

  $filter = "Mispriming" unless $tlxls->[0]->{Rname} eq "Breaksite";
  $filter = "Mispriming" if $tlxls->[0]->{R1_Rend} < $priming_threshold;

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

      if ($cutter =~ /\S/) {
        if (uc($tlxl->{tlx}->{JuncSeq}) =~ $cutter || substr($tlxl->{tlx}->{Qseq},0,$tlxl->{tlx}->{Qstart}-1) =~ $cutter) {
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



      if ($i = 1 && $tlxl->{tlx}->{Rname} eq "Breaksite") {
        $filter = "Breaksite";
        $tlxl->{tlx}->{Filter} = $filter;
        next;
      }

      $outside_breaksite++;
    }
  }

  return $outside_breaksite;
}

sub filter_split_junctions ($) {
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
      $filter = "SplitJunction";


    }
  }

  return $primary_junction;



}


1;