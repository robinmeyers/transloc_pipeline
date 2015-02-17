use strict;
use warnings;
use POSIX qw(floor ceil);
use Math::BigFloat;




sub is_a_junction ($) {
  my $tlx = shift;
  if (! defined $tlx->{Rname} || $tlx->{Rname} eq "" || $tlx->{Rname} eq "Adapter") {
    return(0);
  } else {
    return(1);
  }
}

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

sub filter_unaligned ($) {
  my $read_obj = shift;

  my $tlxs = $read_obj->{tlxs};

  if (! defined $tlxs->[0]->{B_Rname}) {
    # my $junctions = filter_entire_read($tlxs,"unaligned");
    filter_entire_read($tlxs,"unaligned");
    # return(1,$junctions);
  } else {
    # return(0,0);
  }

}


sub filter_baitonly ($) {
  my $read_obj = shift;

  my $tlxs = $read_obj->{tlxs};

  if (defined $tlxs->[0]->{B_Rname} && ! is_a_junction($tlxs->[0])) {
    # my $junctions = filter_entire_read($tlxs,"baitonly");
    filter_entire_read($tlxs,"baitonly");
    # return(1,$junctions)
  } else {
    # return(0,0);
  }
}

sub filter_isjunction ($) {
  my $read_obj = shift;

  my $tlxs = $read_obj->{tlxs};

  foreach my $tlx (@$tlxs) {
    $tlx->{filters}->{isjunction} = 1 if is_a_junction($tlx);
  }
}

sub filter_uncut ($) {
  my $read_obj = shift;

  my $params = $main::params;

  my $tlxs = $read_obj->{tlxs};
  my $tlx = $tlxs->[0];

  if (defined $tlx->{B_Rname}) {
    if ($tlx->{B_Strand} == 1) {
      if ($tlx->{B_Rend} > $params->{brksite}->{breakcoord}) {
        # my $junctions = filter_entire_read($tlxs,"uncut");
        filter_entire_read($tlxs,"uncut",$tlx->{B_Rend} - $params->{brksite}->{breakcoord});
        # return(1,$junctions);
      }
    } else {
      if ($tlx->{B_Rstart} < $params->{brksite}->{breakcoord}) {

        # my $junctions = filter_entire_read($tlxs,"uncut");
        filter_entire_read($tlxs,"uncut",$params->{brksite}->{breakcoord} - $tlx->{B_Rstart});
        # return(1,$junctions);
      }
    }
  }

  # return(0,0);
}

sub filter_misprimed ($) {
  my $read_obj = shift;

  my $params = $main::params;

  my $tlxs = $read_obj->{tlxs};
  my $tlx = $tlxs->[0];

  if (defined $tlx->{B_Strand}) {
    if ($tlx->{B_Strand} == 1) {
      filter_entire_read($tlxs,"misprimed",
        $tlx->{B_Rend} - ($params->{brksite}->{start} + $params->{brksite}->{primer}->length - 1));
      # if ($tlx->{B_Rend} < $params->{brksite}->{misprimed_threshold}) {
      #   my $junctions = filter_entire_read($tlxs,"misprimed");
      #   return(1,$junctions);
      # }
    } else {
      filter_entire_read($tlxs,"misprimed",
        ($params->{brksite}->{end} - $params->{brksite}->{primer}->length) - $tlx->{B_Rstart});
      # if ($tlx->{B_Rstart} > $params->{brksite}->{misprimed_threshold}) {
      #   my $junctions = filter_entire_read($tlxs,"misprimed");
      #   return(1,$junctions);
      # }
    }
  }

  # return(0,0);
}


sub filter_freqcut ($) {
  my $read_obj = shift;
  
  my $params = $main::params;

  my $tlxs = $read_obj->{tlxs};

  my $i = 0;

  if ($params->{cutter} =~ /\S/) {
    foreach my $tlx (@$tlxs) {
      if (uc($tlx->{J_Seq}) =~ $params->{cutter} ||
          uc(substr($tlx->{Seq},0,$tlx->{Qstart}+4)) =~ $params->{cutter}) {
        filter_remainder_of_read($tlx,"freqcut",$i);
      }
      $i++;
    }
  }
}

sub filter_largegap ($) {
  my $read_obj = shift;

  my $tlxs = $read_obj->{tlxs};

  my $i = 0;


  foreach my $tlx (@$tlxs) {
    if (is_a_junction($tlx)) {
      $tlx->{filters}->{largegap} = max(0,$tlx->{Qstart} - $tlx->{B_Qend} - 1);
    }
    
  }

}


sub filter_mapqual ($) {

  my $read_obj = shift;

  my $params = $main::params;

  my $tlxs = $read_obj->{tlxs};

  my $R1_alns = $read_obj->{R1_alns};
  my $R2_alns = $read_obj->{R2_alns};

  my $i = 0;
  foreach my $tlx (@$tlxs) {
    unless (is_a_junction($tlx)) {
      $tlx->{filters}->{mapqual} = 255;
      next;
    }

    my $tlx_R1_aln = $R1_alns->{$tlx->{R1_ID}} if defined $tlx->{R1_ID};
    my $tlx_R2_aln = $R2_alns->{$tlx->{R2_ID}} if defined $tlx->{R2_ID};

    if (defined $tlx_R1_aln && defined $tlx_R2_aln) {

      $tlx_R1_aln->{Primary} = 1;
      $tlx_R1_aln->{Read} = "R1";
      $tlx_R1_aln->{Overlap} = "";

      $params->{mapqfh}->print(join("\t", @{$tlx_R1_aln}{@{$params->{mapq_header}}})."\n");
      
      foreach my $R1_aln_ID (keys $R1_alns) {
        next if $R1_aln_ID eq $tlx->{R1_ID};

        foreach my $R2_aln_ID (keys $R2_alns) {

          my $R1_aln = $R1_alns->{$R1_aln_ID};
          my $R2_aln = $R2_alns->{$R2_aln_ID};
          next if ($R1_aln->{Rname} eq "Adapter");
          next unless pair_is_proper($R1_aln,$R2_aln);

          $R1_aln->{Primary} = 0;
          $R1_aln->{Read} = "R1";

          $R1_aln->{Overlap} = calculate_fraction_overlap($tlx_R1_aln,$R1_aln);

          $params->{mapqfh}->print(join("\t", @{$R1_aln}{@{$params->{mapq_header}}})."\n");
          last;
        }
      }

      $tlx_R2_aln->{Read} = "R2";
      $tlx_R2_aln->{Primary} = 1;
      $tlx_R2_aln->{Overlap} = "";
      $params->{mapqfh}->print(join("\t", @{$tlx_R2_aln}{@{$params->{mapq_header}}})."\n");

      foreach my $R2_aln_ID (keys $R2_alns) {
        next if $R2_aln_ID eq $tlx->{R2_ID};

        foreach my $R1_aln_ID (keys $R1_alns) {

          my $R1_aln = $R1_alns->{$R1_aln_ID};
          my $R2_aln = $R2_alns->{$R2_aln_ID};
          next if ($R2_aln->{Rname} eq "Adapter");

          next unless pair_is_proper($R1_aln,$R2_aln);

          $R2_aln->{Primary} = 0;
          $R2_aln->{Read} = "R2";

          $R2_aln->{Overlap} = calculate_fraction_overlap($tlx_R2_aln,$R2_aln);

          $params->{mapqfh}->print(join("\t", @{$R2_aln}{@{$params->{mapq_header}}})."\n");
          last;
        }
      }



    } elsif (defined $tlx_R1_aln) {
      $tlx_R1_aln->{Primary} = 1;
      $tlx_R1_aln->{Read} = "R1";
      $tlx_R1_aln->{Overlap} = "";

      $params->{mapqfh}->print(join("\t", @{$tlx_R1_aln}{@{$params->{mapq_header}}})."\n");

      # only consider R1 alignments
      foreach my $R1_aln_ID (keys $R1_alns) {
        next if $R1_aln_ID eq $tlx->{R1_ID};

        my $R1_aln = $R1_alns->{$R1_aln_ID};
        next if ($R1_aln->{Rname} eq "Adapter");

        $R1_aln->{Primary} = 0;
        $R1_aln->{Read} = "R1";

        $R1_aln->{Overlap} = calculate_fraction_overlap($tlx_R1_aln,$R1_aln);

        $params->{mapqfh}->print(join("\t", @{$R1_aln}{@{$params->{mapq_header}}})."\n");
      }
    } else {
      $tlx_R2_aln->{Primary} = 1;
      $tlx_R2_aln->{Read} = "R2";
      $tlx_R2_aln->{Overlap} = "";

      $params->{mapqfh}->print(join("\t", @{$tlx_R2_aln}{@{$params->{mapq_header}}})."\n");

      # only consider R2 alignments
      foreach my $R2_aln_ID (keys $R2_alns) {
        next if $R2_aln_ID eq $tlx->{R2_ID};

        my $R2_aln = $R2_alns->{$R2_aln_ID};
        next if ($R2s_aln->{Rname} eq "Adapter");

        $R2_aln->{Primary} = 0;
        $R2_aln->{Read} = "R2";

        $R2_aln->{Overlap} = calculate_fraction_overlap($tlx_R2_aln,$R2_aln);

        $params->{mapqfh}->print(join("\t", @{$R2_aln}{@{$params->{mapq_header}}})."\n");
      }

    }

    # my $tlx_p = Math::BigFloat->new(10)->bpow(-$tlx_sum_base_Q/10);
    # my $competing_p = sum(map { Math::BigFloat->new(10)->bpow(-$_/10) } @competing_sum_base_Q);

    my $map_qual_score;
    # # is there a better way of doing this?
    # # if ($competing_p == 0) {
    #   # $map_qual_score = 0;
    # # } else {

    # my $mapq_incorrect_p = $tlx_p->copy()->bdiv($competing_p)->bmul(-1)->badd(1);
    # my $map_qual_score = $mapq_incorrect_p->is_zero() ?
    #                         Math::BigFloat->new(255) :
    #                         $mapq_incorrect_p->copy()->blog(10)->bmul(-10)->bfloor();
    #                       # floor(-10 * log10(1 - $tlx_p/$competing_p)) 
    # # }

    $map_qual_score = 255;

    # $map_qual_score = Math::BigFloat->new(255) if $map_qual_score->bcmp(255) > 0;

    $tlx->{filters}->{mapqual} = $map_qual_score;

    $i++;
  }

}

sub filter_breaksite ($) {
  my $read_obj = shift;
  
  my $tlxs = $read_obj->{tlxs};
  my $i = 0;
  foreach my $tlx (@$tlxs) {
    if (defined $tlx->{Rname} && $tlx->{Rname} eq "Breaksite") {
      $tlx->{filters}->{breaksite} = 1;
    }
    $i++;
  }
}

sub filter_sequential ($) {
  my $read_obj = shift;

  my $tlxs = $read_obj->{tlxs};

  if (defined $tlxs->[1] && is_a_junction($tlxs->[1])) {
    filter_remainder_of_read($tlxs,"sequential",1);
  }
}


1;