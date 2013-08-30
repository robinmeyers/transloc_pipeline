use strict;
use warnings;
use Switch;
use Time::HiRes qw(time);
use List::Util qw(min max);
use List::MoreUtils qw(firstidx);
use Data::GUID;

sub prepare_reference_genomes ($) {


  my $meta_ref = shift;
  my @masks = ();

  my $GENOME_DB = $ENV{'GENOME_DB'};
  my $BOWTIE2_INDEXES = $ENV{'BOWTIE2_INDEXES'};

  foreach my $expt_id (sort keys %$meta_ref) {
    my $assembly = $meta_ref->{$expt_id}->{assembly};
    my $assemdir = "$GENOME_DB/$assembly";
    croak "Error: could not find reference genome directory $GENOME_DB/$assembly" unless (-d $assemdir);

    my $cleanFa = "$assemdir/$assembly.fa";

    my $mask = $meta_ref->{$expt_id}->{mask};
    next unless ($mask =~ /\S/);
    (my $mask_stub = $assembly."_mask_".$mask) =~ s/[:,\s]+/_/g;
    my $maskdir = "$GENOME_DB/$mask_stub";
    unless (-d $maskdir) {
      mkdir $maskdir or croak "Error: could not create directory for masked genome";
    }
    my $maskFa = "$maskdir/$mask_stub.fa";
    $meta_ref->{$expt_id}->{mask_assembly} = $mask_stub;

    my $bedfile = "$maskdir/$mask_stub.bed"; 

      
    my $bed_fh = IO::File->new(">$bedfile");
    my @loci = split( /\s*,\s*/ , $mask );
    foreach my $locus (@loci) {
      (my $chr, my $start, my $end) = ($locus =~ /(chr\w+):(\d+)-(\d+)/);
      $bed_fh->print(join("\t",$chr,$start,$end)."\n");
    }
    $bed_fh->close;

    my $maskFaSize = -s $maskFa;
    my $cleanFaSize = -s $cleanFa;

    unless (-r $maskFa && $maskFaSize > 0.99 * $cleanFaSize) {
      System("maskFastaFromBed -fi $cleanFa -fo $maskFa -bed $bedfile");
    }

    unless (-r "$BOWTIE2_INDEXES/$mask_stub.1.bt2" && -r "$BOWTIE2_INDEXES/$mask_stub.rev.2.bt2") {     
      System("bowtie2-build $maskFa $BOWTIE2_INDEXES/$mask_stub");
    }
  }

}

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

sub tlx_header {
  my @tlx_header = ( qw(Qname Rname Junction Strand Rstart Rend),
                      qw(B_Rname B_Rstart B_Rend B_Strand B_Qstart B_Qend),
                      qw(Qstart Qend Qlen Seq JuncSeq) );
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
                      qw(R1_Rgap R2_Rgap R1_Cigar R2_Cigar) );
  return (@tlxl_header);
}

sub write_entry ($$$) {

  my $fh = shift;
  my $entry = shift;
  my $header = shift;


  $fh->print(join("\t",map(check_undef($_,""),@{$entry}{@$header}))."\n");
  
}


sub print_aln ($) {
  my $aln = shift;
  print $aln->{Qstart}."-".$aln->{Qend}." ".$aln->{Rname}.":".$aln->{Rstart}."-".$aln->{Rend}.":".$aln->{Strand}."\n";
}

sub find_genomic_distance ($$$) {
  my $aln1 = shift;
  my $aln2 = shift;
  my $brk = shift;

  my $chr1 = $aln1->{Rname};
  my $chr2 = $aln2->{Rname};
  my $strand1 = $aln1->{Strand};
  my $strand2 = $aln2->{Strand};
  my $junc1 = $strand1 == 1 ? $aln1->{Rend} : $aln1->{Rstart};
  my $junc2 = $strand2 == 1 ? $aln2->{Rstart} : $aln2->{Rend};


  my $qend1 = $aln1->{Qend};
  my $qstart2 = $aln2->{Qstart};

  my $g_dist;
  my $q_dist;
  my $r_dist;
  my $g_dist_to_start_of_brk;
  my $g_dist_to_end_of_brk;


  if ($aln1->{Rname} eq $aln2->{Rname}) {

    $q_dist = $strand1 == $strand2 ? $qstart2 - $qend1 - 1 : 0;
    $r_dist = $junc2 - $junc1 - 1;
    $g_dist = $strand1 == 1 ? $r_dist - $q_dist : $r_dist + $q_dist;
  
  } elsif ($chr1 eq "Breaksite" && $chr2 eq $brk->{chr}) {

    if ($brk->{strand} eq "+") {

      if ($junc2 >= $brk->{end}) {
        $junc1 = $brk->{end} - 1 - ($brk->{len} - $junc1);
      } elsif ($junc2 < $brk->{start}) {
        $junc1 = $brk->{start} + $junc1;
      } elsif ($strand2 == 1) {
        $junc1 = $brk->{end} - 1 - ($brk->{len} - $junc1);
      } else {
        $junc1 = $brk->{start} + $junc1;
      }

    } else {

      if ($junc2 >= $brk->{end}) {
        $junc1 = $brk->{end} - $junc1;
      } elsif ($junc2 < $brk->{start}) {
        $junc1 = $brk->{start} + ($brk->{len} - $junc1);
      } elsif ($strand2 == 1) {
        $junc1 = $brk->{end} - $junc1;
      } else {
        $junc1 = $brk->{start} + ($brk->{len} - $junc1);
      }

      $strand1 = -1 * $strand1;

    }

    $q_dist = $strand1 == $strand2 ? $qstart2 - $qend1 - 1 : 0;
    $r_dist = $junc2 - $junc1;
    $g_dist = $strand1 == 1 ? $r_dist - $q_dist : $r_dist + $q_dist;

  } elsif ($chr2 eq "Breaksite" && $chr1 eq $brk->{chr}) {

    if ($brk->{strand} eq "+") {

      if ($junc1 >= $brk->{end}) {
        $junc2 = $brk->{end} - 1 - ($brk->{len} - $junc2);
      } elsif ($junc1 < $brk->{start}) {
        $junc2 = $brk->{start} + $junc2;
      } elsif ($strand1 == 1) {
        $junc2 = $brk->{end} - 1 - ($brk->{len} - $junc2);
      } else {
        $junc2 = $brk->{start} + $junc2;
      }

    } else {

      if ($junc1 >= $brk->{end}) {
        $junc2 = $brk->{end} - $junc2;
      } elsif ($junc1 < $brk->{start}) {
        $junc2 = $brk->{start} + ($brk->{len} - $junc2);
      } elsif ($strand1 == 1) {
        $junc2 = $brk->{end} - $junc2;
      } else {
        $junc2 = $brk->{start} + ($brk->{len} - $junc2);
      }

      $strand2 = -1 * $strand2;

    }

    $q_dist = $strand1 == $strand2 ? $qstart2 - $qend1 - 1 : 0;
    $r_dist = $junc2 - $junc1;
    $g_dist = $strand1 == 1 ? $r_dist - $q_dist : $r_dist + $q_dist;

  } else {
    return undef;
  }

  $g_dist = max(1,abs($g_dist));

  $g_dist = $strand1 == $strand2 ? $g_dist : -$g_dist;

  return $g_dist;
}

sub wrap_alignment ($$) {

  my $pe = shift;

  my $aln = shift;

  return undef unless defined $aln;
  croak "Error: first argument to wrap_alignment must be either \"R1\" or \"R2\"" unless ($pe eq "R1" || $pe eq "R2");
  # my %wrapper :shared;
  my $wrapper = {Qname => $aln->qname};
  if ($aln->unmapped) {
    $wrapper->{Unmapped} = 1;
    $wrapper->{Seq} = $pe eq "R1" ? $aln->query->dna : reverseComplement($aln->query->dna);
    $wrapper->{Qual} = $pe eq "R1" ? $aln->query->qscore : [reverse @{$aln->query->qscore}];
    $wrapper->{Qstart} = 0;
    return $wrapper;
  }

  $wrapper->{Rname} = $aln->seq_id;
  $wrapper->{Rstart} = $aln->start;
  $wrapper->{Rend} = $aln->end;
  $wrapper->{AS} = $aln->aux =~ /AS:i:(\d+)/ ? $1 : 0;
  $wrapper->{Unmapped} = 0;

  if ($pe eq "R1") {
    $wrapper->{Strand} = $aln->strand;
  } else {
    $wrapper->{Strand} = -1 * $aln->strand;
  }


  if ($wrapper->{Strand} == 1) {
    $wrapper->{Seq} = $aln->query->dna;
    $wrapper->{Qual} = $aln->query->qscore;
    $wrapper->{Qstart} = $aln->query->start;
    $wrapper->{Qend} = $aln->query->end;
    $wrapper->{CigarA} = $aln->cigar_array;
  } else {
    $wrapper->{Seq} = reverseComplement($aln->query->dna);
    $wrapper->{Qual} = [reverse @{$aln->query->qscore}];
    $wrapper->{Qstart} = $aln->l_qseq - $aln->query->end + 1;
    $wrapper->{Qend} = $aln->l_qseq - $aln->query->start + 1;
    $wrapper->{CigarA} = [reverse @{$aln->cigar_array}];
  }

  $wrapper->{Cigar} = join("",map {join("",@$_)} @{$wrapper->{CigarA}});
  $wrapper->{Qlen} = length($wrapper->{Seq});
  $wrapper->{ID} = Data::GUID->new->as_string;


  return $wrapper;

}

sub pair_is_proper ($$$) {
  my $R1 = shift;
  my $R2 = shift;
  my $max_frag_len = shift;
  my $max_dovetail = 10;

  return 0 unless $R1->{Rname} eq $R2->{Rname};
  return 0 unless $R1->{Strand} == $R2->{Strand};


  if ($R1->{Strand} == 1) {
    return 0 unless $R1->{Rstart} < $R2->{Rend};
    return 0 if $R2->{Rend} - $R1->{Rstart} + 1 > $max_frag_len;
    return 0 if $R1->{Rstart} > $R2->{Rstart} + $max_dovetail;
    return 0 if $R1->{Rend} > $R2->{Rend} + $max_dovetail;
  } else {
    return 0 unless $R2->{Rstart} < $R1->{Rend};
    return 0 if $R1->{Rend} - $R2->{Rstart} + 1 > $max_frag_len;
    return 0 if $R2->{Rstart} > $R1->{Rstart} + $max_dovetail;
    return 0 if $R2->{Rend} > $R1->{Rend} + $max_dovetail;
  }

  return 1;

}

sub create_tlxl_entries ($) {

  my $OCS_ref = shift;
  my @OCS = @$OCS_ref;

  my @tlxls = ();


  if ( $OCS[0]->{R1}->{Unmapped} ) {
    my $tlxl = { Qname => $OCS[0]->{R1}->{Qname},
                 R1_Seq => $OCS[0]->{R1}->{Seq},
                 R2_Seq => $OCS[0]->{R2}->{Seq},
                 Unmapped => 1 };
    return [$tlxl];
  }


  foreach my $Qseg (@OCS[0..$#OCS]) {
    my $tlxl = {};

    if (defined $Qseg->{R1}) {
      $tlxl->{Qname} = $Qseg->{R1}->{Qname};
      $tlxl->{OCS_score} = $Qseg->{score};
      $tlxl->{R1_Qstart} = $Qseg->{R1}->{Qstart};
      $tlxl->{R1_Qend} = $Qseg->{R1}->{Qend};
      $tlxl->{R1_Qlen} = $Qseg->{R1}->{Qlen};
      $tlxl->{Rname} = $Qseg->{R1}->{Rname};
      $tlxl->{Strand} = $Qseg->{R1}->{Strand};
      $tlxl->{R1_Rstart} = $Qseg->{R1}->{Rstart};
      $tlxl->{R1_Rend} = $Qseg->{R1}->{Rend};
      $tlxl->{R1_Rgap} = $Qseg->{R1_Rgap};
      $tlxl->{R1_Cigar} = $Qseg->{R1}->{Cigar};
      $tlxl->{R1_CigarA} = $Qseg->{R1}->{CigarA};
      $tlxl->{R1_Seq} = $Qseg->{R1}->{Seq};
      $tlxl->{R1_Qual} = $Qseg->{R1}->{Qual};
      $tlxl->{R1_ID} = $Qseg->{R1}->{ID};
      $tlxl->{R1_AS} = $Qseg->{R1}->{AS};
    }

    if (defined $Qseg->{R2}) {

      unless (defined $Qseg->{R1}) {
        $tlxl->{Qname} = $Qseg->{R2}->{Qname};
        $tlxl->{OCS_score} = $Qseg->{score};
        $tlxl->{Rname} = $Qseg->{R2}->{Rname};
        $tlxl->{Strand} = $Qseg->{R2}->{Strand};        
      }
      $tlxl->{R2_Qstart} = $Qseg->{R2}->{Qstart};
      $tlxl->{R2_Qend} = $Qseg->{R2}->{Qend};
      $tlxl->{R2_Qlen} = $Qseg->{R2}->{Qlen};
      $tlxl->{R2_Rstart} = $Qseg->{R2}->{Rstart};
      $tlxl->{R2_Rend} = $Qseg->{R2}->{Rend};
      $tlxl->{R2_Rgap} = $Qseg->{R2_Rgap};
      $tlxl->{R2_Cigar} = $Qseg->{R2}->{Cigar};
      $tlxl->{R2_CigarA} = $Qseg->{R2}->{CigarA};
      $tlxl->{R2_Seq} = $Qseg->{R2}->{Seq};
      $tlxl->{R2_Qual} = $Qseg->{R2}->{Qual};
      $tlxl->{R2_ID} = $Qseg->{R2}->{ID};
      $tlxl->{R2_AS} = $Qseg->{R2}->{AS};

    }

    push(@tlxls,$tlxl);

  }

  return \@tlxls;

}

sub create_tlx_entries ($$) {

  my $tlxls = shift;
  my $refs = shift;
  my $mar = 2;

  if ($tlxls->[0]->{Unmapped}) {
    my $tlx = {};
    $tlx->{Qname} = $tlxls->[0]->{Qname};
    $tlx->{Filter} = "Unaligned";
    $tlxls->[0]->{tlx} = $tlx;
    return;
  }

  my @Qseq = ();
  my @R1_map = ();
  my @R2_map = ();

  my $R1_Qpos = 0;
  my $R2_Qpos = 0;

  my @R1_Qseq = split("",$tlxls->[0]->{R1_Seq});
  my $R2_idx = firstidx {defined $_->{R2_Seq}} @$tlxls;
  my @R2_Qseq = split("",$tlxls->[$R2_idx]->{R2_Seq}) if $R2_idx >= 0;

  foreach my $i (0..$#{$tlxls}) {

    my $tlxl = $tlxls->[$i];

    # print "\nsegment $i\n";

    if ((defined $tlxls->[$i+1] && defined $tlxls->[$i+1]->{R1_ID}) || ! defined $tlxl->{R2_ID}) {
      # print "R1 only\n";
      if ($R1_Qpos < $tlxl->{R1_Qend}) {
        my $old_Qpos = scalar @Qseq;
        push(@Qseq,@R1_Qseq[$R1_Qpos..($tlxl->{R1_Qend}-1)]);
        my $new_Qpos = scalar @Qseq;
        @R1_map[$R1_Qpos..($tlxl->{R1_Qend}-1)] = (($old_Qpos+1)..$new_Qpos);
        $R1_Qpos = $tlxl->{R1_Qend};
      }

    } else {
      # print "merging R1 and R2\n";
    
      my $Rname = $tlxl->{Rname};
      my $Strand = $tlxl->{Strand};

      my @R1_Qual = @{$tlxl->{R1_Qual}};
      @R1_Qual = reverse @R1_Qual unless $Strand == 1;
      my @R2_Qual = @{$tlxl->{R2_Qual}};
      @R2_Qual = reverse @R2_Qual if $Strand == 1;
      my @R1_cigar = ();
      foreach my $i (@{$tlxl->{R1_CigarA}}) {
        next if $i->[0] eq "S";
        push(@R1_cigar,split("",($i->[0] x $i->[1])));
      }
      @R1_cigar = reverse @R1_cigar unless $Strand == 1;
      my @R2_cigar = ();
      foreach my $i (@{$tlxl->{R2_CigarA}}) {
        next if $i->[0] eq "S";
        push(@R2_cigar,split("",($i->[0] x $i->[1])));
      }
      @R2_cigar = reverse @R2_cigar unless $Strand == 1;

      my $R1_Rstart = $tlxl->{R1_Rstart};
      my $R1_Rend = $tlxl->{R1_Rend};
      my $R2_Rstart = $tlxl->{R2_Rstart};
      my $R2_Rend = $tlxl->{R2_Rend};

      my $Rstart = $tlxl->{Strand} == 1 ? $tlxl->{R1_Rstart} : $tlxl->{R2_Rstart};
      my $Rend = $tlxl->{Strand} == 1 ? $tlxl->{R2_Rend} : $tlxl->{R1_Rend};
      my $ref;
      switch ($Rname) {
        case "Breaksite" { $ref = $refs->{brk}; }
        case "Adapter" { $ref = $refs->{adpt}; }
        else { $ref = $refs->{genome}; }
      }

      # Retrieve referense sequence, reverse complement if necessary
      my $Rseq = $ref->seq($Rname,$Rstart,$Rend);
      my @Rseq = $Strand == 1 ? split("",$Rseq) : split("",reverseComplement($Rseq));

      # Artifically adjust Rstarts and Rends so that we can move incrementally up on both Qpos and Rpos
      ($R1_Rstart,$R1_Rend) = $Strand == 1 ? ($R1_Rstart,$R1_Rend) : (-$R1_Rend,-$R1_Rstart);
      ($R2_Rstart,$R2_Rend) = $Strand == 1 ? ($R2_Rstart,$R2_Rend) : (-$R2_Rend,-$R2_Rstart);

      # print "R1: $R1_Rstart-$R1_Rend\n";
      # print join(" ",@R1_cigar)."\n";
      # print "R2: $R2_Rstart-$R2_Rend\n";
      # print join(" ",@R2_cigar)."\n";

      $R2_Qpos = $tlxl->{R2_Qstart} - 1;
      my $Rpos = $R1_Rstart;

      # Push forward on R1 if behind start of alignment (gap between this and last alignment)
      if ($R1_Qpos < $tlxl->{R1_Qstart} - 1) {
        my $old_Qpos = scalar @Qseq;
        push(@Qseq,@R1_Qseq[$R1_Qpos..($tlxl->{R1_Qstart}-2)]);
        my $new_Qpos = scalar @Qseq;
        @R1_map[$R1_Qpos..($tlxl->{R1_Qstart}-2)] = (($old_Qpos+1)..$new_Qpos);
        $R1_Qpos = $tlxl->{R1_Qstart} - 1;
      }

      # Push forward on R1 alignment if start is behind current R1_Qpos (overlap between this and last alignment)
      if ($tlxl->{R1_Qstart} - 1 < $R1_Qpos) {
        my $pos1_tmp = $tlxl->{R1_Qstart} - 1;
        while ($pos1_tmp < $R1_Qpos) {
          my $c1 = shift @R1_cigar;
          switch ($c1) {
            case "M" { $pos1_tmp++; $Rpos++; }
            case "D" { $Rpos++; }
            case "I" { $pos1_tmp++; }
            else { $pos1_tmp++; $Rpos++; }
          }
        }
      }

      # Push forward on R2 if Rstart is behind current Rpos
      if ($R2_Rstart < $Rpos) {
        my $Rpos_tmp = $R2_Rstart;
        while ($Rpos_tmp < $Rpos) {
          my $c2 = shift @R2_cigar;
          switch ($c2) {
            case "M" { $R2_Qpos++; $Rpos_tmp++; }
            case "D" { $Rpos_tmp++; }
            case "I" { $R2_Qpos++; }
            else { $R2_Qpos++; $Rpos_tmp++; }
          }
        }
      }

      # Add R1 alignment up to the start of R2 alignment
      while ($Rpos <= $R1_Rend && $Rpos < $R2_Rstart) {
        my $c1 = shift @R1_cigar;
        switch ($c1) {
          case "M" { 
            push(@Qseq,$R1_Qseq[$R1_Qpos]);
            $R1_map[$R1_Qpos] = scalar @Qseq;
            $R1_Qpos++;
            $Rpos++; }
          case "D" { $Rpos++; }
          case "I" {
            push(@Qseq,$R1_Qseq[$R1_Qpos]);
            $R1_map[$R1_Qpos] = scalar @Qseq;
            $R1_Qpos++; }
        }
      }

      # Fill in gap between alignment if exists
      while ($Rpos > $R1_Rend && $Rpos <= $R2_Rstart) {
        push(@Qseq,lc($Rseq[$Rpos - $R1_Rstart]));
        $Rpos++;
      }

      # Negotiate the R1 and R2 alignments
      while ($Rpos <= $R1_Rend && $Rpos >= $R2_Rstart && $Rpos <= $R2_Rend) {
        my $c1 = shift @R1_cigar;
        my $c2 = shift @R2_cigar;
        my $q1 = $R1_Qseq[$R1_Qpos];
        my $q2 = $R2_Qseq[$R2_Qpos];
        my $r = $Rseq[$Rpos - $R1_Rstart];
        # print "R1: $R1_Qpos $c1 / R2: $R2_Qpos $c2\n";

        switch ($c1) {
          case "M" {
            switch ($c2) {
              case "M" {
                if ($q1 eq $r || $q2 eq $r) {
                  push(@Qseq,$r);
                } elsif ($q1 eq $q2) {
                  push(@Qseq,$q1);
                } elsif (nearby_quals(\@R1_Qual,$R1_Qpos,$mar) >= nearby_quals(\@R2_Qual,$R2_Qpos,$mar)) {
                  push(@Qseq,$q1);
                } else {
                  push(@Qseq,$q2);
                }
                $R1_map[$R1_Qpos] = scalar @Qseq;
                $R2_map[$R2_Qpos] = scalar @Qseq;
                $R1_Qpos++;
                $R2_Qpos++;
                $Rpos++;
              }
              case "D" {
                if ($q1 eq $r || nearby_quals(\@R1_Qual,$R1_Qpos,$mar) >= nearby_quals(\@R2_Qual,$R2_Qpos,$mar)) {
                  push(@Qseq,$q1);
                }
                $R1_map[$R1_Qpos] = scalar @Qseq;
                $R1_Qpos++;
                $Rpos++;
              }
              case "I" {
                if ($q1 ne $r && nearby_quals(\@R1_Qual,$R1_Qpos,$mar) < nearby_quals(\@R2_Qual,$R2_Qpos,$mar)) {
                  push(@Qseq,$q2);
                }
                $R2_map[$R2_Qpos] = scalar @Qseq;
                $R2_Qpos++;
                unshift(@R1_cigar,$c1);
              }
            }
          }
          case "I" {
            switch ($c2) {
              case "M" {
                if ($q2 ne $r && nearby_quals(\@R1_Qual,$R1_Qpos,$mar) >= nearby_quals(\@R2_Qual,$R2_Qpos,$mar)) {
                  push(@Qseq,$q1);
                }
                $R1_map[$R1_Qpos] = scalar @Qseq;
                $R1_Qpos++;
                unshift(@R2_cigar,$c2);
              }
              case "D" {
                if (nearby_quals(\@R1_Qual,$R1_Qpos,$mar) >= nearby_quals(\@R2_Qual,$R2_Qpos,$mar)) {
                  push(@Qseq,$q1);
                }
                $R1_map[$R1_Qpos] = scalar @Qseq;
                $R1_Qpos++;
                unshift(@R2_cigar,$c2);
              }
              case "I" {
                if ($q1 eq $q2 || nearby_quals(\@R1_Qual,$R1_Qpos,$mar) >= nearby_quals(\@R2_Qual,$R2_Qpos,$mar)) {
                  push(@Qseq,$q1);
                } else {
                  push(@Qseq,$q2);
                }
                $R1_map[$R1_Qpos] = scalar @Qseq;
                $R2_map[$R2_Qpos] = scalar @Qseq;
                $R1_Qpos++;
                $R2_Qpos++;
              }
            }
          }
          case "D" {
            switch ($c2) {
              case "M" {
                if ($q2 eq $r || nearby_quals(\@R1_Qual,$R1_Qpos,$mar) < nearby_quals(\@R2_Qual,$R2_Qpos,$mar)) {
                  push(@Qseq,$q2);
                }
                $R2_map[$R2_Qpos] = scalar @Qseq;
                $R2_Qpos++;
                $Rpos++;
              }
              case "D" {
                $Rpos++;
              }
              case "I" {
                if (nearby_quals(\@R1_Qual,$R1_Qpos,$mar) < nearby_quals(\@R2_Qual,$R2_Qpos,$mar)) {
                  push(@Qseq,$q2);
                }
                $R2_map[$R2_Qpos] = scalar @Qseq;
                $R2_Qpos++;
                unshift(@R1_cigar,$c1);
              }
            }
          }
        }
      }

      # Extend to end of R2 alignment
      while ($Rpos >= $R2_Rstart && $Rpos <= $R2_Rend ) {
        my $c2 = shift @R2_cigar;
        switch ($c2) {
          case "M" { 
            push(@Qseq,$R2_Qseq[$R2_Qpos]);
            $R2_map[$R2_Qpos] = scalar @Qseq;
            $R2_Qpos++;
            $Rpos++; }
          case "D" { $Rpos++; }
          case "I" {
            push(@Qseq,$R2_Qseq[$R2_Qpos]);
            $R2_map[$R2_Qpos] = scalar @Qseq;
            $R2_Qpos++; }
        }
      }

      # Add rest of R2
      if ($R2_Qpos < scalar @R2_Qseq) {
        my $old_Qpos = scalar @Qseq;
        push(@Qseq,@R2_Qseq[$R2_Qpos..$#R2_Qseq]);
        my $new_Qpos = scalar @Qseq;
        @R2_map[$R2_Qpos..$#R2_Qseq] = (($old_Qpos+1)..$new_Qpos);
        $R2_Qpos = scalar @R2_Qseq;
      }

      last;

    }
  }

  my $Qseq = join("",@Qseq);
  my $Qlen = length($Qseq) if defined $tlxls->[$#{$tlxls}]->{R2_ID} || $tlxls->[$#{$tlxls}]->{Rname} eq "Adapter";


  TLXL: foreach my $i (0..$#{$tlxls}) {

    my $B_tlxl;
    my $tlxl;

    my $tlx = {};

    switch ($#{$tlxls}) {
      case 0 {
        $B_tlxl = $tlxls->[0];
        $B_tlxl->{tlx} = $tlx;
      }
      case 1 {
        next TLXL if $i == 0;
        $B_tlxl = $tlxls->[0];
        $tlxl = $tlxls->[1];
        next TLXL if (defined $tlxl->{R1_Rgap} && $tlxl->{R1_Rgap} >=0 && $tlxl->{R1_Rgap} < 10) || 
          (defined $tlxl->{R2_Rgap} && $tlxl->{R2_Rgap} >=0 && $tlxl->{R2_Rgap} < 10 );
        $tlxl->{tlx} = $tlx;
      }
      else {
        next TLXL if $i == 0;
        $B_tlxl = $tlxls->[$i-1];
        $tlxl = $tlxls->[$i];
        last TLXL if $tlxl->{Rname} eq "Adapter";
        next TLXL if defined $tlxl->{R1_Rgap} && $tlxl->{R1_Rgap} >=0 && $tlxl->{R1_Rgap} < 10 ||
          (defined $tlxl->{R2_Rgap} && $tlxl->{R2_Rgap} >=0 && $tlxl->{R2_Rgap} < 10 );
        $tlxl->{tlx} = $tlx;
      }
    }

    $tlx->{Qname} = $B_tlxl->{Qname};
    $tlx->{Seq} = $Qseq;
    $tlx->{Qlen} = $Qlen;

    $tlx->{B_Rname} = $B_tlxl->{Rname};
    $tlx->{B_Strand} = $B_tlxl->{Strand};

    my $B_R1_Qstart;
    my $B_R2_Qstart;
    if (defined $B_tlxl->{R1_Qstart} && defined $R1_map[$B_tlxl->{R1_Qstart}-1]) {
      $B_R1_Qstart = $R1_map[$B_tlxl->{R1_Qstart}-1];
    }
    if (defined $B_tlxl->{R2_Qstart} && defined $R2_map[$B_tlxl->{R2_Qstart}-1]) {
      $B_R2_Qstart = $R2_map[$B_tlxl->{R2_Qstart}-1];
    }
    if (defined $B_R1_Qstart && defined $B_R2_Qstart) {
      if ($B_R1_Qstart <= $B_R2_Qstart) {
        $B_R2_Qstart = undef;
      } else {
        $B_R1_Qstart = undef;
      }
    }
    if (defined $B_R1_Qstart) {
      $tlx->{B_Qstart} = $R1_map[$B_tlxl->{R1_Qstart}-1];
      if ($tlx->{B_Strand} == 1) {
        $tlx->{B_Rstart} = $B_tlxl->{R1_Rstart};
      } else {
        $tlx->{B_Rend} = $B_tlxl->{R1_Rend};
      }
    } else {
      $tlx->{B_Qstart} = $R2_map[$B_tlxl->{R2_Qstart}-1];
      if ($tlx->{B_Strand} == 1) {
        $tlx->{B_Rstart} = $B_tlxl->{R2_Rstart};
      } else {
        $tlx->{B_Rend} = $B_tlxl->{R2_Rend};
      }
    }

    my $B_R1_Qend;
    my $B_R2_Qend;
    if (defined $B_tlxl->{R1_Qend} && defined $R1_map[$B_tlxl->{R1_Qend}-1]) {
      $B_R1_Qend = $R1_map[$B_tlxl->{R1_Qend}-1];
    }
    if (defined $B_tlxl->{R2_Qend} && defined $R2_map[$B_tlxl->{R2_Qend}-1]) {
      $B_R2_Qend = $R2_map[$B_tlxl->{R2_Qend}-1];
    }
    if (defined $B_R1_Qend && defined $B_R2_Qend) {
      if ($B_R1_Qend >= $B_R2_Qend) {
        $B_R2_Qend = undef;
      } else {
        $B_R1_Qend = undef;
      }
    }
    if (defined $B_R1_Qend) {
      $tlx->{B_Qend} = $R1_map[$B_tlxl->{R1_Qend}-1];
      if ($tlx->{B_Strand} == 1) {
        $tlx->{B_Rend} = $B_tlxl->{R1_Rend};
      } else {
        $tlx->{B_Rstart} = $B_tlxl->{R1_Rstart};
      }
    } else {
      $tlx->{B_Qend} = $R2_map[$B_tlxl->{R2_Qend}-1];
      if ($tlx->{B_Strand} == 1) {
        $tlx->{B_Rend} = $B_tlxl->{R2_Rend};
      } else {
        $tlx->{B_Rstart} = $B_tlxl->{R2_Rstart};
      }
    }


    next TLXL unless defined $tlxl;

    $tlx->{Rname} = $tlxl->{Rname};
    $tlx->{Strand} = $tlxl->{Strand};

    my $R1_Qstart;
    my $R2_Qstart;
    if (defined $tlxl->{R1_Qstart} && defined $R1_map[$tlxl->{R1_Qstart}-1]) {
      $R1_Qstart = $R1_map[$tlxl->{R1_Qstart}-1];
    }
    if (defined $tlxl->{R2_Qstart} && defined $R2_map[$tlxl->{R2_Qstart}-1]) {
      $R2_Qstart = $R2_map[$tlxl->{R2_Qstart}-1];
    }
    if (defined $R1_Qstart && defined $R2_Qstart) {
      if ($R1_Qstart <= $R2_Qstart) {
        $R2_Qstart = undef;
      } else {
        $R1_Qstart = undef;
      }
    }
    if (defined $R1_Qstart) {
      $tlx->{Qstart} = $R1_map[$tlxl->{R1_Qstart}-1];
      if ($tlx->{Strand} == 1) {
        $tlx->{Rstart} = $tlxl->{R1_Rstart};
      } else {
        $tlx->{Rend} = $tlxl->{R1_Rend};
      }
    } else {
      $tlx->{Qstart} = $R2_map[$tlxl->{R2_Qstart}-1];
      if ($tlx->{Strand} == 1) {
        $tlx->{Rstart} = $tlxl->{R2_Rstart};
      } else {
        $tlx->{Rend} = $tlxl->{R2_Rend};
      }
    }

    my $R1_Qend;
    my $R2_Qend;
    if (defined $tlxl->{R1_Qend} && defined $R1_map[$tlxl->{R1_Qend}-1]) {
      $R1_Qend = $R1_map[$tlxl->{R1_Qend}-1];
    }
    if (defined $tlxl->{R2_Qend} && defined $R2_map[$tlxl->{R2_Qend}-1]) {
      $R2_Qend = $R2_map[$tlxl->{R2_Qend}-1];
    }
    if (defined $R1_Qend && defined $R2_Qend) {
      if ($R1_Qend >= $R2_Qend) {
        $R2_Qend = undef;
      } else {
        $R1_Qend = undef;
      }
    }
    if (defined $R1_Qend) {
      $tlx->{Qend} = $R1_map[$tlxl->{R1_Qend}-1];
      if ($tlx->{Strand} == 1) {
        $tlx->{Rend} = $tlxl->{R1_Rend};
      } else {
        $tlx->{Rstart} = $tlxl->{R1_Rstart};
      }
    } else {
      $tlx->{Qend} = $R2_map[$tlxl->{R2_Qend}-1];
      if ($tlx->{Strand} == 1) {
        $tlx->{Rend} = $tlxl->{R2_Rend};
      } else {
        $tlx->{Rstart} = $tlxl->{R2_Rstart};
      }
    }

    $tlx->{Junction} = $tlx->{Strand} == 1 ? $tlx->{Rstart} : $tlx->{Rend};
    my $ref;
    switch ($tlx->{Rname}) {
      case "Breaksite" { $ref = $refs->{brk}; }
      case "Adapter" { $ref = $refs->{adpt}; }
      else { $ref = $refs->{genome}; }
    }
    $tlx->{JuncSeq} = $tlx->{Strand} == 1 ? $ref->seq($tlx->{Rname},$tlx->{Rstart}-10,$tlx->{Rstart}+9) :
                                            $ref->seq($tlx->{Rname},$tlx->{Rend}-9,$tlx->{Rend}+10);



  }

}


sub nearby_quals ($$$) {
  my $qual_ref = shift;
  my $pos = shift;
  my $mar = shift;

  return mean(@$qual_ref[max(0,$pos-$mar)..min($pos+$mar,$#{$qual_ref})]);
}




1;
