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
                      qw(R1_Rgap R2_Rgap R1_Cigar R2_Cigar) );
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
  $wrapper->{ID} = Data::GUID->new->as_string;

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

  $wrapper->{Cigar} = cigar_array_to_string($wrapper->{CigarA});
  $wrapper->{Qlen} = length($wrapper->{Seq});


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

sub find_random_barcode ($$$$) {
  my $read_obj = shift;
  my $barcode_length = shift;

  my $tlxs = $read_obj->{tlxs};
  my $R1_alns = $read_obj->{R1_alns};
  my $R2_alns = $read_obj->{R2_alns};


  my $barcode = "";


  if (defined $barcode_length && $barcode_length > 0) {
    

    # Search through OCS first

    if (defined $tlxs->[$#$tlxs]->{Rname} && $tlxs->[$#$tlxs]->{Rname} eq "Adapter") {
      my $adapter_aln = $tlxs->[$#$tlxs];
      $barcode = substr($adapter_aln->{Seq},$adapter_aln->{Qstart} - $barcode_length - 1,$barcode_length);
    } else {
      
      my $adapter_aln;

      foreach my $R2_aln (values %$R2_alns) {
        if ($R2_aln->{Unmapped} == 0 && $R2_aln->{Rname} eq "Adapter" && $R2_aln->{Strand} == 1) {
          $adapter_aln = $R2_aln if ! defined $adapter_aln ||
                                    $R2_aln->{Qend} > $adapter_aln->{Qend};
        }
      }

      foreach my $R1_aln (values %$R1_alns) {
        if ($R1_aln->{Unmapped} == 0 && $R1_aln->{Rname} eq "Adapter" && $R1_aln->{Strand} == 1) {
          $adapter_aln = $R1_aln if ! defined $adapter_aln ||
                                    $R1_aln->{Qend} > $adapter_aln->{Qend};
        }
      }

      $barcode = substr($adapter_aln->{Seq},$adapter_aln->{Qstart} - $barcode_length - 1,$barcode_length) if defined $adapter_aln;

    }

  }

  foreach my $tlx (@$tlxs) {
    $tlx->{Barcode} = $barcode;
  }

}


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

  # print("remapping a cigar\n$rseq\n$qseq\n@old_cigar\n");

  while ($Qpos <= @Qseq) {
    # print "$Qpos\n";
    my $c = shift(@old_cigar);
    switch($c) {
      case 'S' {
        push(@new_cigar,"S");
        $Qpos++;
      }
      case 'M' {
        if ($Qseq[$Qpos-1] eq $Rseq[$Rpos-1]) {
          push(@new_cigar,"M");
        } else {
          push(@new_cigar,"X")
        }
        $Qpos++;
        $Rpos++;
      }
      case 'D' {
        push(@new_cigar,"D");
        $Rpos++;
      }
      case 'I' {
        push(@new_cigar,"I");
        $Qpos++;
      }
    }
  }
  return(compact_cigar_array(\@new_cigar));
}


sub create_tlx_entries2 ($$) {

  my $tlxls = shift;
  my $refs = shift;

  my @tlxs;

  # Unmapped reads are easy
  # TLX will be one element long with simply the Qname and nothing else
  if ($tlxls->[0]->{Unmapped}) {
    my $tlx = {};
    $tlx->{Qname} = $tlxls->[0]->{Qname};
    push(@tlxs,$tlx);
    return(\@tlxs);
  }

  # We need to merge the alignments if there is an R1/R2 match
  # And return a query coordinate map for R1 -> merged seq
  # and for R2 -> merged seq
  my $Qseq;
  my $R1_map;
  my $R2_map;
  my $Cigar_array_ref;

  # First decide if we even need a merge
  my $R2_idx = firstidx {defined $_->{R2_Seq}} @$tlxls;

  if ($R2_idx < 0) {
    # Set R1 as Qseq and just a one-to-one map for R1 if no merge
    $Qseq = $tlxls->[0]->{R1_Seq};
    $R1_map = [1..length($Qseq)];
    $R2_map = [];
    $Cigar_array_ref = [];
  } else {
    # print "merging alignments on $R2_idx\n";
    ($Qseq,$R1_map,$R2_map,$Cigar_array_ref) = merge_alignments($tlxls,$refs) ;
  }


  foreach my $i (0..$#$tlxls) {
    # print "TLX for ".$i."th segment\n";


    my $tlx = {};
    my $b_tlxl = $tlxls->[$i];

    $tlx->{Qname} = $b_tlxl->{Qname};
    $tlx->{B_Rname} = $b_tlxl->{Rname};
    $tlx->{B_Strand} = $b_tlxl->{Strand};
    $tlx->{Seq} = $Qseq;
    $tlx->{Qlen} = $R2_idx < 0 ? "" : length($Qseq);

  
    # print "B-R1: ".$b_tlxl->{R1_Qstart}."-".$b_tlxl->{R1_Qend}."\n" if defined $b_tlxl->{R1_Qstart};
    # print "B-R2: ".$b_tlxl->{R2_Qstart}."-".$b_tlxl->{R2_Qend}."\n" if defined $b_tlxl->{R2_Qstart};


    my $b_ref;
    switch ($tlx->{B_Rname}) {
      case "Breaksite" { $b_ref = $refs->{brk}; }
      case "Adapter" { $b_ref = $refs->{adpt}; }
      else { $b_ref = $refs->{genome}; }
    }



    if ($i < $R2_idx || $R2_idx < 0) {
      # Segment earlier than the merge or no merge
      # Use R1
      $tlx->{B_Rstart} = $b_tlxl->{R1_Rstart};
      $tlx->{B_Rend} = $b_tlxl->{R1_Rend};
      $tlx->{B_Qstart} = $R1_map->[$b_tlxl->{R1_Qstart}-1];
      $tlx->{B_Qend} = $R1_map->[$b_tlxl->{R1_Qend}-1];

      $tlx->{B_R1_ID} = $b_tlxl->{R1_ID};

      my $Rseq = $b_ref->seq($tlx->{B_Rname},$tlx->{B_Rstart},$tlx->{B_Rend});
      $Rseq = reverseComplement($Rseq) unless $tlx->{B_Strand} == 1;
      # print "retrieved rseq ".$tlx->{B_Rname}." ".$tlx->{B_Rstart}." ".$tlx->{B_Rend}."\n$Rseq\n";

      $tlx->{B_CigarA} = soft_clip_cigar(remap_cigar($b_tlxl->{R1_CigarA},$b_tlxl->{R1_Seq},$Rseq));

    } elsif ($i == $R2_idx) {
      # This is the segment we merged on
      if ($tlx->{B_Strand} == 1) {
        # "+" strand alignment
        $tlx->{B_Rstart} = $b_tlxl->{R1_Rstart};
        $tlx->{B_Rend} = $b_tlxl->{R2_Rend};

      } else {
        # "-" strand alignment
        $tlx->{B_Rstart} = $b_tlxl->{R2_Rstart};
        $tlx->{B_Rend} = $b_tlxl->{R1_Rend};
      }
      $tlx->{B_Qstart} = $R1_map->[$b_tlxl->{R1_Qstart}-1];
      $tlx->{B_Qend} = $R2_map->[$b_tlxl->{R2_Qend}-1];
      $tlx->{B_R1_ID} = $b_tlxl->{R1_ID};
      $tlx->{B_R2_ID} = $b_tlxl->{R2_ID};
      $tlx->{B_CigarA} = soft_clip_cigar($Cigar_array_ref);

    } else {
      # Segment post merge
      # Use R2
      $tlx->{B_Rstart} = $b_tlxl->{R2_Rstart};
      $tlx->{B_Rend} = $b_tlxl->{R2_Rend};
      $tlx->{B_Qstart} = $R2_map->[$b_tlxl->{R2_Qstart}-1];
      $tlx->{B_Qend} = $R2_map->[$b_tlxl->{R2_Qend}-1];
      $tlx->{B_R2_ID} = $b_tlxl->{R2_ID};

      my $Rseq = $b_ref->seq($tlx->{B_Rname},$tlx->{B_Rstart},$tlx->{B_Rend});
      $Rseq = reverseComplement($Rseq) unless $tlx->{B_Strand} == 1;
      # print "retrieved rseq ".$tlx->{B_Rname}." ".$tlx->{B_Rstart}." ".$tlx->{B_Rend}."\n$Rseq\n";

      $tlx->{B_CigarA} = soft_clip_cigar(remap_cigar($b_tlxl->{R2_CigarA},$b_tlxl->{R2_Seq},$Rseq));

    }

    $tlx->{B_Cigar} = cigar_array_to_string($tlx->{B_CigarA});

    
    if ($i+1 <= $#$tlxls) {
      my $tlxl = $tlxls->[$i+1];

      # print "R1: ".$tlxl->{R1_Qstart}."-".$tlxl->{R1_Qend}."\n" if defined $tlxl->{R1_Qstart};
      # print "R2: ".$tlxl->{R2_Qstart}."-".$tlxl->{R2_Qend}."\n" if defined $tlxl->{R2_Qstart};

      $tlx->{Rname} = $tlxl->{Rname};
      $tlx->{Strand} = $tlxl->{Strand};


      my $ref;
      switch ($tlx->{Rname}) {
        case "Breaksite" { $ref = $refs->{brk}; }
        case "Adapter" { $ref = $refs->{adpt}; }
        else { $ref = $refs->{genome}; }
      }


      if ($i+1 < $R2_idx || $R2_idx < 0) {
        $tlx->{Rstart} = $tlxl->{R1_Rstart};
        $tlx->{Rend} = $tlxl->{R1_Rend};
        $tlx->{Qstart} = $R1_map->[$tlxl->{R1_Qstart}-1];
        $tlx->{Qend} = $R1_map->[$tlxl->{R1_Qend}-1];

        my $Rseq = $ref->seq($tlx->{Rname},$tlx->{Rstart},$tlx->{Rend});
        $Rseq = reverseComplement($Rseq) unless $tlx->{Strand} == 1;
        # print "retrieved rseq ".$tlx->{Rname}." ".$tlx->{Rstart}." ".$tlx->{Rend}."\n$Rseq\n";

        $tlx->{CigarA} = soft_clip_cigar(remap_cigar($tlxl->{R1_CigarA},$tlxl->{R1_Seq},$Rseq));

        $tlx->{R1_ID} = $tlxl->{R1_ID};

      } elsif ($i+1 == $R2_idx) {
        if ($tlx->{Strand} == 1) {
          $tlx->{Rstart} = $tlxl->{R1_Rstart};
          $tlx->{Rend} = $tlxl->{R2_Rend};
        } else{
          $tlx->{Rstart} = $tlxl->{R2_Rstart};
          $tlx->{Rend} = $tlxl->{R1_Rend};
        }
        $tlx->{Qstart} = $R1_map->[$tlxl->{R1_Qstart}-1];
        $tlx->{Qend} = $R2_map->[$tlxl->{R2_Qend}-1];

        $tlx->{CigarA} = soft_clip_cigar($Cigar_array_ref);

        $tlx->{R1_ID} = $tlxl->{R1_ID};
        $tlx->{R2_ID} = $tlxl->{R2_ID};

      } else {
        $tlx->{Rstart} = $tlxl->{R2_Rstart};
        $tlx->{Rend} = $tlxl->{R2_Rend};
        $tlx->{Qstart} = $R2_map->[$tlxl->{R2_Qstart}-1];
        $tlx->{Qend} = $R2_map->[$tlxl->{R2_Qend}-1];


        my $Rseq = $ref->seq($tlx->{Rname},$tlx->{Rstart},$tlx->{Rend});
        $Rseq = reverseComplement($Rseq) unless $tlx->{Strand} == 1;
        # print "retrieved rseq ".$tlx->{Rname}." ".$tlx->{Rstart}." ".$tlx->{Rend}."\n$Rseq\n";

        $tlx->{CigarA} = soft_clip_cigar(remap_cigar($tlxl->{R2_CigarA},$tlxl->{R2_Seq},$Rseq));


        $tlx->{R2_ID} = $tlxl->{R2_ID};

      }

      $tlx->{Cigar} = cigar_array_to_string($tlx->{CigarA});

      if ($tlx->{Strand} == 1) {
        $tlx->{Junction} = $tlx->{Rstart};
      } else {
        $tlx->{Junction} = $tlx->{Rend};
      }
      
      $tlx->{J_Seq} = "";

    }

    push(@tlxs,$tlx);

    last if $i+1 >= $#$tlxls;
  }

  return(\@tlxs);

}

sub merge_alignments ($$) {

  my $tlxls = shift;
  my $refs = shift;

  # First decide if we even need a merge
  my $R2_idx = firstidx {defined $_->{R2_Seq}} @$tlxls;

  if ($R2_idx < 0) {
    # Return a blank seq and some empty arraryrefs if not
    return("",[],[]);
  }

  # Pull the element to merge on out of the tlxl array
  my $tlxl = $tlxls->[$R2_idx];


  # We're gonna spend a while getting setup to do the merge
  # Here's the deal, in order to treat "+" and "-" strands the same
  # we have to take advantage of the symmetry
  # We do this by swapping reads 1 and 2 if its a "-" strand alignment
  # More work now and at the end when we have to flip "-" strands,
  # but way less work when we are actually ticking through and merging the reads

  my $Strand = $tlxl->{Strand};
  my $Rname = $tlxl->{Rname};
  my $Rstart1;
  my $Rend1;
  my $Rstart2;
  my $Rend2;

  my @Qseq1;
  my @Qseq2;
  my @Qual1;
  my @Qual2;
  my @Cigar1 = ();
  my @Cigar2 = ();
  my $Qstart1;
  my $Qend1;
  my $Qstart2;
  my $Qend2;

  # Need to know which sam object to use
  my $ref;
  switch ($Rname) {
    case "Breaksite" { $ref = $refs->{brk}; }
    case "Adapter" { $ref = $refs->{adpt}; }
    else { $ref = $refs->{genome}; }
  }



  if ($Strand == 1) {
    # print "+ strand alignment\n";
    # Query sequence and quality info for each pair split into arrays
    @Qseq1 = split("",$tlxl->{R1_Seq});
    @Qseq2 = split("",$tlxl->{R2_Seq});
    @Qual1 = @{$tlxl->{R1_Qual}};
    @Qual2 = @{$tlxl->{R2_Qual}};

    # Create cigar arrays

    @Cigar1 = @{ expand_cigar_array(soft_clip_cigar($tlxl->{R1_CigarA})) };
    @Cigar2 = @{ expand_cigar_array(soft_clip_cigar($tlxl->{R2_CigarA})) };

    # foreach my $i (@{$tlxl->{R1_CigarA}}) {
    #   next if $i->[0] eq "S";
    #   push(@Cigar1, ($i->[0]) x $i->[1]);
    # }
    # foreach my $i (@{$tlxl->{R2_CigarA}}) {
    #   next if $i->[0] eq "S";
    #   push(@Cigar2, ($i->[0]) x $i->[1]);
    # }

    # Pretty straight forward this time through
    $Qstart1 =  $tlxl->{R1_Qstart};
    $Qend1 = $tlxl->{R1_Qend};
    $Qstart2 = $tlxl->{R2_Qstart};
    $Qend2 = $tlxl->{R2_Qend};

    $Rstart1 =  $tlxl->{R1_Rstart};
    $Rend1 = $tlxl->{R1_Rend};
    $Rstart2 = $tlxl->{R2_Rstart};
    $Rend2 = $tlxl->{R2_Rend};


  } else {
    # print "- strand alignment\n";
    # Query info again, this time flipping R1/R2, reversing (and complementing)
    @Qseq1 = split("",reverseComplement($tlxl->{R2_Seq}));
    @Qseq2 = split("",reverseComplement($tlxl->{R1_Seq}));
    @Qual1 = reverse @{$tlxl->{R2_Qual}};
    @Qual2 = reverse @{$tlxl->{R1_Qual}};

    # Create cigar arrays, same as before, but reversed

    @Cigar1 = reverse @{ expand_cigar_array(soft_clip_cigar($tlxl->{R2_CigarA})) };
    @Cigar2 = reverse @{ expand_cigar_array(soft_clip_cigar($tlxl->{R1_CigarA})) };

    # foreach my $i (reverse @{$tlxl->{R2_CigarA}}) {
    #   next if $i->[0] eq "S";
    #   push(@Cigar1, ($i->[0]) x $i->[1]);
    # }
    # foreach my $i (reverse @{$tlxl->{R1_CigarA}}) {
    #   next if $i->[0] eq "S";
    #   push(@Cigar2, ($i->[0]) x $i->[1]);
    # }

    # This time flip all R1 with R2
    # And since we're reverse complimenting
    $Qstart1 = $tlxl->{R2_Qlen} - $tlxl->{R2_Qend} + 1;
    $Qend1 = $tlxl->{R2_Qlen} - $tlxl->{R2_Qstart} + 1;
    $Qstart2 = $tlxl->{R1_Qlen} - $tlxl->{R1_Qend} + 1;
    $Qend2 = $tlxl->{R1_Qlen} - $tlxl->{R1_Qstart} + 1;

    $Rstart1 =  $tlxl->{R2_Rstart};
    $Rend1 = $tlxl->{R2_Rend};
    $Rstart2 = $tlxl->{R1_Rstart};
    $Rend2 = $tlxl->{R1_Rend};
  
  }

  # Retrieve reference sequence
  my @Rseq = split("",uc($ref->seq($Rname,$Rstart1,$Rend2)));


  # I think we're ready to merge!

  # Here's where we'll save the returned variables
  my @Qseq;
  my @Qmap1 = (0) x @Qseq1;
  my @Qmap2 = (0) x @Qseq2;
  my @Cigar;


  # We will be interating through these
  my $Qpos1 = $Qstart1;
  my $Qpos2 = $Qstart2;
  my $Rpos = $Rstart1;


  # Handle everything on R1 up to Qstart
  # Note!!! For "-" alignments, R1 is actually temporarily R2
  # (but we can forget that because symmetry)
  if ($Qstart1 > 1) {
    # print "premerge\n";

    my @R1_tail = 0..($Qstart1-2);
    push(@Qseq, @Qseq1[@R1_tail]);
    push(@Cigar,("S") x @R1_tail);
    @Qmap1[@R1_tail] = 1..@R1_tail;
  }


  # Slightly tricky... if there is a dovetail at the near end 
  # R1    |------------->
  # R2  <------------------------|
  # We have to push forward on the R2 query so it's caught up
  # Basically this means shifting off the front of the cigar array
  # and incrementing $Qpos2 until we're at the same Ref position as R1
  if ($Rstart2 < $Rstart1) {
    # print "pushing R2 premerge\n";
    my $Rpos_tmp = $Rstart2;
    while ($Rpos_tmp < $Rstart1) {
      my $c2 = shift(@Cigar2);
      switch ($c2) {
        case 'M' {
          $Qpos2++;
          $Rpos_tmp++;
        }
        case 'D' {
          $Rpos_tmp++;
        }
        case 'I' {
          $Qpos2++;
        }
      }
    }
  }

  # Push ahead on R1 where there is no overlap
  # So either until the end of R1 or the start of R2
  # print "pushing R1 at $Rpos\n";
  while ($Rpos <= $Rend1 && $Rpos < $Rstart2) {
    my $c1 = shift(@Cigar1);
    switch ($c1) {
      case 'M' {
        push(@Qseq,$Qseq1[$Qpos1-1]);
        if ($Qseq1[$Qpos1-1] eq $Rseq[$Rpos-$Rstart1]) {
          push(@Cigar,"M");
        } else {
          push(@Cigar,"X");
        }
        $Qmap1[$Qpos1-1] = scalar @Qseq;
        $Qpos1++;
        $Rpos++;
      }
      case 'D' {
        push(@Cigar,"D");
        $Rpos++;
      }
      case 'I' {
        push(@Qseq,$Qseq1[$Qpos1-1]);
        push(@Cigar,"I");
        $Qmap1[$Qpos1-1] = scalar @Qseq;
        $Qpos1++;
      }
    }
  }


  # Handles the gap situation
  # R1   |------------>
  # R2                      <--------------|
  if ($Rstart2 - $Rend1 > 1) {
    # print "handling gap at $Rpos\n";
    my $Gap_start = $Rend1 - $Rstart1 + 1;
    my $Gap_end = $Rstart2 - $Rstart1 - 1;
    push(@Qseq, map {lc $_} @Rseq[$Gap_start..$Gap_end]);
    push(@Cigar,("N") x ($Gap_end-$Gap_start+1));
    $Rpos = $Rstart2;
  }

  # Handles the overlap situation
  # R1  |------------->
  # R2       <--------------|
  if ($Rend1 >= $Rstart2) {
    # print "handling overlap at $Rpos\n";
    while ($Rpos <= $Rend1 && $Rpos <= $Rend2) {
      my $c1 = shift(@Cigar1);
      my $c2 = shift(@Cigar2);
      my $q1 = $Qseq1[$Qpos1-1];
      my $q2 = $Qseq2[$Qpos2-1];
      my $r = $Rseq[$Rpos-$Rstart1];


      switch($c1) {
        case 'M' {
          switch($c2) {
            case 'M' {
              # Match - Match
              # Take base if they match each other,
              # otherwise take the one that matches the reference
              # otherwise take the one with the higher qual score
              # Update maps and push forward on Q1 Q2 and R
              my $n;
              if ($q1 eq $q2) {
                $n = $q1;
              } elsif ($q1 eq $r || $q2 eq $r) {
                $n = $r;
              } else {
                $n = $Qual1[$Qpos1-1] > $Qual2[$Qpos2-1] ? $q1 : $q2;
              }
              push(@Qseq,$n);
              if ($n eq $r) {
                push(@Cigar,"M");
              } else {
                # print "$Rpos, $r, $Qpos1, $q1, $Qpos2, $q2\n";
                push(@Cigar,"X");
              }
              $Qmap1[$Qpos1-1] = scalar @Qseq;
              $Qmap2[$Qpos2-1] = scalar @Qseq;
              $Qpos1++;
              $Qpos2++;
              $Rpos++;
            }
            case 'D' {
              # Match - Del
              # Take matched base if it matches reference
              # Update maps, and push forward on Q1 and R
              if ($q1 eq $r) {
                push(@Qseq,$q1);
                push(@Cigar,"M");
              } else {
                push(@Cigar,"D");
              }
              $Qmap1[$Qpos1-1] = scalar @Qseq;
              $Qpos1++;
              $Rpos++;
            }
            case 'I' {
              # Match - Ins
              # Put cigar back for C1
              # Update maps and push forward on Q2
              unshift(@Cigar1,$c1);
              $Qmap2[$Qpos2-1] = scalar @Qseq;
              $Qpos2++;
            }
          }
        }
        case 'D' {
          switch($c2) {
            case 'M' {
              # Del - Match
              # Take matched base if it matches reference
              # Update maps, and push forward on Q2 and R
              if ($q2 eq $r) {
                push(@Qseq,$q2);
                push(@Cigar,"M");
              } else {
                push(@Cigar,"D");
              }
              $Qmap2[$Qpos2-1] = scalar @Qseq;
              $Qpos2++;
              $Rpos++;
            }
            case 'D' {
              # Del - Del
              # Push forward on R
              push(@Cigar,"D");
              $Rpos++;
            }
            case 'I' {
              # Del - Ins
              # Put cigar back for C1
              # Update maps and push forward on Q2
              unshift(@Cigar1,$c1);
              $Qmap2[$Qpos2-1] = scalar @Qseq;
              $Qpos2++;
            }
          }
        }
        case 'I' {
          switch($c2) {
            case 'M' {
              # Match - Ins
              # Put cigar back for C2
              # Update maps and push forward on Q1
              unshift(@Cigar2,$c2);
              $Qmap1[$Qpos1-1] = scalar @Qseq;
              $Qpos1++;
            }
            case 'D' {
              # Ins - Del
              # Put cigar back for C2
              # Update maps and push forward on Q1
              unshift(@Cigar2,$c2);
              $Qmap1[$Qpos1-1] = scalar @Qseq;
              $Qpos1++;
            }
            case 'I' {
              # Ins - Ins
              # Take higher quality base
              # Update maps and push forward on Q1 and Q2
              push(@Qseq, $Qual1[$Qpos1-1] > $Qual2[$Qpos2-1] ? $q1 : $q2);
              push(@Cigar,"I");
              $Qmap1[$Qpos1-1] = scalar @Qseq;
              $Qmap2[$Qpos2-1] = scalar @Qseq;
              $Qpos1++;
              $Qpos2++;
            }
          }
        }
      }
    }
  }


  # Push ahead on R2
  # print "pushing R2 at $Rpos\n";
  while ($Rpos <= $Rend2) {
    my $c2 = shift(@Cigar2);
    switch ($c2) {
      case 'M' {
        push(@Qseq,$Qseq2[$Qpos2-1]);
        if ($Qseq2[$Qpos2-1] eq $Rseq[$Rpos-$Rstart1]) {
          push(@Cigar,"M");
        } else {
          push(@Cigar,"X");
        }
        $Qmap2[$Qpos2-1] = scalar @Qseq;
        $Qpos2++;
        $Rpos++;
      }
      case 'D' {
        push(@Cigar,"D");
        $Rpos++;
      }
      case 'I' {
        push(@Qseq,$Qseq2[$Qpos2-1]);
        push(@Cigar,"I");
        $Qmap2[$Qpos2-1] = scalar @Qseq;
        $Qpos2++;
      }
    }
  }

  # Same as before, if there's a dovetail at the far end
  # R1  |-------------------------->
  # R2       <-------------------|
  # Just clean up R1 query
  if ($Rend1 > $Rend2) {
    # I think actually we don't have to worry about this end
    # since we don't care about $Qpos1 and @Cigar1 by this point
  }

  # Handle everything on R2 after Qend
  if (scalar @Qseq2 > $Qend2) {
    # print "postmerge\n";
    my @R2_tail = ($Qend2)..$#Qseq2;
    push(@Qseq, @Qseq2[@R2_tail]);
    push(@Cigar,("S") x @R2_tail);
    @Qmap2[@R2_tail] = (scalar @Qseq - scalar @R2_tail - 1)..(scalar @Qseq);
  }

  my $Qseq = join("",@Qseq);


  # print "R1\n";
  # print join("",map {chr($_+33)} @Qual1)."\n";
  # print join("",@Qseq1)."\n";
  # print join(" ",@Qmap1)."\n";
  
  # print "R2\n";
  # print join("",map {chr($_+33)} @Qual2)."\n";
  # print join("",@Qseq2)."\n";
  # print join(" ",@Qmap2)."\n";  

  # print "merge\n";
  # print "$Qseq\n";


  # Flip things around if it was - strand alignment
  if ($Strand == -1) {
    $Qseq = reverseComplement($Qseq);
    @Cigar = reverse @Cigar;
    my @Qmap1_tmp = @Qmap1;
    my @Qmap2_tmp = @Qmap2;
    @Qmap1 = map {$_ == 0 ? 0 : (scalar @Qseq - $_ + 1) } reverse @Qmap2_tmp;
    @Qmap2 = map {$_ == 0 ? 0 : (scalar @Qseq - $_ + 1) } reverse @Qmap1_tmp;

    # print "reversed\n";
    # print join(" ",@Qmap1)."\n";
    # print join(" ",@Qmap2)."\n";  
    # print "$Qseq\n";
  }

  return($Qseq,\@Qmap1,\@Qmap2,compact_cigar_array(\@Cigar));

}


# sub create_tlx_entries ($$) {

#   my $tlxls = shift;
#   my $refs = shift;

#   my $mar = 2;
#   # my $n_juncs;

#   my $tlxs = [];

#   if ($tlxls->[0]->{Unmapped}) {
#     my $tlx = {};
#     $tlx->{Qname} = $tlxls->[0]->{Qname};
#     $tlx->{Filter} = "Unaligned";
#     $tlxls->[0]->{tlx} = $tlx;
#     return;
#   }

#   my @Qseq = ();
#   my @R1_map = ();
#   my @R2_map = ();

#   my $R1_Qpos = 0;
#   my $R2_Qpos = 0;

#   my @R1_Qseq = split("",$tlxls->[0]->{R1_Seq});
#   my $R2_idx = firstidx {defined $_->{R2_Seq}} @$tlxls;
#   my @R2_Qseq = split("",$tlxls->[$R2_idx]->{R2_Seq}) if $R2_idx >= 0;

#   foreach my $i (0..$#{$tlxls}) {

#     my $tlxl = $tlxls->[$i];

#     # print "\nsegment $i\n";

#     unless ( defined $tlxls->[$i]->{R2_ID} ) { #(defined $tlxls->[$i+1] && defined $tlxls->[$i+1]->{R1_ID}) || ! defined $tlxl->{R2_ID}) {
#       # print "R1 only\n";
#       if ($R1_Qpos < $tlxl->{R1_Qend}) {
#         my $old_Qpos = scalar @Qseq;
#         push(@Qseq,@R1_Qseq[$R1_Qpos..($tlxl->{R1_Qend}-1)]);
#         my $new_Qpos = scalar @Qseq;
#         @R1_map[$R1_Qpos..($tlxl->{R1_Qend}-1)] = (($old_Qpos+1)..$new_Qpos);
#         $R1_Qpos = $tlxl->{R1_Qend};
#       }

#     } else {
#       # print "merging R1 and R2\n";
    
#       my $Rname = $tlxl->{Rname};
#       my $Strand = $tlxl->{Strand};

#       my @R1_Qual = @{$tlxl->{R1_Qual}};
#       @R1_Qual = reverse @R1_Qual unless $Strand == 1;
#       my @R2_Qual = @{$tlxl->{R2_Qual}};
#       @R2_Qual = reverse @R2_Qual if $Strand == 1;
#       my @R1_cigar = ();
#       foreach my $i (@{$tlxl->{R1_CigarA}}) {
#         next if $i->[0] eq "S";
#         push(@R1_cigar,split("",($i->[0] x $i->[1])));
#       }
#       @R1_cigar = reverse @R1_cigar unless $Strand == 1;
#       my @R2_cigar = ();
#       foreach my $i (@{$tlxl->{R2_CigarA}}) {
#         next if $i->[0] eq "S";
#         push(@R2_cigar,split("",($i->[0] x $i->[1])));
#       }
#       @R2_cigar = reverse @R2_cigar unless $Strand == 1;

#       my $R1_Rstart = $tlxl->{R1_Rstart};
#       my $R1_Rend = $tlxl->{R1_Rend};
#       my $R2_Rstart = $tlxl->{R2_Rstart};
#       my $R2_Rend = $tlxl->{R2_Rend};

#       my $Rstart = $tlxl->{Strand} == 1 ? $tlxl->{R1_Rstart} : $tlxl->{R2_Rstart};
#       my $Rend = $tlxl->{Strand} == 1 ? $tlxl->{R2_Rend} : $tlxl->{R1_Rend};
#       my $ref;
#       switch ($Rname) {
#         case "Breaksite" { $ref = $refs->{brk}; }
#         case "Adapter" { $ref = $refs->{adpt}; }
#         else { $ref = $refs->{genome}; }
#       }

#       # Retrieve referense sequence, reverse complement if necessary

#       my $Rseq = $ref->seq($Rname,$Rstart,$Rend);
#       # print join("\t",$tlxl->{Qname},$Rname,$Rstart,$Rend,$tlxl->{Strand})."\n$Rseq\n";

#       my @Rseq = $Strand == 1 ? split("",$Rseq) : split("",reverseComplement($Rseq));

#       # Artifically adjust Rstarts and Rends so that we can move incrementally up on both Qpos and Rpos
#       ($R1_Rstart,$R1_Rend) = $Strand == 1 ? ($R1_Rstart,$R1_Rend) : (-$R1_Rend,-$R1_Rstart);
#       ($R2_Rstart,$R2_Rend) = $Strand == 1 ? ($R2_Rstart,$R2_Rend) : (-$R2_Rend,-$R2_Rstart);

#       # print "R1: $R1_Rstart-$R1_Rend\n";
#       # print join(" ",@R1_cigar)."\n";
#       # print "R2: $R2_Rstart-$R2_Rend\n";
#       # print join(" ",@R2_cigar)."\n";

#       $R2_Qpos = $tlxl->{R2_Qstart} - 1;
#       my $Rpos = $R1_Rstart;

#       # Push forward on R1 if behind start of alignment (gap between this and last alignment)
#       if ($R1_Qpos < $tlxl->{R1_Qstart} - 1) {
#         my $old_Qpos = scalar @Qseq;
#         push(@Qseq,@R1_Qseq[$R1_Qpos..($tlxl->{R1_Qstart}-2)]);
#         my $new_Qpos = scalar @Qseq;
#         @R1_map[$R1_Qpos..($tlxl->{R1_Qstart}-2)] = (($old_Qpos+1)..$new_Qpos);
#         $R1_Qpos = $tlxl->{R1_Qstart} - 1;
#       }

#       # Push forward on R1 alignment if start is behind current R1_Qpos (overlap between this and last alignment)
#       if ($tlxl->{R1_Qstart} - 1 < $R1_Qpos) {
#         my $pos1_tmp = $tlxl->{R1_Qstart} - 1;
#         while ($pos1_tmp < $R1_Qpos) {
#           my $c1 = shift @R1_cigar;
#           switch ($c1) {
#             case "M" { $pos1_tmp++; $Rpos++; }
#             case "D" { $Rpos++; }
#             case "I" { $pos1_tmp++; }
#             else { $pos1_tmp++; $Rpos++; }
#           }
#         }
#       }

#       # Push forward on R2 if Rstart is behind current Rpos
#       if ($R2_Rstart < $Rpos) {
#         my $Rpos_tmp = $R2_Rstart;
#         while ($Rpos_tmp < $Rpos) {
#           my $c2 = shift @R2_cigar;
#           switch ($c2) {
#             case "M" { $R2_Qpos++; $Rpos_tmp++; }
#             case "D" { $Rpos_tmp++; }
#             case "I" { $R2_Qpos++; }
#             else { $R2_Qpos++; $Rpos_tmp++; }
#           }
#         }
#       }

#       # Add R1 alignment up to the start of R2 alignment
#       while ($Rpos <= $R1_Rend && $Rpos < $R2_Rstart) {
#         my $c1 = shift @R1_cigar;
#         switch ($c1) {
#           case "M" { 
#             push(@Qseq,$R1_Qseq[$R1_Qpos]);
#             $R1_map[$R1_Qpos] = scalar @Qseq;
#             $R1_Qpos++;
#             $Rpos++; }
#           case "D" { $Rpos++; }
#           case "I" {
#             push(@Qseq,$R1_Qseq[$R1_Qpos]);
#             $R1_map[$R1_Qpos] = scalar @Qseq;
#             $R1_Qpos++; }
#         }
#       }

#       # Fill in gap between alignment if exists
#       while ($Rpos > $R1_Rend && $Rpos <= $R2_Rstart) {
#         push(@Qseq,lc($Rseq[$Rpos - $R1_Rstart]));
#         $Rpos++;
#       }

#       # Negotiate the R1 and R2 alignments
#       while ($Rpos <= $R1_Rend && $Rpos >= $R2_Rstart && $Rpos <= $R2_Rend) {
#         my $c1 = shift @R1_cigar;
#         my $c2 = shift @R2_cigar;
#         my $q1 = $R1_Qseq[$R1_Qpos];
#         my $q2 = $R2_Qseq[$R2_Qpos];
#         my $r = $Rseq[$Rpos - $R1_Rstart];
#         # print "R1: $R1_Qpos $c1 / R2: $R2_Qpos $c2\n";

#         switch ($c1) {
#           case "M" {
#             switch ($c2) {
#               case "M" {
#                 if ($q1 eq $r || $q2 eq $r) {
#                   push(@Qseq,$r);
#                 } elsif ($q1 eq $q2) {
#                   push(@Qseq,$q1);
#                 } elsif (nearby_quals(\@R1_Qual,$R1_Qpos,$mar) >= nearby_quals(\@R2_Qual,$R2_Qpos,$mar)) {
#                   push(@Qseq,$q1);
#                 } else {
#                   push(@Qseq,$q2);
#                 }
#                 $R1_map[$R1_Qpos] = scalar @Qseq;
#                 $R2_map[$R2_Qpos] = scalar @Qseq;
#                 $R1_Qpos++;
#                 $R2_Qpos++;
#                 $Rpos++;
#               }
#               case "D" {
#                 if ($q1 eq $r || nearby_quals(\@R1_Qual,$R1_Qpos,$mar) >= nearby_quals(\@R2_Qual,$R2_Qpos,$mar)) {
#                   push(@Qseq,$q1);
#                 }
#                 $R1_map[$R1_Qpos] = scalar @Qseq;
#                 $R1_Qpos++;
#                 $Rpos++;
#               }
#               case "I" {
#                 if ($q1 ne $r && nearby_quals(\@R1_Qual,$R1_Qpos,$mar) < nearby_quals(\@R2_Qual,$R2_Qpos,$mar)) {
#                   push(@Qseq,$q2);
#                 }
#                 $R2_map[$R2_Qpos] = scalar @Qseq;
#                 $R2_Qpos++;
#                 unshift(@R1_cigar,$c1);
#               }
#             }
#           }
#           case "I" {
#             switch ($c2) {
#               case "M" {
#                 if ($q2 ne $r && nearby_quals(\@R1_Qual,$R1_Qpos,$mar) >= nearby_quals(\@R2_Qual,$R2_Qpos,$mar)) {
#                   push(@Qseq,$q1);
#                 }
#                 $R1_map[$R1_Qpos] = scalar @Qseq;
#                 $R1_Qpos++;
#                 unshift(@R2_cigar,$c2);
#               }
#               case "D" {
#                 if (nearby_quals(\@R1_Qual,$R1_Qpos,$mar) >= nearby_quals(\@R2_Qual,$R2_Qpos,$mar)) {
#                   push(@Qseq,$q1);
#                 }
#                 $R1_map[$R1_Qpos] = scalar @Qseq;
#                 $R1_Qpos++;
#                 unshift(@R2_cigar,$c2);
#               }
#               case "I" {
#                 if ($q1 eq $q2 || nearby_quals(\@R1_Qual,$R1_Qpos,$mar) >= nearby_quals(\@R2_Qual,$R2_Qpos,$mar)) {
#                   push(@Qseq,$q1);
#                 } else {
#                   push(@Qseq,$q2);
#                 }
#                 $R1_map[$R1_Qpos] = scalar @Qseq;
#                 $R2_map[$R2_Qpos] = scalar @Qseq;
#                 $R1_Qpos++;
#                 $R2_Qpos++;
#               }
#             }
#           }
#           case "D" {
#             switch ($c2) {
#               case "M" {
#                 if ($q2 eq $r || nearby_quals(\@R1_Qual,$R1_Qpos,$mar) < nearby_quals(\@R2_Qual,$R2_Qpos,$mar)) {
#                   push(@Qseq,$q2);
#                 }
#                 $R2_map[$R2_Qpos] = scalar @Qseq;
#                 $R2_Qpos++;
#                 $Rpos++;
#               }
#               case "D" {
#                 $Rpos++;
#               }
#               case "I" {
#                 if (nearby_quals(\@R1_Qual,$R1_Qpos,$mar) < nearby_quals(\@R2_Qual,$R2_Qpos,$mar)) {
#                   push(@Qseq,$q2);
#                 }
#                 $R2_map[$R2_Qpos] = scalar @Qseq;
#                 $R2_Qpos++;
#                 unshift(@R1_cigar,$c1);
#               }
#             }
#           }
#         }
#       }

#       # Extend to end of R2 alignment
#       while ($Rpos >= $R2_Rstart && $Rpos <= $R2_Rend ) {
#         my $c2 = shift @R2_cigar;
#         switch ($c2) {
#           case "M" { 
#             push(@Qseq,$R2_Qseq[$R2_Qpos]);
#             $R2_map[$R2_Qpos] = scalar @Qseq;
#             $R2_Qpos++;
#             $Rpos++; }
#           case "D" { $Rpos++; }
#           case "I" {
#             push(@Qseq,$R2_Qseq[$R2_Qpos]);
#             $R2_map[$R2_Qpos] = scalar @Qseq;
#             $R2_Qpos++; }
#         }
#       }

#       # Add rest of R2
#       if ($R2_Qpos < scalar @R2_Qseq) {
#         my $old_Qpos = scalar @Qseq;
#         push(@Qseq,@R2_Qseq[$R2_Qpos..$#R2_Qseq]);
#         my $new_Qpos = scalar @Qseq;
#         @R2_map[$R2_Qpos..$#R2_Qseq] = (($old_Qpos+1)..$new_Qpos);
#         $R2_Qpos = scalar @R2_Qseq;
#       }

#       last;

#     }
#   }

#   my $Qseq = join("",@Qseq);
#   my $Qlen = length($Qseq) if defined $tlxls->[$#{$tlxls}]->{R2_ID} || $tlxls->[$#{$tlxls}]->{Rname} eq "Adapter";


#   TLXL: foreach my $i (0..$#{$tlxls}) {

#     my $B_tlxl;
#     my $tlxl;

#     my $tlx = {};

#     switch ($#{$tlxls}) {
#       case 0 {
#         $B_tlxl = $tlxls->[0];
#         $B_tlxl->{tlx} = $tlx;
#       }
#       case 1 {
#         next TLXL if $i == 0;
#         $B_tlxl = $tlxls->[0];
#         $tlxl = $tlxls->[1];
#         # next TLXL if (defined $tlxl->{R1_Rgap} && $tlxl->{R1_Rgap} >=0 && $tlxl->{R1_Rgap} < 10) || 
#         #   (defined $tlxl->{R2_Rgap} && $tlxl->{R2_Rgap} >=0 && $tlxl->{R2_Rgap} < 10 );
#         $tlxl->{tlx} = $tlx;
#       }
#       else {
#         next TLXL if $i == 0;
#         $B_tlxl = $tlxls->[$i-1];
#         $tlxl = $tlxls->[$i];
#         last TLXL if $tlxl->{Rname} eq "Adapter";
#         # next TLXL if defined $tlxl->{R1_Rgap} && $tlxl->{R1_Rgap} >=0 && $tlxl->{R1_Rgap} < 10 ||
#         #   (defined $tlxl->{R2_Rgap} && $tlxl->{R2_Rgap} >=0 && $tlxl->{R2_Rgap} < 10 );
#         $tlxl->{tlx} = $tlx;
#       }
#     }

#     $tlx->{Qname} = $B_tlxl->{Qname};
#     $tlx->{Seq} = $Qseq;
#     $tlx->{Qlen} = $Qlen;

#     $tlx->{B_Rname} = $B_tlxl->{Rname};
#     $tlx->{B_Strand} = $B_tlxl->{Strand};

#     my $B_R1_Qstart;
#     my $B_R2_Qstart;
#     if (defined $B_tlxl->{R1_Qstart} && defined $R1_map[$B_tlxl->{R1_Qstart}-1]) {
#       $B_R1_Qstart = $R1_map[$B_tlxl->{R1_Qstart}-1];
#     }
#     if (defined $B_tlxl->{R2_Qstart} && defined $R2_map[$B_tlxl->{R2_Qstart}-1]) {
#       $B_R2_Qstart = $R2_map[$B_tlxl->{R2_Qstart}-1];
#     }
#     if (defined $B_R1_Qstart && defined $B_R2_Qstart) {
#       if ($B_R1_Qstart <= $B_R2_Qstart) {
#         $B_R2_Qstart = undef;
#       } else {
#         $B_R1_Qstart = undef;
#       }
#     }
#     if (defined $B_R1_Qstart) {
#       $tlx->{B_Qstart} = $R1_map[$B_tlxl->{R1_Qstart}-1];
#       if ($tlx->{B_Strand} == 1) {
#         $tlx->{B_Rstart} = $B_tlxl->{R1_Rstart};
#       } else {
#         $tlx->{B_Rend} = $B_tlxl->{R1_Rend};
#       }
#     } else {
#       $tlx->{B_Qstart} = $R2_map[$B_tlxl->{R2_Qstart}-1];
#       if ($tlx->{B_Strand} == 1) {
#         $tlx->{B_Rstart} = $B_tlxl->{R2_Rstart};
#       } else {
#         $tlx->{B_Rend} = $B_tlxl->{R2_Rend};
#       }
#     }

#     my $B_R1_Qend;
#     my $B_R2_Qend;
#     if (defined $B_tlxl->{R1_Qend} && defined $R1_map[$B_tlxl->{R1_Qend}-1]) {
#       $B_R1_Qend = $R1_map[$B_tlxl->{R1_Qend}-1];
#     }
#     if (defined $B_tlxl->{R2_Qend} && defined $R2_map[$B_tlxl->{R2_Qend}-1]) {
#       $B_R2_Qend = $R2_map[$B_tlxl->{R2_Qend}-1];
#     }
#     if (defined $B_R1_Qend && defined $B_R2_Qend) {
#       if ($B_R1_Qend >= $B_R2_Qend) {
#         $B_R2_Qend = undef;
#       } else {
#         $B_R1_Qend = undef;
#       }
#     }
#     if (defined $B_R1_Qend) {
#       $tlx->{B_Qend} = $R1_map[$B_tlxl->{R1_Qend}-1];
#       if ($tlx->{B_Strand} == 1) {
#         $tlx->{B_Rend} = $B_tlxl->{R1_Rend};
#       } else {
#         $tlx->{B_Rstart} = $B_tlxl->{R1_Rstart};
#       }
#     } else {
#       $tlx->{B_Qend} = $R2_map[$B_tlxl->{R2_Qend}-1];
#       if ($tlx->{B_Strand} == 1) {
#         $tlx->{B_Rend} = $B_tlxl->{R2_Rend};
#       } else {
#         $tlx->{B_Rstart} = $B_tlxl->{R2_Rstart};
#       }
#     }


#     next TLXL unless defined $tlxl;

#     $tlx->{Rname} = $tlxl->{Rname};
#     $tlx->{Strand} = $tlxl->{Strand};

#     my $R1_Qstart;
#     my $R2_Qstart;
#     if (defined $tlxl->{R1_Qstart} && defined $R1_map[$tlxl->{R1_Qstart}-1]) {
#       $R1_Qstart = $R1_map[$tlxl->{R1_Qstart}-1];
#     }
#     if (defined $tlxl->{R2_Qstart} && defined $R2_map[$tlxl->{R2_Qstart}-1]) {
#       $R2_Qstart = $R2_map[$tlxl->{R2_Qstart}-1];
#     }
#     if (defined $R1_Qstart && defined $R2_Qstart) {
#       if ($R1_Qstart <= $R2_Qstart) {
#         $R2_Qstart = undef;
#       } else {
#         $R1_Qstart = undef;
#       }
#     }
#     if (defined $R1_Qstart) {
#       $tlx->{Qstart} = $R1_map[$tlxl->{R1_Qstart}-1];
#       if ($tlx->{Strand} == 1) {
#         $tlx->{Rstart} = $tlxl->{R1_Rstart};
#       } else {
#         $tlx->{Rend} = $tlxl->{R1_Rend};
#       }
#     } else {
#       $tlx->{Qstart} = $R2_map[$tlxl->{R2_Qstart}-1];
#       if ($tlx->{Strand} == 1) {
#         $tlx->{Rstart} = $tlxl->{R2_Rstart};
#       } else {
#         $tlx->{Rend} = $tlxl->{R2_Rend};
#       }
#     }

#     my $R1_Qend;
#     my $R2_Qend;
#     if (defined $tlxl->{R1_Qend} && defined $R1_map[$tlxl->{R1_Qend}-1]) {
#       $R1_Qend = $R1_map[$tlxl->{R1_Qend}-1];
#     }
#     if (defined $tlxl->{R2_Qend} && defined $R2_map[$tlxl->{R2_Qend}-1]) {
#       $R2_Qend = $R2_map[$tlxl->{R2_Qend}-1];
#     }
#     if (defined $R1_Qend && defined $R2_Qend) {
#       if ($R1_Qend >= $R2_Qend) {
#         $R2_Qend = undef;
#       } else {
#         $R1_Qend = undef;
#       }
#     }
#     if (defined $R1_Qend) {
#       $tlx->{Qend} = $R1_map[$tlxl->{R1_Qend}-1];
#       if ($tlx->{Strand} == 1) {
#         $tlx->{Rend} = $tlxl->{R1_Rend};
#       } else {
#         $tlx->{Rstart} = $tlxl->{R1_Rstart};
#       }
#     } else {
#       $tlx->{Qend} = $R2_map[$tlxl->{R2_Qend}-1];
#       if ($tlx->{Strand} == 1) {
#         $tlx->{Rend} = $tlxl->{R2_Rend};
#       } else {
#         $tlx->{Rstart} = $tlxl->{R2_Rstart};
#       }
#     }

#     $tlx->{Junction} = $tlx->{Strand} == 1 ? $tlx->{Rstart} : $tlx->{Rend};
#     my $ref;
#     switch ($tlx->{Rname}) {
#       case "Breaksite" { $ref = $refs->{brk}; }
#       case "Adapter" { $ref = $refs->{adpt}; }
#       else { $ref = $refs->{genome}; }
#     }
#     $tlx->{J_Seq} = $tlx->{Strand} == 1 ? $ref->seq($tlx->{Rname},$tlx->{Rstart}-10,$tlx->{Rstart}+9) :
#                                             $ref->seq($tlx->{Rname},$tlx->{Rend}-9,$tlx->{Rend}+10);



#   }

#   return(1,$n_juncs);

# }


sub nearby_quals ($$$) {
  my $qual_ref = shift;
  my $pos = shift;
  my $mar = shift;

  return mean(@$qual_ref[max(0,$pos-$mar)..min($pos+$mar,$#{$qual_ref})]);
}




1;
