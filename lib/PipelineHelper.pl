use strict;
use warnings;



sub manage_blast_options ($$) {

  my $defaultoptstr = shift;
  my $useroptstr = shift;

  return $defaultoptstr unless $useroptstr =~ /\S/;
  return $useroptstr unless $defaultoptstr =~ /\S/;

  my @defaultopt = split(/\s+-/,$defaultoptstr);
  my @useropt = split(/\s+-/,$useroptstr);

  my %opt_hash;
  my $optkey;
  my $optval;

  my @return_opt = ();

  my $i = 1;
  while ($i < scalar @defaultopt) {
    if ($defaultopt[$i] =~ /^\d/) {
      $defaultopt[$i-1] .= " -".$defaultopt[$i];
      splice(@defaultopt,$i,1);
      next;
    }
    $i++;
  }
  $i = 1;
  while ($i < scalar @useropt) {
    if (defined $useropt[$i] && $useropt[$i] =~ /^\d/) {
      $useropt[$i-1] .= " -".$useropt[$i];
      splice(@useropt,$i,1);
      next;
    }
    $i++;
  }

  foreach (@defaultopt) {
    if ( /^-?(\S+)\s+(\S+)$/ ) {
      $opt_hash{$1} = $2;
    } elsif ( /^-?(\S+)$/ ) {
      $opt_hash{$1} = "";
    } else {
      croak "Error: incorrect default program option string $defaultoptstr";
    }
  }

  foreach (@useropt) {
    if ( /^-?(\S+)\s+(\S+)$/ ) {
      $opt_hash{$1} = $2;
      delete $opt_hash{$1} if $2 eq "off";
    } elsif ( /^-?(\S+)$/ ) {
      $opt_hash{$1} = "";
    } else {
      croak "Error: incorrect program option string $useroptstr";
    }
  }

  foreach $optkey (sort keys %opt_hash) {
    $optval = $opt_hash{$optkey};
    if ( $optval eq "" ) {
      push(@return_opt,"-$optkey");
    } else {
      push(@return_opt,"-$optkey $optval");
    }
  }

  return(join(" ",@return_opt));


}



sub manage_blat_options ($$) {

  my $defaultoptstr = shift;
  my $useroptstr = shift;

  return $defaultoptstr unless $useroptstr =~ /\S/;
  return $useroptstr unless $defaultoptstr =~ /\S/;

  my @defaultopt = split(/\s+-/,$defaultoptstr);
  my @useropt = split(/\s+-/,$useroptstr);

  my %opt_hash;
  my $optkey;
  my $optval;

  my @return_opt = ();

  my $i = 1;
  while ($i < scalar @defaultopt) {
    if ($defaultopt[$i] =~ /^\d/) {
      $defaultopt[$i-1] .= " -".$defaultopt[$i];
      splice(@defaultopt,$i,1);
      next;
    }
    $i++;
  }
  $i = 1;
  while ($i < scalar @useropt) {
    if (defined $useropt[$i] && $useropt[$i] =~ /^\d/) {
      $useropt[$i-1] .= " -".$useropt[$i];
      splice(@useropt,$i,1);
      next;
    }
    $i++;
  }

  foreach (@defaultopt) {
    if ( /^-?(\S+)=(\S+)$/ ) {
      $opt_hash{$1} = $2;
    } elsif ( /^-?(\S+)$/ ) {
      $opt_hash{$1} = "";
    } else {
      croak "Error: incorrect default program option string $defaultoptstr";
    }
  }

  foreach (@useropt) {
    if ( /^-?(\S+)=(\S+)$/ ) {
      $opt_hash{$1} = $2;
      delete $opt_hash{$1} if $2 eq "off";
    } elsif ( /^-?(\S+)$/ ) {
      $opt_hash{$1} = "";
    } else {
      croak "Error: incorrect program option string $useroptstr";
    }
  }

  foreach $optkey (sort keys %opt_hash) {
    $optval = $opt_hash{$optkey};
    if ( $optval eq "" ) {
      push(@return_opt,"-$optkey");
    } else {
      push(@return_opt,"-$optkey $optval");
    }
  }

  return(join(" ",@return_opt));


}

sub prepare_reference_genomes ($$) {
  # First argument is local genome database file path
  # which should contain subdirectories for each assembly
  # Second argument is a reference to a hash of hashes
  # $meta_ref->{expt_id}->{assembly}
  #                     ->{mask}

  my $dir = shift;
  my $meta_ref = shift;
  my @masks = ();


  foreach my $expt_id (sort keys %$meta_ref) {
    my $assembly = $meta_ref->{$expt_id}->{assembly};
    my $assemdir = "$dir/$assembly";
    croak "Error: could not find reference genome directory $dir/$assembly" unless (-d $assemdir);
    my $clean2bit= "$assemdir/$assembly.2bit";
    my $cleanFa = "$assemdir/$assembly.fa";

    my $mask = $meta_ref->{$expt_id}->{mask};
    unless ($mask =~ /\S/) {
      $meta_ref->{$expt_id}->{refdb} = $clean2bit;
      next;
    }
    (my $mask_stub = $assembly."_mask_".$mask) =~ s/[:,\s]+/_/g;
    my $maskFa = "$assemdir/$mask_stub.fa";
    my $mask2bit = "$assemdir/$mask_stub.2bit";
    $meta_ref->{$expt_id}->{refdb} = $mask2bit;

    next if (-r $mask2bit);

    my $beddir = "$assemdir/maskBED";
    my $bedfile = "$beddir/$mask_stub.bed"; 

    unless (-r $clean2bit) {
      if (-r $cleanFa) {
        System("faToTwoBit $cleanFa $clean2bit");
      } else {
        System("faToTwoBit $assemdir/chr*.fa $clean2bit");
      }
    }

    

    unless (-r $bedfile) {
      unless (-d $beddir) {
        mkdir $beddir or croak "Error: cannot create maskBED directory";
      }
      my $bed_fh = IO::File->new(">$bedfile");
      my @loci = split( /\s*,\s*/ , $mask );
      foreach my $locus (@loci) {
        (my $chr, my $start, my $end) = ($locus =~ /(chr\w+):(\d+)-(\d+)/);
        $bed_fh->print(join("\t",$chr,$start,$end)."\n");
      }
      $bed_fh->close;
    }



    unless (-r $mask2bit) {
      unless (-r $cleanFa) {
        System("twoBitToFa $clean2bit $cleanFa");
      }
      System("maskFastaFromBed -fi $cleanFa -bed $bedfile -fo $maskFa");
      System("faToTwoBit $maskFa $mask2bit");
      System("rm $maskFa");
    }
  }

}

# Generates red and blu primer+continued fasta files in the expt/sequences/ directory
# Updates meta data hash to contain file locations
sub create_sequence_files ($$) {
  # First argument is experiment ID
  # Second argument contains a hash reference with the following keys:
  # exptdir (where files are to be output), mid (mid sequence),
  # redp (red primer sequence), redc (red cont sequence), blup, bluc

  my $expt_id = shift;
  my $expt_hash = shift;

  my $seqdir = $expt_hash->{exptdir} . "/sequences";

  unless (-d $seqdir) {
    mkdir $seqdir or croak "Error: could not create sequences directory for $expt_id";
  } 
  $expt_hash->{redfa} = "$seqdir/red.fa";
  $expt_hash->{redpfa} = "$seqdir/redp.fa";
  $expt_hash->{blufa} = "$seqdir/blu.fa";
  $expt_hash->{blupfa} = "$seqdir/blup.fa";

  my $redio = IO::File->new(">".$expt_hash->{redfa}) or croak "Error: could not write to red primer fasta file";
  $redio->print("> $expt_id MID + Red Primer + Red Continued\n");
  $redio->print(uc( $expt_hash->{mid} . $expt_hash->{redp} . $expt_hash->{redc} ));
  $redio->close;

  my $redpio = IO::File->new(">".$expt_hash->{redpfa}) or croak "Error: could write to red primer fasta file";
  print $redpio "> $expt_id Red Primer\n";
  print $redpio uc($expt_hash->{redp});
  $redpio->close;

  my $bluio = IO::File->new(">".$expt_hash->{blufa}) or croak "Error: could not write to blu primer fasta file";
  $bluio->print("> $expt_id Blu Continued + Blu Primer\n");
  $bluio->print(uc( $expt_hash->{bluc} . $expt_hash->{blup} ));
  $bluio->close;

  my $blupio = IO::File->new(">".$expt_hash->{blupfa}) or croak "Error: could write to blu primer fasta file";
  print $blupio "> $expt_id Blu Primer\n";
  print $blupio uc($expt_hash->{blup});
  $blupio->close;

}

sub blast_to_sequence_files ($$$) {
  my $expt_id = shift;
  my $expt_hash = shift;
  my $blastopt = shift;


  my $aligndir = $expt_hash->{exptdir} . "/alignments";
  unless (-d $aligndir) {
    mkdir $aligndir or croak "Error: could not create alignments directory for $expt_id";
  }

  $expt_hash->{redblast} = "$aligndir/red_blast.txt";
  $expt_hash->{redpblast} = "$aligndir/redp_blast.txt";
  $expt_hash->{blublast} = "$aligndir/blu_blast.txt";
  $expt_hash->{blupblast} = "$aligndir/blup_blast.txt";

  System(join(" ","blastn","-subject",$expt_hash->{redfa},"-query",$expt_hash->{raw},"-out",$expt_hash->{redblast},$blastopt))
                  or croak "Error: $expt_id failed during blast to red sequence";

  System(join(" ","blastn","-subject",$expt_hash->{blufa},"-query",$expt_hash->{raw},"-out",$expt_hash->{blublast},$blastopt))
                  or croak "Error: $expt_id failed during blast to blu sequence";

  System(join(" ","blastn","-subject",$expt_hash->{redpfa},"-query",$expt_hash->{raw},"-out",$expt_hash->{redpblast},$blastopt))
                  or croak "Error: $expt_id failed during blast to red sequence";

  System(join(" ","blastn","-subject",$expt_hash->{blupfa},"-query",$expt_hash->{raw},"-out",$expt_hash->{blupblast},$blastopt))
                  or croak "Error: $expt_id failed during blast to blu sequence";

}



sub align_to_genome ($$$) {

  my $expt_id = shift;
  my $expt_hash = shift;
  my $blatopt = shift;

  my $aligndir = $expt_hash->{exptdir} . "/alignments";

  unless (-d $aligndir) {
    mkdir $aligndir or croak "Error: could not create alignments directory for $expt_id";
  }

  $expt_hash->{psl} = "$aligndir/$expt_id.psl";

  System(join(" ","blat",$expt_hash->{refdb},$expt_hash->{raw},$expt_hash->{psl},$blatopt))
                  or croak "Error: $expt_id failed during blat to genome";
  compressHeader($expt_hash->{psl},1);


}


sub make_tlxl ($$) {

  my $expt_id = shift;
  my $expt_hash = shift;
  

  my %query;



  # Build raw sequence index
  # Open raw file
  my $rawio = IO::File->new("<".$expt_hash->{raw});
  # Create index
  $expt_hash->{rawidx} = $expt_hash->{raw} . ".idx";
  my $rawidxio = IO::File->new("+>".$expt_hash->{rawidx});
  build_index($rawio,$rawidxio);
  # Then read line numbers for each read into the query hash
  seek($rawio,0,0);
  my $line_no = 1;
  while ( (my $seq_id, my $seq) = read_fasta($rawio) ) {
    (my $qname = $seq_id) =~ s/^\s*(\S+)\s.*$/$1/;
    croak "Error: multiple raw sequences with same ID in $expt_id" if defined $query{$qname};
    $query{$qname}->{rawidx} = $line_no + 1;
    $line_no += 2;
  }

  # Build red sequence alignment index
  # Open the red alignment file
  my $redio = IO::File->new("<".$expt_hash->{redblast});
  # Create index file
  $expt_hash->{redidx} = $expt_hash->{redblast} . ".idx";
  my $redidxio = IO::File->new("+>".$expt_hash->{redidx});
  build_index($redio,$redidxio);
  # Then read line numbers for each read into the query hash
  seek($redio,0,0);
  my $redcsv = Text::CSV->new({sep_char => "\t"});
  $redcsv->column_names(blast_header());
  $line_no = 0;
  while ( my $row = $redcsv->getline_hr($redio) ) {
    $line_no += 1;
    my $qname = $row->{qseqid};
    unless (defined $query{$qname}->{redidx}) {
      $query{$qname}->{redidx} = $line_no;
    } else {
      $query{$qname}->{redidx} .= ",$line_no";
    }
  }

  # Build red primer alignment index
  # Open red primer alignment file
  my $redpio = IO::File->new("<".$expt_hash->{redpblast});
  # Create index file
  $expt_hash->{redpidx} = $expt_hash->{redpblast} . ".idx";
  my $redpidxio = IO::File->new("+>".$expt_hash->{redpidx});
  build_index($redpio,$redpidxio);
  # Then read line numbers for each read into the query hash
  seek($redpio,0,0);
  my $redpcsv = Text::CSV->new({sep_char => "\t"});
  $redpcsv->column_names(blast_header());
  $line_no = 0;
  while ( my $row = $redpcsv->getline_hr($redpio) ) {
    $line_no += 1;
    my $qname = $row->{qseqid};
    unless (defined $query{$qname}->{redpidx}) {
      $query{$qname}->{redpidx} = $line_no;
    } else {
      $query{$qname}->{redpidx} .= ",$line_no";
    }
  }

  # Build blu sequence alignment index
  # Open blu alignment file
  my $bluio = IO::File->new("<".$expt_hash->{blublast});
  # Create index file
  $expt_hash->{bluidx} = $expt_hash->{blublast} . ".idx";
  my $bluidxio = IO::File->new("+>".$expt_hash->{bluidx});
  build_index($bluio,$bluidxio);
  # Then read line numbers for each read into the query hash
  seek($bluio,0,0);
  my $blucsv = Text::CSV->new({sep_char => "\t"});
  $blucsv->column_names(blast_header());
  $line_no = 0;
  while ( my $row = $blucsv->getline_hr($bluio) ) {
    $line_no += 1;
    my $qname = $row->{qseqid};
    unless (defined $query{$qname}->{bluidx}) {
      $query{$qname}->{bluidx} = $line_no;
    } else {
      $query{$qname}->{bluidx} .= ",$line_no";
    }
  }

  # Build blu primer alignment index
  # Open blu primer file
  my $blupio = IO::File->new("<".$expt_hash->{blupblast});
  # Create index
  $expt_hash->{blupidx} = $expt_hash->{blupblast} . ".idx";
  my $blupidxio = IO::File->new("+>".$expt_hash->{blupidx});
  build_index($blupio,$blupidxio);
  # Then read line numbers for each read into the query hash
  seek($blupio,0,0);
  my $blupcsv = Text::CSV->new({sep_char => "\t"});
  $blupcsv->column_names(blast_header());
  $line_no = 0;
  while ( my $row = $blupcsv->getline_hr($blupio) ) {
    $line_no += 1;
    my $qname = $row->{qseqid};
    unless (defined $query{$qname}->{blupidx}) {
      $query{$qname}->{blupidx} = $line_no;
    } else {
      $query{$qname}->{blupidx} .= ",$line_no";
    }
  }

  # Create bedfile to extract reference sequence around junctions
  # Open alignment file
  my $pslio = IO::File->new("<".$expt_hash->{psl});
  my $pslcsv = Text::CSV->new({sep_char => "\t"});
  my $header = $pslcsv->getline($pslio);
  $pslcsv->column_names(@$header);
  # Create bedfile
  (my $bedfile = $expt_hash->{psl}) =~ s/\.psl$/_marg.bed/;
  my $bedio = IO::File->new(">$bedfile");
  # Read through pslx line by line writing margin coordinates and Qname to the bedfile
  while (my $psl = $pslcsv->getline_hr($pslio)) {
    my $qname = $psl->{Qname};
    my $tname = $psl->{Tname};
    my $strand = $psl->{strand};
    my $tstart = $psl->{Tstart};
    my $tend = $psl->{Tend};
    my $junction = $strand eq "+" ? $tstart : $tend;
    my $rev_marg = 10;
    my $for_marg = 6;
    if ($strand eq "+") {
      $bedio->print(join("\t",$tname,$junction-$rev_marg,$junction+$for_marg,"$qname:$tname:$tstart-$tend")."\n");
    } else {
      $bedio->print(join("\t",$tname,$junction-$for_marg,$junction+$rev_marg,"$qname:$tname:$tstart-$tend")."\n");
    }
  }
  # Extract sequences using twoBitToFa
  # When extracting large number of sequences,
  # twoBitToFa is much faster when called once from a bedfile
  # instead of once for each sequence
  ($expt_hash->{margfa} = $bedfile) =~ s/\.bed$/.fa/;
  System("twoBitToFa ".$expt_hash->{refdb}." ".$expt_hash->{margfa}." -bed=".$bedfile );
  ## Build margin index
  ## Open margin file
  #my $margio = IO::File->new("<".$expt_hash->{margfa});
  ## Create index
  #$expt_hash->{margidx} = $expt_hash->{margfa} . ".idx";
  #my $margidxio = IO::File->new(">".$expt_hash->{margidx});
  #build_index($margio,$margidxio);
  #$margidxio->close;
  ## Then read line numbers for each read into the query hash
  #seek($margio,0,0);
  #my $line_no = 1;
  #while ( (my $seqid, my $seq) = read_fasta($margio) ) {
  # my ($qname, $tname, $tstart, $tend) = ($seqid =~ /(\S+):(\w+):(\d+)-(\d+)/)
  # unless (defined $query{$qname}->{margidx}) {
  #   $query{$qname}->{margidx} = $line_no + 1;
  # } else {
  #   $query{$qname}->{margidx} .= ",".$line_no+1;
  # }
  # $line_no += 2;
  #}





  # Prepare tlxl file
  $expt_hash->{tlxl} = $expt_hash->{exptdir} . "/" . $expt_id . ".tlxl";
  my $tlxlio = IO::File->new(">".$expt_hash->{tlxl});
  my @tlxlheader = qw( Qname Tname Tstart Tend Strand Qstart Qend Qsize Match Mismatch
    blockSizes tStarts qStarts RedQstart RedQend RedpQstart RedpQend
    BluQstart BluQend BlupQstart BlupQend Margin Raw);
  $tlxlio->print(join("\t",@tlxlheader)."\n");

  # Finally read in genome alignments line by line
  seek($pslio,0,0);
  $header = $pslcsv->getline($pslio);

  my $margio = IO::File->new("<".$expt_hash->{margfa});


  while (my $psl = $pslcsv->getline_hr($pslio)) {
    my $qname = $psl->{Qname};
    

    # Retrieve sequence
    my $raw = line_with_index($rawio,$rawidxio,$query{$qname}->{rawidx});

    # Retrieve reference sequence around junction margin
    my ($seqid,$margin) = read_fasta($margio);
    my ($red_aln,$redp_aln,$blu_aln,$blup_aln);
    # Retrieve red sequence alignments
    if (defined $query{$qname}->{redidx}) {
      my @redidxs = split(",",$query{$qname}->{redidx});
      foreach my $line (@redidxs) {
        my $red_curr = line_with_index_hr($redio,$redcsv,$redidxio,$line);
        unless (defined $red_aln) {
          $red_aln = $red_curr;
        } else {
          $red_aln = $red_curr if $red_curr->{qstart} < $red_aln->{qstart};
        }
      }
    }
    # Retrieve red primer alignments
    if (defined $query{$qname}->{redpidx}) {
      my @redpidxs = split(",",$query{$qname}->{redpidx});
      foreach my $line (@redpidxs) {
        my $redp_curr = line_with_index_hr($redpio,$redpcsv,$redpidxio,$line);
        unless (defined $redp_aln) {
          $redp_aln = $redp_curr;
        } else {
          $redp_aln = $redp_curr if $redp_curr->{qstart} < $redp_aln->{qstart};
        }
      }
    }
    # Retrieve blu sequence alignments
    if (defined $query{$qname}->{bluidx}) {
      my @bluidxs = split(",",$query{$qname}->{bluidx});
      foreach my $line (@bluidxs) {
        my $blu_curr = line_with_index_hr($bluio,$blucsv,$bluidxio,$line);
        unless (defined $blu_aln) {
          $blu_aln = $blu_curr;
        } else {
          $blu_aln = $blu_curr if $blu_curr->{qend} > $blu_aln->{qend};
        }
      }
    }
    # Retrieve blu primer alignments
    if (defined $query{$qname}->{blupidx}) {
      my @blupidxs = split(",",$query{$qname}->{blupidx});
      foreach my $line (@blupidxs) {
        my $blup_curr = line_with_index_hr($blupio,$blupcsv,$blupidxio,$line);
        unless (defined $blup_aln) {
          $blup_aln = $blup_curr;
        } else {
          $blup_aln = $blup_curr if $blup_curr->{qend} > $blup_aln->{qend};
        }
      }
    }

    unless (defined $red_aln) {
      $red_aln = {qstart => "", qend => ""};
    }
    unless (defined $redp_aln) {
      $redp_aln = {qstart => "", qend => ""};
    }
    unless (defined $blu_aln) {
      $blu_aln = {qstart => "", qend => ""};
    }
    unless (defined $blup_aln) {
      $blup_aln = {qstart => "", qend => ""};
    }

    $tlxlio->print(join("\t", $psl->{Qname} ,
                              $psl->{Tname} ,
                              $psl->{Tstart} ,
                              $psl->{Tend} ,
                              $psl->{strand} ,
                              $psl->{Qstart} ,
                              $psl->{Qend} ,
                              $psl->{Qsize} ,
                              $psl->{match} ,
                              $psl->{'mis-match'} ,
                              $psl->{blockSizes} ,
                              $psl->{tStarts} ,
                              $psl->{qStarts} ,
                              $red_aln->{qstart} ,
                              $red_aln->{qend} ,
                              $redp_aln->{qstart} ,
                              $redp_aln->{qend} ,
                              $blu_aln->{qstart} ,
                              $blu_aln->{qend} ,
                              $blup_aln->{qstart} ,
                              $blup_aln->{qend} ,
                              $margin,
                              $raw ));
  }
}




1;
