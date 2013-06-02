#!/usr/bin/perl


use strict;
use warnings;
use Getopt::Long;
use Carp;
use IO::Handle;
use IO::File;
use Text::CSV;
use File::Basename;
use File::Which;
use Bio::PrimarySeq;
use Bio::Factory::EMBOSS;
use Bio::Tools::Run::StandAloneBlastPlus;
use Bio::AlignIO;
use Bio::SeqIO;
use List::Util qw(min max);
use Interpolation 'arg:@->$' => \&argument;
use Time::HiRes qw(gettimeofday tv_interval);

use Cwd qw(abs_path);
use FindBin;
use lib abs_path("$FindBin::Bin/../lib");

require "PerlSub.pl";
require "TranslocHelper.pl";


# Flush output after every write
select( (select(STDOUT), $| = 1 )[0] );

##
## This program
## 
## run with "--help" for usage information
##
## Robin Meyers

# Forward declarations
sub parse_command_line;
sub read_in_translocations;
sub align_to_primers;
sub print_html_header ($);
sub print_html_footer ($);
sub write_marked_read ($$);


# Global flags and arguments, 
# Set by command line arguments
my $tlxfile;
my $htmlfile;
my $for_primer = "";
my $rev_primer = "";


# Global variables
my %tlx;
my @qids;
my @qseqs;

#
# Start of Program
#

parse_command_line;

my $t0 = [gettimeofday];

read_in_translocations;

align_to_primers;


my $htmlfh = IO::File->new(">$htmlfile") or croak "Error: cannot write to html file $htmlfile";
print_html_header($htmlfh);

foreach my $qid (@qids) {

  write_marked_read($htmlfh,$tlx{$qid});

}



print_html_footer($htmlfh);
$htmlfh->close;



my $t1 = tv_interval($t0);

printf("\nFinished all processes in %.2f seconds.\n", $t1);



#
# End of program
#

sub write_marked_read ($$) {
  my $fh = shift;
  my $tl = shift;

  my @seq = split("",$tl->{Seq});
  my %marks;
  my $marktypes = {alignment => {start => "Qstart", end => "Qend"},
                   breaksite => {start => "BrkQstart", end => "BrkQend"},
                   forward_primer => {start => "ForQstart", end => "ForQend"},
                   reverse_primer => {start => "RevQstart", end => "RevQend"}};

  foreach my $marktype (sort keys %$marktypes) {

    if (exists $tl->{$marktypes->{$marktype}->{start}}) {

      my $start = $tl->{$marktypes->{$marktype}->{start}} - 1;
      my $end = $tl->{$marktypes->{$marktype}->{end}} - 1;

      if (exists $marks{$start}) {
        push(@{$marks{$start}->{start}},$marktype);
      } else {
        $marks{$start}->{start} = [$marktype];
      }

      if (exists $marks{$end}) {
        push(@{$marks{$end}->{end}},$marktype);
      } else {
        $marks{$end}->{end} = [$marktype];
      }
    }
  }

  my %open_marks;
  my $marked_read = "<div class=\"sequence\">";

  foreach my $i (0..$#seq) {

    if (exists $marks{$i}->{start}) {

      my @marks = @{$marks{$i}->{start}};

      $marked_read .= "</span>" if keys %open_marks > 0;
      $marked_read .= "<span class=\"" . join(" ",@marks) . 
                    (keys %open_marks > 0 ? " ".join(" ",keys %open_marks) : "") . "\">";

      foreach my $mark (@marks) {
        $open_marks{$mark} = 1;
      }
    }

    $marked_read .= $seq[$i];

    if (exists $marks{$i}->{end}) {
      my @marks = @{$marks{$i}->{end}};

      $marked_read .= "</span>";

      foreach my $mark (@marks) {
        delete $open_marks{$mark};
      }

      $marked_read .= "<span class=\"" . join(" ",keys %open_marks) . "\">" if keys %open_marks > 0;

    }

  }

  $marked_read .= "</span>" if keys %open_marks > 0;

  $marked_read .= "</div>";

  my $qname = $tl->{Qname};
  my $rname = $tl->{Rname};
  my $rstart = $tl->{Rstart};
  my $rend = $tl->{Rend};
  my $strand = $tl->{Strand};

  my $marked_id = "<div class=\"seqname\">>$qname $rname:$rstart-$rend $strand</div>";
  $fh->print("$marked_id\n");
  $fh->print("$marked_read\n");


}

sub read_in_translocations {
  my $tlxfh = IO::File->new("<$tlxfile");
  my $csv = Text::CSV->new({sep_char => "\t"});
  my $header = $csv->getline($tlxfh);
  $csv->column_names(@$header);

  my $ct = 0;

  while (my $tl = $csv->getline_hr($tlxfh)) {

    (my $qname = $tl->{Qname}) =~ s/\W/_/g;
    my $qid = "q_".++$ct;

    $tlx{$qid} = $tl;

    push(@qids,$qid);

    push(@qseqs,Bio::PrimarySeq->new( -id => $qid,
                                     -seq => $tl->{Seq} ));

  }

  $tlxfh->close;

}

sub align_to_primers {

  my $f = Bio::Factory::EMBOSS -> new();
  my $water = $f->program('water');




  unless ($for_primer eq "") {
    my $for_seq = Bio::PrimarySeq->new( -id => "forward_primer", -seq => $for_primer);

    # my $fac = Bio::Tools::Run::StandAloneBlastPlus->new( -subject => $for_seq);

    # my $result = $fac->bl2seq( -method => 'blastn',
    #                            -query => \@qseqs,
    #                            -subject => $for_seq,
    #                            -method_args => [ -word_size => 10, -max_target_seqs => 1 ] );

    # while ( my $hit)

    (my $for_water = $tlxfile) =~ s/\.tlx/_for.water/;
    $water->run({ -asequence => $for_seq,
                  -bsequence    => \@qseqs,
                  -gapopen   => '10.0',
                  -gapextend => '0.5',
                  -outfile   => $for_water});
    my $alnio = Bio::AlignIO->new( -format => 'emboss',
                                   -file   => $for_water);
    while (my $for_aln = $alnio->next_aln) {

      next unless $for_aln->percentage_identity > 90 && $for_aln->length >= length($for_primer) - 3;

      my $qseq = $for_aln->get_seq_by_pos(2);
      my $qid = $qseq->id;

      next if exists $tlx{$qid}->{ForQstart};

      $tlx{$qid}->{ForQstart} = $qseq->start;
      $tlx{$qid}->{ForQend} = $qseq->end;
    }


    # (my $for_water = $tlxfile) =~ s/\.tlx/_for.water/;
    # $water->run({ -asequence => $for_seq,
    #               -bsequence    => \@qseqs,
    #               -gapopen   => '10.0',
    #               -gapextend => '0.5',
    #               -outfile   => $for_water,
    #               -aformat3 => 'fasta'});

    

    # my $for_aln = Bio::SeqIO->new( -format => 'fasta',
    #                                -file   => $for_water);
    # while ( my $seq = $for_aln->next_seq ) {
    #   my $qname = $seq->id;
    #   next if $qname =~ /primer/;
    #   (my $subseq = $seq->seq) =~ s/-//g;
      
    #   next if $subseq ne $for_seq->seq && levenshtein($subseq,$for_seq->seq) < 3;

    #   my $start = index($tlx{$qname}->{Seq},$subseq);
    #   $tlx{$qname}->{ForQstart} = $start + 1;
    #   $tlx{$qname}->{ForQend} = $start + length($subseq);
    # }

    unlink $for_water;

  }

  unless ($rev_primer eq "") {
    my $rev_seq = Bio::PrimarySeq->new( -id => "reverse_primer", -seq => $rev_primer);
    # (my $rev_water = $tlxfile) =~ s/\.tlx/_rev.water/;
    # $water->run({ -asequence => $rev_seq,
    #               -bsequence    => \@qseqs,
    #               -gapopen   => '10.0',
    #               -gapextend => '0.5',
    #               -outfile   => $rev_water,
    #               -aformat3 => 'fasta'});
    # my $rev_aln = Bio::SeqIO->new( -format => 'fasta',
    #                                -file   => $rev_water);
    # while ( my $seq = $rev_aln->next_seq ) {
    #   my $qname = $seq->id;
    #   next if $qname =~ /primer/;
    #   (my $subseq = $seq->seq) =~ s/-//g;

    #   next if $subseq ne $rev_seq->seq && levenshtein($subseq,$rev_seq->seq) < 3;
      
    #   my $start = index($tlx{$qname}->{Seq},$subseq);
    #   $tlx{$qname}->{RevQstart} = $start + 1;
    #   $tlx{$qname}->{RevQend} = $start + length($subseq);
    # }


    (my $rev_water = $tlxfile) =~ s/\.tlx/_rev.water/;
    $water->run({ -asequence => $rev_seq,
                  -bsequence    => \@qseqs,
                  -gapopen   => '10.0',
                  -gapextend => '0.5',
                  -outfile   => $rev_water});
    my $alnio = Bio::AlignIO->new( -format => 'emboss',
                                   -file   => $rev_water);
    while (my $rev_aln = $alnio->next_aln) {

      next unless $rev_aln->percentage_identity > 90 && $rev_aln->length >= length($rev_primer) - 3;

      my $qseq = $rev_aln->get_seq_by_pos(2);
      my $qid = $qseq->id;

      next if exists $tlx{$qid}->{RevQstart};

      $tlx{$qid}->{RevQstart} = $qseq->start;
      $tlx{$qid}->{RevQend} = $qseq->end;
    }

  unlink $rev_water;

  }



}

sub parse_command_line {
  my $help;

  usage() if (scalar @ARGV == 0);

  my $result = GetOptions ( "for=s" => \$for_primer,
                            "rev=s" => \$rev_primer,
                            "help" => \$help
                          ) ;
  
  usage() if ($help);

  $tlxfile = shift(@ARGV);
  $htmlfile = shift(@ARGV);

  #Check options

  croak "Error: cannot read tlx file $tlxfile" unless -r $tlxfile;
  
  exit unless $result;
}

sub print_html_header ($) {
  my $fh = shift;
  print $fh <<EOF
<!DOCTYPE html>
<html>
<head>
  <style type="text/css">
      #sequences {
        word-wrap:break-word;
      }
      .sequence {
        letter-spacing:-1px;
        font-family:courier;
      }
      .forward_primer {
        color:red;
      }
      .reverse_primer {
        color:blue;
      }
      .breaksite {
        text-decoration:underline;
      }
      .alignment {
        background-color:yellow;
      }
  </style>
</head>

<body>
<div id="sequences">
EOF
}

sub print_html_footer ($) {
  my $fh = shift;
  print $fh <<EOF
</div>
</body>
</html>
EOF
}


sub usage()
{
print<<EOF;
Title, by Robin Meyers, ddmonthyyyy

This program .


Usage: $0 arg1 arg2 arg3 ...
        [--option VAL] [--flag] [--help]

Arguments (defaults in parentheses):

$arg{"tlxfile"," "}
$arg{"htmlfile"," "}
$arg{"--for"," ",$for_primer}
$arg{"--rev"," ",$rev_primer}
$arg{"--help","This helpful help screen."}


EOF

exit 1;
}
