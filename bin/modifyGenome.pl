#!/usr/bin/env perl

use strict;
use warnings;
use Getopt::Long;
use Carp;
use File::Basename;
use File::Copy;
use File::List;
use File::Spec;
use IO::Handle;
use IO::File;
use Bio::SeqIO;
use Bio::DB::Fasta;
use Text::CSV;
use Interpolation 'arg:@->$' => \&argument;
use Time::HiRes qw(gettimeofday tv_interval);

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
sub load_genome;
sub read_in_annotation_file;
sub modify_genome;
sub build_bowtie2_index;


# Global variables set by command line arguments
my $genome;
my $annotfile;
my $mod_dir;
my $output;
my $force_single_file;
my $force_multi_file;
my $multi_file;
my $use_links;
my $bt2_index;

# More global variables
my $genome_obj;
my %modifications;
my @additions;
my $outfh;
my $output_tag;

#
# Start of Program
#

parse_command_line;

if ((-d $genome && ! $force_single_file) || $force_multi_file) {
  $multi_file = 1;
  print "Writing modified genome to individual chromosome files...\n";
} else {
  $multi_file = 0;
  print "Writing modified genome to a single file...\n";
}

load_genome;

read_in_annotation_file;

mkdir $output;

$output_tag = basename($output);



modify_genome;

build_bowtie2_index if $bt2_index;

#
# End of program
#

sub modify_genome {
  print "Making modifications and writing to file...\n";

  $outfh = Bio::SeqIO->new(-file => ">$output/$output_tag.fa", -format => 'fasta') unless $multi_file;

  foreach my $chr ($genome_obj->ids) {
    print "Looking for modifications on $chr...\n";


    

    if (exists $modifications{$chr}) {

      $outfh = Bio::SeqIO->new(-file => ">$output/$chr.fa", -format => 'fasta') if $multi_file;
      
      my $seq = $genome_obj->get_Seq_by_id($chr);

      my $temp_seq = $seq->seq;
      my $newseq;


      foreach my $mod (@{$modifications{$chr}}) {
        my $insert;
        if (defined $mod->{file}) {
          my $fh = Bio::SeqIO->new(-file=>"<".$mod_dir."/".$mod->{file},-format=>'fasta');
          my $insert_obj = $fh->next_seq;
          $insert = $insert_obj->seq;
          print join(" ","Replacing",$mod->{chr},$mod->{start},$mod->{end},"with",$insert_obj->id,"\n");
          print join(" ","This operation replaces a",$mod->{end}-$mod->{start},"bp sequence with a",$insert_obj->length,"bp sequence\n");
        } else {
          $insert = "N" x ($mod->{end} - $mod->{start});
          print join(" ","Masking",$mod->{end}-$mod->{start},"bps at",$mod->{chr},$mod->{start},$mod->{end},"\n");
        }
        substr($temp_seq,$mod->{start}-1,$mod->{end}-$mod->{start},$insert);
        $newseq = Bio::Seq->new(-id => $seq->id, -seq => $temp_seq);
      }

      print "made ". scalar @{$modifications{$chr}} ." modifications, writing to output...\n";

      $outfh->write_seq($newseq);
      $outfh->close if $multi_file;

    } else {

      print "No modifications...";

      if ($multi_file) {
        if (-d $genome) {
          if ($use_links) {
            my $relpath = File::Spec->abs2rel($genome,$output);
            print "creating symlink...\n";
            symlink "$relpath/$chr.fa" , "$output/$chr.fa";
            next;
          } else {
            print "copying chromosome file...\n";
            copy "$genome/$chr.fa" , "$output/$chr.fa";
            next;
          }
        } else {
          $outfh = Bio::SeqIO->new(-file => ">$output/$chr.fa", -format => 'fasta');
        }
      }

      my $seq = $genome_obj->get_Seq_by_id($chr);
      print "writing to output file...\n";
      $outfh->write_seq($seq);
      $outfh->close if $multi_file;

    }





  }

  foreach my $add (@additions) {

    my $fh = Bio::SeqIO->new(-file=>"<".$add->{file},-format => 'fasta');
    my $insert = $fh->next_seq;

    $outfh = Bio::SeqIO->new(-file => ">$output/".$insert->id.".fa", -format => 'fasta') if $multi_file;

    print join(" ","Adding",$insert->id,"to output\n");

    $outfh->write_seq($insert);

    $outfh->close if $multi_file;
  }

  $outfh->close unless $multi_file;


}

sub build_bowtie2_index {

  my $index_dir = $ENV{"BOWTIE2_INDEXES"};
  unless (defined $index_dir && $index_dir =~ /\S/) {
    carp "Warning: could not find BOWTIE2_INDEXES directory, writing index to output directory";
    $index_dir = $output;
  }

  my @files;

  if ($multi_file) {
    my $output_search = File::List->new($output);
    @files = @{ $output_search->find(qr/\.fa$/) };

  } else {
    @files = ("$output/$output_tag.fa");
  }

  my $bt2_build_cmd = join(" ", "bowtie2-build",
                                join(",",@files),
                                "$index_dir/$output_tag");
  print $bt2_build_cmd . "\n";
  system($bt2_build_cmd);

}

sub read_in_annotation_file {
  print "Reading in annotation file...\n";

  my $annotfh = IO::File->new("<$annotfile");
  my $csv = Text::CSV->new({sep_char => "\t"});

  while (my $annot = $csv->getline($annotfh)) {
    next unless $annot->[0] =~ /\S/;
    my $mod = {};
    if (grep { $_ eq $annot->[0] } $genome_obj->ids) {
      
      $mod->{chr} = $annot->[0];

      croak "Error: start coordinate is outside range of chromosome"
        unless ($annot->[1] > 0 && $annot->[1] <= $genome_obj->length($mod->{chr}));
      $mod->{start} = $annot->[1];
      
      croak "Error: end coordinate is outside range of chromosome"
        unless ($annot->[2] > 0 && $annot->[1] <= $genome_obj->length($mod->{chr}));
      $mod->{end} = $annot->[2];

      if (defined $annot->[3]) {
        croak "Error: cannot find or read new sequence file"
          unless (-r "$mod_dir/".$annot->[3]);
        $mod->{file} = $annot->[3];
      }

    $modifications{$mod->{chr}} = [] unless exists $modifications{$mod->{chr}};

    push(@{$modifications{$mod->{chr}}},$mod);

    } elsif (-r "$mod_dir/".$annot->[0]) {
      
      $mod->{file} = $annot->[0];

      push(@additions,$mod);

    } else {
      croak "Error: first column neither an existing chromosome nor a new sequence file";
    }


  }

}

sub load_genome {
  print "Loading original genome...\n";
  $genome_obj = Bio::DB::Fasta->new($genome);
  my @genome_ids = $genome_obj->ids;
  print join(" ","Genome composed of",@genome_ids,"\n");
}


sub parse_command_line {
  my $help;

  usage() if (scalar @ARGV==0);

  my $result = GetOptions (
        "genome=s" => \$genome,
        "annotation=s" => \$annotfile,
        "output=s" => \$output,
        "mod-dir=s" => \$mod_dir,
        "force-single-file" => \$force_single_file,
        "force-multi-file" => \$force_multi_file,
        "use-links" => \$use_links,
        "bt2-index" => \$bt2_index,
        "help" => \$help
      ) ;
  
  usage() if ($help);
 
  #Check options

  unless (-e $genome) {
    print "Error: must specify input genome\n";
    usage();
  }

  unless (defined $annotfile && -r $annotfile) {
    print "Error: must specify a readable annotation file\n";
    usage();
  }

  unless (defined $mod_dir) {
    $mod_dir = dirname($annotfile);
  }

  unless (defined $output) {
    print "Error: must specify output genome\n";
    usage();
  }

  if (-d $output) {
    my $output_search = File::List->new($output);
    croak "Error: output directory already contains .fa files" if @{$output_search->find(qr/\.fa$/)} > 0;
  }

  exit unless $result;
}

sub argument {
  my $var = shift;
  my $description = shift;
  my $default = shift;

  return sprintf("  \%\-16s - %s",
    (defined $default ? sprintf("%s (%s)",$var,$default) : sprintf("%s",$var)),
    $description);
}

sub usage()
{
print<<EOF;
Title, by Robin Meyers, ddmonthyyyy

This program .


Usage: $0 --genome (FILE|DIR) --annotation FILE --output PATH/NAME [--opts]

Arguments (defaults in parentheses):

$arg{"--genome","Input genome - either the file path of a single fasta, or a directory of several fasta"}
$arg{"--annotation","Annotation file of modifications"}
$arg{"--output","Directory to write the modified genome - the genome will be named according to the base folder"}
$arg{"--mod-dir","Directory containing sequences to be included in modifications"}
$arg{"--force-single-file","When reading in multiple fastas, force output to single file"}
$arg{"--force-multi-file","When reading in a single fasta, force output of multiple files"}
$arg{"--use-links","When in multi-file mode, create symbolic link to unmodified chromosomes instead of copying them"}
$arg{"--bt2-index","Create bowtie2 index when finished"}
$arg{"--help","This helpful help screen."}


EOF

exit 1;
}