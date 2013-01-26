#!/usr/bin/perl

use strict;
use warnings;
use Carp;
use File::Basename;
use IO::Handle;
use IO::File;
use CGI;
use List::Util qw(sum);


sub mean {
    return sum(@_)/@_;
}

sub argument {
    my $var = shift;
    my $description = shift;
    my $default = shift;

    return sprintf("  \%\-16s - %s\n",
        (defined $default ? sprintf("%s (%s)",$var,$default) : sprintf("%s",$var)),
        $description);
}

sub read_fastq {

    my $fh = shift;
    croak "Error: invalid file handle" unless $fh->opened;
    my $seq_name = $fh->getline();
    return () unless defined $seq_name;
    chomp($seq_name);
    croak "Error: unexpected format - missing '\@' symbol in sequence header $seq_name" unless $seq_name =~ /^\@.*$/;
    my $seq_bases = $fh->getline();
    croak "Error: bad input file, expecting line with sequences" unless defined $seq_bases;
    croak "Error: unexpected base in sequence $seq_bases" unless $seq_bases =~ /^[AGCTagctNn]*$/;
    chomp($seq_bases);
    my $seq_name2 = $fh->getline();
    croak "Error: bad input file, expecting line with sequence name2" unless defined $seq_name2;
    croak "Error: unexpected format - missing '\+' symbol in sequence name2 $seq_name2" unless $seq_name2 =~ /^\+.*$/;
    chomp($seq_name2);
    my $seq_quals = $fh->getline();
    croak "Error: bad input file, expecting line with quality scores" unless defined $seq_quals;
    chomp($seq_quals);
    croak "Error: sequence and quality must be same length" unless length($seq_bases) == length($seq_quals);

    return ($seq_name,$seq_bases,$seq_name2,$seq_quals);
}

sub write_fastq {
    croak "Error: too few arguments to 'write_fastq' method" if scalar @_ < 5;

    my $fh = shift;
    croak "Error: invalid file handle" unless $fh->opened;
    my $seq_name = shift;
    my $seq_bases = shift;
    my $seq_name2 = shift;
    my $seq_quals = shift;

    $seq_name = "@" . $seq_name unless $seq_name =~ /^\@.*$/;
    croak "Error: sequence contains unexpected base" unless $seq_bases =~ /^[AGCTagctNn]+$/;
    $seq_name2 = "+" . $seq_name2 unless $seq_name2 =~ /^\+.*$/;
    croak "Error: sequence and quality must be same length" unless length($seq_bases) == length($seq_quals);

    print $fh $seq_name,"\n";
    print $fh $seq_bases,"\n";
    print $fh $seq_name2,"\n";
    print $fh $seq_quals,"\n";

}

sub quals_to_ascii (@) {
    my @quals = @_;
    return join("",map( chr($_ + 33), @quals ));
}

sub quals_from_ascii ($) {
    my $quals = shift;
    return map( ord($_) - 33, split("", $quals));
}

#Quickly calculate hamming distance between two strings
#
#NOTE: Strings must be same length.
#      returns number of different characters.
#see  http://www.perlmonks.org/?node_id=500235
sub mismatch_count($$) { length( $_[ 0 ] ) - ( ( $_[ 0 ] ^ $_[ 1 ] ) =~ tr[\0][\0] ) }

sub reverseComplement {
    my $seq = shift;
    (my $rc = reverse($seq)) =~ tr/ACGTacgtNn/TGCAtgcaNn/;
    return $rc;
}

sub listFilesInDir 
{
    my $path = shift;
    croak "Error: path $path not found" unless -d $path;
    my ($record, @data, @filenames);
    my @list = `ls -l $path`;

    foreach $record (@list) {
        next if $record =~ /^d/;
        @data = split(/\s+/,$record);
        next unless defined($data[8]);
        push(@filenames, join("/",$path,$data[8]));
    }
    
    carp "Warning: nothing found in path $path" unless @filenames > 0;

    return @filenames;
}

sub parseFilename {
    my $fullname = shift;
    my ($name, $path, $ext) = fileparse($fullname, qr/\.\w{2,5}$/);
    return ($path, $name, $ext);
}

sub readChromsizeFile {

    my $chrsizefile = shift;
    my %chrsize;

    # Open Chrsize file to get the total length of each Chr
    open (CHRSZ, "<", $chrsizefile) or croak "Error: cannot open $chrsizefile for reading";
    while (<CHRSZ>) {
        chomp;
        my ($chrnum, $len) = split(/\t/);
        next unless $chrnum =~ /[Cc]hr/;
        $chrsize{$chrnum} = $len;
    }

    return %chrsize;
}

#invert exit status of system call
sub System {
    my $cmd = shift;
    print "$cmd\n";
    my $status = system($cmd);
    return !$status;
}

sub dateFileFormat
{
  my ($sec, $min, $hour, $day,$month,$year) = localtime time;
  return( sprintf( "%04d%02d%02d_%02d%02d%02d",$year+1900,$month+1,$day,$hour,$min,$sec) );
}

sub cgiUpload {
    my $cgi = shift;
    my $param = shift;
    my $outdir = shift;
    my $safe_chars =  "A-Za-z0-9_.-";
    my $file = $cgi->param($param);
    my ($path,$name,$ext) = parseFilename($file);
    $file = $name.$ext;
    $file =~ s/ /_/g;
    $file =~ s/[^$safe_chars]//g;
    $file = $outdir.$file;
    print $cgi->p("Uploading goodfile to server at $file...");

    my $fh  = $cgi->upload($param);
    if (defined $fh) {
        # Upgrade the handle to one compatible with IO::Handle:
        my $io_handle = $fh->handle;
        open OUTFILE,">",$file or croak "Error: could not open $file for writing";
        while (my $bytesread = $io_handle->read(my $buffer,1024)) {
            print OUTFILE $buffer;
        }
    }
    System("perl -pi -e 's/\r/\n/g' $file");

    return $file;
}

1;
