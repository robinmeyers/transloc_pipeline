#!/usr/bin/perl

use strict;
use warnings;
use Clone 'clone';
use Cwd qw(abs_path);
use FindBin;
use lib abs_path("$FindBin::Bin/../lib");

require "PerlSub.pl";
require "TranslocHelper.pl";
require "TranslocFilters.pl";


# Global variables
my $tests_passed = 0;
my $tests_failed = 0;



sub create_random_alignment {

  my $strand = shift;

  $strand = random 1 or -1 unless defined $strand;

}

sub test_pair_is_proper {

  # Create random plus strand alignment
  my $p_aln1 = create_random_alignment(1);

  # Create alignment that should succeed
  my $aln2 = clone $aln1;

  $aln2->{Rstart} = $aln1->{Rend} + $aln1->{Strand} * $params->{max_pe_gap};
  $aln2->{Rend} = $aln1->{Rend} + $params->{max_pe_gap} - 1;

  # Create random minus strand alignment
  my $m_aln1 = create_random_alignment(-1);



}