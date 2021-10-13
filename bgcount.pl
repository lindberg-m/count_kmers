#!/usr/bin/env perl

use strict;
use warnings;

my $usage = << "EOF";
usage:
  bgcount.pl [-h,-v,-p,-i,-u,-l,-b bedfile, -m bedfile] SIZE < fasta_file.fa

  Optional Arguments
    -h, --help        Show this message and exit.
    -v, --verbose     Increase terminal output.
    -p, --pyrimidine  When kmers are of odd length, use pyrimidine based
                      contexts.
    -i, --ignore-amb  Ignore ambiguous contexts [default=true], setting this
                      param turn this feature off.
    -u, --upper       Turn all sequences uppercase before counting.
    -l, --ignore-low  Ignore contexts with lower case characters in them.
    -b, --bed  FILE   Bedfile to use for subsetting file.
    -m, --mask FILE   Bedfile to use for masking fasta. Positions overlapping
                      mask regions will be transformed to 'n' characters.

  Positional arguments:
    SIZE              Size of kmers to calculate bgfreq on.

This program reads a fasta formatted file from standard input and
calculate background frequencies of specified length. Optionally,
subset fasta based on a bedfile first.

EOF

my %PARAMS = (
  VERBOSE    => 0,
  PYRIMIDINE => 0,
  IGNORE_AMB => 1,
  MK_UPPER   => 0,
  IGNORE_LOW => 0,
  BEDFILE    => '',
  MASK       => '',
);

sub main {
  my $kmer_size = parse_args(\%PARAMS);
  my %bg_counts;
  my $regions = $PARAMS{BEDFILE} ? parse_bed($PARAMS{BEDFILE}) : {};
  my $mask    = $PARAMS{MASK}    ? parse_bed($PARAMS{MASK})    : {};

  # Parse fasta from STDIN and perform subsetting and
  # counting
  my $last_chrom = '';
  my $seq        = '';
  while (<STDIN>) {
    chomp;

    if ( /^>/ ) {
      s/^>//; print STDERR "$_\n" if ($PARAMS{VERBOSE});

      if ($seq) {
        my $seqparts = subset_regions($seq, $last_chrom, $regions);
        update_counts($seqparts, \%bg_counts, $kmer_size, \%PARAMS);
      }
      $last_chrom = $_;
      $seq = '';
    } else {
      $seq .= $_;
    }
  }
  if ($seq) {
    my $seqparts = subset_regions($seq, $last_chrom, $regions);
    update_counts($seqparts, \%bg_counts, $kmer_size, \%PARAMS);
  }

  # Print results
  if ($PARAMS{PYRIMIDINE}) {
    # Add combine pyrimidine based kmers with its reverse
    # complement counterparts
    my %bg_counts_pyr;
    my $mid_point = ($kmer_size - 1) / 2;
    my $pyrimidines = 'ctCT';
    for my $k (keys %bg_counts) {
      if (substr($k, $mid_point, 1) =~ /[$pyrimidines]/) {
        $bg_counts_pyr{$k} += $bg_counts{$k};
      } else {
        my $krc = revcomp($k);
        $bg_counts_pyr{$krc} += $bg_counts{$krc};
      }
    }

    # Then print them
    for my $k (sort keys %bg_counts_pyr) {
      print "$k\t$bg_counts_pyr{$k}\n";
    }
  } else {
    for my $k (sort keys %bg_counts) {
      print "$k\t$bg_counts{$k}\n";
    }
  }
}

sub subset_regions {
  # Accept a sequence with its name
  # and look up the start and stop positions
  # for that sequence in a regions hash, parsed
  # from a BED-file. Returns all parts of the sequence
  # covered by the specified regions
  my $seq       = shift; # Arref,  complete DNA sequences for ...
  my $seqname   = shift; # String, ... this region/chromosome
  my $regions_r = shift; # Hashref, each key stores an array of (start, end) arrays
  my @seqparts;

  my %regions   = %{$regions_r};
  my @chroms    = keys %regions;


  if (@chroms) {
    my ($start, $offset, $seqpart);
    for my $posref (@{$regions{$seqname}}) {
      $start   = $posref->[0] - 1;
      $offset  = $posref->[1] - $start;
      $seqpart = substr($seq, $start, $offset);
      push @seqparts, $seqpart;
    }
  } else {
    @seqparts = ( $seq );
  }
  return \@seqparts;
}

sub mask_sequence {
  my $seq  = shift;
  my $mask = shift;

  my @mask_regions = sort {$a->[0] <=> $b->[0]} @{$mask};
  if (scalar @mask_regions > 0) {
    my $last_end = 0;
    my $seqlen;
    for my $r (@mask_regions) {
      $seqlen = $r->[0] - $last_end - 1;
      substr($seq, $last_end, $seqlen) = 'n' x $seqlen;
      $last_end = $r->[1];
    }
  }
  return $seq;
}

sub update_counts {
  my $sequence_r = shift; # Arref,   DNA sequences
  my $counts     = shift; # Hashref, kmer counts per chrom
  my $ks         = shift; # Integer, size of kmers
  my $params     = shift; # Hashref, parameters

  my @seqparts = @{$sequence_r};

   # Make all parts uppercase if necessary
  if ($params->{MK_UPPER}) {
    @seqparts = map (uc, @{$sequence_r}) 
  }
  
  # Remove lowercase characters (if necessary) simply
  # by splitting them away. 
  if ($params->{IGNORE_LOW}) {
    my @aux;
    for my $sp (@seqparts) {
      for my $sp2 (split /[a-z]+/, $sp){
        push @aux, $sp2
      }
    }
    @seqparts = @aux;
  }

  # Remove ambiguous characters (if necessary) simply
  # by splitting them away. 
  if ($params->{IGNORE_AMB}) {
    my @aux;
    for my $sp (@seqparts) {
      for my $sp2 (split /[nxywrNXYWR]+/, $sp) {
        push @aux, $sp2;
      }
    }
    @seqparts = @aux
  }

  for my $seqpart (@seqparts) {
    my $end = $ks;
    my $seqlen = length($seqpart);
    for (my $i = 0; $end < $seqlen; $i++) {
      my $ctx = substr($seqpart, $i, $ks);
      if ($ctx =~ /[0-9]+/) {
      }
      $counts->{$ctx}++;
      $end++
    }
  }
  return 0;
}

sub parse_bed {
  my $bedfile = shift;
  my %regions;
  open BED, '<', $bedfile or die "Cannot open bedfile: $bedfile\n";
  while (<BED>) {
    chomp;
    my ($chrom, $start, $stop) = split;
    push @{$regions{$chrom}}, [ $start - 1, $stop ];
  }
  return \%regions;
}

sub parse_args {
  my $params = shift;
  my $j = 0;
  my $size;
  for (my $i = 0; $i<@ARGV; $i++) {
    if ($ARGV[$i] =~ /^-/) {
      if ($ARGV[$i] eq '-h' || $ARGV[$i] eq '--help') {
        die $usage;
      } elsif ($ARGV[$i] eq '-v' || $ARGV[$i] eq '--verbose') {
        $params->{VERBOSE} = 1;
      } elsif ($ARGV[$i] eq '-p' || $ARGV[$i] eq '--pyrimidine') {
        $params->{PYRIMIDINE} = 1;
      } elsif ($ARGV[$i] eq '-i' || $ARGV[$i] eq '--ignore-amb'){
        $params->{IGNORE_AMB} = 0;
      } elsif ($ARGV[$i] eq '-u' || $ARGV[$i] eq '--upper'){
        $params->{MK_UPPER} = 1;
      } elsif ($ARGV[$i] eq '-l' || $ARGV[$i] eq '--ignore-low'){
        $params->{IGNORE_LOW} = 1;
      } elsif ($ARGV[$i] eq '-b' || $ARGV[$i] eq '--bed') {
        $i++; $params->{BEDFILE} = $ARGV[$i];
      } elsif ($ARGV[$i] eq '-m' || $ARGV[$i] eq '--mask') {
        $i++; $params->{MASK} = $ARGV[$i];
      } else {
        die "Unrecongnized argument: $ARGV[$i]\n";
      }
    } else {
      $j++;
      $size = $ARGV[$i];
    }
  }

  die "Need 1 positional argument but got $j\n" . $usage unless ($j == 1);
  if ($size % 2 == 0 && $params->{PYRIMIDINE}) {
    die "Cannot use pyrimidine based calculation on an even-sized k-mers\n";
  }
  return $size;
}

sub revcomp {
  my $sequence = shift;
  $sequence =~ tr/actgACTG/tgacTGAC/;
  return reverse $sequence;
}

main()
