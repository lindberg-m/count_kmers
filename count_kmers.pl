#!/usr/bin/perl

use strict;
use warnings;

my $usage = << "EOF";
usage:
  $0 [-h,-v,-p,-i,-u,-l,-b bedfile, -m bedfile] [SIZE] < fasta_file.fa

  Optional Arguments
    -h, --help        Show this message and exit.
    -v, --verbose     Increase terminal output.
    -p, --pyrimidine  When kmers are of odd length, use pyrimidine based
                      contexts with the middle nucleotide as a reference point
    -i, --ignore-amb  Ignore ambiguous contexts [default=true], setting this
                      param turn this feature off.
    -u, --upper       Turn all sequences uppercase before counting.
    -l, --ignore-low  Ignore contexts with lower case characters in them.
    -b, --bed  FILE   Bedfile to use for subsetting file.
    -m, --mask FILE   Bedfile to use for masking fasta. Positions overlapping
                      mask regions will be transformed to 'n' characters.

  Positional arguments:
    SIZE              Size of kmers to calculate bgfreq on. Defaults to 3

This program reads a fasta formatted file from standard input and
calculate k-mer frequencies of specified length. Optionally,
subset or mask sequences based on bedfile(s) first.

EOF

my %PARAMS = (
  SIZE       => 3,  # K-mer size
  VERBOSE    => 0,  # Terminal Verbosity
  PYRIMIDINE => 0,  # Report pyrimidine based kmers only
  IGNORE_AMB => 1,  # Remove ambiguous nucleotides
  MK_UPPER   => 0,  # Make sequences uppercase before any processing
  IGNORE_LOW => 0,  # Remove lower case nucleotides
  BEDFILE    => '', # Path to bed-file for regions to do analysis on
  MASK       => '', # Path to bed-file to use for masking the genome
);

sub main {
  parse_args(\%PARAMS);
  my $bg_counts = {};
  my $regions = $PARAMS{BEDFILE} ? parse_bed($PARAMS{BEDFILE}) : {};
  my $mask    = $PARAMS{MASK}    ? parse_bed($PARAMS{MASK})    : {};

  # Parse fasta from STDIN and perform subsetting and counting
  my $last_chrom = '';
  my $seq        = '';
  while (<STDIN>) {
    chomp;
    if ( /^>/ ) {
      s/^>//; print STDERR "$_\n" if ($PARAMS{VERBOSE});
      count_kmers(\$seq, $mask->{$last_chrom}, $regions->{$last_chrom}, $bg_counts, \%PARAMS) if ($seq);
      $last_chrom = $_;
      $seq = '';
    } else {
      $seq .= $_;
    }
  }

  # Last sequence isn't handeled in the while-loop, hence this:
  count_kmers(\$seq, $mask->{$last_chrom}, $regions->{$last_chrom}, $bg_counts, \%PARAMS) if ($seq);

  # Print results
  if ($PARAMS{PYRIMIDINE}) {
    # Combine pyrimidine based kmers with its reverse
    # complement counterparts
    my %bg_counts_pyr;
    my $mid_point = ($PARAMS{SIZE} - 1) / 2;
    my $pyrimidines = 'ctCT';
    for my $k (keys %{$bg_counts}) {
      if (substr($k, $mid_point, 1) =~ /[$pyrimidines]/) {
        $bg_counts_pyr{$k} += $bg_counts->{$k};
      } else {
        my $krc = revcomp($k);
        $bg_counts_pyr{$krc} += $bg_counts->{$krc};
      }
    }

    # Then print them
    for my $k (sort keys %bg_counts_pyr) {
      print "$k\t$bg_counts_pyr{$k}\n";
    }
  } else {
    for my $k (sort keys $bg_counts->%*) {
      print "$k\t$bg_counts->{$k}\n";
    }
  }
}

sub count_kmers {
  my $seq = shift; 
  my $mask = shift;
  my $regions = shift;
  my $counts = shift;
  my $params = shift;

  mask_sequence($seq, $mask) if ($params->{MASK});
  $$seq = uc $$seq           if ($params->{MK_UPPER});

  if ($params->{BEDFILE}) {
    foreach my $region (@{$regions}) {
      my $subseq = subset_region($$seq, $region);
      my $seqparts = remove_ambig_and_softmasked($subseq, $params->{IGNORE_LOW}, $params->{IGNORE_AMB});
      map { update_counts($_, $counts, $params->{SIZE}) } @{$seqparts};
    }
  } else {
    my $seqparts = remove_ambig_and_softmasked($seq, $params->{IGNORE_LOW}, $params->{IGNORE_AMB});
    map { update_counts($_, $counts, $params->{SIZE}) } @{$seqparts};
  }
}

sub remove_ambig_and_softmasked {
  my $seq          = shift; # Sequence
  my $rm_lowercase = shift; # Bool
  my $rm_ambig     = shift; # Bool

  my (@temp, @res);
  if ($rm_lowercase) {
    foreach (split /[a-z]+/, $seq) {
      push @temp, $_;
    }
  } else {
    push @temp, $seq;
  }

  if ($rm_ambig) {
    foreach my $sp (@temp) {
      foreach (split /[nxywrNXYWR]+/, $sp) {
        push @res, $_;
      }
    }
  }
  return \@res;
}

sub subset_region {
  my $seq = shift;
  my $region = shift;

  my $start = $region->[0];
  my $end   = $region->[1];
  my $len   = $end - $start;
  return substr($seq, $start, $len);
}

sub mask_sequence {
  my $seqr = shift; # Scalar Ref to sequence
  my $mask = shift; # Arref to (start, end) pairs

  my @mask_regions = sort {$a->[0] <=> $b->[0]} @{$mask};

  if (scalar @mask_regions > 0) {
    my $last_end = 0;
    for my $r (@mask_regions) {
      my $start  = $r->[0];
      my $end    = $r->[1];
      my $seqlen = $end - $start;

      substr($$seqr, $start, $seqlen) = 'n' x $seqlen;
    }
  }
}

sub update_counts {
  my $seq       = shift; # Arref,   DNA sequences
  my $counts    = shift; # Hashref, kmer counts per chrom
  my $kmer_size = shift; # Int,     size of kmers to count

  my $ctx;
  my $i = 0;
  my $seqlen = length($seq);
  for (my $end = $kmer_size; $end < $seqlen; $end++ ) {
    $ctx = substr($seq, $i++, $kmer_size);
    $counts->{$ctx}++;
  }
}

sub parse_bed {
  my $bedfile = shift;
  my %regions;
  open BED, '<', $bedfile or die "Cannot open bedfile: $bedfile\n";
  while (<BED>) {
    chomp;
    my ($chrom, $start, $stop) = split;
    push @{$regions{$chrom}}, [ $start, $stop ];
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
                $params->{BEDFILE} = $ARGV[++$i];
            } elsif ($ARGV[$i] eq '-m' || $ARGV[$i] eq '--mask') {
                $params->{MASK} = $ARGV[++$i]
            } else {
                die "Unrecongnized argument: $ARGV[$i]\n";
            }
        } else {
            $j++;
            if ($i == scalar @ARGV) {
                $PARAMS{SIZE} = $ARGV[$i];
            }
        }
    }

    die "Need at most 1 positional argument but got $j\n" . $usage if ($j > 1);
    if ($PARAMS{SIZE} % 2 == 0 && $params->{PYRIMIDINE}) {
        die "Cannot use pyrimidine based calculation on an even-sized k-mers\n";
    }
}

sub revcomp {
  my $sequence = shift;
  $sequence =~ tr/actgACTG/tgacTGAC/;
  return reverse $sequence;
}

main()
