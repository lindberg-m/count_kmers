#!/usr/bin/env perl

use strict;
use warnings;

my $usage = << "EOF";
usage:
  bgcount.pl [-h] [-v] [-p] [-b bedfile] SIZE < fasta_file.fa

  Optional Arguments
    -h, --help        Show this message and exit
    -v, --verbose     Increase terminal output
    -p, --pyrimidine  When kmers are of odd length, use pyrimidine based
                      contexts
    -i, --ignore-amb  Ignore ambiguous contexts [default=true], setting this param
                      turn this feature off.
    -u, --upper       Turn all sequences uppercase before counting
    -l, --ignore-low  Ignore contexts with lower case characters in them
    -b, --bed         Bedfile to use for subsetting file

  Positional arguments:
    SIZE              Size of kmers to calculate bgfreq on

This program reads a fasta formatted file from standard input and
calculate background frequencies of specified length. Optionally,
subset fasta based on a bedfile first.

EOF

my %PARAMS = (
  VERBOSE => 0,
  PYRIMIDINE => 0,
  IGNORE_AMB => 1,
  MK_UPPER => 0,
  IGNORE_LOW => 0,
  BEDFILE => ''
);

sub main {
  my $kmer_size = parse_args(\%PARAMS);
  my %bg_counts;
  my $regions = $PARAMS{BEDFILE} ? parse_bed($PARAMS{BEDFILE}) : {};

  my $seq_name = <STDIN>; # Assumes fasta file start with a sequence identifier
  $seq_name =~ s/^>//;
  my $old_seq_name = '';
  if ($PARAMS{BEDFILE}) {
    my $chroms = keys %{$regions};
  } else {
    my ($seq, $fasta_eof);

    $old_seq_name = $seq_name; chomp $old_seq_name;
    ($seq, $seq_name, $fasta_eof) = read_fasta();
    while ($seq) {
      update_counts($seq, \%bg_counts, $kmer_size, \%PARAMS);
      
      #DEBUG
      #      {
      #        my @seqparts = split /[nN]+/, $seq;
      #        my $minseq = $seq =~ s/[nN]+//gr;
      #        my $a = length($seq);
      #        my $b = length($minseq);
      #        my $c = scalar @seqparts;
      #        my $d = length (join '', @seqparts);
      #
      #        print "$old_seq_name\t$a\t$b\t$c\t$d\n";
      #      }
      #END DEBUG
      
      last if ($fasta_eof);
      ($seq, $seq_name, $fasta_eof) = read_fasta();
    }
    for my $k (sort keys %bg_counts) {
      print "$k\t$bg_counts{$k}\n";
    }
  }
}

sub update_counts {
  my $sequence = shift;
  my $counts   = shift;
  my $ks       = shift;
  my $params   = shift;

  $sequence = uc $sequence if ($params->{MK_UPPER});
  my @seqparts = $params->{IGNORE_LOW} ? split /[a-z]+/, $sequence : ( $sequence );

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
      $counts->{$ctx}++;
      $end++
    }
  }
}

sub read_fasta {
  my $seq           = '';
  my $EOF           = 1;
  my $next_seq_name = '';
  while (<STDIN>) {
    chomp;
    if (/^>/) {
      $next_seq_name = $_; $EOF = 0; last;
    } else {
      $seq .= $_;
    }
  }
  return ($seq, $next_seq_name, $EOF);
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
        $i++;
        $params->{BEDFILE} = $ARGV[$i];
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
