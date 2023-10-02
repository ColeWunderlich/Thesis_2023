#!/bin/env perl

use strict;
use warnings;
use File::Basename;
use Getopt::Std;

my %opts;
getopts('u:c:', \%opts);
my ($UMItag,$BCtag);
(defined $opts{'u'}) ? ($UMItag = $opts{'u'}) : ($UMItag = "UR");
(defined $opts{'c'}) ? ($BCtag = $opts{'c'}) : ($BCtag = "CB");
#(defined $opts{'1'}) ? ($baseZero = "N") : ($baseZero = "Y");

usage() if scalar @ARGV < 1;
my $BAM = shift @ARGV;

open my $ifh, "-|", "samtools", "view", "-h", $BAM or die $!;

my %seenID;

while(my $line = <$ifh>){
  if($line =~ /^\@/){
    next;
  }
  next if $line !~ /$UMItag:Z:/;
  next if $line !~ /$BCtag:Z:/;
  chomp $line;
  my @tmp = split "\t", $line;
  my $readID = $tmp[0];
  next if exists $seenID{$readID};
  my ($BC,$UMI);
  foreach my $field (@tmp){
    if($field =~ /^$UMItag:Z:/){
      $UMI = $field;
      $UMI =~ s/^$UMItag:Z://;
    }
    if($field =~ /^$BCtag:Z:/){
      $BC = $field;
      $BC =~ s/^$BCtag:Z://;
    }
    else{
      next;
    }
  }
  next unless (defined $UMI && defined $BC);
  next if $BC eq "-";  # Skip lines if barcode not correctable
  $seenID{$readID} = 1;
  print "$readID\t$BC\t$UMI\n";
}

sub usage{
  print <<EOF

  Usage: perl extract_BC_UMI.pl (options) [BAM]

  Options:
  -u [TAG]    Specify UMI tag to use as read name (default UR)
  -c [TAG]    Specify cell barcode tag to use as read name (default CB)

EOF
    ;
  exit 1;
}
