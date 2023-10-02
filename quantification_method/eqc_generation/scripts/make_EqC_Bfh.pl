#!/bin/env perl

use strict;
use warnings;
use File::Basename;

die "Usage: make_EqC_Bfh.pl [BCUMI annotation]\n" if scalar @ARGV < 1;
my $infile = shift;
my $outbase = basename($infile);
$outbase =~ s/_BCUMIannot\.txt//;
my $eqcOutput = $outbase . "_eqc.txt";
open my $ifh, "<", $infile or die $!;
open my $ofh, ">", $eqcOutput or die $!;

my $annotCounter = 0;
my %annot;
my $BC_counter = 0;
my %BC;
my $EqC_counter = 0;
my $currEqC;
my %EqC;

while(my $line = <$ifh>){
  chomp $line;
  my ($completeAnn, $bc, $umi, $reads) = split "\t", $line;
  $currEqC = $completeAnn if ! defined $currEqC;
  if($currEqC ne $completeAnn){
    output_EqC();
    $EqC_counter ++;
    %EqC = ();
    $currEqC = $completeAnn;
  }
  my @anns = split ",", $completeAnn;
  @anns = get_uniq(@anns);
  foreach my $ann (@anns){
    if(! exists $annot{$ann}){
      $annot{$ann} = $annotCounter;
      $annotCounter ++;
    }
  }
  if(! exists $BC{$bc}){
    $BC{$bc} = $BC_counter;
    $BC_counter ++;
  }
  $EqC{$completeAnn}{$bc}{$umi} = $reads;
}
close $ifh or die $!;
output_EqC();
$EqC_counter ++;
close $ofh or die $!;

my $headerFile = $outbase . "_header.txt";
open $ofh, ">", $headerFile or die $!;
print {$ofh} scalar keys %annot, "\n";
print {$ofh} scalar keys %BC, "\n";
print {$ofh} $EqC_counter, "\n";
foreach my $out (sort { $annot{$a} <=> $annot{$b} } keys %annot){
  print {$ofh} $out, "\n";
}
foreach my $out (sort { $BC{$a} <=> $BC{$b} } keys %BC){
  print {$ofh} $out, "\n";
}
close $ofh or die $!;

sub output_EqC{
  my @annIndex;
  my @BCUMI_info;
  my $totalReads = 0;
  foreach my $eqc (keys %EqC){
    my @anns = split ",", $eqc;
    foreach my $ann (@anns){
      push @annIndex, $annot{$ann};
    }
    push @BCUMI_info, scalar keys %{$EqC{$eqc}};
    foreach my $bc (keys %{$EqC{$eqc}}){
      push @BCUMI_info, $BC{$bc};
      push @BCUMI_info, scalar keys %{$EqC{$eqc}{$bc}};
      foreach my $umi (keys %{$EqC{$eqc}{$bc}}){
	push @BCUMI_info, $umi;
	my @reads = split ",", $EqC{$eqc}{$bc}{$umi};
	$totalReads += scalar @reads;
	push @BCUMI_info, scalar @reads;
      }
    }
  }
  print {$ofh} scalar @annIndex, "\t", join("\t",@annIndex), "\t", $totalReads, "\t", join("\t",@BCUMI_info), "\n";
}

sub get_uniq{
  my @array = @_;
  my %seen;
  my @unique;
  foreach my $value (@array){
    if (! $seen{$value}++){
      push @unique, $value;
    }
  }
  return @unique;
}
