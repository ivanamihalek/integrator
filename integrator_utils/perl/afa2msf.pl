#! /usr/bin/perl -w
use strict;
use warnings;
use IO::Handle;         #autoflush
# FH -> autoflush(1);

# if fasta already aligned, convert to msf

defined $ARGV[0]  ||
    die "Usage: afa2msf.pl  <afa_file>.\n"; 

my $fasta =  $ARGV[0];

my @names = ();
open ( FASTA, "<$fasta") ||
    die "Cno $fasta: $!\n";
my %sequence = ();
my $name = "";
TOP: while ( <FASTA> ) {
    next if ( !/\S/);
    chomp;
    if (/^>\s*(\S+)/ ) {
		$name = $1;
		push @names,$name;
		$sequence{$name} = "";
    } else  {
		s/\./\-/g;
		s/\#/\-/g;
		s/\s//g;
		#s/x/\./gi;
	    $sequence{$name} .= $_;
    } 
}
close FASTA;

my $longest_name = -1;
foreach $name (@names) {
    if ( length $name > $longest_name ) {
	$longest_name = length $name 
    }
}
$longest_name ++;
( $longest_name < 20 ) && ($longest_name=20);

my $seqlen = length $sequence{$name};
print "PileUp\n\n";
print "            GapWeight: 30\n";
print "            GapLengthWeight: 1\n\n\n";
printf ("  MSF: %d  Type: N    Check:  9554   .. \n\n",$seqlen) ;

my $format = " Name: %-$longest_name"."s   Len: %5d   Check: 9554   Weight: 1.00\n";

foreach $name ( @names  ) {
    printf ( $format, $name, $seqlen);
}
printf "\n//\n\n\n\n";

$format = "%-$longest_name"."s";
for (my $j=0; $j  < $seqlen; $j += 50) {
    foreach $name ( @names  ) {
	printf $format, $name;
	for (my $k = 0; $k < 5; $k++ ) {
	    if ( $j+$k*10+10 >= $seqlen ) {
		printf ("%-10s ",   substr ($sequence{$name}, $j+$k*10 ));
		last;
	    } else {
		printf ("%-10s ",   substr ($sequence{$name}, $j+$k*10, 10));
	    }
	}
	printf "\n";
    } 
    printf "\n";
}
