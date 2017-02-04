#! /usr/bin/perl

use strict;
use warnings;
my (@aux, @new, $blah);
while ( <>) {
    chomp;
    @aux = split "\t";
    @new = ();
    foreach (@aux) {
	$blah = $_;
	$blah =~ s/\s//g;
	($blah eq '-') && ($blah='');
	if (! $blah) {
	    push @new, '\N';
	} else {
	    # somebody in BioGRID thought it would be cute to use backslash
	    # as quotes
	    $blah = $_;
	    $blah =~ s/\\//g;
	    push @new, $blah;
	}
    }

    print join ("\t", @new), "\n";
}
