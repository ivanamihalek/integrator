#!/usr/bin/perl -w 
# shorten the names staring with gi (like in entres\z return)


use strict;
use warnings;
while (<STDIN>) {
    next if ( !/\S/);
    if (0) {
    } elsif ( />\s*tr\|(\w+?)\|/ ) { # for trembl  names
		print ">$1\n";
	} elsif ( />\s*(\S+)/ ) { # for uniprot
		print ">$1\n";
    } else {
		print;
    }
 }


