#!/usr/bin/perl

# tx idsv... though in this particuar file (proteins.dat) it seems
# we do not have the info stored under groups, only individua species
# human 9606
# mammals 40674
# vertebrates 7742
# metazoans 33208

# plants 33090
# archaea 2157
# fungi 4751


use strict;
use warnings FATAL => 'all';

@ARGV==1 || die "Usage: $0 <datfile name>\n";
my $datfile = $ARGV[0];

open (INFILE,  "<$datfile") || die "error opeining $datfile: $!\n";

# first find attributes -- nice idea, except that some datfiles
# do not have the attributes listed in the header

my %attribute = ();
while (<INFILE>)  {
    if (/^([A-Z\-]+)\s\-\s\S+/) {
        $1 eq "CREDITS"  && next;
        $1 eq "DBLINKS"  && next;
        defined $attribute{$1} && next;
        $attribute{$1} = 1;
    }
}
my @attributes = sort { $a cmp $b }(keys %attribute);
@attributes || die "no attributes found\n";

seek INFILE, 0, 0; # rewind
my %values = ();
my $current_key = 0;
while (<INFILE>) {
    if (/^\#/) {
        next;
    } elsif (/^\/\//) {
        for my $key (@attributes) {
            defined $values{$key} || next;
            print "-------------------------\n";
            print "$key\n";
            print "$values{$key}\n";
        }
        print "\n\n";
        %values = ();

    } elsif (/^\//) { # this is continuation:
        $values{$current_key} .= substr $_, 1; # get rid of the slash

    } elsif (/^([A-Z\-]+)\s\-\s(.+)/) {
        $current_key = $1;
        if ( defined $values{$current_key}){
            $values{$current_key} .= "; ".$2;
        } else {
            $values{$current_key} = $2;
        }
    }
}






1;