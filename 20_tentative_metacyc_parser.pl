#!/usr/bin/perl
use strict;
use warnings FATAL => 'all';
my $datfile = "/databases/meta_cyc/data/datfiles/proteins.dat";

open (INFILE,  "<$datfile") || die "error opeining $datfile: $!\n";

my $record = "";
my $human  = 0;
my $polypep = 0;
my $human_ct = 0;
my $polypeptide_ct = 0;

while (<INFILE>)  {
        if (/^\/\//) {
            $_ =~ s/\^//g; # not sure what's the dire need to use this shit
            # the end of record
            $human_ct += $human;
            if ($human) {
                $polypeptide_ct += $polypep;
                print ("$record\n\n");
            }
            $record = $_;
            $human  = 0;
            $polypep = 0;

        } else {

            $record .= $_;
            if (/^SPECIES/ ) {
                my ($fldnm, $taxid) = split " - ";
                $taxid =~ s/\s//g;
                if ($taxid eq "TAX-9606") {
                    $human = 1;
                }
            } elsif (/^TYPES/ ) {
                my ($fldnm, $type) = split " - ";
                $type =~ s/\s//g;
                $type = lc $type;
                if ($type =~ "polypeptide") {
                    $polypep =1;
                }
            }
        }
}

printf "human records: $human_ct\n";
printf "      polypeptides: $polypeptide_ct\n";

1;