#!/usr/bin/perl -w
use strict;
use warnings;


my $infile = "/databases/mesh/mesh2017.txt";
# https://www.nlm.nih.gov/mesh/2015/mesh_trees/C16.html
my $inborn_diseases_of_metabolism_code = "16.320.565";
sub process(@);

open (INF, "<$infile") || die "Cno $infile: $!\n";

my $current = "";
while (<INF>){
     if ($current && !/\S/ ) {
          process($current);
          $current = "";
     } else {
          $current .= $_;
     }

} ;


######################################
sub process (@) {

    my $entry = $_[0];
    $entry =~ $inborn_diseases_of_metabolism_code && $entry=~/OMIM/ && print $entry."------------------\n";

}