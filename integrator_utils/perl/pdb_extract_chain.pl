#! /usr/bin/perl -w

use strict;
use warnings;
@ARGV || die "Usage: $0   <pdb_file>  [-c<chain_name> | -i | -p].\n".
            "No space between -c and the chain name.\n".
            "Use -i to invert the selection and -p to extract peptide only.\n".
            "This script does not know how to hadle modified atoms (HETATM) in the peptide chain.\n";
my $pdbfile = shift @ARGV;
my $query_chain_name = "";
my $peptide_only = 0;
my $invert_selection = 0;
foreach my $arg (@ARGV) {
    substr($arg, 0, 2)eq'-c' && ($query_chain_name=substr($arg,2,1));
    substr($arg, 0, 2)eq'-p' && ($peptide_only=1);
    substr($arg, 0, 2)eq'-i' && ($invert_selection = 1);
}


open ( IF, "<$pdbfile") || die "Cno $pdbfile: $!.\n";

my %seen = ();
while ( <IF> ) {

    last if ( /^ENDMDL/);

    $peptide_only && !/^ATOM/  && next;
    next if ( !/^ATOM/ && !/^HETATM/ );

    my $chain_name = substr ($_,  21, 1) ;
    if ($query_chain_name) {
        if ($invert_selection){
            ($chain_name eq $query_chain_name) && next;
        } else {
            ($chain_name ne $query_chain_name) && next;
        }
    }

    my $res_name  = substr $_, 17, 3; $res_name =~ s/\s//g;
    #try to handle insertion code cases:
    my $res_seq   = substr $_, 22, 5;  $res_seq=~ s/\s//g;
    my $atom_name = substr $_, 12, 4; $atom_name=~ s/\s//g;

    if ( ! defined $seen{"$chain_name $res_seq $atom_name" } ){
        $seen{"$chain_name $res_seq $atom_name"} = 1;
        print $_;
    }
    
}

close IF;
