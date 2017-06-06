#! /usr/bin/perl -w 

use strict;
use warnings;
my $CUTOFF_R = 5.0;

(defined $ARGV[0] && defined $ARGV[1]  ) ||
    die "usage: $0  <protein_pdb_file> <ligand_pdb_file> [<cutoff_distance>].\n";

my $pdbfile1 = $ARGV[0];
my $pdbfile2 = $ARGV[1];
if ( defined $ARGV[2] ) {
    $CUTOFF_R =  $ARGV[2];
}
my $inverse = 0;
if ( defined $ARGV[3] ) {
    $inverse = 1;
}

open (PDBFILE1, "<$pdbfile1") ||
    die "could not open $pdbfile1 $!\n";

open (PDBFILE2, "<$pdbfile2") ||
    die "could not open $pdbfile2: $!\n";


my $ctr = 0;
my %atom_1_ctr = ();
my %atom_1     = ();
my @labels_1   = ();
my @res_name_1 = ();
my $res_label      = "";
my $prev_res_label = "";
while ( <PDBFILE1> ) {
    /^ATOM/  || /^HETATM/ || next;
    chomp;
    my $serial = substr $_, 6, 5;  $serial =~ s/\s//g;
    my $res_seq   = substr $_, 22, 4; $res_seq=~ s/\s//g;
    my $chain  = substr $_, 21, 1;
    my $x = substr $_,30, 8;  $x=~ s/\s//g;
    my $y = substr $_,38, 8;  $y=~ s/\s//g;
    my $z = substr $_, 46, 8; $z=~ s/\s//g;
    $res_label= $chain.$res_seq;
    if ( $res_label ne $prev_res_label ) {
        if ($prev_res_label) {
            $atom_1_ctr{$prev_res_label} = $ctr;
            push @labels_1, $prev_res_label;
        } 
        
        $ctr = 0;
        $prev_res_label = $res_label;
    }
    $atom_1{$res_label}[$ctr][0] = $serial;
    $atom_1{$res_label}[$ctr][1] = $x;
    $atom_1{$res_label}[$ctr][2] = $y;
    $atom_1{$res_label}[$ctr][3] = $z;
    $res_name_1[$res_seq] = substr $_,  17, 3; $res_name_1 [$res_seq]=~ s/\s//g;
    $ctr++;
}
if ($prev_res_label) {
    $atom_1_ctr{$prev_res_label} = $ctr;
    push @labels_1, $prev_res_label;
}


$ctr = 0;
my @atom_2 = ();
while ( <PDBFILE2> ) {
    /^ATOM/  || /^HETATM/ || next;
    chomp;
    my $serial = substr $_, 6, 5;  $serial =~ s/\s//g;
    my $res_seq  = substr $_, 22, 4;  $res_seq=~ s/\s//g;
    my $x = substr $_,30, 8;  $x=~ s/\s//g;
    my $y = substr $_,38, 8;  $y=~ s/\s//g;
    my $z = substr $_, 46, 8; $z=~ s/\s//g;
    $atom_2[$ctr][0] = $serial;
    $atom_2[$ctr][1] = $x;
    $atom_2[$ctr][2] = $y;
    $atom_2[$ctr][3] = $z;
    $ctr++;
}
my $atom_2_ctr = $ctr;

#****************************************************************
foreach my $aa1 (@labels_1) {
    defined $atom_1{$aa1}  || next;
    my $nearest = 1000;
    my ($x1, $y1, $z1, $x2, $y2, $z2);

        for (my $ctr1=0; $ctr1<$atom_1_ctr{$aa1}; $ctr1++) {
            $x1 = $atom_1{$aa1}[$ctr1][1];
            $y1 = $atom_1{$aa1}[$ctr1][2];
            $z1 = $atom_1{$aa1}[$ctr1][3];

            for (my $ctr2=0; $ctr2<$atom_2_ctr; $ctr2++) {
                $x2 = $atom_2[$ctr2][1];
                $y2 = $atom_2[$ctr2][2];
                $z2 = $atom_2[$ctr2][3];

                my  $r = ($x1-$x2)*($x1-$x2) + ($y1-$y2)*($y1-$y2) + ($z1-$z2)*($z1-$z2);
                $r = sqrt ($r);
                if ($r <= $CUTOFF_R && $nearest > $r) {
                    $nearest = $r;
                }
             }
        }
    

    if ( $nearest <  1000 ) {
        printf ("%s   %8.1f  \n", $aa1,  $nearest);
    }
}
