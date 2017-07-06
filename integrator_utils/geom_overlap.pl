#! /usr/bin/perl -w 
# return fraction of the first molecule (in terms of fraction of the number of residues) that are closer than the cutoff to the second
use strict;
use warnings;

(defined $ARGV[0] && defined $ARGV[1]  ) ||
    die "usage: $0  <protein_pdb_file> <ligand_pdb_file>.\n";

my $pdbfile1 = $ARGV[0];
my $pdbfile2 = $ARGV[1];


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
    my $z = substr $_,46, 8; $z=~ s/\s//g;
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
    $res_name_1[$res_seq] = substr $_,  17, 3; $res_name_1[$res_seq]=~ s/\s//g;
    $ctr++;
}
if ($prev_res_label) {
    $atom_1_ctr{$prev_res_label} = $ctr;
    push @labels_1, $prev_res_label;
}

$ctr = 0;
my %atom_2_ctr = ();
my %atom_2     = ();
my @labels_2   = ();
my @res_name_2 = ();
$res_label      = "";
$prev_res_label = "";
while ( <PDBFILE2> ) {
    /^ATOM/  || /^HETATM/ || next;
    chomp;
    my $serial = substr $_, 6, 5;  $serial =~ s/\s//g;
    my $res_seq   = substr $_, 22, 4; $res_seq=~ s/\s//g;
    my $chain  = substr $_, 21, 1;
    my $x = substr $_,30, 8;  $x=~ s/\s//g;
    my $y = substr $_,38, 8;  $y=~ s/\s//g;
    my $z = substr $_,46, 8; $z=~ s/\s//g;
    $res_label= $chain.$res_seq;
    if ( $res_label ne $prev_res_label ) {
        if ($prev_res_label) {
            $atom_2_ctr{$prev_res_label} = $ctr;
            push @labels_2, $prev_res_label;
        }
        $ctr = 0;
        $prev_res_label = $res_label;
    }
    $atom_2{$res_label}[$ctr][0] = $serial;
    $atom_2{$res_label}[$ctr][1] = $x;
    $atom_2{$res_label}[$ctr][2] = $y;
    $atom_2{$res_label}[$ctr][3] = $z;
    $res_name_2[$res_seq] = substr $_,  17, 3; $res_name_2[$res_seq]=~ s/\s//g;
    $ctr++;
}
if ($prev_res_label) {
    $atom_2_ctr{$prev_res_label} = $ctr;
    push @labels_2, $prev_res_label;
}


#****************************************************************
my %geom_ctr_1 = ();
foreach my $aa1 (@labels_1) {
    defined $atom_1{$aa1}  || next;
    @{$geom_ctr_1{$aa1}} = (0,0,0);
    for (my $ctr1=0; $ctr1<$atom_1_ctr{$aa1}; $ctr1++) {
        for (my $i=0; $i<3; $i++) {$geom_ctr_1{$aa1}[$i] += $atom_1{$aa1}[$ctr1][$i+1]};
    }
    for (my $i=0; $i<3; $i++) {$geom_ctr_1{$aa1}[$i] /= $atom_1_ctr{$aa1}};
}

my $covered = 0;
foreach my $aa1 (@labels_1) {
    foreach my $aa2 (@labels_2) {
        my @min = (1000,1000,1000);
        my @max = (-1000,-1000,-1000);
        for (my $ctr2=0; $ctr2<$atom_2_ctr{$aa2}; $ctr2++) {
            for (my $i = 0; $i < 3; $i++) {
                $min[$i] = $atom_2{$aa2}[$ctr2][$i + 1] < $min[$i] ? $atom_2{$aa2}[$ctr2][$i + 1] : $min[$i];
                $max[$i] = $atom_2{$aa2}[$ctr2][$i + 1] > $max[$i] ? $atom_2{$aa2}[$ctr2][$i + 1] : $max[$i];
            }
        }
        my $inside = 1;
        for (my $i=0; $i<3; $i++) {
            if ( $geom_ctr_1{$aa1}[$i]< $min[$i] || $geom_ctr_1{$aa1}[$i]> $max[$i]) {
                $inside = 0;
                last;
            }
        }

        if ($inside) {
            $covered += 1;
            last;
        }
    }
}

printf "%.2f\n", $covered/@labels_1;