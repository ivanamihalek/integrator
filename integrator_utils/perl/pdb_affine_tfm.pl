#! /usr/bin/perl -w 



defined ( $ARGV[0])  && defined ( $ARGV[1] ) ||
    die "Usage: pdb_affine_tfm.pl   <pdb_file>  <tfm_file>".
    " [-trnsl_first|-tr_rot_tr].\n";

$pdbfile = $ARGV[0];
$tfmfile = $ARGV[1];

$translate_first = 0;
$tr_rot_tr = 0;

if ( defined $ARGV[2] ) {
    if ( $ARGV[2] eq "-trnsl_first"  ) { 
	$translate_first = 1;
    } elsif ( $ARGV[2] eq "-tr_rot_tr"  ) { 
	# translate, rotate & translate back
	$tr_rot_tr = 1;
    } else {
	die "argument  $ARGV[2] not recognized.\n";
    }
}

open ( TFM, "<$tfmfile") ||
    die "Cno $tfmfile: $!.\n";
# expects input in the format
# 
#$A[0][0]     $A[0][1]   $A[0][2]   $x0
#$A[1][0]     $A[1][1]   $A[1][2]   $y0
#$A[2][0]     $A[2][1]   $A[2][2]   $z0

$i = 0;
while ( <TFM> ) {
    @aux = split;
    for $j ( 0 .. 2) {
	$A[$i][$j] = $aux[$j];
    }
    $t[$i] = $aux[3];
    $i++;
}

close TFM;

#print "$A[0][0]     $A[0][1]   $A[0][2]  \n"; 
#print "$A[1][0]     $A[1][1]   $A[1][2]  \n"; 
#print "$A[2][0]     $A[2][1]   $A[2][2] \n"; 

open ( PDB, "<$pdbfile") ||
    die "Cno $pdbfile: $!.\n";

while ( <PDB> ) {
    if ( ! ( /^ATOM/ || /^HETATM/) ){
	print;
	next;
    }
    chomp;
    $crap = substr ($_, 0, 30);
    $crap2 = substr ($_, 54);
    $x = substr $_, 30, 8;  $x=~ s/\s//g;
    $y = substr $_, 38, 8;  $y=~ s/\s//g;
    $z = substr $_, 46, 8;  $z=~ s/\s//g;

    $xtr  = $x;
    $ytr  = $y;
    $ztr  = $z;

    # translate
    if ( $translate_first ) {
	$xtr  +=  $t[0];
	$ytr  +=  $t[1];
	$ztr  +=  $t[2];
    } elsif ( $tr_rot_tr ) {
	$xtr  -=  $t[0];
	$ytr  -=  $t[1];
	$ztr  -=  $t[2];
	
    }
    # rotate
    $xnew = $A[0][0]*$xtr +   $A[0][1]*$ytr  +  $A[0][2]*$ztr;
    $ynew = $A[1][0]*$xtr +   $A[1][1]*$ytr  +  $A[1][2]*$ztr;
    $znew = $A[2][0]*$xtr +   $A[2][1]*$ytr  +  $A[2][2]*$ztr;
    if ( !$translate_first || $tr_rot_tr ) {
	# translate
	$xnew += $t[0];
	$ynew += $t[1];
	$znew += $t[2];
    }
    printf  "%30s%8.3f%8.3f%8.3f%s \n",
	   $crap,  $xnew, $ynew, $znew, $crap2;
   
}
close PDB;


