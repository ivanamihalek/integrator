#! /usr/bin/perl -w
# read in PDB file and a seqeunce
# delete all the sidechains, keeping Cbravo's only
# where possible and appropriate, and replace the name by the name from the 
# sequence

use IO::Handle;         #autoflush
# FH -> autoflush(1);

(defined $ARGV[1] ) ||
    die "Usage: switch_backbone_identity.pl  <pdbfile> ".
    "<new seq (fasta)> [<pdbid_start>  <pdbid_end>].\n";

$pdbfile =  $ARGV[0];
$seqfile =  $ARGV[1];
$pdb_start = -1;
( defined  $ARGV[2] ) &&  ($pdb_start =  $ARGV[2]);
$pdb_end =  10000;
( defined  $ARGV[3] ) &&  ($pdb_end =  $ARGV[3]);


# extension: replace only selected  piece of pdb

# check the number of aa's in the pdb
$filename = $pdbfile;
open ( IF, "<$filename" ) || die "Cno $filename: $!.\n";
$res_ctr = 0;
$old_res_seq = -100;
$old_res_name  ="";
$frag_length = 0;
while ( <IF> ) {
    next if ( ! /^ATOM/);
    $res_seq  = substr $_, 22, 4;  $res_seq=~ s/\s//g;
    $res_name = substr $_,  17, 3; $res_name=~ s/\s//g;
    if ( $res_seq != $old_res_seq  ||  ! ($res_name eq $old_res_name) ){
	$old_res_seq =  $res_seq;
	$old_res_name =  $res_name;
	$res_ctr++;
	if ( $pdb_start <=  $res_seq &&  $res_seq <= $pdb_end ) {
	    $frag_length ++;
	} 
    }
    
}
close IF;
$pdblen = $frag_length;
# read in the sequence
$filename = $seqfile;
open ( IF, "<$filename" ) || die "Cno $filename: $!.\n";
$seq = "";
while ( <IF> ) {
    next if ( /^>/ ) ;
    chomp;
    $seq .= $_;
}
close IF;
$seq =~ s/\s//g;
$seqlen = length $seq;

# check that the lengths of the two match
if  ( $seqlen !=  $frag_length )  {
    print "no of res in pdb: $frag_length.\n";
    print "no of res in seq: $seqlen.\n";
    die "The two seq lengths do not match.";
}


%letter_code = ( 'GLY', 'G', 'ALA', 'A',  'VAL', 'V', 'LEU','L', 'ILE','I',
           'MET', 'M', 'PRO', 'P',  'TRP', 'W', 'PHE','F', 'SER','S',
           'CYS', 'C', 'THR', 'T',  'ASN', 'N', 'GLN','Q', 'TYR','Y',
               'LYS', 'K', 'ARG', 'R',  'HIS', 'H', 'ASP','D', 'GLU','E');
foreach  $tri ( keys %letter_code  ) {
    $letter2three{ $letter_code{$tri} } = $tri;
}

@backbone = ("N", "CA", "C", "O", "CB");
@newres = split '', $seq;

# output the pdb:
# if new and old type are the same, keep everything
#  deleting everything but cbravo
# if new identity GLY, skip Cbravo too
# use new id as type ...
$filename = $pdbfile;
open ( IF, "<$filename" ) || die "Cno $filename: $!.\n";
$res_ctr = -1;
$old_res_seq = -100;
$old_res_name  ="";



while ( <IF> ) {

    if ( ! /^ATOM/ ) {
	print ;
	next;
    }
    $name = substr $_,  12, 4 ;  $name =~ s/\s//g; 
    $name =~ s/\*//g; 
    $alt_loc = substr $_,16, 1 ;  $alt_loc =~ s/\s//g;
    $res_seq  = substr $_, 22, 4;  $res_seq=~ s/\s//g;
    $res_name = substr $_,  17, 3; $res_name=~ s/\s//g;
    if ($res_seq < $pdb_start || $res_seq > $pdb_end ) {
	print ;
	next;
    }
    next if ( $alt_loc =~ "B" );
    $newline = $_;
    substr ($newline,16, 1 ) = " "; # alt loc

    if ( $res_seq != $old_res_seq  ||  ! ($res_name eq $old_res_name) ){
	$old_res_seq =  $res_seq;
	$old_res_name =  $res_name;
	$res_ctr++;
    }
    if ( ($newres[$res_ctr] eq ".") ||  ($res_name eq   $letter2three{ $newres[$res_ctr]}) ) {
	if ( ! defined $newline ) {
	    print "blah";
	}
	print $newline;
    } else {
	$newtype = $newres[$res_ctr];
	foreach $bb_atom ( @backbone ) {
	    next if ( $newtype eq "G" &&  $bb_atom eq"CB");
	    if ( $name eq  $bb_atom  ) {
		substr( $newline,  17, 3) = $letter2three{ $newtype };
		print $newline;
		last;
	    }
	}
    }
    
    
}
close IF;
