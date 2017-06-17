#! /usr/bin/perl -w
# read in PDB file and a seqeunce
# delete all the sidechains, keeping Cbravo's only
# where possible and appropriate, and replace the name by the name from the 
# sequence

use IO::Handle;         #autoflush
# FH -> autoflush(1);

(defined $ARGV[0] ) ||
    die "Usage: pdb2seq.pl  <pdbfile> [<chain>].\n";
$pdbfile =  $ARGV[0];
$chain = "";

(defined $ARGV[1] ) && ($chain = $ARGV[1]);
# check that the lengths of the two match
%letter_code = ( 'GLY', 'G', 'ALA', 'A',  'VAL', 'V', 'LEU','L', 'ILE','I',
                 'MET', 'M', 'PRO', 'P',  'TRP', 'W', 'PHE','F', 'SER','S', 'SCY', 'C',
                 'CYS', 'C', 'THR', 'T',  'ASN', 'N', 'GLN','Q', 'TYR','Y',
                 'LYS', 'K', 'ARG', 'R',  'HIS', 'H', 'ASP','D', 'GLU','E', 'PTR', 'Y',
                 'MSE', 'M' ); 
@modifications = ("ACE", "OXT");

$res_ctr = 0;
$old_res_seq = -100;
$old_res_name  ="";
$filename = $pdbfile;
@aux = split '\/', $pdbfile;
$rootname = pop @aux;
$rootname =~ s/\.pdb$//;
#($rootname, $blah) = split '\.', $rootname;
$rootname .= $chain;
#print ">$rootname \n";

open ( IF, "<$filename" ) || die "Cno $filename: $!.\n";
while ( <IF> ) {

    last if (  /^ENDMDL/  );
    if ( ! /^ATOM/ && ! /^HETATM/  ) {
	next;
    }

    if ( $chain ) {
	$chain_name = substr ( $_,  21, 1) ; $chain_name=~ s/\s//g;
	next if ( $chain_name ne $chain );
    }

    $name = substr $_,  12, 4 ;  $name =~ s/\s//g; 
    $name =~ s/\*//g; 
    $alt_loc  = substr $_, 16, 1;  $alt_loc =~ s/\s//g;
    $res_seq  = substr $_, 22, 5;  $res_seq =~ s/\s//g;
    $res_name = substr $_, 17, 4; $res_name =~ s/\s//g;
    defined $letter_code{$res_name} || next;
    grep ( /$res_name/, @modifications) && next;
    next if ( $alt_loc =~ "B" );
    $newline = $_;
    substr ($newline,16, 1 ) = " "; # alt loc

    if ( $res_seq ne $old_res_seq  ||  ! ($res_name eq $old_res_name) ){
	$old_res_seq =  $res_seq;
	$old_res_name =  $res_name;
	$res_ctr++;
	$letter = $letter_code{$res_name};
	if ( defined $letter ) {
	    print $letter;
	    if ( ! ($res_ctr %50 ) ){
		print "\n";
	    }
	} else {
	    $res_name =~ s/\s//g;
	    if ( $res_name =~ /[catugCATUG]/ ) {
		print "$rootname contains nucleotides:\n$_\n";
		exit;
	    } else {
		die  "res type $res_name not recognized\n";
	    }
	}
    }
}
close IF;
print "\n";
