#! /usr/bin/perl -w 
use IO::Handle;         #autoflush
use File::Copy;     # copy a file (no kidding)
# FH -> autoflush(1);
defined $ARGV[0]  ||
    die "Usage: $0  <msffile> \n"; 


$home = `pwd`;
chomp $home;
$name = $ARGV[0] ;

open ( MSF, "<$name" ) ||
    die "Cno: $name  $!\n";
	

while ( <MSF>) {
    last if ( /\/\// );
}

@names = ();
while ( <MSF>) {
    next if ( ! (/\w/) );
    chomp;
    @aux = split;
    $seq_name = $aux[0];
    if ( defined $sequence{$seq_name} ){
	$sequence{$seq_name} .= join ('', @aux[1 .. $#aux]);
    } else { 
	$sequence{$seq_name}  = join ('', @aux[1 .. $#aux]);
	push @names, $seq_name ;
    }
}

close MSF;


# turn the msf into a table (first index= sequence, 2bd index= position
$seq = 0;
foreach $seq_name ( keys %sequence ) {
    @aux = split '', $sequence{$seq_name};
    foreach $pos ( 0 .. $#aux ) {
	$array[$seq][$pos] = $aux[$pos];
    }
    $seq++;
    
}

$max_seq = $seq-1; # max index a seq can have
$max_pos = $#aux;  # max index a position can have


for $pos (0 .. $max_pos) {
    $all_gap = 1;
    for $i (0 .. $max_seq) {
	if ( $array[$i][$pos] !~ /[\.\-]/) {
	    $all_gap = 0;
	}
    }
    if ( $all_gap ) {
	$skip{$pos} = 1;
    }

}

# get rid of all-gap positions
for ($pos=$max_pos; $pos >=0; $pos--) {

    next if ( ! defined $skip{$pos} );
    foreach $seq_name ( @names  ) {
	$sequence{$seq_name} = substr ($sequence{$seq_name},0, $pos).
	    substr ($sequence{$seq_name},$pos+1);
    }
}



	
$seqlen = length $sequence{$seq_name};
print "PileUp\n\n";
print "            GapWeight: 30\n";
print "            GapLengthWeight: 1\n\n\n";
printf ("  MSF: %d  Type: P    Check:  9554   .. \n\n",$seqlen) ;
foreach $seq_name ( @names  ) {
    printf (" Name: %-20s   Len: %5d   Check: 9554   Weight: 1.00\n", $seq_name, $seqlen);
}
printf "\n//\n\n\n\n";

for ($j=0; $j  < $seqlen; $j += 50) {
    foreach $seq_name ( @names  ) {
	printf "%-20s", $seq_name;
	for ( $k = 0; $k < 5; $k++ ) {
	    if ( $j+$k*10+10 >= $seqlen ) {
		printf ("%-10s ",   substr ($sequence{$seq_name}, $j+$k*10 ));
		last;
	    } else {
		printf ("%-10s ",   substr ($sequence{$seq_name}, $j+$k*10, 10));
	    }
	}
	printf "\n";
    } 
    printf "\n";
}
