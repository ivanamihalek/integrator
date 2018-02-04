#! /usr/bin/perl -w


$MAX_GAPS_FRACTION  = 0.34;

sub replace (@);
(@ARGV>=3 ) ||
    die "Usage: $0  <msf_file_name>  <output afa  name> <threshold fraction> [<protected seq name>].\n"; 


($msf, $afa, $threshold) = @ARGV;

$protected = "";

if( defined $ARGV[3] ) {
    $protected = $ARGV[3];
}

##########################################
# read in the msf file:
open ( MSF, "<$msf") ||
    die "Cno $msf: $!\n";

while ( <MSF> ) {
    last if ( /\/\//);
}
@name = ();
%sequence = ();
do {
    if ( /\w/ ) {
	@aux = split;
	$name = $aux[0];
	$aux_str = join ('', @aux[1 .. $#aux] );
	if ( defined $sequence{$name} ) {
	    $sequence{$name} .= $aux_str;
	} else {
	    push @name, $name;
	    $sequence{$name}  = $aux_str;
	}
		
    } 
} while ( <MSF>);
close MSF;


###########################################
# turn the msf into a table (first index= sequence, 2nd index= position)
$seq = 0;
foreach $name ( @name ) {
    @aux = split '', $sequence{$name};
    foreach $pos ( 0 .. $#aux ) {
	$array[$seq][$pos] = $aux[$pos];
	$new_array[$seq][$pos] = $aux[$pos];
    }
    $seq++;   
}
$no_seqs = $seq;   # number of seqs
$max_seq = $seq-1; # max index a seq can have
$max_pos = $#aux;  # max index a position can have
    
# sanity check:
$no_seqs || die "Error: no seqs found.\n"; 

###########################################
# tag gapped positions
$seq = 0;

foreach $pos ( 0 .. $max_pos ) {
    $gap_ct = 0;
    foreach $seq(0 .. $max_seq  ) {
	if ( $array[$seq][$pos] =~ /[\.\-X]/ ) {
	    $gap_ct ++;
	}
    }
    if ( $gap_ct/$no_seqs > $MAX_GAPS_FRACTION ) {
	$gapped[$pos] = 1;
    } else {
	$gapped[$pos] = 0;
    }
}

##########################################
# calculate similarity table
for $seq1 ( 0 .. $max_seq) {


    $len1 = 0;
    for $pos ( 0 .. $max_pos) {
        ($array[$seq1][$pos] =~ /[\.\-X]/ ) || $len1++;
    }
    $len1 ||  die "Sequence of length zero (?).\n";

    $similarity[$seq1][$seq1] = -1;

    for $seq2 ( $seq1+1 .. $max_seq) {
 
       $len2 = 0;
       
       for $pos ( 0 .. $max_pos) {
	   ($array[$seq2][$pos] =~ /[\.\-X]/ ) || $len2++;
       }
       $len2 ||  die "Sequence of length zero (?).\n";
    
 	
       $common = 0;
       $common_length = 0;
       for $pos ( 0 .. $max_pos) {
	   if ( $array[$seq1][$pos] !~ /[\.\-X]/  && 
		$array[$seq2][$pos] !~ /[\.\-X]/) {
	       $common_length ++;
	       if ( $array[$seq1][$pos] eq $array[$seq2][$pos])    {
		   $common ++;
	       }
	   }
       } 
       #print "$name[$seq1] $name[$seq2]  $common   $common_length  \n";
       if ( $common_length ) {
	   $id =  $common/$common_length;
       } else {
	   $id = 0;
       }
       $similarity[$seq1][$seq2] = $id;
       $similarity[$seq2][$seq1] = $id;
    }
}


for $seq1 ( 0 .. $max_seq) {
    next if ($name[$seq1] eq $protected );
    next if ( $sequence{$name[$seq1]} !~  /[\.\-X]/);
    @replacementinorder = ();
    foreach $seq2 ( sort {$similarity[$seq1][$b]<=>$similarity[$seq1][$a]} ( 0 .. $max_seq) ){
        push @replacementinorder, $seq2; 
    }
    print "nearest replacement for $name[$seq1] is $name[$replacementinorder[0]].\n";
    printf "(similarity %4.2f)\n", $similarity[$seq1][$replacementinorder[0]];
    if (  $similarity[$seq1][$replacementinorder[0]]  < $threshold ) {
	printf "\t no replacement done\n";
    } else {
	replace ($seq1, @replacementinorder);
    }
    printf "\n";

}



open (AFA, ">$afa") || die "Cno $afa: $!.\n";
for $seq1 ( 0 .. $max_seq) {
    print AFA  ">$name[$seq1]\n";
    for $pos ( 0 .. $max_pos) {
	print AFA $new_array[$seq1][$pos]  ;
	( ($pos+1) % 50)  || print AFA "\n";
    }
    print AFA "\n";
}
close AFA;

########################################################
sub replace (@) {
    my $seq1 = shift @_;
    my @replacementinorder = @_;
    my $pos;
    my $seq2;

    for $pos ( 0 .. $max_pos) {

	next if ( $gapped[$pos] );
	next if ($array[$seq1][$pos] !~ /[\.\-X]/ );

	for $seq2 (  @replacementinorder ) {
	    last if ($similarity[$seq1][$seq2] < $threshold);
	    if ( $array[$seq2][$pos] !~ /[\.\-X]/ ) {
		print "\t replacing pos $pos with $name[$seq2]\n";
		$new_array[$seq1][$pos] = $array[$seq2][$pos];
		last;
	    }
	}
    }
       
}
