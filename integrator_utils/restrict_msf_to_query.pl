#! /usr/bin/perl -w
use IO::Handle;         #autoflush
# FH -> autoflush(1);

 
defined $ARGV[1]  ||
    die "Usage: restrict_msf_to_query.pl  <msf_file_name>  <protected sequence> [ ...<more protected seqs>... ]\n"; 

$msf             = shift @ARGV;
@protected_names = @ARGV;



open ( MSF, "<$msf") ||
    die "Cno $msf: $!\n";

# read in the msf file:
while ( <MSF> ) {
    last if ( /\/\//);
}

%sequence = ();
@names = ();
do {
    if ( /\w/ ) {
	@aux = split;
	$name = $aux[0];
	$aux_str = join ('', @aux[1 .. $#aux] );
	if ( defined $sequence{$name} ) {
	    $sequence{$name} .= $aux_str;
	} else {
	    $sequence{$name}  = $aux_str;
	    push @names, $name;
	}
		
    } 
} while ( <MSF>);

# sanity check:
@names || die "Error in $0: no seqs found.\n"; 


# turn the msf into a table (first index= sequence, 2nd index= position
$max_pos = -1;
foreach $name ( @names ) {
    @{$array{$name}} = split '', $sequence{$name};
    if ( $max_pos < length($sequence{$name}) ) {
	$max_pos = length($sequence{$name}); 
    }
}

$max_pos--;

foreach $pos ( 0 .. $max_pos ) {
    
    #if ( $array[$query_seq][$pos] =~ '\.' ) {
    ###########################################################################################################################
    ##since in afa2msf.pl, zonghong has changed all '.' to '-' in order for seaview to count properly, here we have to change
    ##to detect '-'  
    ##########################################################################################################################

    $delete[$pos] = 1;
    for $query_seq (@protected_names) {
	if (! defined $array{$query_seq} ) {
	    next;
	}
	if ( $array{$query_seq}[$pos] !~ /[\-\.]/) { #|| is added by zong hong since some seq has \. as gaps
	    $delete[$pos] = 0;
	    last;
	}
    }
}



foreach $name ( @names) {
    $seq_new{$name} = "";
    foreach $pos ( 0 .. $max_pos ) {
	if ( ! $delete[$pos] ) {
	    $seq_new{$name} .= $array{$name}[$pos];
	}
    }
}


$deleted = 0;
foreach $pos ( 0 .. $max_pos ) {
    $deleted += $delete[$pos];
}

$new_length = length $seq_new{$name};

@aux = split '\.', $msf;

$longest_name = -1;
foreach $name ( @names ) {
    if ( length $name > $longest_name ) {
	$longest_name = length $name;
    }
}
$longest_name ++;
( $longest_name < 20 ) && ($longest_name=20);
$format = "%-$longest_name"."s";


print  "PileUp\n\n";
print  "            GapWeight: 30\n";
print  "            GapLengthWeight: 1\n\n\n";
printf   ("  MSF: %d  Type: P    Check:  9554   .. \n\n",$new_length) ;
foreach $name ( @names ) {
    printf  (" Name: ".$format."   Len: %5d   Check: 9554   Weight: 1.00\n", $name, $new_length);
}
print  "\n//\n\n\n\n";

for ($j=0; $j  < $new_length; $j += 50) {
    $seq = 0;
    foreach $name ( @names ) {
	printf  $format, $name;
	for ( $k = 0; $k < 5; $k++ ) {
	    if ( $j+$k*10+10 >= $new_length ) {
		printf  "%-10s ",   substr ($seq_new{$name}, $j+$k*10 );
		last;
	    } else {
		printf  "%-10s ",   substr ($seq_new{$name}, $j+$k*10, 10);
	    }
	}
	print  "\n";
	$seq++;
    } 
    print  "\n";
}

$max_pos ++;
#printf "removed $deleted columns   (out of $max_pos)\n";


exit(0);
