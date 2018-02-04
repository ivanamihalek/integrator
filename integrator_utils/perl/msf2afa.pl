#! /usr/bin/perl -w 
use IO::Handle;         #autoflush
use File::Copy;     # copy a file (no kidding)
# FH -> autoflush(1);
defined $ARGV[0]  ||
    die "Usage: msf2afa.pl <msffile> \n"; 


$home = `pwd`;
chomp $home;
$name = $ARGV[0] ;

open ( MSF, "<$name" ) ||
    die "Cno: $name  $!\n";
	

while ( <MSF>) {
    last if ( /\/\// );
    last if ( /CLUSTAL FORMAT for T-COFFEE/ );
}
@names = ();
while ( <MSF>) {
    next if ( ! (/\w/) );
    chomp;
    @aux = split;
    $seq_name = $aux[0];
    if ( defined $seqs{$seq_name} ){
	$seqs{$seq_name} .= join ('', @aux[1 .. $#aux]);
    } else { 
	$seqs{$seq_name}  = join ('', @aux[1 .. $#aux]);
	push @names, $seq_name;
    }
}

close MSF;


	

foreach $seq_name ( @names ) {
    $seqs{$seq_name} =~ s/\./-/g;
    @seq = split ('', $seqs{$seq_name});
    print  ">$seq_name \n";
    $ctr = 0;
    for $i ( 0 .. $#seq ) {
	print   $seq[$i];
	$ctr++;
	if ( ! ($ctr % 50) ) {
	    print  "\n";
	}
    }
    print  "\n";
}

 


