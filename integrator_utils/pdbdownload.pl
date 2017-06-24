#!/usr/bin/perl

defined ( $ARGV[0] ) ||
    die "Usage: pdbdownload.pl <pdbname>.\n";

my @pdbnames = ();
if ( -e $ARGV[0]  ) {
    @pdbnames = split "\n", `cat $ARGV[0]`;

} else {
    push @pdbnames, $ARGV[0];
}

use strict;
use warnings;
use Net::FTP;

my $PDB_REPOSITORY = "/databases/pdb/structures";
(-e $PDB_REPOSITORY) || die "pdb repository not found.\n";


foreach my $pdbname (@pdbnames) {

    $pdbname =~ s/\s//g;
    my @aux  = split ('\.', $pdbname); # get rid of extension
    $pdbname =  lc substr ($aux[0], 0, 4);
    if (  -e "$PDB_REPOSITORY/$pdbname.pdb" ) {
	print "$pdbname.pdb found in $PDB_REPOSITORY\n";
	next;
    }
    print $pdbname, " \n"; 
    
    my $ftp = Net::FTP->new("ftp.wwpdb.org", Debug => 0, Passive=> 1)
	or die "Cannot connect to ftp.wwpdb.org: $@";

    $ftp->login("anonymous",'-anonymous@')
	or die "Cannot login ", $ftp->message;

    $ftp->cwd("/pub/pdb/data/structures/all/pdb")
	or die "Cannot change working directory ", $ftp->message;
    $ftp->binary;
    $ftp->get("pdb$pdbname.ent.gz")
        or die "get failed ", $ftp->message;

    system ( "gunzip pdb$pdbname.ent.gz" ) && 
	die "error uncompressing pdb$pdbname.ent.gz.\n";
	
    `mv  pdb$pdbname.ent $PDB_REPOSITORY/$pdbname.pdb`;

    print "\t downloaded to $PDB_REPOSITORY \n"; 
}
