#! /usr/bin/perl -w
use strict;
use warnings;
use IO::Handle;         #autoflush
# FH -> autoflush(1);

defined ( $ARGV[0]  ) ||
    die "Usage: $0   <pdb_file> ".
    "  [<old_chain_name> <new_chain_name> <res from> <res to>].\n";

my $pdbfile = $ARGV[0];
my $old_name = defined $ARGV[1]? $ARGV[1]:"";
my $new_name = defined $ARGV[2]? $ARGV[2]:"A";
my $from     = defined $ARGV[3]? $ARGV[3]:-100;
my $to       = defined $ARGV[4]? $ARGV[4]:100000;

open ( IF, "<$pdbfile") ||
    die "Cno $pdbfile: $!.\n";

while ( <IF> ) {

    if ( ! /^ATOM/ && ! /^HETATM/ ) {
        print  ;
        next;
    }
    my $res_seq  = substr $_, 22, 4;  $res_seq=~ s/\s//g;
    my $chain_name = substr ( $_,  21, 1);

    if ($old_name && $chain_name ne $old_name) {
        print;
        next;
    }
    if ( $res_seq < $from  ||  $res_seq > $to) {
        print;
        next;
    }
    my $line = $_;
    substr ( $line,  21, 1) = $new_name;
    print $line;
}

close IF;
