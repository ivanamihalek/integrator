#!/usr/bin/perl

# the tables I am currently using
# metacyc_classes.csv    metacyc_enzrxns.csv  metacyc_pathways.csv
# metacyc_proteins.csv		metacyc_pubs.csv       metacyc_regulations.csv
# metacyc_compounds.csv  metacyc_genes.csv    metacyc_protein_features.csv
# metacyc_protligandcplxes.csv	metacyc_reactions.csv  metacyc_rnas.csv


use strict;
use warnings FATAL => 'all';

@ARGV==1 || die "Usage: $0 <datfile name>\n";
my $datfile = $ARGV[0];

open (INFILE,  "<$datfile") || die "error opeining $datfile: $!\n";

# first find attributes -- nice idea, except that some datfiles
# do not have the attributes listed in the header

my %attribute = ();
while (<INFILE>)  {
    if (/^([A-Z\-]+)\s\-\s\S+/) {
        $1 eq "CREDITS"  && next;
        $1 eq "DBLINKS"  && next;
        defined $attribute{$1} && next;
        $attribute{$1} = 1;
    }
}
my @attributes = sort { $a cmp $b }(keys %attribute);
@attributes || die "$datfile: no attributes found\n";


my $datfile_root = $datfile;
$datfile_root =~ s/\.dat$//;
$datfile_root =~ s/\-/_/g;
my $tablename = "metacyc_$datfile_root";
#make the table name end in s to comply with RoR expectation
$tablename =~ /s$/ || ($tablename .= "s");

open (OUTF, ">$tablename.csv") || die "Cno $tablename: $!\n";

seek INFILE, 0, 0; # rewind
my %values = ();
my %max_length = ();
my $current_key = 0;
my $count = 0;
while (<INFILE>) {
    chomp;
    if (/^\#/) {
        next;
    } elsif (/^\/\//) {
        $count += 1;
         my $outstr = $count;
         for my $key (@attributes) {
            $outstr .= "\t";
            defined $values{$key} || next;
            #something is wrong here - some fields are humongous; looks like I have newlines
            $outstr .= $values{$key};
            #print "-------------------------\n$key\n$values{$key}\n";
            my $current_length = length($values{$key});
            if (! defined  $max_length{$key} || $current_length>$max_length{$key}) {
                $max_length{$key} = $current_length;
            }
        }
        #print "\n\n";
        print OUTF $outstr."\n";
        %values = ();

    } elsif (/^\//) { # this is continuation:
        $values{$current_key} .= substr $_, 1; # get rid of the slash

    } elsif (/^([A-Z\-]+)\s\-\s(.+)/) {
        $current_key = $1;
        if ( defined $values{$current_key}){
            $values{$current_key} .= "; ".$2;
        } else {
            $values{$current_key} = $2;
        }
    }
}
close OUTF;


open (OUTF, ">$tablename.sql") || die "Cno $tablename.sql: $!\n";
print OUTF "CREATE TABLE `$tablename` (\n";
print OUTF "     id INT(11) NOT NULL,\n"; # int(11) is the type that rails uses for id

for my $key (@attributes) {
    defined $max_length{$key} || next;
    #print "$key    max entry length: $max_length{$key}\n";
    my $field_name = lc $key;
    # make sure that the field names are not reserved words:
    for ("left", "right", "signal") {
        $field_name eq $_ && ($field_name = "rxn_".$field_name );
    }
    $field_name =~ s/\-/_/g;
    # if length < 100, take twice the value and fix the length (CHAR)
    my $length = $max_length{$key};
    my $field_type = "LONTEXT";
    if ($length < 100) {
        $field_type = sprintf "CHAR(%d)", $length*2;
    } elsif ($length < 1000) { # if length <1000, take twice the value and store as VARCHAR
        $field_type = sprintf "VARCHAR(%d)", $length*2;
    } elsif ($max_length{$key} < 65000) {# if lnegth > 1000  and < 65,000 store as text
        $field_type = "TEXT";
    }  # and finally if > 65,000 keep LONGTEXT
    if ($field_name eq "unique_id") {
        print OUTF "     $field_name  $field_type   NOT NULL,\n";

    } else {
        print OUTF"     $field_name  $field_type   DEFAULT NULL,\n";
    }
}
print OUTF "     PRIMARY KEY (`id`) \n";
print OUTF ") ENGINE=MyISAM;\n";

close OUTF;
print "\n===============================\n";
print "to create the table: \n";
print "mysql -u root -p blimps_environment < $tablename.sql\n";
print "\n";
print "to load data, use mysqlimport:\n";
print "mysqlimport -u root -p --local  blimps_development $tablename.csv \n";
print "(mysqlimport expects the name of the table to be the name of the input file, minus extension)";
print "===============================\n\n";

1;