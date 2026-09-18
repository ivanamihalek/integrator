#! /usr/bin/perl -w

# store as uniprot_basic_infos.csv
# and load with
#  mysqlimport --local -u root -p blimps_development  uniprot_basic_infos.csv
use strict;
use warnings;
sub parse();


@ARGV ||
	die "Usage:  $0  <file name> \n";

my $filename = $ARGV[0];
open (IF, "<$filename" )
	|| die "Cno $filename: $!.\n";


my $entry = "";
while ( <IF> ) {
	if (/^\/\//) {
		parse();
		$entry = $_;
	} else {
		$entry .= $_;
	}
}
parse();
close IF;

sub parse() {
	$entry || return;
	$entry =~ /\nOS   Homo sapiens/ || return;
	#print $entry;
	my ($uniprot_ids, $hgnc_id, $approved_symbol) = ("", "", "");
	foreach (split "\n", $entry) {
		if  (/^ID/) {
		} elsif (/^AC/) {
			# it looks like the first iD is current,
			# and the rest are the older ids
			# # AC can also stretch through multiple lines
			$_ =~ s/^AC//;
			$uniprot_ids .= $_;

		} elsif (/^DR   HGNC; /) {
			$_  =~ s/^DR   HGNC; //;
			$_  =~ s/[\s\.]//g;
			($hgnc_id, $approved_symbol)  = split ";";
		}
	}
	$uniprot_ids =~ s/\s//g;
	my @aux = split (";", $uniprot_ids);
	my $uniprot_id = shift @aux;
	my $old_uniprot_ids  = join ";", @aux;
	print "$uniprot_id\t$approved_symbol\t$old_uniprot_ids\n";

}

