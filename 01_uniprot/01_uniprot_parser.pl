#! /usr/bin/perl -w

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
	if (/^ID/) {
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
	my ($uniprot_id, $full_name,$gene_name, $tissue, $fn) = ("", "", "", "", "");
	my $reading_function  = 0;
	foreach (split "\n", $entry) {
		if  (/^ID/) {
		} elsif (/^AC/) {
			# it looks like the first iD is current,
			# and the rest are the older ids
			$_  =~ s/^AC//;
			($uniprot_id)  = split ";";
			$uniprot_id =~ s/\s//g;
		} elsif (/^DE   RecName: /) {
			$_  =~ s/^DE   RecName: //;
			($full_name)  = split ";";
			$full_name  =~ s/Full=//g;
		} elsif (/^GN   Name=/) {
			$_ =~ s/^GN   Name=//;
			($gene_name) = split /[;\s]/; # it ocurred to some idiot to put cross0references here
		} elsif (/^RC   TISSUE/) {
			$_ =~ s/^RC   TISSUE=//;
			my $tiss = $_;
			# their attempt to xlnk are just making things messy; sometimes we can forget to close the bracket too
			$tiss =~ s/\{.*?\}//g;
			$tiss =~ s/\{.*?$//g;
			(substr $tiss, -1) eq ',' && chop $tiss;
			(substr $tiss, -1) eq ';' || ($tiss.=";");
			($tissue =~ $tiss) ||  ($tissue .= $tiss);
		} elsif (/^CC/) {
			if (/\-!\-/) {
				$reading_function  && last;

				if (/FUNCTION/) {
					$reading_function  = 1;
					$_ =~ s/CC   -!- FUNCTION: //;
					$fn .= $_." ";
				}
			} elsif ( $reading_function )  {
				$_ =~ s/CC      //;
				$fn .= $_." ";
			}
		}
	}
	print "$uniprot_id\t$gene_name\t$full_name\t$tissue\t";
	print "$fn\n";
	#exit;

}
