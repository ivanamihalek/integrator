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


my %manual_fix_for_ensembl = ();
#not sure what is going on here:
#  O43826 is annotated an manually reviewed entry for SLC37A4, associated with ENSG00000281500 --> patch in ensembl GRCh38
# U3KQS2 is unreviewd, un-annotated, also says SLC37A4, and associated with ENSG00000137700 ensembl GRCh37
$manual_fix_for_ensembl{'O43826'} = 'ENSG00000137700';

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
	my ($uniprot_ids, $full_name,$gene_name, $tissue, $loc, $fn, $ec_number, $ensembl_id,$aa_lengths) =
		("", "", "", "", "", "", "", "");
	my $reading_function  = 0;
	my $reading_cofactors = 0;
	my $reading_location  = 0;
	my @cofs = ();
	foreach (split "\n", $entry) {
		if  (/^ID/) {
		} elsif (/^AC/) {
			# it looks like the first iD is current,
			# and the rest are the older ids
			# # AC can also stretch through multiple lines
			$_ =~ s/^AC//;
			$uniprot_ids .= $_;

		} elsif (/^DE   RecName: /) {
			$_  =~ s/^DE   RecName: //;
			$_  =~ s/\{.+?\}//;
			($full_name)  = split ";";
			$full_name  =~ s/Full=//g;
		} elsif (/^DE\s*EC\=([\d\.]+)[\s\;\n]+/){
			$ec_number  = $1;
		} elsif (/^DR   Ensembl.*(ENSG\d{11}).*/){
			if ($ensembl_id !~ $1) {
				$ensembl_id && ($ensembl_id .= ";");
				$ensembl_id .= $1;
			}
		} elsif (/^GN   Name=/) {
			$_ =~ s/^GN   Name=//;
			($gene_name) = split /[;\s]/; # it ocurred to some idiot to put cross0references here
		} elsif (/^RC   TISSUE/) {
			$_ =~ s/^RC   TISSUE=//;
			my $tiss = $_;
			# their attempts to xlnk are just making things messy; sometimes we can forget to close the bracket too
			$tiss =~ s/\{.*?\}//g;
			$tiss =~ s/\{.*?$//g;
			$tiss =~ s/^.*?\}//g;
			(substr $tiss, -1) eq ',' && chop $tiss;
			(substr $tiss, -1) eq ';' || ($tiss.=";");
			($tissue =~ $tiss) ||  ($tissue .= $tiss);
		} elsif (/^CC/) {
			if (/\-!\-/) {
				if ($reading_function) {
					$fn =~ s/\{.*?\}//g;
					$reading_function  = 0;
				}

				if (/COFACTOR/) {
					$reading_cofactors  = 1;
					@cofs = ();
				} elsif (/FUNCTION/) {
					$reading_function  = 1;
					$_ =~ s/CC   -!- FUNCTION: //;
					$fn .= $_." ";
				} elsif (/SUBCELLULAR LOCATION/) {
					$reading_location  = 1;
					$_ =~ s/CC   -!- SUBCELLULAR LOCATION: //;
					$loc .= $_." ";
				}
			} elsif ( $reading_function )  {
				$_ =~ s/CC      //;
				$fn .= $_." ";
			} elsif ( $reading_cofactors )  {
				/Name\=(\S+?)\;/;
				defined $1 && push @cofs, $1;
				$reading_cofactors = 0; # one at a time
			}
		} elsif (/^SQ\s+SEQUENCE\s+(\d+)\s*AA/) {
			#$aa_lengths && ($aa_lengths.=";");
			# as of this writing, there is singe SQ field, presumably corresponding to the canonical sequence
			$aa_lengths || ($aa_lengths .= $1);
		}
	}
	$uniprot_ids =~ s/\s//g;
	my @aux = split (";", $uniprot_ids);
	my $uniprot_id = shift @aux;
	defined $manual_fix_for_ensembl{$uniprot_id} && ($ensembl_id = $manual_fix_for_ensembl{$uniprot_id});
	my $old_uniprot_ids  = join ";", @aux;
	my $cofactors = "";
	@cofs && ($cofactors = join ";", @cofs);
	print "$uniprot_id\t$gene_name\t$ensembl_id\t$ec_number\t$cofactors\t$aa_lengths\t$full_name\t$tissue\t$loc\t";
	print "$fn\t$old_uniprot_ids\n";


}
