#! /usr/bin/perl -w

# store as uniprot_basic_infos.csv
# and load with
#  mysqlimport --local -u root -p monogenic_development  uniprot_seqs.csv
use strict;
use warnings;
sub parse();


my $filename = "/databases/uniprot/uniprot_sprot.dat";
my $exon_coords_file = "/databases/ucsc/ensGene.txt";
my $seqfile = "/databases/uniprot/blast/uniprot_sprot.fasta";
my $blastextract = "/usr/local/bin/blastdbcmd";
foreach ($filename, $exon_coords_file, $seqfile, $blastextract) {
	-e $_ || die "$_ not found\n";
}

open (IF, "<$filename" )
	|| die "Cno $filename: $!.\n";


my %manual_fix_for_ensembl = ();
#not sure what is going on here:
#  O43826 is annotated an manually reviewed entry for SLC37A4, associated with ENSG00000281500 --> patch in ensembl GRCh38
# U3KQS2 is unreviewd, un-annotated, also says SLC37A4, and associated with ENSG00000137700 ensembl GRCh37
$manual_fix_for_ensembl{'O43826'} = 'ENST00000545985;ENST00000357590;ENST00000330775;ENST00000538950';

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
	my ($uniprot_ids, $ensembl_transcript_ids) = ("","","");
	my $reading_function  = 0;
	foreach (split "\n", $entry) {
		if  (/^ID/) {
		} elsif (/^AC/) {
			# it looks like the first iD is current,
			# and the rest are the older ids
			# # AC can also stretch through multiple lines
			$_ =~ s/^AC//;
			$uniprot_ids .= $_;
		} elsif (/^DR   Ensembl.*(ENST\d{11}).*/){
			my $enst = $1;
			/(ENSP\d{11}).*/;
			my $ensp  = $1;
			my $combo = $enst."=".$ensp;
			if ($ensembl_transcript_ids !~ $combo) {
				$ensembl_transcript_ids && ($ensembl_transcript_ids .= ";");
				$ensembl_transcript_ids .= $combo;
			}
		}
	}
	$uniprot_ids =~ s/\s//g;
	my @aux = split (";", $uniprot_ids);
	my $uniprot_id = shift @aux;

	defined $manual_fix_for_ensembl{$uniprot_id} && ($ensembl_transcript_ids = $manual_fix_for_ensembl{$uniprot_id});
    foreach my $combo (split ";",$ensembl_transcript_ids) {
		my ($enst, $ensp) = split ("=", $combo);
		$ensp && $ensp ne"" || next;
		my $ret = `grep $enst $exon_coords_file`;
		my @field = split "\t", $ret;
		my ($chrom, $strand, $cds_start, $cds_end, $exon_starts, $exon_ends) =
			($field[2],  $field[3] ,$field[6],  $field[7], $field[9],  $field[10]);
		$exon_starts && $exon_starts ne"" || next;
		$exon_ends && $exon_ends ne"" || next;
		my $seq = `$blastextract -db $seqfile -entry $uniprot_id | grep -v $uniprot_id`;
		$seq eq ""  && next;
		$seq =~ s/\n//g;
		(substr $exon_starts, -1)eq',' && chop($exon_starts);
		(substr $exon_ends, -1)eq',' && chop($exon_ends);
		$chrom =~ s/chr//;

		print join "\t", ($uniprot_id,  $enst,  $ensp, $chrom, $strand,
				$exon_starts, $exon_ends, $cds_start, $cds_end, $seq,"","",""); # I have some extr fields in the database table
		print "\n";
	}

}
