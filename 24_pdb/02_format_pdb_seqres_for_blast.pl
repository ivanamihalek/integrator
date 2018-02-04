#!/usr/bin/perl
use strict;
use warnings;
use feature qw(say);

@ARGV==1 ||  die "Usage: $0 <pdbseqres file, full path>\n";

my @file_path = split ('/', $ARGV[0]);
my $filename = pop @file_path ;
my $dirpath = join ('/', @file_path);


print "  $dirpath   $filename \n";

# check filename has extension txt
my @name_components = split '\.',$filename;
my $extension = pop @name_components;
# otherwise require so
$extension eq 'txt' || die "I expected extension 'txt' in pdbseqres file; instead got '$extension'.\n";

# make new filename
push (@name_components, 'fasta');
my $newname = join ('.', @name_components);
my $new_filename_full = $dirpath."/".$newname;

# make filename for the original headers
pop @name_components;
push (@name_components, 'orig_headers');
my $new_headers_full = $dirpath."/".join ('.', @name_components);

# make filename for the name resolution table
pop @name_components;
push (@name_components, 'name_resolution');
my $name_resolution =  $dirpath."/".join ('.', @name_components);

#################
open (INF, "<$ARGV[0]")  || die "Cno $ARGV[0]: $!\n";
open (OUTF,">$new_filename_full")  || die "Cno $new_filename_full: $!\n";
open (OUT_HDRS,">$new_headers_full")  || die "Cno $new_headers_full: $!\n";
open (OUT_NAMERES,">$name_resolution")  || die "Cno $name_resolution: $!\n";

my %seen = (); # check for names that would be duplicates if the case is ignored
while (<INF>) {
    if (/^>/) {
        if (/^>([\d\w]{4})_([\d\w]*)[\s\n]/) {
            print OUT_HDRS $_;
            # $1 is pdbname;
            # $2 is chain;
            # change format to pdb\|entry\|chain 	(example: pdb\|1I4L\|D)
            # as described here: https://ncbi.github.io/cxx-toolkit/pages/ch_demo#ch_demo.T5
            # the format is described as requirement here https://www.ncbi.nlm.nih.gov/books/NBK279688/
            # print OUTF  ">pdb|$1|$2\n";
            #
            # the problem (Jan 2018): if formated as PDB, makeblastdb rejects two-character chain names
            # if formated as local (lcl) makeblastdb ignores case, and complains that, for example,
            # 3j3q_ga, and 3j3q_gA are duplicated (both 3J3Q_GA in makeblastdb book)
            #
            # as some makeshift solution, rename seqs and make name resolution file
            my $orig_name = "$1\_$2";
            my $new_name = uc "$1$2";
            if (defined $seen{$new_name}) {
                $seen{$new_name} += 1;
            } else {
                $seen{$new_name} = 1;
            }
            $new_name .= "_".$seen{$new_name};

           print OUT_NAMERES "$orig_name   $new_name\n";
            print OUTF ">$new_name\n";
        } else {
            die "Unexpected format for the sequence name:\n$_\n";
        }
    } else {
        print OUTF $_;
    }
}

close INF;
close OUTF;
close OUT_HDRS;
close OUT_NAMERES;

`/usr/bin/makeblastdb -in $new_filename_full -parse_seqids -dbtype prot`;