#!/usr/bin/perl
use strict;
use warnings;

my @tables = split "\n", `ls /home/ivana/scratch/gnomad_freqs_tmp`;


for my $table (@tables) {
    my @fields_tsv = split "\t", `tail -n1 /home/ivana/scratch/gnomad_freqs_tmp/$table`;
    my $last_pass_position_in_tsv = $fields_tsv[1];

    my $table_orig = $table;
    $table_orig =~ s/freqs_chr_/gnomad.exomes.r2.1.1.sites./;
    $table_orig =~ s/.tab/.liftover_grch38.vcf/;

    my $tail = 500;
    my $cmd = "tail -n $tail /storage/databases/gnomad/$table_orig | cut -f 1-7 | awk '\$7==\"PASS\" {print \$2}' | tail -n1";
    my $ret = `$cmd`;
    if ($ret) {
        my $last_pass_position_in_orig = $ret;
        chomp $last_pass_position_in_orig;
        print("$table ");
        if ($last_pass_position_in_tsv == $last_pass_position_in_orig) {
            print "$last_pass_position_in_tsv   $last_pass_position_in_orig  OK \n";
        } else {
            print("last pass pos not the same: $last_pass_position_in_tsv $last_pass_position_in_orig\n");
        }

    } else {
         print("$last_pass_position_in_tsv   no pass in the last $tail lines of the orig file\n");
    }
}