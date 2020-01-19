#!/usr/bin/perl
use strict;
use warnings;

sub annovar (@) {

    for ("/home/ivana/third/annovar/table_annovar.pl", "/home/ivana/third/annovar/humandb/") {
        -e $_ || die "$_ not found.\n"
    }

    my $filename = $_[0];

    my $tmp_filename = "annov$$"; # $$ is the process id
    my $cmd = "/home/ivana/third/annovar/table_annovar.pl $filename  /home/ivana/third/annovar/humandb/ ";
    $cmd .= " -buildver hg19 -out  $tmp_filename -protocol refGene  -operation g  -nastring . -vcfinput ";
   
    system ($cmd) && die "error running $cmd: $!\n";
    my $new_filename = $filename;
    $new_filename=~ s/\.vcf/.annov.vcf/; # annovar will add .vcf
    my $annovar_filename =  `ls $tmp_filename*.vcf`; chomp $annovar_filename;

    $annovar_filename || die "no annovar output found (expected $tmp_filename*.vcf)\n";
    `mv $annovar_filename $new_filename`;
    `rm $tmp_filename*`;

    return $new_filename;

}

@ARGV >= 1 || die "Usage $0 <vcf filename>.\n";
my $vcf = $ARGV[0];
(-e $vcf) || die "$vcf not found.\n";
(-z $vcf) && die "$vcf seems to be empty.\n";

annovar($vcf);
