#!/usr/bin/perl
use strict;
use warnings FATAL => 'all';

my $repo            = "/home/ivana/monogenic/public/alignments";
my $afa2msf         = "/home/ivana/pypeworks/integrator/integrator_utils/afa2msf.pl";
my $msf2afa         = "/home/ivana/pypeworks/integrator/integrator_utils/msf2afa.pl";
my $restrict_to_qry = "/home/ivana/pypeworks/integrator/integrator_utils/restrict_msf_to_query.pl";

for ($repo,$afa2msf,$msf2afa,$restrict_to_qry ) {
    -e $_ || die "$_ not found.\n";
}

chdir $repo;

for my $uni_fa ( split "\n", `ls *.fa` ) {
    my $uni = $uni_fa;
    $uni =~ s/\.fa//;
    my $afa = "$uni.patched.afa";
    my $afa_with_uni = "$uni.with_uni.afa";
    `muscle -profile -in1 $uni_fa -in2 $afa -out $afa_with_uni`;
    my $msf_with_uni = "$uni.with_uni.msf";
    `$afa2msf $afa_with_uni > $msf_with_uni`;
    my $uni_restr = "$uni.restr.msf";
    `$restrict_to_qry  $msf_with_uni $uni > $uni_restr`;
    my $uni_restr_afa = "$uni.restr.afa";
    `$msf2afa $uni_restr > $uni_restr_afa`;
    `rm -f $afa_with_uni $msf_with_uni $uni_restr`
}

1;