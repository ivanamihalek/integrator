#! /usr/bin/perl
use strict;
use warnings FATAL => 'all';

# Therrel tabels give number of cases and total number of cases screened
# for each state (see the original paper)

my $table_dir = "/databases/therrell";

my %incidence = ();
for my $table (split ("\n", `ls $table_dir`)) {

    open(INF, "<$table_dir/$table") || die "Cno $table_dir/$table: $!";
    while (<INF>) {
        /\S/ || next;
        /^PMC full text/ && next;
        my @fields = split "\t";
        my $disease = shift @fields;
        defined $incidence{$disease} || (@{$incidence{$disease}}=());
        for my $i (map { 2 * $_ } 0 .. $#fields/2) {
            my $nr_cases = $fields[$i];
            $nr_cases =~ s/\s//g;
            $nr_cases =~ /\D/ && next;
            my $nr_screens = $fields[$i + 1];
            $nr_screens =~ s/[\s\,]//g;
            $nr_screens =~ /\D/ && next;
           push @{$incidence{$disease}}, $nr_cases."_".$nr_screens;
        }
    }
    close INF;

}

for my $disease (keys %incidence) {

    #print "$disease\n";
    my $tot_cases   = 0;
    my $tot_screens = 0;
    my $avg = 0;
    my $avg_sq = 0;
    my $min = 5000;
    my $max = -1;
    my $nr_states = 0;
    for my $store_string (@{$incidence{$disease}}) {
        my ($nr_cases, $nr_screens) = split "_", $store_string;
        $tot_cases   += $nr_cases;
        $tot_screens += $nr_screens;
        $avg         += $nr_cases/ $nr_screens;
        $avg_sq      += ($nr_cases/$nr_screens)**2;
        #printf "\t  %s  %s  \n", $nr_cases, $nr_screens; #,  $nr_cases/$nr_screens;
        $nr_cases==0 && ($nr_cases=1);
        my $permil = $nr_cases/$nr_screens*1.e6;
        ($min>$permil) &&  ($min=$permil);
        ($max<$permil) &&  ($max=$permil);
        $nr_states += 1;
    }
    $avg /= $nr_states;
    $avg_sq /= $nr_states;

    my $stdev = sqrt($avg_sq - $avg*$avg);
    #printf "** %15s  %5d  %10d          %4.1f        %4.1f   %4.1f        %4.1f  %4.1f \n",
    #    $disease , $tot_cases,$tot_screens,  $tot_cases/$tot_screens*1.e6,  $avg*1.e6, $stdev*1.e6, $min, $max;
    printf  "update therrell_incidences set cases=%d, screens=%d, ",  $tot_cases, $tot_screens;
    printf  "avg_incidence_per_mil_us_states=%4.1f, ", $avg*1.e6;
    printf  "stdev_incidence_per_mil_us_states=%4.1f  where name_short='%s';\n", $stdev*1.e6, $disease;
}

1;