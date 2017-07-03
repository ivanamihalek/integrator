#!/usr/bin/perl
use strict;
use warnings FATAL => 'all';

## ===========================================================================
## adapted from web_blast.pl, NCBI's own script
## could not get this to take in HITLIST_SIZE and format parameters

use URI::Escape;
use LWP::UserAgent;
use HTTP::Request::Common qw(POST);

my $ua = LWP::UserAgent->new;

@ARGV>=1 || die "Usage: $0 <fasta file>.\n";

my $program  = "blastp";
my $database = "refseq_protein";
my $query_file = $ARGV[0];
my $evalue = 1.0e-20;
my $max_num_seqs = 1000;
# read and encode the query
my $encoded_query = "";
open(QUERY,"<$query_file") || die "Cno $query_file: $!.\n";
while(<QUERY>){
    $encoded_query = $encoded_query . uri_escape($_);
}

# build the request
my $args = "CMD=Put&PROGRAM=$program&DATABASE=$database&FORMAT_TYPE=Tabular&";
$args .= "HITLIST_SIZE=$max_num_seqs&EXPECT=$evalue&";
print "$args\n"; exit;
$args .= "QUERY=" . $encoded_query;

my $req = new HTTP::Request POST => 'https://blast.ncbi.nlm.nih.gov/blast/Blast.cgi';
$req->content_type('application/x-www-form-urlencoded');
$req->content($args);

# get the response
my $response = $ua->request($req);

# parse out the request id
$response->content =~ /^    RID = (.*$)/m;
my $rid=$1;

# parse out the estimated time to completion
$response->content =~ /^    RTOE = (.*$)/m;
my $rtoe=$1;
print "estimated time to completion: $rtoe secs\n";
# wait for search to complete
sleep $rtoe;
print "polling for results ... \n";
# poll for results
while (1) {
    sleep 1;

    $req = new HTTP::Request GET => "https://blast.ncbi.nlm.nih.gov/blast/Blast.cgi?CMD=Get&FORMAT_OBJECT=SearchInfo&RID=$rid";
    $response = $ua->request($req);

    if ($response->content =~ /\s+Status=WAITING/m)
        {
        # print STDERR "Searching...\n";
        next;
        }

    if ($response->content =~ /\s+Status=FAILED/m)
        {
        print STDERR "Search $rid failed; please report to blast-help\@ncbi.nlm.nih.gov.\n";
        exit 4;
        }

    if ($response->content =~ /\s+Status=UNKNOWN/m)
        {
        print STDERR "Search $rid expired.\n";
        exit 3;
        }

    if ($response->content =~ /\s+Status=READY/m)
        {
        if ($response->content =~ /\s+ThereAreHits=yes/m)
            {
            #  print STDERR "Search complete, retrieving results...\n";
            last;
            }
        else
            {
            print STDERR "No hits found.\n";
            exit 2;
            }
        }

    # if we get here, something unexpected happened.
    exit 5;
    } # end poll loop

# retrieve and display results
$req = new HTTP::Request GET => "https://blast.ncbi.nlm.nih.gov/blast/Blast.cgi?CMD=Get&FORMAT_TYPE=Text&RID=$rid";
$response = $ua->request($req);

print $response->content;
exit 0;
1;