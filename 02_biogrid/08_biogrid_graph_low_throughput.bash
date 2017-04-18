#!/bin/bash

# not that for daddy, this is replaced my a more careful consideration using blimps's own
# rake db:populate:mark_problematic_biogrid
# rake db:output:biogrid_edges
awk -F '\t' '$18~"Low" {if ($2 < $3) print $2 " " $3; else print $3 " " $2 }' BIOGRID-ORGANISM-Homo_sapiens-3.4.142.tab2.clean.txt  | sort | uniq > biogrid.low_thruput_graph.txt
