#!/bin/bash
awk -F '\t' '{if ($2 < $3) print $2 " " $3; else print $3 " " $2 }' BIOGRID-ORGANISM-Homo_sapiens-3.4.140.tab2.clean.txt  | sort | uniq > biogrid.whole.graph.txt
