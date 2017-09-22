#!/bin/bash

esearch -db nucleotide -query "$1 [Organism]AND $2 [All Fields] AND Plasmid[filter]" | efetch --format fasta > output.fasta
