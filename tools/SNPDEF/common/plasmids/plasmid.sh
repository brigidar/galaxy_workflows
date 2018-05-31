#!/bin/bash

esearch -db nucleotide -query "$1 [Organism] AND Plasmid[filter]" | efetch --format fasta > output.fasta
