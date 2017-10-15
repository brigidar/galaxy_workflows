#!/usr/bin/env python

#########################################################################################
#                                            #
# Name          :    rename.py                                #
# Version     : 0.1                                    #
# Project     :  SNPDEF                          #
# Description :  rename leafs       #
# Author      : Brigida Rusconi                                #
# Date        : Sept 26th, 2017                            #
#                                            #
#########################################################################################

import argparse, os, sys, csv,pdb
import Bio
from Bio import Phylo
import pandas
from pandas import *

parser = argparse.ArgumentParser()
parser.add_argument('-s', '--snp_tree', help="snp table to sort")
parser.add_argument('-n', '--names', help="list of names")
args = parser.parse_args()
input_file = args.snp_tree
names1=args.names

tree=Phylo.read(input_file,'newick')
nam=read_csv(names1,sep='\t',dtype=object,header=None)
nam1=dict(zip(nam[0],nam[1]))

for clade in tree.find_clades():
    if clade.name is None :
        continue
    else:
        clade.name=nam1[clade.name]

with open('tree.nwk','w') as output:
    Phylo.write(tree,output,'newick')
