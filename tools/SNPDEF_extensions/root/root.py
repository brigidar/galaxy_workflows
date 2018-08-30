#!/usr/bin/env python

#########################################################################################
#                                            #
# Name          :    root.py                                #
# Version     : 0.3                                    #
# Project     :  SNPDEF                          #
# Description :  change location of bootstrap and root      #
# Author      : Brigida Rusconi                                #
# Date        : Aug 30th, 2018                            #
#                                            #
#########################################################################################

import argparse, os, sys, csv,pdb
import Bio
from Bio import Phylo

parser = argparse.ArgumentParser()
parser.add_argument('-t', '--tree', help="snp table to sort")
parser.add_argument('-s', '--start', help="first leaf name")
parser.add_argument('-f', '--stop', help="second leaf name")
args = parser.parse_args()
input_file = args.tree
start=args.start
stop=args.stop

tree=Phylo.read(input_file,'newick')
if start="mid":
    tree.root_at_midpoint()
else:
#pdb.set_trace()
    tree.root_with_outgroup(start,stop)
pdb.set_trace()
with open('tree1.nwk','w') as output:
    Phylo.write(tree,output,'newick',format_branch_length='%1.10f')


