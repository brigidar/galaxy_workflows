#!/usr/bin/env python

#########################################################################################
#                                            #
# Name          :    rename_panel.py                                #
# Version     : 0.1                                    #
# Project     :  SNPDEF                          #
# Description :  rename snp panel       #
# Author      : Brigida Rusconi                                #
# Date        : Nov 5gh, 2017                            #
#                                            #
#########################################################################################

import argparse, os, sys, csv,pdb
import Bio
from Bio import Phylo
import pandas
from pandas import *

parser = argparse.ArgumentParser()
parser.add_argument('-s', '--snp_panel', help="snp table to sort")
parser.add_argument('-n', '--names', help="list of names to replace")
args = parser.parse_args()
input_file = args.snp_panel
names1=args.names

panel=read_csv(input_file,dtype=object,sep='\t')
names=read_csv(names1,dtype=object,sep='\t',header=None)
nam1=dict(zip(names[0],names[1]))

count_qbase=list(panel.columns.values)
qindexes=[]
qnames=[]
for i, v in enumerate(count_qbase):
    if 'qbase:' in v:
        qindexes.append(i)
        qnames.append(v.split(':')[1])

for i,v in enumerate(qnames):
    count_qbase[qindexes[i]]=':'.join(['qbase',nam1[v]])

panel.columns=count_qbase


with open('output.txt','w') as output2:
    panel.to_csv(output2, sep='\t',index=False)
