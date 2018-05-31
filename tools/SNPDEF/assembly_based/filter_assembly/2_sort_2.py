#!/usr/bin/env python

#########################################################################################
#											#
# Name	      :	2_sort_2.py								#
# Version     : 0.7									#
# Project     : sort & merge SNP tables							#
# Description : Script to sort out no hits, indels, identical lines and double hits		#
# Author      : Brigida Rusconi								#
# Date        : May 31st, 2018							#
#											#
#########################################################################################
#for replacement of a given value with NaN
#http://stackoverflow.com/questions/18172851/deleting-dataframe-row-in-pandas-based-on-column-value

# to remove any symbol and replace it with nan
#http://stackoverflow.com/questions/875968/how-to-remove-symbols-from-a-string-with-python

# for isin information
#http://pandas.pydata.org/pandas-docs/stable/indexing.html

# for selecting rows that have an indel:
#http://stackoverflow.com/questions/14247586/python-pandas-how-to-select-rows-with-one-or-more-nulls-from-a-dataframe-without



#------------------------------------------------------------------------------------------
import argparse, os, sys, csv,pdb
#import pdb
import numpy
from numpy import *
import Bio
from pandas import *
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
from Bio import SeqIO
#------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------
# add transition and transversion information function developed by Mando Rodriguez

def get_trans_state(base, hit):
    
    if base == hit:
        return "--"
    elif base == 'A' and hit == 'G':
        return "transition"
    elif base == 'G' and hit == 'A':
        return "transition"
    elif base == 'C' and hit == 'T':
        return "transition"
    elif base == 'T' and hit == 'C':
        return "transition"
    
    elif base == 'A' and hit == 'C':
        return "transversion"
    elif base == 'C' and hit == 'A':
        return "transversion"
    
    elif base == 'A' and hit == 'T':
        return "transversion"
    elif base == 'T' and hit == 'A':
        return "transversion"
    
    elif base == 'C' and hit == 'G':
        return "transversion"
    elif base == 'G' and hit == 'C':
        return "transversion"
    
    elif base == 'G' and hit == 'T':
        return "transversion"
    elif base == 'T' and hit == 'G':
        return "transversion"
    
    else:
        return "--"
#------------------------------------------------------------------------------------------
def get_snp(table):
    snp_n=[]
    ident=[]
    for i,item in enumerate(table.index):
        #only append snp
        snp_u=[n for n in unique(table.iloc[i,:])[:] if n!=table.refbase[i]]
        if len(snp_u)>0:
            snp_n.append(snp_u)
        else:
            snp_n.append([n for n in unique(table.iloc[i,:])[:]])
            ident.append(i)
    return snp_n, ident

#------------------------------------------------------------------------------------------
def missing_char(str, pos,n):
    item=list(str)
    item[pos]=n
    str="".join(item)
    return str
#------------------------------------------------------------------------------------------

#output and input file name to give with the script
parser = argparse.ArgumentParser()
parser.add_argument('-o', '--output', help="fasta snps")
parser.add_argument('-s', '--snp_table', help="snp table to sort")
parser.add_argument('-t', '--total', help="inverted table to look at", default= "snp_filtered_table.txt")
parser.add_argument('-r','--remove', help="remove non-canonical nucleotides", default= "True")



args = parser.parse_args()
output_file = args.output
input_file = args.snp_table
output2_file = args.total
remove=args.remove

#------------------------------------------------------------------------------------------
#read in file as dataframe
df =read_csv(input_file,sep='\t', dtype=object)
df=df.set_index(['molecule','refpos']).fillna('--')

#------------------------------------------------------------------------------------------
# only columns with qbase and refbase in table
count_qbase=list(df.columns.values)

qindexes=[]
for i, v in enumerate(count_qbase):
    if 'qbase:' in v:
        qindexes.append(i)
df1=df.iloc[:,qindexes]

#------------------------------------------------------------------------------------------
#replaces lines with "No Hits" with NaN and removes lines with NaN in qbase columns

df=df.mask(df=='No Hit').dropna()

#print " No Hits and were replaced with N"

#removes lines that had no hit because alignment was too short
#df1=df1.mask(df1=='--').dropna()
ref=df['refbase']
df=df[df.index.isin(df1.index)]

#------------------------------------------------------------------------------------------

#longest alignment for each position
mindexes=[]
for i, v in enumerate(count_qbase):
    if 'maxlen:' in v:
        mindexes.append(i)
df_m=df.iloc[:,mindexes]


#multiple hits check for alignment length create dataframe with information
bindexes=[]
for i, v in enumerate(count_qbase):
    if 'blengths:' in v:
        bindexes.append(i)

df_l=df.iloc[:,bindexes]
df_l=df_l.replace({'No Hit':'0'},regex=True)

ind=[]
#check how many full alignments are present and drop position if more than one.
for i in range(0, df_l.index.size):
    #check if there is more than one hit
    if any(df_l.iloc[i,:].str.contains('/')):
        ind2=[]
        for v,x in enumerate(df_l.iloc[i,:]):
            #if there is more than one hit is it bad quality or good?
            if '/' in x:
                #count the number of times that a maximum length alignment is present.
                if x.split('/').count(df_m.iloc[i,v])!=1:
                    continue
                else:
                    ps=x.split('/').index(df_m.iloc[i,v])
                    #replace nucleotide with only valid nucleotide
                    df.iloc[i,v+2]=df.iloc[i,v+2].split('/')[ps]
                    ind2.append(i)
            else:
                ind2.append(i)
        if len(ind2)!=df_l.columns.size:
            #if not all positions remain with only one hit remove the location
            continue
        else:
            ind.append(i)

    else:
        ind.append(i)

#drop all positions that do not have an unique high quality blast hit.
df.reset_index(inplace=True)
df2=df[df.index.isin(ind)]
#remove blength and maxlength
df_l.reset_index(inplace=True)
df_l1=df_l[df_l.index.isin(ind)]
df2.set_index(['molecule','refpos'],inplace=True)
df2.drop(df2.columns[bindexes],axis=1,inplace=True)
f2.drop(df2.columns[mindexes],axis=1,inplace=True)

#------------------------------------------------------------------------------------------
# remove identical line

bases=['A','C','G','T','-']
df3=df2.iloc[:,qindexes].join(df2.iloc[:,1])
cols=df3.columns

# remove non-canononical snps (optional)
if remove=="True":
    slash3=df3[~df3[cols].isin(bases).all(axis=1)].dropna(how='all')
    df3=df3[~df3.index.isin(slash3.index)]
    for i in bases:
        id=df1
        id=id.replace({'N':i},regex=True)
        id=id[id !=i].dropna(how='all')
        df1=df1[df1.index.isin(id.index)]
else:
    for i in bases:
        slash2=slash2[slash2 !=i].dropna(how='all').fillna(i)
        for i in bases:
            id=df1
            id=id.replace({'N':i},regex=True)
            id=id[id !=i].dropna(how='all')
            df1=df1[df1.index.isin(id.index)]
    
#print "identical lines removed: SNP left " + str( slash2.index.size)

#------------------------------------------------------------------------------------------
#replaces lines with indel

indel=slash2.mask(slash2=='indel').dropna()


#------------------------------------------------------------------------------------------
#print "removed lines with short alignment %s lines left" % (str(indel.index.size))
final2 =indel.reset_index(drop=True).T #drops indexes of molecule and refpos
final2.reset_index(inplace=True)

#-------------------------------------------------------------------------------------
#save file to fasta

tab2=[]
for i in range(0,final2.index.size):
    tab=[]
    if ':' in final2.iloc[i,0]:
        tab.append(final2.iloc[i,0].split(':')[1])
        tab.append(''.join(final2.iloc[i,1:]))
    else:
        tab.append(final2.iloc[i,0])
        tab.append(''.join(final2.iloc[i,1:]))
    tab2.append(tab)

with open('table', 'w') as t:
    for i in range(0,len(tab2)):
        t.write("\t".join(tab2[i])+"\n")

with open('table','rU') as input:
    with open(output_file,'w') as output:
        sequences = SeqIO.parse(input, "tab")
        count = SeqIO.write(sequences, output, "fasta")

#------------------------------------------------------------------------------------------
#recalculate snps/gene gene length and dn/ds and transition/transversion

df=df[df.index.isin(indel.index)]

df.drop(['snps_per_gene','snps/gene_length'],axis=1,inplace=True)
temp=df.mask(df['gene_name']=='intergenic').dropna()
snps_gene=temp.groupby('gene_name').size().reset_index()

snps_gene2=dict(zip(snps_gene['gene_name'].tolist(),snps_gene[0].tolist()))
df['snps_per_gene']=df['gene_name'].map(snps_gene2)
df['snps/gene_length']=to_numeric(df['snps_per_gene'],errors='coerce').divide(to_numeric(df['gene_length'],errors='coerce'),axis=0)


pdb.set_trace()
#------------------------------------------------------------------------------------------
cod=df.dropna(subset=['gene_name'])
cod=cod.mask(cod['gene_name']=='intergenic').dropna()
cod.reset_index(inplace=True)

pdb.set_trace()
count_qbase2=list(cod.columns.values)
qindexes2=[]
for i, v in enumerate(count_qbase2):
    if 'qbase:' in v:
        qindexes2.append(i)

# get query base information
df3=cod.iloc[:,qindexes2].join(cod.refbase)
#position in codon
def mod(x):
    t= x % 3
        #last position in codon becomes 3 instead of 0
    if t==0:
        t=3
    #move back to 0-2 range to modify string in codon
    return t-1

cod['pos1']=cod.pos_in_gene.astype(int).apply(mod)
pos1=cod.pos1.tolist()
ref_codon=cod.ref_codon.astype(str).tolist()
#get allele for each position


snp_nb, idn =get_snp(df3)


query_codon=[]
for i,v in enumerate(snp_nb):
    if ref_codon[i]=='nan':
        query_codon.append('nan')
    else:
        ts=list()
        if len(v)==1 and v!='N':
            query_codon.append(missing_char(str(ref_codon[i]),pos1[i],v[0]))
        else:
            #multiallelic position gets codon for each
            for n in v:
                if n != 'N':
                    ts.append(missing_char(str(ref_codon[i]),pos1[i],n[0]))
            query_codon.append('/'.join(ts))


query_aa=[]
for i,v in enumerate(query_codon):
    qq=[]
    if v=='nan':
        query_aa.append('nan')
    else:
        if '/' in v:
            gl=v.split('/')
            for n in gl:
                qq.append(str(Seq(n, generic_dna).translate(table=11)))
                query_aa.append('/'.join(qq))
        else:
            query_aa.append(str(Seq(v, generic_dna).translate(table=11)))
cod.drop(['query_codon','query_aa','transition/transversion'],axis=1,inplace=True)
pdb.set_trace()
cod.insert(cod.columns.size,'query_codon',query_codon)
cod.insert(cod.columns.size,'query_aa',query_aa)

#print "Read query codons and aa"
#-------------------------------transition/transversion -------------------------------
#
df4=df.iloc[:,qindexes].join(df.refbase)


snp_nb2, ident=get_snp(df4)
ts_tv=[]
for i,v in enumerate(snp_nb2):
    ts=[]
    if len(v)==1 and v[0]!= 'N':
        ts_tv.append(get_trans_state(df.refbase[i],v[0]))
    else:
        for n in v:
            if n!='N':
                ts.append(get_trans_state(df.refbase[i],n))
        ts_tv.append('/'.join(ts))
df.drop(['query_codon','query_aa','transition/transversion','syn/nsyn/intergenic'],axis=1,inplace=True)
df.insert(df.columns.size,'transition/transversion',ts_tv)
#print "Read transition/transversion"


cod.set_index(['molecule','refpos'],inplace=True)
fin=df.join(cod.loc[:,['query_codon','query_aa']])
# -------------------------------synonymous nonsynonymous-------------------------------
query_aa=fin.query_aa.astype(str).tolist()
ref_aa=fin.ref_aa.astype(str).tolist()
syn=[]
for i,item in enumerate(query_aa):
    if item=='nan':
        if len(snp_nb2[i])==1:
            syn.append('intergenic')
        else:
            g=[n for n in snp_nb2[i] if n!='No Hit']
            syn.append('/'.join(repeat('intergenic',len(g))))
    elif '/' in item:
        mult=list()
        for n in item.split('/'):
            if n==ref_aa[i]:
                mult.append('SYN')
            else:
                mult.append('NSYN')
        syn.append('/'.join(mult))
    else:
        if item==ref_aa[i]:
            syn.append('SYN')
        else:
            syn.append('NSYN')

fin.insert(0,'syn?',syn)

#genes that are not CDS are tagged as genic
genic=dict()
genic['No CDS']='genic'
fin.reset_index(inplace=True)
fin.set_index('product',inplace=True)
fin.replace({'syn?':genic},inplace=True)
fin.reset_index(inplace=True)
fin.set_index(['molecule','refpos'],inpace=True)
#SNPs that are actually identical are replaced with No SNP
fin['syn?'][ident]='No SNP'
fin.reset_index(inplace=True)
#-------------------------------dn/ds-------------------------------
dn=fin.groupby(['gene_name','syn?']).size().reset_index()
dn_2=dn[dn['syn?'].str.contains('SYN')]
dn_2.rename(columns={0:'count'},inplace=True)
dn_ds=dict()

for i,v in enumerate(dn_2['gene_name']):
    t=dn_2[dn_2['gene_name']==v]
    if any(t['syn?'].str.contains('/')):
        ns=dict()
        try:
            ns['NSYN']=t[t['syn?']=='NSYN']['count'].values[0]
        except IndexError:
            ns['NSYN']=0
        try:
            ns['SYN']=t[t['syn?']=='SYN']['count'].values[0]
        except IndexError:
            ns['SYN']=0
        bl=t[t['syn?'].str.contains('/')]
        if bl.index.size==1:
            c=bl.iloc[0,1].split('/')
            for d in c:
                if d=='NSYN':
                    ns['NSYN']=(ns['NSYN']+bl.iloc[0,2])
                else:
                    ns['SYN']=(ns['SYN']+bl.iloc[0,2])
        else:
            for g in range(0,bl.index.size):
                c=bl.iloc[g,1].split('/')
                for d in c:
                    if d=='NSYN':
                        ns['NSYN']=(ns['NSYN']+bl.iloc[g,2])
                    else:
                        ns['SYN']=(ns['SYN']+bl.iloc[g,2])
        dn_ds[v]=ns['NSYN']/ns['SYN']
    else:
        try:
            dn_ds[v]=t[t['syn?']=='NSYN']['count'].values[0] / t[t['syn?']=='SYN']['count'].values[0]
        except IndexError:
            dn_ds[v]=0

fin['dn_ds']=fin['gene_name'].map(dn_ds)
fin1=fin.iloc[:,1:(max(qindexes)+3)]
fin1.set_index(['molecule','refpos'],inplace=True)
fin2=fin.reindex_axis(['molecule','refpos','gene_name','gene_start','gene_end','gene_length','pos_in_gene','ref_codon','ref_aa','query_codon','query_aa','product','transition/transversion','snps_per_gene','snps/gene_length','dn_ds'],axis=1)#'dn/ds'
fin2.set_index(['molecule','refpos'],inplace=True)
final=fin1.join(fin2)
final.reset_index(inplace=True)
final.sort_values(by=['molecule','refpos'],inplace=True)

#------------------------------------------------------------------------------------------
#save total file for plotting -t option
with open(output2_file,'w') as output2:
    final.to_csv(output2, sep='\t', index=False)










