#!/usr/bin/env python
#adapted from Lee Bergstrand Seqextract.py script
#https://github.com/LeeBergstrand/Genbank-Downloaders/blob/master/SeqExtract.py

#AVRS01000000

import Bio
from Bio import Entrez
from Bio import SeqIO
import argparse, os, sys
import re


#import pdb

parser = argparse.ArgumentParser()
parser.add_argument('-g', '--genomes', help="ncbi id")
parser.add_argument('-e','--emails', help='email address')

args = parser.parse_args()
genomes=args.genomes
emails=args.emails

#-------------------------------------------------
def isSSProject(sequence):
    is_ssp = False
    
    if sequence.annotations.get('wgs'):
        is_ssp = True
    
    return is_ssp
#-------------------------------------------------

def getSeqRecords(seqList):
    try:
        handle = Entrez.efetch(db="nucleotide", id=seqList, rettype="gb",retmode="genbank")  # Gets records and stores them.
        SeqRecords = list(SeqIO.parse(handle, "genbank"))  # Creates a list of SeqRecord objects from genbank files e.
        handle.close()  # Closes handle since it is no longer needed.
    except IOError:
        #print("Failed to connect to NCBI server. ")
        sys.exit(1)
    return SeqRecords
#-------------------------------------------------


AccessionBaseRegex = re.compile("^[a-zA-Z]{4}\d{2}")
WGSSProjectRegex = re.compile("[a-zA-Z_]{4,7}\d{8,10}")

Entrez.email= emails

records= getSeqRecords(genomes)

contigList = []
for item in records:
    if isSSProject(item):
        itemrange = item.annotations["wgs"]
        if len(itemrange) == 2:
    # Extracts the accession base from first contig accession number.
            AccessionBase = (AccessionBaseRegex.findall(itemrange[0])[0])
        
        # Takes both the the min and max accession and slices off (using python's string slice syntax s[start:end:step])
        # the accession base code leaving the numerical difference between the contigs.
        # Converts these differences to integers.
            itemrangeMin = int(itemrange[0][6:])
            itemrangeMax = int(itemrange[1][6:])
# WGSS accession number length actually varies. Its normally 12 characters but I have seen 13 before.
# The code block below accounts for this.
            zeroOffset = 6
            accessionLength = len(itemrange[0])
            if accessionLength != 12:
                zeroOffset = accessionLength - 6  # 6 is the length of the standard accession base.

            # Creates accession list
            for x in range(itemrangeMin, (itemrangeMax + 1)):
                contigAccession = AccessionBase
                contigAccession += ("{0:0" + str(zeroOffset) + "d}").format(x)  # Uses zero offset to make accessions proper length.
                contigList.append(contigAccession)
        else:
            contigList.append(itemrange[0])# If one contig, simply append it to the list.

    else:
        with open("output.gb",'w') as output2:
            output2.write(item.format("genbank"))
        with open("output.fasta",'w') as output2:
            output2.write(item.format("fasta"))

if len(contigList)>0:
    contigrecords=getSeqRecords(contigList)
    with open("output.gb",'w') as output2:
        for item in contigrecords:
            output2.write(item.format("genbank"))

    with open("output.fasta",'w') as output2:
        for item in contigrecords:
            output2.write(item.format("fasta"))

