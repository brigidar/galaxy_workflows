# Galaxy Workflows

## Assembly and Annotation (workflow A)
The first workflow uses SPAdes v.3.7 for de novo assembly of reads and PROKKA for annotation. 
  * Please install the modified xml files for SPADes and PROKKA to allow for renaming of files in the batch submission mode. 
  * Reads should be imported as fastqsanger format and grouped into a paired dataset list.
  * If you are assembling longer reads please change the k-mer size in the SPADES option
     * Reads 250 bp or more use k-mer 21,33,55,77,99,127
  * Make a tab delimited file to import the annotation options for PROKKA.
     * You only need to provide a value of 1 for kingdom if you are annotating viruses.

### Example annotation file  
    strain	 locustag	 centre	 genus	      species	plasmid	kingdom
    strainA	 xxx	     centre	 Escherichia	coli		
    strainB	 yyy	     centre  Escherichia	coli	  p90
  



## Reads-based SNP Discovery (workflow B)
This repository contains multiple workflows for Galaxy to discover and quality control SNPs based on reads (Illumina). 
 * Install the adapted xml file for Bowtie2

The SNP discovery is reference based and you will need to define the mobilome within the reference genome to avoid SNP miscalling.
For the selected reference you will need a fasta and an annotated genbank file. For the query genomes you will need reads and 
assembled contigs.

### B1 Bowtie2 Mapping & Freebayes discovery
If you are using MiSeq reads lower the coverage to 10 in the Freebayes options.

### B2 Excluded regions
If you have a close reference with known mobilome regions this step can be skipped. Provide a bed formatted file with the coordinates
of the rgions to exclude. 
If you are working with a draft genome: 
 * Identify IS elements with ISFinder and create a multifasta file
 * Find phage regions with PHAST and create multifasta file
 * Create plasmid multifasta to remove contigs that might match plasmids
 * Provide the fasta file of the reference
If your species does not have some of these mobile elements you can modify the workflow to remove them.

### B3 Verify SNPs
You will need a fasta and genbank file of the reference, the excluded regions, and the predicted SNPs.
for the species of interest. 
