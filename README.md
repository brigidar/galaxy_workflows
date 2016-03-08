# Galaxy Workflows

## Assembly and Annotation (workflow A)
The first workflow uses SPAdes v.3.7 for *de novo* assembly of reads and PROKKA for annotation. 
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
  



## Reads-based SNPDV (workflow B)
3 workflows in Galaxy to discover and quality control SNPs based on reads (Illumina). 
 * Install the adapted xml file for Bowtie2

The SNP discovery is reference based and you will need to define the mobilome within the reference genome to avoid SNP miscalling.
For the selected reference you will need a fasta and an annotated genbank file. For the query genomes you will need reads and 
assembled contigs.

### B1 Bowtie2 Mapping & Freebayes Discovery
If you are using MiSeq reads lower the coverage to 10 in the Freebayes options.

### B2 Excluded Regions 
If you have a close reference with known mobilome and repeated regions this step can be skipped. Provide a bed formatted file with the coordinates of the regions to exclude. If you don't know the mobilome use PHAST for phage region prediction and save regions as bed formatted file. Use ISFinder to predict IS elements. If your species carries any other mobile genetic elements please add the regions to the bed file (e.g. resistance cassettes). Once you have them use the **closed genomes excluded regions**.


#### Example bed file
NC_011353.1	273971	274038
NC_011353.1	275213	276349
NC_011353.1	302573	314525

If you are working with a draft genome: 
 * Identify IS elements with ISFinder and create a multifasta file
 * Find phage regions with PHAST and create multifasta file
 * Create plasmid multifasta to remove contigs that might match plasmids
 * Provide the fasta file of the reference
 * Take predicted SNPs from Freebayes of reference genome against itself
 
If your species does not have some of these mobile elements you can modify the workflow to remove them.

### B3 Verify SNPs
You will need a fasta and genbank file of the reference, the excluded regions, and the predicted SNPs
for the species of interest. You will need the assemblies of the query genomes as separate fasta files and concatenated.
Increase the threads to the total amount of query genomes to speed up the curation. The output is a filtered table and a multifasta with the curated SNPs for each query genome.

### B4 Post-processing (optional)
The SNPs multifasta is converted to phylip to run on PhyML. The filtered table is processed to provide genotypes and summary
of SNP characteristics.

## Contig-based SNPDV
You will need a genbank and fasta file for the reference genome and fasta files of the query genomes.

### C1 SNP discovery
SNP discovery is based on NUCmer with delta-filter and show-snps. Provide a list of query genomes (fasta files) to predict SNPs against reference genome (fasta file).


### C2 Excluded Regions 
If you have a close reference with known mobilome and repeated regions this step can be skipped. Provide a bed formatted file with the coordinates of the regions to exclude. If you don't know the mobilome use PHAST for phage region prediction and save regions as bed formatted file. Use ISFinder to predict IS elements. If your species carries any other mobile genetic elements please add the regions to the bed file (e.g. resistance cassettes). Once you have them use the **closed genomes excluded regions**.

#### Example bed file
NC_011353.1	273971	274038
NC_011353.1	275213	276349
NC_011353.1	302573	314525

If you are working with a draft genome: 
 * Identify IS elements with ISFinder and create a multifasta file
 * Find phage regions with PHAST and create multifasta file
 * Create plasmid multifasta to remove contigs that might match plasmids
 * Provide the fasta file of the reference
If your species does not have some of these mobile elements you can modify the workflow to remove them.

### C3 SNP Curation
You will need a fasta and genbank file of the reference, the excluded regions, and the predicted SNPs
for the species of interest. You will need the assemblies of the query genomes as separate fasta files and concatenated.
If you choose the threaded version you can increase the threads to the total amount of query genomes. The output is a 
filtered table and a multifasta with the curated SNPs for each query genome.
  * If you have both reads and draft genomes in our queries you can combine them at this step to get all the SNPs in one single merged table.
 The output is a filtered table and a multifasta with the curated SNPs for each query genome.

### C4 Post-processing (optional)
The SNPs multifasta is converted to phylip to run on PhyML. The filtered table is processed to provide genotypes and summary
of SNP characteristics.



