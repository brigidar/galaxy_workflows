# SNPDV

## Galaxy Setup

### System Requirements
Linux or Unix based only. Some of the tools used in these workflows were developed by others and do not come in binaries for Windows.

Before installing Galaxy please check your python version (2.7 or higher). You will also need the following packages: Numpy, pandas, matplotlib, and biopython. We recommend installing [anaconda python] (https://docs.continuum.io/anaconda/install), which comes with all the packages except biopython that needs to be installed with the command "conda install biopython". 

Galaxy is an open, web-based platform for accessible, reproducible, and transparent computational biomedical research. 
To install Galaxy please refer to the [Galaxy Wiki] (https://wiki.galaxyproject.org/FrontPage).
Galaxy can be run locally, on a server (cluster or single server), and on the cloud. If you are planning on running Galaxy on the cloud, we recommend you select a linux-based operating system.

You can retrieve all the tools and workflows required from the Galaxy toolshed under the repository SNPDV. This can be done from the browser interface directly. For more information on how to install tools, please refer to the [Galaxy Toolshed Wiki] (https://wiki.galaxyproject.org/ToolShed).

## Workflows

### Assembly and Annotation (workflows A)
The first workflow uses SPAdes v.3.7 for *de novo* assembly of reads and PROKKA for annotation. 
  * Please install the modified xml files for SPADes and PROKKA to allow for renaming of files in the batch submission mode. 
  * Reads should be imported as fastqsanger format and grouped into a paired dataset list.
  * If you are assembling longer reads please change the k-mer size in the SPADES option
     * Reads 250 bp or more use k-mer 21,33,55,77,99,127
   * For more information on SPAdes please refer to the [manual] (http://spades.bioinf.spbau.ru/release3.7.0/manual.html)
   
  * Make a tab delimited file to import the annotation options for PROKKA.
     * You only need to provide a value of 1 for kingdom if you are annotating viruses.
   * For more information on PROKKA please refer to the [github repository] (https://github.com/tseemann/prokka)
   * We recommend building a genus database of curated annotations to install in your Galaxy instance. This will improve the annotation by preferentially using the genus database during annotation.

##### Example annotation file  
    strain	 locustag	 centre	 genus	      species	plasmid	kingdom
    strainA	 xxx	     centre	 Escherichia	coli		
    strainB	 yyy	     centre  Escherichia	coli	   p90
  

## Excluding Mobilome Regions in SNPDV

###**Excluding with a closed reference:**
If you have a close reference with known mobilome and repeated regions this step can be skipped. Provide a bed formatted file with the coordinates of the regions to exclude. If you don't know the mobilome use PHAST for phage regions prediction and save regions as bed formatted file. Use ISFinder to predict IS elements. If your species carries any other mobile genetic elements please add the regions to the bed file (e.g. resistance cassettes). Once you have them use the **closed genomes excluded regions**.

##### Example bed file
    NC_011353.1	273971	274038
    NC_011353.1	275213	276349
    NC_011353.1	302573	314525
    
### Reads-based discovery SNPDV (workflows B)
3 workflows in Galaxy to discover and quality control SNPs based on reads (Illumina). 
 * Install the adapted xml file for Bowtie2

The SNP discovery is reference based and you will need to define the mobilome within the reference genome to avoid SNP miscalling.
For the selected reference you will need a fasta and an annotated genbank file. For the query genomes you will need reads and 
assembled contigs.

### B1 Bowtie2 Mapping & Freebayes Discovery
If you are using MiSeq reads lower the coverage to 10 in the Freebayes options.






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

## Contig-based SNPDV (wrokflows C)
You will need a genbank and fasta file for the reference genome and fasta files of the query genomes.

### C1 SNP discovery
SNP discovery is based on NUCmer with delta-filter and show-snps. Provide a list of query genomes (fasta files) to predict SNPs against reference genome (fasta file).


### C2 Excluded Regions 
If you have a close reference with known mobilome and repeated regions this step can be skipped. Provide a bed formatted file with the coordinates of the regions to exclude. If you don't know the mobilome use PHAST for phage region prediction and save regions as bed formatted file. Use ISFinder to predict IS elements. If your species carries any other mobile genetic elements please add the regions to the bed file (e.g. resistance cassettes). Once you have them use the **closed genomes excluded regions workflow**.

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

#### Example filtered SNP table
    molecule	    refpos	 syn?	 refbase	 qbase:abht	 gene name	          gene start	 gene end	 gene length	 snps per gene	 pos in gene/	
    NZ CM000662	    101797	 NSYN	 A	         T	         ESCCO14588_RS26480	  101814	     101644	     171	         1	             18//	
    
    //ref codon  ref aa      query codon	query aa	 transition/transversion	//		maxlen:abht	 blengths:abht	 product
      TTA        L           TTT	        F	         transversion	           // 	    40	         40	             hypothetical

The table not only provides the SNP, but also the corresponding annotation from the reference genome. It provides the length of the
blastn hit(maxlen), as well if there are additional hits of lower quality (blenghts).

### C4 Post-processing (optional)
The SNPs multifasta is converted to phylip to run on PhyML. The filtered table is processed to provide genotypes and summary
of SNP characteristics.



