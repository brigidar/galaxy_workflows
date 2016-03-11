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

###### Example annotation file  
    strain	 locustag	 centre	 genus	      species	plasmid	kingdom
    strainA	 xxx	     centre	 Escherichia	coli		
    strainB	 yyy	     centre  Escherichia	coli	   p90
  

### SNPDV (workflows B)

#### Excluding Mobilome Regions in SNPDV

#####Excluding with a closed reference:
If you have a close reference with known mobilome and repeated regions this step can be skipped. Provide a interval formatted file with the coordinates of the regions to exclude. If you don't know the mobilome use [PHAST] (http://phast.wishartlab.com/) for phage regions prediction and save regions as interval formatted file. Use [ISFinder] (https://www-is.biotoul.fr/) to predict IS elements for your species of interest, this needs to be done only once per species. If your species carries any other mobile genetic elements please add the regions to the interval file (e.g. resistance cassettes). Use the workflow **"Excluded Regions (closed reference)"**.

#####Excluding with a draft genome reference (no reads):
If you have a draft genome reference without reads you can predict the mobilome with [PHAST] (http://phast.wishartlab.com/), but you will need to save the predicted regions as multifasta, since the PHAST intervals correspond to the concatenated draft genome and not the intervals on each contig. For IS elements you can use  [ISFinder] (https://www-is.biotoul.fr/).If your species carries any other mobile genetic elements please add the regions to the interval file (e.g. resistance cassettes). Since it is a draft genome, any plasmids present will be mixed in with the contigs. If there are known plasmids to be present in all strains  provide a multifasta to exclude them. Use the workflow **"Excluded Regions (draft genome)"**.

#####Excluding with a draft genome reference (reads):
If you have a draft genome reference with reads you can predict the mobilome with [PHAST] (http://phast.wishartlab.com/), but you will need to save the predicted regions as multifasta, since the PHAST intervals correspond to the concatenated draft genome and not the intervals on each contig. For IS elements you can use  [ISFinder] (https://www-is.biotoul.fr/).If your species carries any other mobile genetic elements please add the regions to the bed file (e.g. resistance cassettes). In addition you will add the Freebayes output for the reference reads against itself to exclude any base miscalling of the assembly. Since it is a draft genome, any plasmids present will be mixed in with the contigs. If there are known plasmids to be present in all strains  provide a multifasta to exclude them Use the workflow **"Excluded Regions (reads)"**.


####### Example interval file
    NC_011353.1	273971	274038
    NC_011353.1	275213	276349
    NC_011353.1	302573	314525
    
#### Reads-based discovery SNPDV (workflow B.1.1)

For the selected reference you will need a fasta and an annotated genbank file. For the query genomes you will need reads and assembled contigs. SNPs will be predicted by Bowtie2 mapping & Freebayes. If you are using MiSeq reads lower the coverage to 10 in the Freebayes options.

#### Contig-based Discovery SNPDV (workflow B.1.2)
If no reads are available for query genomes, SNP discovery is based on NUCmer with delta-filter and show-snps. Provide a list of query genomes (fasta files) to predict SNPs against reference genome (fasta file).

#### SNP validation (workflow B.2)
You will need a fasta and genbank file of the reference, the excluded regions (prepared according to reference genome), and the predicted SNPs for the species of interest. You will need the assemblies of the query genomes as separate fasta files and concatenated.
Increase the threads to the total amount of query genomes to speed up the curation.  If you have both reads and draft genomes in your queries you can combine them at this step to get all the SNPs in one single merged table. The output is a filtered table and a multifasta with the curated SNPs for each query genome.


#### Example filtered SNP table
    molecule	    refpos	 syn?	 refbase	 qbase:abht	 gene name	          gene start	 gene end	 gene length	 snps per gene	 pos in gene/	
    NZ CM000662	    101797	 NSYN	 A	         T	         ESCCO14588_RS26480	  101814	     101644	     171	         1	             18//	
    
    //ref codon  ref aa      query codon	query aa	 transition/transversion	//		maxlen:abht	 blengths:abht	 product
      TTA        L           TTT	        F	         transversion	           // 	    40	         40	             hypothetical

The table not only provides the SNP, but also the corresponding annotation from the reference genome. It provides the length of the
blastn hit(maxlen), as well if there are additional hits of lower quality (blenghts).


### Post-processing (workflow B.3 optional)
The SNPs multifasta is converted to phylip to run on PhyML. The filtered table is processed to provide genotypes and summary
of SNP characteristics. Please provide the name of the reference genome in the genotyper summary tab.





