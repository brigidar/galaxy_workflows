# SNPDV

## General Considerations before starting

This pipeline is thought for outbreak inclusion and exclusion analyses based on single nucleotide polymorphism, as well as phylogenetic inference of closely related bacterial species or serotypes. For more distantly related bacteria the amount of homoplastic events will preclude an accurate reconstruction of the evolutionary relationship of the species of interest and other methods should be considered. We also recommend to have some prior knowledge of possible mobile genetic elements, such as resistance cassettes, phages and insertion elements, of the species of interest, as they will interfere with the SNP analysis if not correctly excluded (see exclusion criteria).

## Files required
 * Reference Genome (fasta and annotated genbank, see A workflows to generate a reference from in-house sequencing or SRA files)
 * IS elements of reference genome (see exclusion criteria)
 * Plasmids database of species of interest (only required if reference is a draft genome or assembled from in-house sequencing)
 * Reads or assemblies from genomes of interest

## Galaxy Setup

### System Requirements
Linux or Unix based only. Some of the tools used in these workflows were developed by others and do not come in binaries for Windows.

Before installing Galaxy please check your python version (2.7 or higher). You will also need the following packages: Numpy, Pandas, Matplotlib, and Biopython. We recommend installing [anaconda python] (https://docs.continuum.io/anaconda/install), which comes with all the packages except biopython that needs to be installed with the command "conda install biopython". 

Galaxy is an open, web-based platform for accessible, reproducible, and transparent computational biomedical research. 
Galaxy can be found on [Github] (https://github.com/galaxyproject) and the manual in the [Galaxy Wiki] (https://docs.galaxyproject.org/en/latest/index.html).
Galaxy can be run locally, on a server (cluster or single server), and on the cloud. If you are planning on running Galaxy on the cloud, we recommend you select a linux-based operating system.

All tools and workflows required can be retrieved from the Galaxy toolshed under the repository SNPDV. This can be done from the browser interface directly. For more information on how to install tools, please refer to the [Galaxy Toolshed Wiki] (https://docs.galaxyproject.org/en/latest/ts_api_doc.html).

Once installed you can upload your files through the browser or retrieve reads from NCBI with 

## Workflows

### Assembly and Annotation (workflows A)
**A.1 Retrieve SRR reads**

The workflow fetches SRR reads based on Bioproject number (NCBI) and saves it as a paired list.

**A.2 Assembly**

The workflow uses SPAdes v.3.7.1 for *de novo* assembly of reads. 
  * Reads should be imported as fastqsanger format and grouped into a paired dataset list.
  * If you are assembling longer reads please change the k-mer size in the SPADES option
     * Reads 250 bp or more use k-mer 21,33,55,77,99,127
     ![kmer option] (https://github.com/brigidar/galaxy_workflows/blob/master/kmer.png)

  * Scaffolds are filtered according to size and coverage
     * Coverage or length settings should be changed according to sequencing technology (MiSeq, NextSeq, HiSeq)  
     ![coverage option] (https://github.com/brigidar/galaxy_workflows/blob/master/coverage.png)

  * For more information on SPAdes please refer to the [manual] (http://spades.bioinf.spbau.ru/release3.7.0/manual.html)
  * An assembled genome is required as a reference for read-based SNP calling.
   
**A.3 Annotation**

Annotation provided by PROKKA
  * Make a tab delimited file to import the annotation options for PROKKA.
     * Use 1 in kingdom column if you are annotating viruses, otherwise leave blank.
  * For more information on PROKKA please refer to the [github repository] (https://github.com/tseemann/prokka)
  * We recommend building a genus database of curated annotations to install in your Galaxy instance. This will improve the annotation by preferentially using the genus database during annotation.
  * Annotated genbank file is required for reference sequence only in SNPDV.

###### Example annotation file  
    strain	 locustag	 centre	 genus	      species	plasmid	kingdom
    strainA	 xxx	     centre	 Escherichia	coli		
    strainB	 yyy	     centre  Escherichia	coli	   p90
  

### Excluding regions (workflows B)

#### Excluding Mobilome Regions in SNPDV

#####Excluding with a closed reference:
Use the workflow **"Excluded Regions (closed reference)"**.
If you have a closed reference with known mobilome and repeated regions this step can be skipped. Provide a interval formatted file with the coordinates of the regions to exclude. If you don't know the mobilome use [PHASTER] (http://phaster.ca/) for phage regions prediction and save regions as an interval formatted file. Use [ISFinder] (https://www-is.biotoul.fr/) to predict IS elements for your species of interest. This needs to be done only once per species. If your species carries any other mobile genetic elements within the chromosome please add the regions to the interval file (e.g. resistance cassettes). 

#####Excluding with a draft genome reference (no reads):
Use the workflow **"Excluded Regions (draft genome)"**.
If you have a draft genome reference without reads you can predict the mobilome with [PHAST] (http://phast.wishartlab.com/), but you will need to save the predicted regions as multifasta, since the PHAST intervals correspond to the concatenated draft genome and not the intervals on each contig. For IS elements you can use  [ISFinder] (https://www-is.biotoul.fr/). If your species carries any other mobile genetic elements please add the regions to the interval file (e.g. resistance cassettes). Since it is a draft genome, any plasmids present will be mixed in with the chromosomal contigs, therefore provide a multifasta to exclude them. 

#####Excluding with a draft genome reference (reads):
Use the workflow **"Excluded Regions (reads)"**.
If you have a draft genome reference with reads you can predict the mobilome with [PHAST] (http://phast.wishartlab.com/), but you will need to save the predicted regions as multifasta, since the PHAST intervals correspond to the concatenated draft genome and not the intervals on each contig. For IS elements you can use [ISFinder] (https://www-is.biotoul.fr/). If your species carries any other mobile genetic elements please add the regions to the bed file (e.g. resistance cassettes). In addition you will add the Freebayes output for the reference reads against itself to exclude any base miscalling of the assembly. Since it is a draft genome, any plasmids present will be mixed in with the contigs. therefore provide a multifasta to exclude them. 

####### Example interval file
    NC_011353.1	273971	274038
    NC_011353.1	275213	276349
    NC_011353.1	302573	314525
    
### SNPDV (workflows C)

#### Reads-based Discovery (workflow C.1.1)

For the selected reference you will need a fasta and an annotated genbank file. For the query genomes you will need reads (fastqsanger). The reads should be provided as a paired list if not downloaded from SRR. SNPs will be predicted by Bowtie2 mapping & Freebayes variant calling. If you are using MiSeq reads lower the coverage to 10 in the Freebayes options.
![coverage] (https://github.com/brigidar/galaxy_workflows/blob/master/coverage_fb.png)

#### Count predicted SNPs (optional)
This script can help to predict an outlier before validation. It will count the amount of variants predicted for each genome in the output list of workflow B.1.1. Outlier genomes can be removed from the analysis before proceeding to validation. Keeping outlier genomes can cause the loss of valid SNP locations due to lack of matching sequences (no hits).

#### SNP Validation reads (workflow C.1.2)
A reference genome (fasta & genbank) is required as well as the excluded regions and predicted SNP list. SNPs are compared to excluded regions for filtering. Remaining SNPs are combined into a single table and no hits and positions with ambigous nucleotides are removed (optional). The output is a filtered table and a multifasta with the curated SNPs for each query genome.

#### Contig-based Discovery (workflow C.2.1)
If no reads are available for query genomes, SNP discovery is based on NUCmer with delta-filter and show-snps. Provide a list of query genomes (fasta files) to predict SNPs against reference genome (fasta file).

#### SNP Validation contigs (workflow C.2.2)
You will need a fasta and genbank file of the reference, the excluded regions (prepared according to reference genome format), and the predicted SNPs for the query genomes. You will need the assemblies of the query genomes as separate fasta files and multifasta. The multifasta file can be generated in Galaxy with the concatenate fasta tool. If you have both read and contig-based SNP discovery outputs you can combine them at this step to get all the SNPs in one single merged table.
If your computational settings allow it the verification step can be threaded. Number of threads will dependon on number of SNPs (nb of threads = 30% of predicted SNPs or max amount of available threads). To change to a threaded version open the workflow and modify the SNP verify block by selecting threaded.  
![change to threaded] (https://github.com/brigidar/galaxy_workflows/blob/master/threaded.png)

The output is a filtered table and a multifasta with the curated SNPs for each query genome. 



##### Example filtered SNP table
    molecule	    refpos	 syn?	 refbase	 qbase:abht	 gene name	          gene start	 gene end	 gene length	 snps per gene	 pos in gene/	
    NZ CM000662	    101797	 NSYN	 A	         T	         ESCCO14588_RS26480	  101814	     101644	     171	         1	             18//	
    
    //ref codon  ref aa      query codon	query aa	 transition/transversion	//		maxlen:abht	 blengths:abht	 product
      TTA        L           TTT	        F	         transversion	           // 	    40	         40	             hypothetical

The table not only provides the SNP, but also the corresponding annotation from the reference genome. For contig based discovery it provides the length of the blastn hit(maxlen), as well if there are additional hits of lower quality (blenghts). Positions with multiple high quality blastn hits are excluded, as well as indels, no hits, and invariant sites.


#### Genotyping (workflow C.3 optional)
The SNPs multifasta is converted to phylip to run on PhyML (GTR, gamma, 1000 bootstraps). If you have included the reference genome in the query genomes during the SNP verify step (last option),  provide the identifier in the drop fasta option.
![drop fasta] (https://github.com/brigidar/galaxy_workflows/blob/master/drop_fasta.png)

The filtered table is processed to provide genotypes and summary of SNP characteristics.
![summary genotypes] (https://github.com/brigidar/galaxy_workflows/blob/master/genotyper.png)

If you have included the reference genome in the query genomes during the SNP verify step (last option), provide the column numbers to the left and right of the column to be removed before genotyping. 
![drop column] (https://github.com/brigidar/galaxy_workflows/blob/master/drop_column.png)

Provide the name of the reference genome in the genotyper summary tab. 
![reference genome name] (https://github.com/brigidar/galaxy_workflows/blob/master/genotyper_select.png)

To plot SNPs along the genome povide the length of the genome and the bin size for the SNPs. This option is currently only available for closed genomes.
![snp location] (https://github.com/brigidar/galaxy_workflows/blob/master/figure_1.png)





