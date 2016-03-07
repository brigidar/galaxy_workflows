# galaxy_workflows

'''Assembly and Annotation (workflow A)'''
The first workflow uses SPAdes v.3.7 for de novo assembly of reads and PROKKA for annotation. Please install the modified xml
files for SPADes and PROKKA to allow for renaming of files in the batch submission mode. Reads should be imported as fastqsanger format.
If you are assembling longer reads please change the k-mer size in the SPADES option.Make a tab delimited file to import the 
annotation options for PROKKA. You only need to provide a value of 1 for kingdom if you are annotating viruses.




'''Reads-based SNP Discovery (workflow B)''
This repository contains multiple workflows for Galaxy to discover and quality control SNPs based on reads (Illumina. 
There is a second option that provides SNP discovery based on contigs. 
