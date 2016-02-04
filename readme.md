# Sailfish-cir #

This Saifish-cir is a tool for estimating the relative abundance of Circular RNA transcript. 
It use output of CIRI to specify the list of circular RNA candidates. 
this list will collaborate with genomic sequence file and gene annotation to produce a fasta file containing reference transcipts. 
In the end , Sailfish will estimate the expression of your transcripts , linear and circular, in the help of index built on this fasta file . 

## Things before you use: ##

1. first, we assume that this script is running under unix-like enviroments. 
2. Since Gffread, a part of cufflink package , is used as sequence extractor and Sailfish used as a estimator of RNA expression.  make sure they are properly installed. here , Sailfish version should be 0.7.0 or above. 
3. add path of  binary executive file of gffread and sailfish to your PATH variable. 
4. gffutils , a Python package, is need to build a database for your annotation file.


## Files needed to run this script : ##

1. In order to get reference sequence of both circRNA and linear RNA transcript, fasta/fastq files that contain genomic sequence and .gtf format annotation file are both needed. 
2. a CIRI output will be needed to specify the circular RNA sequence. for more information , refer to https://sourceforge.net/projects/ciri/
3. your reads file in fasta/fastq format. 


## A typical run example: ##

python path/to/your/sailfish_cir.py -g path/to/your/genomic/sequence/foo.fa -a path/to/your/annotation/bar.gtf -1 path/to/your/reads/mate1.fastq -2 /path/to/your/reads/mate2.fastq -o /path/to/where/you/want/your/result -c /path/to/your/CIRI/output/file


    -g  path to genomic sequence fasta file
    -a  path to genomic annotation file, ie, .gtf or .gff files
    -r  path to your single-end sequencing reads file.
    -1  path to your pair-end reads , mate 1
    -2  path to your pair-end reads, mate 2
    -o  output folder that contains your quantification results ,also include index built by sailfish
    -c  path to CIRI output file to specify circular RNA
    -h/--help   print this help message
	
some optional arguments:
    
	-k  Kmer size used by sailfish to built index. default is 19
	 --libtype   format string describing the library type of your reads. default is "IU" 



 python here means Python 2.x, (we test this under Python2.7.x)
 
-k parameter , Since Sailfish  v0.7.0 , this parameter means the minimum match size during the quantification phase of sailfish. and this parameter must be an odd number. 

--libtype parameter defines the characteristic of your reads. for more information , refer to http://sailfish.readthedocs.org/en/master/library_type.html

output file will be found at a subdirectory named "quant_circular" under the path of  your -o parameter .

## Known issues: ##

1. gtf file 
   In order to generate reference sequence of circular RNA transcripts, a data base of your genomic annotation will be build. It is recommended to  create the database file manually with the suffix "db" and put it under the same directory of gtf file. 
Because the format of  gtf file changes after grch37.75.  It will be very slow to build the database using the same option for earlier gtf version. 

2. CIRI file 
   This script only support CIRI input to specify the circular RNA isoform information now {2016/02/04},  using CIRI v1.2. File format may change as CIRI updates. other file format will be supported soon.


   
   
   


