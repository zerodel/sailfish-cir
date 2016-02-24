## Introduction ##

Saifish-cir is a software tool for estimating the relative abundance of Circular RNA transcript from high-throughput RNA-seq data.

It accepts [CIRI](https://sourceforge.net/projects/ciri/) output to specify reference set of circular RNA transcripts in RNA-seq data, transforms circular transcripts to pseudo-linear transcripts, and then estimates the expression of RNA transcripts, both linear and circular, in [Sailfish](http://www.cs.cmu.edu/~ckingsf/software/sailfish/) framework. 


## Prerequisites ##

The following three tools should be installed before running this script.
	
    1. Cufflinks
	   Gffread,a part of Cufflinks Package, is used as sequnce extractor here.
    2. Sailfish
	   Sailfish is used to quantify RNA expression here. Sailfish 0.7.0 or above is recommended.
    3. gffutils
	   This Python package is needed to build a database for your genomic annotation file.
	
We assume that these pre-installed scripts are running under a Unix-like enviroment, and paths of gffread and sailfish should be added to your ``PATH`` variable.

Other than these scripts, you should provide the following data files for Saifish-cir run as well.
	
	1. Your reads file in fasta/fastq format;
	2. A CIRI output file.
	3. Genomic annotation file in .gtf format, and a single or multi-fasta file containing RNA sequences. 


## Usage ##

    ```
	python path/to/your/sailfish_cir.py -g path/to/your/genomic/sequence/foo.fa -a path/to/your/annotation/bar.gtf -1 path/to/your/reads/mate1.fastq -2 /path/to/your/reads/mate2.fastq -o /path/to/where/you/want/your/result -c /path/to/your/CIRI/output/file
    ```
    
    -g  path to genomic sequence fasta file
    -a  path to genomic annotation file, ie, .gtf or .gff files
    -r  path to your single-end sequencing reads file.
    -1  path to your pair-end reads, mate 1
    -2  path to your pair-end reads, mate 2
    -o  output folder that contains your quantification results and the index built by sailfish
    -c  path to CIRI output file to specify circular RNA
    -h/--help	print this help message
	
optional arguments:
    
    -k  k-mer size used by sailfish to built index. default is 19
	--libtype   format string describing the library type of your reads. default is "IU", [read more on libtype of Sailfish](http://sailfish.readthedocs.org/en/master/library_type.html)


## Output file  ##

Sailfish-cir output files can be found at a subdirectory named ``quant_circular`` under the path set by ``-o`` parameter.
	

## Notes ##

1. gtf file 
In order to generate reference sequence of circular RNA transcripts, a database of your genomic annotation will be build. It is recommended to create the database file manually with the suffix "db" and put it under the same directory of gtf file. 
Because the format of gtf file changes after GRCH37.75. It will be very slow to build the database using the same option for earlier gtf version. 

2. CIRI file 
This script only supports CIRI outputs (CIRI v1.2) to specify the circular RNA isoforms for current version {2016/02/04}. File formats may change as CIRI updates, and other file formats will be supported soon.


## Reference ##
Musheng Li, Xueying Xie, Jing Zhou, Tong Zhou* and Wanjun Gu*(2016) Sailfish-cir: a model-based tool to quantify circular RNA expression from high-throughput RNA-seq data. Submitted.


## Contact ##
Musheng li (zerodel@126.com) 




## File DownLoad ##

### Script files ###

You can download the scripts by [click here](https://github.com/zerodel/sailfish-cir/archive/master.zip). 
If you are familiar with git, use ``git@github.com:zerodel/sailfish-cir.git`` .


### Example files ###

1. Genomic Annotation file   [hg19 .gtf file](https://mega.nz/#!spA1BYZS!ab7EEWilWhUsvp6LeAPic1ia32dkO049sN17OB3foww)
    [binary gffutils database file](https://mega.nz/#!5pYmXKwZ!Fxgr5nc2LncyTojDXR1jTxbBs4RyDmDBgglg55udCbM)
2. Genome file [Genome file](https://mega.nz/#!QpBkXArT!HCyijZK6av5MRCwFnPaf7OS0eHC8sRa3szTP5Tt_Qas)
3. RNA-seq data [mate1](https://mega.nz/#!t9o2SLIJ!hlOsdJ2RC6XfvAHY1o6GD-KZ6PzPwGwjGmM7ORYXlSU) [mate2](https://mega.nz/#!Ulx2AQZa!G0Cu01LG3bE6z8WaiSP4gw5ohBEOgzfXtTJDulN0Az8)
4. CIRI output file [CIRI output of RNA-seq data above](https://mega.nz/#!JspGGBzA!316V2Y8OsBFedMctE2LU0RtGAaZdvCWfsF583gocsAI)