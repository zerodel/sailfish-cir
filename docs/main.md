## Introduction ##

Sailfish-cir is a software tool for estimating the relative abundance of Circular RNA transcript from high-throughput RNA-seq data.

It accepts [CIRI](https://sourceforge.net/projects/ciri/) output and BED format file to specify reference set of circular RNA transcripts in RNA-seq data, transforms circular transcripts to pseudo-linear transcripts, and then estimates the expression of RNA transcripts, both linear and circular, in [Sailfish](http://www.cs.cmu.edu/~ckingsf/software/sailfish/) framework. 


## Prerequisites ##

The following three tools should be installed before running this script.
	
    1. Cufflinks
	   Gffread,a part of Cufflinks Package, is used as sequence extractor here.
    2. Sailfish
	   Sailfish is used to quantify RNA expression here. Sailfish 0.7.0 or above is recommended.
    3. gffutils
	   This Python package is needed to build a database for your genomic annotation file.
	
We assume that these pre-installed scripts are running under a Unix-like environment, and paths of gffread and sailfish binary executive should be added to your ``PATH`` variable.

Other than these scripts, you should provide the following data files for Sailfish-cir as well.
	
	1. Your reads file in fasta/fastq format;
	2. A CIRI output file. Or a BED file 
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

    -k  k-mer size used by sailfish to built index. default is 21
	--libtype   format string describing the library type of your reads. default is "IU", [read more on libtype of Sailfish](http://sailfish.readthedocs.org/en/master/library_type.html)
    --bed  path to bed file which contains circular RNA.
    --mll mean library length, this option is to fix up the effective length.
    

    
## Output file  ##

Sailfish-cir output files can be found at a subdirectory named ``quant_circular`` under the path set by ``-o`` parameter.
	

## Notes ##

1. .gtf file 

    In order to generate reference sequence of circular RNA transcripts, a database of your genomic annotation is needed. It is recommended to create the database file manually using the same filename with extension ".db" in the same directory.
Since gffutils don't guess GTF file format ,and GTF format changes after GRCH37.75. It will be time consuming to build the database using the same option for earlier gtf version. 

2. CIRI file , or BED file

    This script now accept CIRI outputs (CIRI v1.2) or BED file to specify the circular RNA locus. It also provides a small utility to convert CIRI output to BED file.
usage as follows:
    ```
    python path/to/your/sailfish_cir.py convert your/CIRI/output/file
    ```
    this will produce a .bed file and a .mapping file (a tab-delimited file contains circular RNA names and its host gene in each line) in the same folder. 


## Reference ##
Musheng Li, Xueying Xie, Jing Zhou, Mengyin Sheng, Eun-A Ko, Tong Zhou* and Wanjun Gu*(2016) Sailfish-cir: a model-based tool to quantify circular RNA expression from high-throughput RNA-seq data. Submitted.


## Contact ##
Musheng li (zerodel@126.com)
 
 
 
## ChangeLog ##

Release v0.11

* BED file format support : now '--bed' can load a bed file as input
* Convert utility: convert CIRI output file into .bed format
* Effective length issue: to fix up the effective length gap between circular RNA and linear RNA, use '--mll' to specify the mean length of fragment. 
* Add some filter to unclutter the GTF entries. 

Release v0.1 

* the first version.  snapshot with documents is [here](https://github.com/zerodel/sailfish-cir/releases/tag/v0.1).


