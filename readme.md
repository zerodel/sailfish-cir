## Introduction ##

Sailfish-cir is a computational tool to estimate the relative abundance of circular RNA transcripts from high-throughput RNA-seq data.

It accepts [CIRI](https://sourceforge.net/projects/ciri/) output or a BED-format file to specify the reference set of circular RNA transcripts in RNA-seq data. Then, it transforms all circular transcripts to pseudo-linear transcripts. Finally, it estimates the expression of both linear and circular transcripts using [Sailfish](http://www.cs.cmu.edu/~ckingsf/software/sailfish/) framework. 


## Prerequisites ##

The following three tools should be installed before running Sailfish-cir.
	
    1. Cufflinks 
	   Gffread, a part of Cufflinks Package (http://cole-trapnell-lab.github.io/cufflinks/), is used to extract sequence.
    2. Sailfish 
	   Sailfish (http://www.cs.cmu.edu/~ckingsf/software/sailfish/) is used to quantify RNA expression. Sailfish 0.7.0 or above is recommended.
    3. gffutils 
	   gffutils (https://pypi.python.org/pypi/gffutils) is used to build a database for reference set of all RNA transcripts.
	
We assume these pre-installed tools are running under a Unix-like environment, and paths of gffread and sailfish binary executive should be added to your ``PATH`` variable.

Other than above third party tools, several data files are required for Sailfish-cir as well.
	
	1. Raw sequencing reads file in fasta/fastq format;
	2. A CIRI output file or a BED file;
	3. Gene annotation file in .gtf format, and a single fasta file or multi-fasta files containing RNA sequences. 


## Usage ##

```
python path/to/your/sailfish_cir.py -g path/to/your/genomic/sequence/foo.fa -a path/to/your/annotation/bar.gtf -1 path/to/your/reads/mate1.fastq -2 /path/to/your/reads/mate2.fastq -o /path/to/where/you/want/your/result -c /path/to/your/CIRI/output/file
```
    

    -g  path to genomic sequence fasta file
    -a  path to gene annotation file, ie, .gtf or .gff files
    -r  path to single-end raw sequencing reads file.
    -1  path to raw pair-end reads, mate 1
    -2  path to raw pair-end reads, mate 2
    
    -c  path to CIRI output file to specify circular RNA
    --bed  path to bed file which contains circular RNA.
    --circRNA_finder path to circRNA_finder output file.
    --KNIFE_report_folder path to KNIFE output folder, make sure it has the subdirectory "circReads"
    
    -o  output folder that contains the index built by sailfish and quantification results
    -k  minimum match size used during sailfish quantification,   default is 21
	--libtype   format string describing the library type of your reads. default is "IU", [read more on libtype of Sailfish](http://sailfish.readthedocs.org/en/master/library_type.html)
    --mll mean library length, this option is to fix up the effective length.
    -h/--help	print this help message

    

    
## Output file  ##

Sailfish-cir expression estimats can be found at a subdirectory named ``quant_circular`` under the path set by ``-o`` parameter.
	

## Notes ##

1. .gtf file 

    In order to generate reference sequences of circular RNA transcripts, a database of your genomic annotation is needed. It is recommended to create the database file manually using the same filename with extension ".db" in the same directory. Since gffutils don't guess GTF file format, and GTF format was changed after GRCH37.75. It will be time-consuming when building the database using the same option for earlier gtf version. 

2. output file from circRNA identification tools

    This script accepts BED-format file or output files of circRNA identification tools (CIRI, circRNA_finder, KNIFE are supported) to specify the circular RNA transcripts. 
    
    We also provide a small utility to convert those outputs to a BED-format file.
usage as follows:
    ```
    python path/to/your/sailfish_cir.py convert_CIRI your/CIRI/output/file
    python path/to/your/sailfish_cir.py convert_KNIFE your/KNIFE/output/folder
    ```
   
    For CIRI, this will create a .bed file and a .mapping file (a tab-delimited file contains circular RNA name and its host gene in each line) in the same folder. 
    For KNIFE, this will create a "summarized_knife_junction.bed" under your KNIFE output folder.
    circRNA_finder output is actually BED-format.
    
    If you want to quantify circRNA from multiple circRNA detection results, convert those output files into BED format and merge them into a single .bed file, then use "--bed" option. 

## ChangeLog ##

Release v0.11

* BED-file format support: use '--bed' option to input a bed file as the reference set of circular RNA transcripts
* Convert utility: convert CIRI outputs into a BED-file format
* Issue of effective length correction: fix the effective length correction when estimating circular RNA expression values in real rRNA-depleted RNA-seq data. Please use '--mll' to specify the estimated effective fragment length of RNA-seq data. 
* Add some filters to unclutter the GTF entries. 

Release v0.10 

* the first version. Snapshot with documents is [here](https://github.com/zerodel/sailfish-cir/releases/tag/v0.10).


## Example files ##

Example files are available in [MEGA](https://mega.nz)
    
1. Human gene Annotation file:   [hg19 gtf file](https://mega.nz/#!coZEBY5D!-w5VbydDbNFW4peA2yK3gYjX0kb7mUBdMlBII6HOtpg)
    OR [binary database file by gffutils](https://mega.nz/#!Z1QFHBYb!2lYvqCDzNXh6X1othSvPwA0NQb1RlhtMoHAqveOxmSM)
2. Genome sequence file: [hg19 fasta file](https://mega.nz/#!40JiUDJK!9oC5PSleQSZjgIlFWUaRODYKh5nYxIW_Lfexwlk9QJc)
3. Sailfish (ver 0.9.0): [Sailfish Linux binary](https://mega.nz/#!hopk3IzA!7b39ya6xy9YlCYmnSDO9I6xXSEw8-PTlTiXxs7CE3UU)
4. RNA-seq simulation scripts [simulation](https://mega.nz/#!NxwniILD!Ysmy4ybcZaQUfx9pe2h6Rsysn5vZDodiVynkONJSgEs) 
5. Simulated RNA-seq datasets [11 simulated datasets (uniform expresion distribution)](https://mega.nz/#!FwhEgSoa!lE-vZ5Hv9Ib3UAEiNhoyUWvZfdgu5Md_OPMoYFDath8) , [one simulated dataset (empirical expression distribution)](https://mega.nz/#!AohTlQZT!69BBJSfze0cmDioRd9gBn0kdG125eivRyZBWMiw1buQ) , [uniform sequencing error with error rate 0.01](https://mega.nz/#!oh4VlTaJ!6Dt3_vENKRbWH2Jc2uIwe3ne0bYDVtrVp3DFThlU44g) , [uniform sequencing error with error rate 0.02](https://mega.nz/#!cwREQRLb!x5dm81OHhkx5CCVyiMeZcXtYE54dyNZHA1D5FU0sVWg), [empirical error model: illumina 4](https://mega.nz/#!kwAw3DJY!jtVKO78aLLKzNne-5FRhTj7fX_sToVflWIgYtaXKlz4), [empirical error model: illumina 5](https://mega.nz/#!pphAECAY!yHMgVI9h3B5pcYd1Aaj3AQLgkFjeWjHdQX6CbA-cdmo), [simulation size two fold](https://mega.nz/#!gkZgHLza!wXOL58ioIODUhHq-sq-5-RXyY-YXzWOsGE61HHP9K4Y), [simulation size 4 fold](https://mega.nz/#!ApIB1DZJ!lMsh7SL6MLKyWph0cAGzhkazwBtNLMQTt2QQwkoNUls), [simulation using whole genome](https://mega.nz/#!d4RCTTII!hhGW_2dky-7sZW4Z8t5330oiKNVeV6D7M8BFEUDaZe8)


## Reference ##
Musheng Li, Xueying Xie, Jing Zhou, Mengyin Sheng, Xiaofeng Yin, Eun-A Ko, Tong Zhou* and Wanjun Gu*(2016) Sailfish-cir: a model-based tool to quantify circular RNA expression from high-throughput RNA-seq data. Submitted.


## Contact ##
Musheng Li (zerodel@126.com)
 
 
 

 
