# BLASTnet
This python project implements a fully automated pipeline that allows to perform a BLAST ALL search, gives as output a csv file with the results of the BLASTALL run and generates network graphs that can be further analyzed with the software package CLANS, GEPHI and CYTOSCAPE.
The network graphs are an intermediate step to a network graph analysis that allows to take the similarities between sequences and cluster them in groups.

## Background:
This script is for bioinformaticians that want to visualize their sequence data by similarities.
BLAST is known as "basic local alignment search tool" it allows to search a database of sequences with a query and outputs all the sequences that are similar to it in terms of nucleotide/aminoacid sequence.

Graphs are mathematical structures used to study pairwise relationships between objects and entities.
Graphs provide a better way of dealing with abstract concepts like relationships and interactions. They also offer an intuitively visual way of thinking about these concepts (extracted from. https://www.analyticsvidhya.com/blog/2018/04/introduction-to-graph-theory-network-analysis-python-codes/ ).

## Conceptual Pipeline:
1. mySequences (collection of my sequences in FASTA format)
2. generation of DB of mySequences = mySequencesDB
3. Run BLASTP, BLASTN, PARNASSUS BLAST or BLASTR(will be removed soon) on mySequencesDB with query mySequences.
4. Convert the output file to a CLANS, GEPHI and CYTOSCAPE readable format.
5. Proceed with the analysis through the software package.


## Installation steps:
**For running this script:**
* Clone this repo: `git@github.com:akaped/blastnet.git`
* Install the right BLAST+ package for your system
  ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/
* Be sure to have Python3.6 or latest and pip installed
* install all the python dependencies python complains about
**For the graphical analysis:**
* CLANS software (https://www.eb.tuebingen.mpg.de/protein-evolution/software/clans/)
* GEPHI (https://gephi.org/)
* CYTOSCAPE (https://cytoscape.org/)

## How to run this tool
The script I’ve made for you (blastnet.py) comprises the full pipeline from *1* to *4*.


`python blastnet.py mySequences.fasta -info`

If the database is present it will use that, otherwise it will generate it.
If the -p (protein) flag is used the script will perform a BLASTP analysis.


`python blastnet.py mySequences.fasta -p`

If the -n (nucleotide) flag is used the script will perform a BLASTN analysis.

`python blastnet.py mySequences.fasta -n`

If -eval (evalue) flag is set the user can specify a value for the blast analysis, otherwise it will use the default value of: 10 .

`python blastnet.py mySequences.fasta -n -eval “1E-10”`

If the -parnassus flag is set then the script will use PARNASSUS blastP binary.
PARNASSUS is a tool we are developing in IIMCB Warszawa.
If you don't know what it is, it's ok :)
This package doesn't contain all the necessary binaries to use this function, released in a private repository.

`python blastnet.py mySequences.fasta -parnassus`


The script will process the data, create a folder with the same name of mySequence.fasta, and place there all the run files + results. The results will be the output of  .gephi, .cytoscape and .clans file.

This is the output structure of the generated folder:
```
mySequences
|  
|__ runfiles
|   |__mySequences.fasta
|   |__mySequences.db*
|
|__results
|   |__ mySequences.gephi
|   |__ mySequences.clans
|   |__ mySequences.cytoscape
|   |__ mySequences.tsv
|
|__report.txt
```
