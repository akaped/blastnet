# BLASTnet
This python project implements a pipeline that allows to perform a BLAST ALL search, gives as output a cst file with the results of the BLASTALL run and generats network graphs that can be further analised with the software package CLANS or GEPHI. 
The network graphs are in intermediate step to a network graph analysis that allows to take be similarities between sequences and cluster them by their similarities in groups. 


## Background: 
This script is for bioinformaticians that want to visualise their sequence data by similarities.
BLAST is known as "basic local alignment search tool" it allows to search a database of sequences with a query and outputs all the sequences that are similar to it in terms of nucleotide/aminoacid sequence.


## Conceptual Pipeline:
1. mySequences (collection of my sequences in FASTA format)
2. generation of DB of mySequences = mySequencesDB
3. Run BLASTP or BLASTN on mySequencesDB with query mySequences 
4. Convert the output file to a CLANS or GEPHI readable format.
5. Proceed with the analysis through the software package. 


## Installation steps:
- For running this script: 
* Clone this repo: `git@github.com:akaped/blastnet.git`
* Install the right BLAST+ package for your system
  ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/
* Be sure to have Python3 and pip installed
* install all the python dependencies python complains about
- For the graphical analysis:
* CLANS software (https://www.eb.tuebingen.mpg.de/protein-evolution/software/clans/)
* GEPHI (https://gephi.org/) 

## How to run this tool
The script I’ve made for you (blastnet.py) comprises the full pipeline from *1* to *4*. 

	
`python blastnet.py mySequences.fasta -info`

If the database is present it will use that, otherwise it will generate it. 
If the -p (protein) flag is used the script will perform a BLASTP analysis.

	
`python blastnet.py mySequences.fasta -p`

If the -n (nucleotide) flag is used the script will perform a BLASTN analysis.

`python blastnet.py mySequences.fasta -n`

If -eval (evalue) flag is set the user can specify a value for the blast analysis, otherwise it will use the default value of: --- . 

`python blastnet.py mySequences.fasta -n -eval “1E-10”` 


The script will process the data, create a folder with the same name of mySequence.fasta, and place there all the run files + results. The results will be the output of  .gephi and .clans file. 

This is the output structure of the generated folder:
`
mySequence 
|     |     | 
|     |     |_ runfiles
|     |         |     |_mySequences.fasta 
|     |         |____mySequences.db*
|     |____results
|_______report.txt
`
