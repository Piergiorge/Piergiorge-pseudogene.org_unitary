This archive consists of a standalone version of pseudopipe pipeline and a simple try-out data. 
Since the pipeline was originally designed and automated to work with ensembl data, so some manual settings are required to run it with other input data.

===
There are three folders within a parent directory "pgenes" after extraction of the archive:
===
- pseudopipe: which contains pipeline code;
- ppipe_input: which contains input data;
- ppipe_output: which contains output data.


===
Input data:
===
You may create a separate folder within the ppipe_input (and ppipe_output) for each species. There need to be three folders for each species genomic input data,
- dna: contains a file named dna_rm.fa, which is  entire repeat masked dna from that species, and a list file for all unmasked dna divided into different chromosomes in FASTA format;
- pep: contains a FASTA file for all the proteins in the species;
- mysql: contains a list of files named as "chr1_exLocs", "chr2_exLocs", etc. to specify exons coordinates, one for each chromosome. Only thing matters for these files are their third and fourth columns, which should be start and end coordinates of exons.


===
Environment setting:
===
You'll need python, blast and tfasty to run the pipeline. Their paths should be indicated at the end of  /pseudopipe/bin/env.sh


===
Run the pipeline
===
First go to the folder pseudopipe/bin, and run with  command line in the form of: ./pseudopipe.sh [output dir] [masked dna dir] [input dna dir] [input pep dir] [exon dir] 0.

An example using the try-out data is as follow:
./pseudopipe.sh ~/pgenes/ppipe_output/caenorhabditis_elegans_62_220a ~/pgenes/ppipe_input/caenorhabditis_elegans_62_220a/dna/dna_rm.fa ~/pgenes/ppipe_input/caenorhabditis_elegans_62_220a/dna/Caenorhabditis_elegans.WS220.62.dna.chromosome.%s.fa ~/pgenes/ppipe_input/caenorhabditis_elegans_62_220a/pep/Caenorhabditis_elegans.WS220.62.pep.fa ~/pgenes/ppipe_input/caenorhabditis_elegans_62_220a/mysql/chr%s_exLocs 0
(This command line assumes you extract the archive in your home directory, i.e., "~/". Please note that the paths in the command line need to be absolute, and chromosome and exon files are specified with wild card "%s".)

The blast step is already included in the pipeline.


===
Output:
===
The output can be found at ppipe_output/caenorhabditis_elegans_62_220a/pgenes/ppipe_output_pgenes.txt, given the above command line.


===
Run time
===
On a single laptop (2.6GHz, 4GB RAM): The most time consuming step is tblastn. It may take around one day to finish an entire genome in a comparable size of C. elegans. The  following steps will finish in a few hours.

We've implemented the pipeline to run parallel in cluster machines. However, the pipeline I sent can only run on a single machine. The parallel implementation is currently hard-coded to our local settings.
