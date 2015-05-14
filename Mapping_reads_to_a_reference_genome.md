# Mapping Reads to a Reference Genome

Depending on your organism of study, there may or may not be a relatively closely related genome sequence to work with.  Depending on your research question, this may or may not be useful.  In our example study on Tonkean macaques, we are ultimately interested in quantifying molecular polymorphism on the X chromosome and comparing it to polymorphism on the autosomes.  For this reason, the genomic location of the data is crucial and we can benefit from the complete genome sequence of a closely related species of macaque monkey, the rhesus macaque.  We will use a program called [`bwa`] (http://sourceforge.net/projects/bio-bwa/files) and also [`samtools`](http://samtools.sourceforge.net/), to map our data to the rhesus macaque genome.

## Preparing your reference genome

Reference genomes for many sequences are available at multiple publicly available databases.  We can download the complete genome sequence for the rhesus macaque from the [USC genome browser](http://hgdownload.cse.ucsc.edu/downloads.html#rhesus).  I did this earlier because it takes a while.  It is a fasta-formatted file, and is located in this directory:

** insert path to rhesus macaque genome here **

Before we map our data to this reference genome, we need to generate some files that will be used in the mapping process.

./bwa index -a bwtsw **path_to_rhesus_genome**/**rhesus_genome_fasta_file**

./samtools faidx /work/ben/XL_unigene/Xl.seq.uniq
