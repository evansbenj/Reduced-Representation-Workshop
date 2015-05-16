# Mapping Reads to a Reference Genome

(or you can go back to the de-multipexing page [here](https://github.com/evansbenj/Reduced-Representation-Workshop/blob/master/De-multiplexing.md)).

Depending on your organism of study, there may or may not be a relatively closely related genome sequence to work with.  Depending on your research question, this may or may not be useful.  In our example study on Tonkean macaques, we are ultimately interested in quantifying molecular polymorphism on the X chromosome and comparing it to polymorphism on the autosomes.  For this reason, the genomic location of the data is crucial and we can benefit from the complete genome sequence of a closely related species of macaque monkey, the rhesus macaque.  We will use a program called [`bwa`] (http://sourceforge.net/projects/bio-bwa/files) and also [`samtools`](http://samtools.sourceforge.net/), to map our data to individual chromosomes of the genome of a rhesus macaque (*Macaca mulatta*).  Normally one would map reads to an entire genome because the data were generated from a complete genome but in our case we are doing only an example analysis and we have to work within time constraints.

## Preparing your reference genome

Reference genomes for many sequences are available at multiple publicly available databases.  We can download the complete genome sequence for the rhesus macaque from the [USC genome browser](http://hgdownload.cse.ucsc.edu/downloads.html#rhesus).  I did this earlier because it takes a while.  It is a fasta-formatted file, and is located in this directory:

** insert path to rhesus macaque genome here **

Please go to this directory using this command:

`cd ** insert path to rhesus macaque genome here **`

Before we map our data to this reference genome, we need to generate some files that will be used in the mapping process.  This can be done in three steps:

1. ./bwa index -a bwtsw **path_to_rhesus_genome**/**rhesus_genome_fasta_file**

  The `./bwa` command tells the computer to execute the bwa program.  The `index` command tells `bwa` to generate index files from the rhesus genome file that is indivated by the `**path_to_rhesus_genome**/**rhesus_genome_fasta_file**`.  The `-a bwtsw` flag specifies the indexing algorithm for `bwa` to use.  This step will a few minutes.

2. We now need to to generate another file using `samtools`.  Please type this:

  ./samtools faidx **path_to_rhesus_genome**/**rhesus_genome_fasta_file**

  Here, the `./samtools` command tells the computer to execute the `samtools` program.  The `faidx` option tells samtools to generate a file called `**rhesus_genome_fasta_file**.fai` in which each line has information for one the contigs within the reference genome including the contig name, size, location and other information.  Our reference genome has a contig for each chromosome.

3.  The third thing we need to do is to generate a `.dict` file with a program called [`piccard`](http://broadinstitute.github.io/picard/).  Please type this command:

  `java -jar picard.jar CreateSequenceDictionary REFERENCE=**rhesus_genome_fasta_file** OUTPUT=**rhesus_genome_fasta_file**.dict`

  This should generate a file called "**rhesus_genome_fasta_file**.dict"

## Mapping the data to the reference genome

Now we can align the data from each individual to the reference genome using `bwa` as follows:

`./bwa aln reference_genome individual_1.fastq > individual_1.sai`

Add a header and generate `.sam` files

`./bwa samse -r "@RG\tID:FLOWCELL1.LANE6\tSM:Individual_1\tPL:illumina" Individual_1.fa Individual_1.sai Individual_1.fastq > Individual_1.sam`

Generate a `.bam` file:

`./samtools view -bt Individual_1.fa -o Individual_1.bam Individual_1.sam`

Sort the `.bam` file:

`./samtools sort Individual_1.bam Individual_1_sorted`

Make a `.bai` file:

`./samtools index Individual_1_sorted.bam`

OK, if this all went smoothly we are now ready to make some genotype calls.
