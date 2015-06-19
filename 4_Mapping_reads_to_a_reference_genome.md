# Mapping Reads to a Reference Genome

(or you can go back to the de-multipexing page [here](https://github.com/evansbenj/Reduced-Representation-Workshop/blob/master/3_De-multiplexing.md)).

Depending on your organism of study, there may or may not be a relatively closely related genome sequence to work with.  Depending on your research question, this may or may not be useful.  In our example study on Tonkean macaques, we are ultimately interested in quantifying molecular polymorphism on the X chromosome and comparing it to polymorphism on the autosomes.  For this reason, the genomic location of the data is crucial and we can benefit from the complete genome sequence of a closely related species of macaque monkey, the rhesus macaque.  We will use a program called [`bwa`] (http://sourceforge.net/projects/bio-bwa/files) and also [`samtools`](http://samtools.sourceforge.net/), to map our data to individual chromosomes of the genome of a rhesus macaque (*Macaca mulatta*).  Normally one would map reads to an entire genome because the data were generated from a complete genome, but in our case we are doing only an example analysis and we have to work within time constraints.  Ben will assign each of you a chromosome to work on.

## Preparing your reference genome

Reference genomes for many sequences are available at multiple publicly available databases.  We can download the complete genome sequence for the rhesus macaque from the [USC genome browser](http://hgdownload.cse.ucsc.edu/downloads.html#rhesus).  I did this earlier because it takes a while.  The whole genome comes as a fasta-formatted file, and I split it up into individual fasta files corresponding with each of the chromosomes.  These are located in this directory:

`/home/datasets/2015_Ben_Evans/rhesus_chromosomes/`

Please go to this directory using this command:

`cd /home/datasets/2015_Ben_Evans/rhesus_chromosomes/`

Now check out what is in this directory by typing this:

`ls`

Before we map our data to this reference genome, we need to generate some files that will be used in the mapping process.  This can be done in three steps:

1. `/apps/bwa/0.7.12/bwa index -a bwtsw /home/datasets/2015_Ben_Evans/rhesus_chromosomes/chrXXX.fa`

  The `/apps/bwa/0.7.12/bwa` command tells the computer to execute the bwa program.  The `index` command tells `bwa` to generate index files from the rhesus genome file that is indicated by the `/home/datasets/2015_Ben_Evans/rhesus_chromosomes/chrXXX.fa`. The `-a bwtsw` flag specifies the indexing algorithm for `bwa` to use.  You will need to change the `chrXXX.fa` to match whatever chromosome Ben tells you to work on.  For example, if you are working on chromosome 9, you should type this:

  `/apps/bwa/0.7.12/bwa index -a bwtsw /home/datasets/2015_Ben_Evans/rhesus_chromosomes/chr9.fa`  
  
  This step will take a few minutes.

2. We now need to to generate another file using `samtools`.  Please type this:

  `/apps/samtools/0.1.19/samtools faidx /home/datasets/2015_Ben_Evans/rhesus_chromosomes/chrXXX.fa`

  Here, the `/apps/samtools/0.1.19/samtools` command tells the computer to execute the `samtools` program.  The `faidx` option tells samtools to generate a file called `chrXXX.fai` in which each line has information for one the contigs within the reference genome, including the contig name, size, location and other information.  Our reference genome has a contig for each chromosome.

3.  The third thing we need to do is to generate a `.dict` file with a program called [`picard`](http://broadinstitute.github.io/picard/).  Please type this command:

  `java -jar /apps/picard-tools/1.131/picard.jar CreateSequenceDictionary REFERENCE=/home/datasets/2015_Ben_Evans/rhesus_chromosomes/chrXXX.fa OUTPUT=/home/datasets/2015_Ben_Evans/rhesus_chromosomes/chrXXX.dict`

  As before, you will need to change the `chrXXX` in this command to match the chromosome you are working with.  This should generate a file called `/home/datasets/2015_Ben_Evans/rhesus_chromosomes/chrXXX.dict`

## Mapping the data to the reference genome

Now we can align the data from each individual to the reference genome using `bwa`.  First let's go back to the directory that has our de-multiplexed data in it like this:

`cd ~/monkey`

You can see the demultiplexed fastq files by typing the `ls` command.

Now let's map data from one individual to the reference genome as follows:

`/apps/bwa/0.7.12/bwa aln reference_genome individual_1.fastq > individual_1.sai`

For example, for the first individual (PF515) we could type this

`/apps/bwa/0.7.12/bwa aln /home/datasets/2015_Ben_Evans/rhesus_chromosomes/chrXXX.fa PF515.fastq > PF515.sai`

(but with the `chrXXX.fa` changed to match the chromosome you are working on.)

This command generates an intermediate file with the `.sai` suffix (which stands for `suffix array index`). Now we need to generate a `.sam` formatted file from our `.sai` files.  A `.sam` file is a tab delimited text file that contains the alignment data.  The format for this command is:

`/apps/bwa/0.7.12/bwa samse reference_genome.fa individual_1.sai individual_1.fastq > individual_1.sam`

We also need to add a header to each `.sam` file, so we can type this command:

`/apps/bwa/0.7.12/bwa samse -r "@RG\tID:FLOWCELL1.LANE6\tSM:Individual_1\tPL:illumina" reference.fa Individual_1.sai Individual_1.fastq > Individual_1.sam`

Now we can generate a `.bam` file.  A `.bam` formatted file is a binary version of the `.sam` file.

`/apps/samtools/0.1.19/samtools view -bt reference_genome -o Individual_1.bam Individual_1.sam`

Sort the `.bam` file:

`/apps/samtools/0.1.19/samtools sort Individual_1.bam Individual_1_sorted`

Make a `.bai` file:

`/apps/samtools/0.1.19/samtools index Individual_1_sorted.bam`

## Problem 3: Assessing coverage

Samtools can provide information on the number of reads for each position of the reference sequence for which there are data.  You can see this information by typing this:

`/apps/samtools/0.1.19/samtools depth XXX_sorted.bam`

Where `XXX` is the sampleid number.  If you want to know the average depth across all sites, you could type this:

`/apps/samtools/0.1.19/samtools depth XXX_sorted.bam | awk '{sum+=$3} END { print "Average = ",sum/NR}'`

Here the vertical bar `|` is a "pipe" that sends the information from the command before it to the command after it.  So the data you generated will be parsed with the unix `awk` command.  This will add the values of the third column `$3` to a variable called `sum` and then at the end (`END`) print out the word `Average` followed by the quotient `sum/NR` where NR is the number of rows.

Now it is your turn.  Using the manual for [samtools](http://www.htslib.org/doc/samtools-0.1.19.html) please quantify how many reads mapped to yoru chromosome for each individual.



## OK, if this all went smoothly we are now ready to make some genotype calls.  Please click [here](https://github.com/evansbenj/Reduced-Representation-Workshop/blob/master/5_Automating_alignment_with_a_bash_script.md) to go to the next page.
