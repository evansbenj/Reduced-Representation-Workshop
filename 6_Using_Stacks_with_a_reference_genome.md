# Using Stacks with a reference genome

(Or you can go back to using a bask script [here](https://github.com/evansbenj/Reduced-Representation-Workshop/blob/master/5_Automating_alignment_with_a_bash_script.md).

[Stacks]{http://creskolab.uoregon.edu/stacks/manual/} us a software suite that compiles reduced representation genome sequencing data first witih each individual in your analysis, and then across individuals.  It can calculate summary statistics from your data and also output your results in useful formats that other programs can take as input.

Stacks actually can do some of the tasks that we have already done in this workshop, such as de-multiplexing Illumina data, trimming off linker sequences, and filtering sequences based on quality (the proportion of Ns in a read).  There is a very nice online manual [here](http://creskolab.uoregon.edu/stacks/manual/) and we will only scratch the surface of what this software can do.

The portion of Stacks that we will use consists of three main steps:
  1. For each individual, sort the data into 'loci' that have one or two alleles (i.e. that are homozygous or heterozygous, respectively).  For analyses with a reference genome, this is accomplished using the program `pstacks` using a `.bam` file as input.  Alternatively if you lack a reference genome, you can use `ustacks` instead of `pstacks`.
  2. Across all individuals, generate a catalog of loci that is a comprehensive list of all genomic regions that have data from at least one individual.  This is accomplished with the `cstacks` program.
  3. Once a catelog of all loci is made, we can compile the data for all individuals to generate a multi-individual genotype for each locus.  This is done with the `sstacks` program

These programs can be run in a batch using a `Perl` script that comes with the program called `refmap.pl`.  This script functions in a similar way to the bash scripts we have used already but it has some added features, such as allowing options to be specified using flags.

