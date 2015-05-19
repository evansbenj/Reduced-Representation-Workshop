# Using Stacks with a reference genome

(Or you can go back to using a bask script [here](https://github.com/evansbenj/Reduced-Representation-Workshop/blob/master/5_Automating_alignment_with_a_bash_script.md)).

## Stacks pipeline and setting up a run

[Stacks](http://creskolab.uoregon.edu/stacks/manual/) us a software suite that compiles reduced representation genome sequencing data first witih each individual in your analysis, and then across individuals.  It can calculate summary statistics from your data and also output your results in useful formats that other programs can take as input.

Stacks actually can do some of the tasks that we have already done in this workshop, such as de-multiplexing Illumina data, trimming off linker sequences, and filtering sequences based on quality (the proportion of Ns in a read).  There is a very nice online manual [here](http://creskolab.uoregon.edu/stacks/manual/) and we will only scratch the surface of what this software can do.

The portion of Stacks that we will use consists of three main steps:
  1. For each individual, sort the data into 'loci' that have one or two alleles (i.e. that are homozygous or heterozygous, respectively).  For analyses with a reference genome, this is accomplished using the program `pstacks` using a `.bam` file as input.  Alternatively if you lack a reference genome, you can use `ustacks` instead of `pstacks`.
  2. Across all individuals, generate a catalog of loci that is a comprehensive list of all genomic regions that have data from at least one individual.  This is accomplished with the `cstacks` program.
  3. Once a catelog of all loci is made, we can compile the data for all individuals to generate a multi-individual genotype for each locus.  This is done with the `sstacks` program

These programs can be run in a batch using a `Perl` script that comes with the program called `refmap.pl`.  This script functions in a similar way to the bash scripts we have used already but it has some added features, such as allowing options to be specified using flags.

To get started, lets first make a directory within the `monkey` directory that has our data called `Stacks_Results`.  To do this, first make sure you are in the `monkey` directory by typing this:

`pwd`

If the path that is shown is not the `monkey` directory, please change to it by typing this:

`cd path_to_monkey_directory`

OK, now please make a new direcotry by typing this:

`mkdir Stacks_Results`

## Analysis of population structure

Let's first examine whether population structure is present within our sample by defining two populations.  This can be done in `Stacks` by making a file called a `population map`.  Use emails to generate a file like this:

`emacs population_map`

and now enter this information:

```
PF515 population_1 
PM561 population_1                                                                                 
PM565 population_1                                                                                    
PM566 population_1                                                                       
PM567 population_1                                                                                    
PM582 population_2                                                                                    
PM584 population_2                                                                                    
PM592 population_2                                                                                    
PM602 population_2
```

This file will be used to tell Stacks that the first five samples are from one population and the last four samples are from another population.

One way to quantify population structure is using the F-statistic (F<sup>ST</sup>).

## Summary Statistics

## Whitelists and Blacklists

## Other information
