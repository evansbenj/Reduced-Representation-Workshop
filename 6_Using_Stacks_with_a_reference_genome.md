# Using Stacks with a reference genome

(Or you can go back to using a bask script [here](https://github.com/evansbenj/Reduced-Representation-Workshop/blob/master/5_Automating_alignment_with_a_bash_script.md)).

## Stacks pipeline and setting up a run

[Stacks](http://creskolab.uoregon.edu/stacks/manual/) us a software suite that compiles reduced representation genome sequencing data first witih each individual in your analysis, and then across individuals.  It can calculate summary statistics from your data and also output your results in useful formats that other programs can take as input.

Stacks actually can do some of the tasks that we have already done in this workshop, such as de-multiplexing Illumina data, trimming off linker sequences, and filtering sequences based on quality (the proportion of Ns in a read).  There is a very nice online manual [here](http://creskolab.uoregon.edu/stacks/manual/) and we will only scratch the surface of what this software can do.

The portion of Stacks that we will use consists of three main steps:
  1. For each individual, sort the data into 'loci' that have one or two alleles (i.e. that are homozygous or heterozygous, respectively).  For analyses with a reference genome, this is accomplished using the program `pstacks` using a `.bam` file as input.  Alternatively if you lack a reference genome, you can use `ustacks` instead of `pstacks`.
  2. Across all individuals, generate a catalog of loci that is a comprehensive list of all genomic regions that have data from at least one individual.  This is accomplished with the `cstacks` program.
  3. Once a catelog of all loci is made, we can compile the data for all individuals to generate a multi-individual genotype for each locus.  This is done with the `sstacks` program

Once loci are compiled within and across individuals, we can use the program `populations` within Stacks to analyze the data. We can also output the data in different formats that can be analyzed with other software such as [`Structure`](http://pritchardlab.stanford.edu/structure.html) and [`Phylip`](http://evolution.genetics.washington.edu/phylip/getme.html).   

The pipeline of programs within `Stacks` can be run in a batch using a `Perl` script that comes with the program called `refmap.pl`.  This script functions in a similar way to the bash scripts we have used already but it has some added features, such as allowing options to be specified using flags.

To get started, lets first make a directory within the `monkey` directory that has our data called `Stacks_Results`.  To do this, first make sure you are in the `monkey` directory by typing this:

`pwd`

If the path that is shown is not the `monkey` directory, please change to it by typing this:

`cd ~/monkey`

OK, now please make a new direcotry by typing this:

`mkdir Stacks_Results`

## Analysis of population structure

Let's first examine whether population structure is present within our sample by defining two populations.  This can be done in `Stacks` by making a file called a `population map`.  Use emails to generate a file like this:

`emacs population_map`

and now enter this (tab delimited) information:

```
PF515_sorted</t>population_1
PM561_sorted/tpopulation_1
PM565_sorted	population_1
PM566_sorted	population_1
PM567_sorted	population_1
PM582_sorted	population_2
PM584_sorted	population_2
PM592_sorted	population_2
PM602_sorted	population_2
```

Note that the term tab-delimited means that there is a tab between the columns of information.  Please type `Ctrl-X` and then `Ctrl-S` to save this file and then `Ctrl-X` and then `Ctrl-C` to exit emacs.  This file will be used to tell Stacks that the first five samples are from one population and the last four samples are from another population.

One way to quantify population structure is using the F-statistic (F<sub>ST</sub>).  F<sub>ST</sub> is an index of population structure that ranges from zero (no population structure) to one (two populations are each fixed for different alleles.  Let's calulate F<sub>ST</sub> between the two populations specified avove using `Stacks`.  To do this, please type:

```
/apps/stacks/1.29/bin/ref_map.pl -S -b 1 -n 0 \
	-O ~/monkey/population_map \
	-o ~/monkey/Stacks_Results \
   	-s ~/monkey/PF515_sorted.bam \
    -s ~/monkey/PM561_sorted.bam \
    -s ~/monkey/PM565_sorted.bam \
    -s ~/monkey/PM566_sorted.bam \
    -s ~/monkey/PM567_sorted.bam \
    -s ~/monkey/PM582_sorted.bam \
    -s ~/monkey/PM584_sorted.bam \
    -s ~/monkey/PM592_sorted.bam \
    -s ~/monkey/PM602_sorted.bam \
   	-e /apps/stacks/1.29/bin -X "populations:--fstat" \
```

In this command, the backslashes `\` just indicate that the command is continued on the next line.  The program we are executing is a perl script caled "ref_map.pl".  Similar to the bash scripts we wrote earlier, this program just executes a bunch of other prorgams.  We can pass some of these programs additional commands using the `-X` flag.  Here we have used this flag at the end to pass the program `populations` a this flag: `--fstat`, which tells the program `populations` to calculate F<sub>ST</sub> using the population map that we specified using the `-O` flag.  We have additionally specified a directory to write our results to using teh `-o` flag.  The `-e` flag tells the computer where the executable files that are referenced by `ref_map.pl` are (these are programs such as `cstacks` and `populations`).

If you now go to the `Stacks_Results` directory (`cd path_to_Stacks_Results_directory`) and list the files in this directory (`ls`) you should see a bunch of compressed files that have a `gz` suffix.  For each sample we have a file whose name includes the word `alleles`, one with `matches`, one with `snps`, and one with `tags`.  Some details of the contents of these files is available in the [Stacks manual](http://creskolab.uoregon.edu/stacks/manual/).  The file called `batch_1.fst_population_1-population_2.tsv` has the results of the F<sub>ST</sub> calculation we requested from the module `populations`.  We can look at this by typing:

`more batch_1.fst_population_1-population_2.tsv`

## Now let's move on to [learn more about Stacks](https://github.com/evansbenj/Reduced-Representation-Workshop/blob/master/7_More_on_Stacks.md).







