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

`cd path_to_monkey_directory`

OK, now please make a new direcotry by typing this:

`mkdir Stacks_Results`

## Analysis of population structure

Let's first examine whether population structure is present within our sample by defining two populations.  This can be done in `Stacks` by making a file called a `population map`.  Use emails to generate a file like this:

`emacs population_map`

and now enter this (tab delimited) information:

```
PF515_sorted	population_1
PM561_sorted	population_1
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
**Path_to_stacks***/stacks-1.30/scripts/ref_map.pl -S -b 1 -n 0 \
	-O **Path_to_data***/population_map \
	-o **Path_to_data***/Stacks_Results \
   	-s **Path_to_data***/PF515_sorted.bam \
    -s **Path_to_data***/PM561_sorted.bam \
    -s **Path_to_data***/PM565_sorted.bam \
    -s **Path_to_data***/PM566_sorted.bam \
    -s **Path_to_data***/PM567_sorted.bam \
    -s **Path_to_data***/PM582_sorted.bam \
    -s **Path_to_data***/PM584_sorted.bam \
    -s **Path_to_data***/PM592_sorted.bam \
    -s **Path_to_data***/PM602_sorted.bam \
   	-e **Path_to_stacks***/stacks-1.30/ -X "populations:--fstat" \
```

In this command, the backslashes `\` just indicate that the command is continued on the next line.  The program we are executing is a perl script caled "ref_map.pl".  Similar to the bash scripts we wrote earlier, this program just executes a bunch of other prorgams.  We can pass some of these programs additional commands using the `-X` flag.  Here we have used this flag at the end to pass the program `populations` a this flag: `--fstat`, which tells the program `populations` to calculate F<sub>ST</sub> using the population map that we specified using the `-O` flag.  We have additionally specified a directory to write our results to using teh `-o` flag.  The `-e` flag tells the computer where the executable files that are referenced by `ref_map.pl` are (these are programs such as `cstacks` and `populations`).

If you now go to the `Stacks_Results` directory (`cd path_to_Stacks_Results_directory`) and list the files in this directory (`ls`) you should see a bunch of compressed files that have a `gz` suffix.  For each sample we have a file whose name includes the word `alleles`, one with `matches`, one with `snps`, and one with `tags`.  Some details of the contents of these files is available in the [Stacks manual](http://creskolab.uoregon.edu/stacks/manual/).  The file called `batch_1.fst_population_1-population_2.tsv` has the results of the F<sub>ST</sub> calculation we requested from the module `populations`.  We can look at this by typing:

`more batch_1.fst_population_1-population_2.tsv`

## Using individual modules within Stacks

Now lets try to do some more stuff using the `populations` module of stacks.  For example, we can output a file that can be analyzed with the program `Structure` to give us an idea about whether or not our sample has population structure.  To accomplish this, lets work with a larger dataset that I made earlier.  Please switch to this directory:

`cd **insert_directory_with_complete_data_here***`

Let's first uncompress one of the results file like this:

`gunzip batch_1.catalog.tags.tsv.gz`

and then look at the compressed file like this:

`more batch_1.catalog.tags.tsv`

The first line begins with a hash (`#`) symbol and is reserved for comments.  The next lines have columns of text.  The 4th column lists the chromosome number and position of each tag.  We are going to generate a file in which we sample random SNPs from 1000 tags but we want to exclude data from the X chromosome because there are differences in copy number between males and females (i.e. two in XX females and one in XY males).  In order to do this, we can create a `blacklist` of tags to exclude, which is just a list of the numbers in the 3rd column that correspond with the `chrX` in the 4th column.  To generate a `blacklist` based on this criterion, please use this `unix` command:

`awk '$4 ~ /chrX/ {print $3}' batch_1.catalog.tags.tsv > chrX_blacklist`

This uses a `unix` function called `awk`.  It says to print the number in column 3 to a file called `chrX_blacklist` whenever the value in column 4 is equal to `chrX`.

We can view this file by typing:

`more chrX_blacklist`

The program `Structure` can not handle all of our data, so let's select 1000 loci randomly to analyze like this:

'shuf -n 1000 batch_1.catalog.tags.tsv | awk '{print $3}'  > 1000_randoms'

This uses another `unix` command called `shuf`.  This says to print a randomly selected value from column 3 1000 times to a file called `1000_randoms`.

Now we are ready to generate an input file for `Structure`.  Please type this:

`/work/ben/workshop_software/stacks-1.30/populations -P /work/ben/2015_workshop/complete_data/monkey/Stacks_Results -b 1 -r 1 -t 36 --structure --write_single_snp -W ./1000_randoms -B chrX_blacklist`

This command directs the `populations` module of `Stacks` to reute results to a directory specified by the `-P` flag.  It tells `populations` to output a single snp (the `--write_single_snp` flag) from tags specified by the `1000_randoms` file (the `-W` tag) but not to include any snps from chromosome X (the `-B` flag).  The `-r` flag says we want to only print data where 100% of the individuals have a genotype.  The `-b`, and `-t` flags specify, respectively, that `populations` should focus on batch_ID number 1 (you can work with special IDs if you need to but this is beyond the scope of this workshop) and that `populations` should use 36 threads to do the calculations.




