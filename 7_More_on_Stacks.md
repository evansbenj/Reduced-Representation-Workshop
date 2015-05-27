# More on Stacks:  Whitelists, blacklists, using individual modules, and summary statistics

(Or you can go back to [Using Stacks with a Reference Genome](https://github.com/evansbenj/Reduced-Representation-Workshop/blob/master/6_Using_Stacks_with_a_reference_genome.md)).

Now lets try to calculate summary statistics from specific genomic regions using the `populations` module of stacks. To accomplish this, lets work with a larger dataset that I made earlier.  Please switch to this directory:

`cd **insert_directory_with_complete_data_here***`

Let's first uncompress one of the results file like this:

`gunzip batch_1.catalog.tags.tsv.gz`

and then look at the compressed file like this:

`more batch_1.catalog.tags.tsv`

The first line begins with a hash (`#`) symbol and is reserved for comments.  The next lines have columns of text.  The 4th column lists the chromosome number and position of each tag.  We are going to generate a file in which we sample random SNPs from 1000 tags but we want to exclude data from the X chromosome because there are differences in copy number between males and females (i.e. two in XX females and one in XY males).  In order to do this, we can create a list of tags to include (a `whitelist`) or exclude (a `blacklist`), which is just a list of the numbers in the 3rd column that correspond with the `chrX` in the 4th column.  To generate a list based on this criterion, please use this `unix` command:

`awk '$4 ~ /chrX/ {print $3}' batch_1.catalog.tags.tsv > chrX_list`

This uses a `unix` function called `awk`.  It says to print the number in column 3 to a file called `chrX_list` whenever the value in column 4 is equal to `chrX`.

We can view this file by typing:

`more chrX_list`

Now, if we want to calculate pairwise nucleotide diversity on the autosomes, we can use the `populations` module of `Stacks` with the `blacklist` flag as follows:

`/work/ben/workshop_software/stacks-1.30/populations -P /work/ben/2015_workshop/complete_data/monkey/Stacks_Results -b 1 -r 1 -t 36 -B chrX_list`

We can view the pairwise nucleotide diversity statistic in a file called `batch_1.sumstats_summary.tsv`.

Similarly, if we want to calculate pairwise nucleotide diversity on only the X chromosome, we can use the `chrX_list` as a `whitelist` as follows:

`/work/ben/workshop_software/stacks-1.30/populations -P /work/ben/2015_workshop/complete_data/monkey/Stacks_Results -b 1 -r 1 -t 36 --structure --write_single_snp -W ./1000_randoms -W chrX_list`

## Running `Structure`

`Structure` is a software that attempts to assign individuals to *k* populations in such a way as to minimize Hardy-Weinberg and linkage disequilibrium.  We run structure by specifying multiple values of *k* and then seeing which value(s) maximuze the likelihood of the data given the model of population structure. We can output a file that can be analyzed with the program `Structure` to give us an idea about whether or not our sample has population structure. 

The program `Structure` can not handle all of our data, so let's select 1000 loci randomly to analyze like this:

'shuf -n 1000 batch_1.catalog.tags.tsv | awk '{print $3}'  > 1000_randoms'

This uses another `unix` command called `shuf`.  This says to print a randomly selected value from column 3 1000 times to a file called `1000_randoms`.

Now we are ready to generate an input file for `Structure`.  Please type this:

`/work/ben/workshop_software/stacks-1.30/populations -P /work/ben/2015_workshop/complete_data/monkey/Stacks_Results -b 1 -r 1 -t 36 --structure --write_single_snp -W ./1000_randoms -B chrX_list`

This command directs the `populations` module of `Stacks` to reute results to a directory specified by the `-P` flag.  It tells `populations` to output a single snp (the `--write_single_snp` flag) from tags specified by the `1000_randoms` file (the `-W` tag) but not to include any snps from chromosome X (the `-B` flag).  The `-r` flag says we want to only print data where 100% of the individuals have a genotype.  The `-b` and `-t` flags specify, respectively, that `populations` should focus on batch_ID number 1 (you can work with special IDs if you need to but this is beyond the scope of this workshop) and that `populations` should use 36 threads to do the calculations.

This will generate a file called `batch_1.structure.tsv` which can be used as an input file for the program `Structure`.
