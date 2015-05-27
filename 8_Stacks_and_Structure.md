# Using `Stacks` to export an input file for `Structure`

Or you can go back to using `Stacks` to calculate summary statistics [here](https://github.com/evansbenj/Reduced-Representation-Workshop/blob/master/7_More_on_Stacks.md)).

`Structure` is a software that attempts to assign individuals to *k* populations in such a way as to minimize Hardy-Weinberg and linkage disequilibrium.  We run structure by specifying multiple values of *k* and then seeing which value(s) maximuze the likelihood of the data given the model of population structure. We can output a file that can be analyzed with the program `Structure` to give us an idea about whether or not our sample has population structure. 

The program `Structure` can not handle all of our data, so let's select 1000 loci randomly to analyze like this:

'shuf -n 1000 batch_1.catalog.tags.tsv | awk '{print $3}'  > 1000_randoms'

This uses another `unix` command called `shuf`.  This says to print a randomly selected value from column 3 1000 times to a file called `1000_randoms`.

Now we are ready to generate an input file for `Structure`.  Please type this:

`/work/ben/workshop_software/stacks-1.30/populations -P /work/ben/2015_workshop/complete_data/monkey/Stacks_Results -b 1 -r 1 -t 36 --structure --write_single_snp -W ./1000_randoms -B chrX_list`

This command directs the `populations` module of `Stacks` to output a single snp (the `--write_single_snp` flag) from tags specified by the `1000_randoms` file (the `-W` tag) but not to include any snps from chromosome X (the `-B` flag).  

This will generate a file called `batch_1.structure.tsv` which can be used as an input file for the program `Structure`.
