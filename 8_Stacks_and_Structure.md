# Using `Stacks` to export an input file for `Structure`

(Or you can go back to using `Stacks` to calculate summary statistics [here](https://github.com/evansbenj/Reduced-Representation-Workshop/blob/master/7_More_on_Stacks.md)).

[`Structure`](http://pritchardlab.stanford.edu/structure_software/release_versions/v2.3.4/html/structure.html) is a software that attempts to assign individuals to *k* populations in such a way as to minimize Hardy-Weinberg and linkage disequilibrium.  We run `Structure` by specifying multiple values of *k* and then seeing which value(s) maximuze the likelihood of the data given the model of population structure. We can generate an input file for this program using `Stacks`. 

The program `Structure` can not handle all of our data, so let's select 1000 loci randomly to analyze like this:

`shuf -n 1000 batch_1.catalog.tags.tsv | awk '{print $3}'  > 1000_randoms`

This uses another `unix` command called `shuf`.  This says to print a randomly selected value from column 3 1000 times to a file called `1000_randoms`.

Now we are ready to generate an input file for `Structure`.  Please type this:

`/work/ben/workshop_software/stacks-1.30/populations -P /work/ben/2015_workshop/complete_data/monkey/Stacks_Results -b 1 -r 1 -t 36 --structure --write_single_snp -W ./1000_randoms -B chrX_list`

This command directs the `populations` module of `Stacks` to output a single snp (the `--write_single_snp` flag) from tags specified by the `1000_randoms` file (the `-W` tag) but not to include any snps from chromosome X (the `-B` flag).  

This will generate a file called `batch_1.structure.tsv` which can be used as an input file for the program `Structure`.

Before we can run `Structure` with out data, we have some clerical issues to take care of.  Firstly, lets delete the first row from our input file  - this row has a comment that we don't need.  Please type this:

`tail -n+2 batch_1.structure.tsv > simple_structure.tsv`

This uses the Unix command `tail` to feed all lines beginning with the second line to a new file called `simple_structure.tsv`.

Because we required that there be no missing data from the SNPs in our analysis, we can get a differing number of SNPs in our input file depending on which sites were randomly chosen. So we need to count how many columns we have in this file.  Please type this:

`head -n1 simple_structure.tsv |  sed 's/\t/\n/g' | wc -l`

This should output the number of tab spaces in the first column of the file `simple_structure.tsv`.  Because the first column begins with a tab space, the actual number of loci is equal to this number **minus one**.

Now we are ready to run `Structure`.  Please type this command:

`/work/ben/workshop_software/console/structure -m /work/ben/workshop_software/console/mainparams -e /work/ben/workshop_software/console/extraparams -K 3 -L 279 -N 9 -i simple_structure.tsv -o output_K_3`

This tells the system to execute the `Structure` program and it specifies the paths to two files (`mainparams` and `extraparams`) that are used in the analysis.  It then has flags for the number of populations (`-K`), the number of loci (`-L`; based on the number we got above from the `head` command), the number of individuals (`-N`; our study has 9 individuals), the input file (`-i`) and an output file (`-o`).

 Assuming the command executes without error, you can check out the results in the file `output_K_3` like this:
 
 `more output_K_3`
 
 You could scroll down to the population assignments, which should look something like this:
 ```
 Inferred ancestry of individuals:
        Label (%Miss) Pop:  Inferred clusters
  1 PF515_sorte    (0)    1 :  0.687 0.150 0.164 
  2 PM561_sorte    (0)    1 :  0.583 0.194 0.222 
  3 PM565_sorte    (0)    1 :  0.635 0.198 0.167 
  4 PM566_sorte    (0)    1 :  0.665 0.171 0.164 
  5 PM567_sorte    (0)    1 :  0.591 0.205 0.204 
  6 PM582_sorte    (0)    1 :  0.674 0.149 0.177 
  7 PM584_sorte    (0)    1 :  0.613 0.184 0.203 
  8 PM592_sorte    (0)    1 :  0.591 0.200 0.209 
  9 PM602_sorte    (0)    1 :  0.521 0.226 0.253 
```

This tells us, for each individual, what the probability that that individual is assigned to each one of *k* populations (three in this case).

## A quick example of plotting with `R`

We can plot this by making a file and pasting these data in this file:

`emacs assignments`

If you now paste the data beginning with the line `1 PF515_sorte...` and save it (`Ctrl-X` and then `Ctrl-S`) and exit (`Ctrl-X` and then `Ctrl-C`) we can easily plot the data with R.  To do this type: `R`.  This should open up the `R` environment.  Now import the data:

`>` `dat<-read.table(file="assignments")`

and make a pdf..

`>`  `pdf("temp.pdf")`

and make a barplot

`>` `barplot(as.matrix(t(dat[,6:8])))`

This command tells R to make a boxplot using columns 6 thru 8 of the table called `dat`.  We have transformed these data for this plot using `t()`.

and now exit the `R` environment:

`>` `q()`


