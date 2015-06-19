## Using `Stacks` to export a file for phylogenetic analysis.

(Or you can go back to Stacks and Structure [here](https://github.com/evansbenj/Reduced-Representation-Workshop/blob/master/8_Stacks_and_Structure.md)).

Much in the same way that we used `Stacks` to export an input file for the program `Structure`, we can also use `Stacks` to export an input file for phylogenetic analysis.  Let's try this now.

Because this (again) takes a little while, Ben did it in advance using the following commands (You do not need to type these).  First he went to the results directory from the complete dataset:

`cd /home/datasets/2015_Ben_Evans/complete_data/monkey/Stacks_Results`

Then he made a population map file called `population_map_all_different` in which all individuals were assigned to a different population.  This tab-delimited file looks like this:

```
PF515_sorted    population_1
PM561_sorted    population_2
PM565_sorted    population_3
PM566_sorted    population_4
PM567_sorted    population_5
PM582_sorted    population_6
PM584_sorted    population_7
PM592_sorted    population_8
PM602_sorted    population_9
```

...and then he typed this command:

`/apps/stacks/1.29/bin/populations -P ./ -b 1 -r 1 -t 36 -M population_map_all_different --phylip --phylip_var`

Most of these flags were discussed previously.  The `--phylip` flag combined with the `--phylip_var` flag tells `Stacks` to output sites that are variable between and within populations in `Phylip` format, which is the format of the input file for Joe Felsenstein's `Phylip` package.  This can be easily modified for other programs, such as the `nexus` format.  The `-M` flag tells `Stacks` to use the new population map file in which each individual is assigned to a different population.

This generated a file called `batch_1.phylip.tsv`.  

Please copy this file to your home directory like this:

`scp /home/datasets/2015_Ben_Evans/complete_data/monkey/Stacks_Results/batch_1.phylip ~/monkey/Stacks_Results`

Let's have a look at this file now.  Please type this:

`more batch_1.phylip`

You should be able to see a Phylip formatted file. You can press the space bar to scroll down the file.  The first line has the number of taxa (9 in this case) followed by the number of characters of data (70484 in this case).  The next line has a taxon name (`1`, which corresponds to the first sample PF515) followed by the sequence data.  Some of these data are regular nucleotides (A, C, G, or T) and others are [IUPAC symbols](http://www.bioinformatics.org/sms/iupac.html) that indicate heterozygous SNPs (Y for C/T, R for A/G, etc).  The lines after this give data for the next samples.

## Making a Quick phylogeny using `Phylip`

Now we can use a program called `Phylip` to make a quick phylogenetic tree.  We will make a neighborjoining tree because it is very quick to do.  For your studies I recommend to instead use maximum likelihood or Bayesian methods to make phylogenetic trees, for example using software such as [MrBayes](http://mrbayes.sourceforge.net/), [BEAST](http://beast.bio.ed.ac.uk/), and [RaxML](http://sco.h-its.org/exelixis/web/software/raxml/index.html).

To make a phylogenetic tree with Phylip, please type this:

`/apps/PHYLIP/3.696/exe/neighbor batch_1.phylip`

This should generate a tree file called `batch_1.phylip.tree` in your home directory.  You can copy it to your local computer by opening up another session and typing this:

`scp username@caf-hpc1.sun.ac.za:~/monkey/Stacks_Results/batch_1.phylip.tree .`

You can now view the tree you made using the [FigTree](http://tree.bio.ed.ac.uk/software/figtree/) program that should be available on your local computer.




