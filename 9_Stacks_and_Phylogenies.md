## Using `Stacks` to export a file for phylogenetic analysis.

(Or you can go back to Stacks and Structure [here](https://github.com/evansbenj/Reduced-Representation-Workshop/blob/master/8_Stacks_and_Structure.md).

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

`/apps/stacks/1.29/bin/populations -P ./ -b 1 -r 1 -t 36 -M population_map_all_different --phylip`

Most of these flags were discussed previously.  The `--phylip` flag tells `Stacks` to output data in `Phylip` format, which is the format of the input file for Joe Felsenstein's `Phylip` package.  This can be easily modified for other programs, such as the `nexus` format.  The `-M` flag tells `Stacks` to use the new population map file in which each individual is assigned to a different population.

This generated a file called `batch_1.phylip.tsv`.  

Please copy this file to your home directory like this:

`scp /home/datasets/2015_Ben_Evans/complete_data/monkey/Stacks_Results/batch_1.phylip.tsv ~/monkey/Stacks_Results`


