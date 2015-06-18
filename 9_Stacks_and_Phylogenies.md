## Using `Stacks` to export a file for phylogenetic analysis.

(Or you can go back to Stacks and Structure [here](https://github.com/evansbenj/Reduced-Representation-Workshop/blob/master/8_Stacks_and_Structure.md).

Much in the same way that we used `Stacks` to export an input file for the program `Structure`, we can also use `Stacks` to export an input file for phylogenetic analysis.  Let's try this now.

Because this (again) takes a little while, Ben did it in advance using the following commands (You do not need to type these).  First he went to the results directory from the complete dataset:

`cd /home/datasets/2015_Ben_Evans/complete_data/monkey/Stacks_Results`

...and then he typed this command:

`/apps/stacks/1.29/bin/populations -P ./ -b 1 -r 1 -t 36 --phylip --phylip_var`

Most of these flags were discussed previously.  The `--phylip` flag tells `Stacks` to output data in `Phylip` format, which is the format of the input file for Joe Felsenstein's `Phylip` package.  This can be easily modified for other programs, such as the `nexus` format.  The `--phylip_var` flag tells `Stacks` to output sites that are variable within populations.  We have only one population specified here because we are not using the population map file we created earlier.  If we did not include the `--phylip_var` flag, no sites would be exported at all because by default `Stacks` only exports sites that vary between populations.

This generated a file called `xxx`.  

Please copy this file to your home directory like this:

`scp /home/datasets/2015_Ben_Evans/complete_data/monkey/Stacks_Results/batch_1.phylip.tsv ~/monkey/Stacks_Results`


