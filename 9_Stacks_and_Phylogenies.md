## Using `Stacks` to export a file for phylogenetic analysis.

(Or you can go back to Stacks and Structure [here](https://github.com/evansbenj/Reduced-Representation-Workshop/blob/master/8_Stacks_and_Structure.md).

Much in the same way that we used `Stacks` to export an input file for the program `Structure`, we can also use `Stacks` to export an input file for phylogenetic analysis.  Let's try this now.

Because this (again) takes a little while, Ben did it in advance using these commands.  First he went to the results directory from the complete dataset:

`cd /home/datasets/2015_Ben_Evans/complete_data/monkey/Stacks_Results`

And then he typed this command:

`/apps/stacks/1.29/bin/populations -P ./ -b 1 -r 1 -t 36 --phylip`


