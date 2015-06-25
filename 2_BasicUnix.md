# Basic Unix use

(Or you can go back to the Readme page [here](https://github.com/evansbenj/Reduced-Representation-Workshop/blob/master/1_README.md)).

## A quick note about genetic samples

Your statistical power, precision, and accuracy will depend on the quantity and quality of your data.  These factors, of course, depend on the statring material you use for sequencing.  Excitingly, a bunch of recent studies suggest that RRGS can be used on sub-optimal samples, for example from museum specimens, fecal extractions, and the like.  Nonetheless, it is in your interest to use as high quality DNA as possible.  For animal tissue, I recommend using ethanol or RNAlater to preserve your tissues and I suggest chilling your samples (at -20 or -80 degrees) as soon as possible after collection.

For DNA extraction, I've had success using Qiagen DNAEasy extraction kits.  These kits can be used in any lab that has a heat block and a centrifuge.  When a small amount of starting material is being used I recommend reducing the volume of elution buffer you use in order to concentrate the DNA.  You can always dilute it later, and it is harder and less efficient to make a sample more concentrated.

If you want to outsource the library preparation (I always do), I recommend using an agarose gel to normalize the concentration of your samples and ensure that the quality is as good as possible.  Different methods have different requirements in terms of the volume and concentration of gDNA, but most expect these parameters to be uniform across the samples you submit.  If you have many samples, I suggest extracting and comparing more than you plan to run so you can choose the best ones. 

## How much does this cost and where can I do it?

I've done RADseq multiple times at [Floragenix](http://www.floragenex.com/), which is located in Oregon, USA.  I've also done Genotype by Sequencing at [Cornell University](http://www.biotech.cornell.edu/brc/genomic-diversity-facility) in New York, USA.  The prices for a 95 RADseq sample run, including library preparation and one lane of single end Illumina sequencing but no bioinformatics was  US$7725.  The price of a 95 GBS run with two single end lanes of Illumina sequencing was $6,080 (considerably less expensive). I think RADseq services are also available at the University of Edinburgh.  I anticipate that the cost of these services will decline considerably over the next few years and that more sequencing centers will offer this service.

## Some Unix Basics

All of the software in this workshop is available for free and runs on a Unix operating system.  In order to use these applications, we need to know/review a few basic commands that will allow us to move through directories, copy and edit files, and see what is in a directory.

To connect to the Unix system we are using, please use thr application PuTTy to open a connection to this server:

`caf-hpc1.sun.ac.za`

You will be prompted for your `username` and password. Once you are connected, to see what directory you are in, type this:

`pwd`

and then hit `return`.  

This command (short for "print working directory") should print out the home directory.  You should see something like this:

`/home/username`

Each forwardslash divides the name of a directory; the directories listed on the right side are inside the ones that are on the left.  

This cluster (like most clusters) consists of a 'head' node that acts as a portal between the cluster and the rest of the world.  Often the head node is not designed for intensive computation - instead it is just used to dispatch jobs to other processors that have more power.  The cluster we are using has a quequing system.  For the sake of simplicity, we have been granted permission by the system administrator to run jobs directly on a powerful node without using the quequing system. To access this node, please type this:

`ssh n01.hpc`

If you type `pwd`, you will see that you are still in the same directory as before, even though we are now working in a different computer.  The files we generate on this computer will be accessible from any node on the cluster.

To move to a different directory you can use the change directory command:

`cd`

followed by the path (the directory structure) that you want to go to.  Let's go to the directory that contains the RADseq dataset.  Please type this:

`cd /home/datasets/2015_Ben_Evans/data` 

Now let see what is in this directory.  Please type this:

`ls`

This command (short for "list") should list all of the files in the directory you are in.  

If you want more information about these files you can type this:

`ls -l`

Now each file gets its own line with the name of the file on the right side. You should see some letters and numbers which specify who can access and manipulate each file, and also the file size and date that the file was last modified.

Sometimes we need to see what is inside a file. For text files you can use the more command to scroll through the file from the begining.  Please type this:

`more ` **filename**.

where **filename** is the full name (including the suffix) of a file you can see in this directory.  For example you could look at the datafile we will work with by typing this:

`more forward_subset.fastq`

You can use the spacebar to scroll down or type `q` to return to the prompt.  Later, if you happen to want to look at a compressed file that has the suffix `.gz`, you can type this:

`zmore ` **filename**.

This current directory doesn't have any compressed files, but we will use this command later.

If you need to rename a file you can type this:

`mv oldfinename newfilename`

Here `mv` is short for "move".

If you would like to copy a file from the server to your local computer, you can use the `scp` command.  I usually open a separate window for this and then type:

`scp username@server_address:path_to_file/file path_to_destination_where_I_want_the_file`

For example, if we wanted to save the datafile we will use from the server to our local machine, we could type this:

`scp username@caf-hpc1.sun.ac.za:/home/datasets/2015_Ben_Evans/data/forward_subset.fastq .`

Note that a colon `:` separates the server address from the path on the server.  Also a space separates the source file name from the destination path.

## Text editing and making a bash script

It is often useful to edit text files on a server and most servers have software installed to facilitate this.  The one I usually use is called [emacs](http://www.gnu.org/software/emacs/); another popular text editor is [vim](http://www.vim.org/about.php).

As an example, lets use `emacs` to make a simple script that we can run on the server.  Please type this:

`emacs fruit_script`

This should open an `emacs` session.  Please now type this:

```
#!/bin/bash

fruits="apples pears grapefruits pineapples"

for type_of_fruit in $fruits
do
	echo ${type_of_fruit}
done
```

Now save this by typing `Ctrl-x` and `Ctrl-s`.  (As a reminder, this just means that you should hold down the `Ctrl` button while you press the `x` button and then hold down the `Ctrl` button while you press the `s` button.  Now you can exit by typing `Ctrl-x` and then `Ctrl-c`.

Using the list command described above (`ls`) you should be able to see your file in the directory you are in.

In order to make this text file into an executeable program, we need to change a characteristic of the program called its permissions.  This can be done like this:

`chmod +x fruit_script`

Here the `+x` means that we want to add the ability to execute the program.  And we can run the program like this:

`./fruit_script`

Here the dot forwardslash (`./`) tells the computer that you want it to look in the directory that you are currently in for the program called `fruit_script`.

If your script lacks errors, you should see the names of the fruit you had in the program.  If it gives you an error, please let Ben know.

## Problem 1: Now try it on your own

As an exercise, please write and execute a bash script that will execute the script we just made (`fruit_script`) one hundred times.  Please use a for loop such as one of the ones described [here](http://www.cyberciti.biz/faq/bash-for-loop/).  When you get this script working, please use the `scp` command to copy your script to your local computer.  


## OK, now lets learn how to de-multiplex Illumina data by clicking [here](https://github.com/evansbenj/Reduced-Representation-Workshop/blob/master/3_De-multiplexing.md).
