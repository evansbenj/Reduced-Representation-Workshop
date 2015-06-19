# Basic Unix use

(Or you can go back to the Readme page [here](https://github.com/evansbenj/Reduced-Representation-Workshop/blob/master/1_README.md)).

All of the software in this workshop is available for free and runs on a Unix operating system.  In order to use these applications, we need to know/review a few basic commands that will allow us to move through directories, copy and edit files, and see what is in a directory.

## Some Unix Basics

To connect to the Unix system we are using, please type this:

`ssh username@caf-hpc1.sun.ac.za`

Where `username` is your username.  You will be prompted for a password. Then, to see what directory you are in, type this:

`pwd`

and then hit `return`.  

This command (short for "print working directory") should print out the home directory.  You should see something like this:

`/home/username`

Each forwardslash divides the name of a directory; the directories listed on the right side are inside the ones that are on the left.  To move to a different directory you can use the change directory command:

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

As an exercise, please write and execute a bash script that will execute the script we just made (`fruit_script`) one hundred times.  Please use a for loop such as one of the ones described [here](http://www.cyberciti.biz/faq/bash-for-loop/).


## OK, now lets learn how to de-multiplex Illumina data by clicking [here](https://github.com/evansbenj/Reduced-Representation-Workshop/blob/master/3_De-multiplexing.md).
