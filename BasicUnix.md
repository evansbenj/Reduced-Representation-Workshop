# Basic Unix use

All of the software in this workshop is available for free and runs on a Unix operating system.  In order to use these applications, we need to know/review a few basic commands that will allow us to move through directories, copy and edit files, and see what is in a directory.

## Some Unix Basics

To connect to the Unix system we are using, do this **(to be edited later)**:

`ssh XXXX`

To see what directory you are in, type this:

`pwd`

and then hit `return`.

This command (short for "print working direcotry" should print out the home directory.  Each forwardslash divides the name of a directory; the directories listed on the right side are inside the ones that are on the left.

To move to a different directory you can use the change directory command:

`cd`

followed by the path (the directory structure) that you want to go to.  Let's go to the directory that contains the RADseq dataset.  Please type this:

`cd` **(add the path here)**

Now let see what is in this directory.  Please type this:

`ls`

This command (short for "list") should list all of the files in the directory you are in.  

If you want more information about these files you can type this:

`ls -l`

Now each file gets its own line with the name of the file on the right side. You should see some letters and numbers which specify who can access and manipulate each file, and also the file size and date that the file was last modified.

Sometimes we need to see what is inside a file. For text files you can use the more command to scroll through the file from the begining.  Please type this:

`more ` **filename**.

where **filename** is the full name (including the suffix) of a file you can see in this directory.  You can use the spacebar to scroll down or type `q` to return to the prompt.

If you need to rename a file you can type this:

`mv oldfinename newfilename`

Here `mv` is short for "move".

## Using `screen` in Unix

The `screen` command allows you to execute a job on a server and then log out.  This is useful for our purposes because RRGS generates large datasets that sometimes take a while to analyze.  We generally do data processing on a server (and not on your desktop computer), and we don't want to have to sit and watch our computer for hours or days while a program runs.  Like all commands in Unix, you can look at the manual for a command by typing `man ' and then the command.  So, please type this:

`man screen`

This should bring up a description of the details of the command.  You can scroll through it by hitting any key or exit by typing `q`.

When I am running multiple screens, I usually name them like this:

`screen -S add_description_here`

This brings you into a screen session called `add_description_here`.  You can exit it by typing `Ctrl-A` and then `Crtl-D`, meaning you hold down the control button and type `A` and then hold down the control button again and type `D`.

You can then see what screens are running by typing this:

`screen -list`

and return to an active screen by typing this:

`screen -r -S add_description_here`

or end an active screen like this:

`screen -X -S add_description_here kill`

Here the `-X` flag tells screen that a command is to be delivered to the screen named `add_description_here`.  The `kill` command at the end tells screen to end this screen.

## Text editing and making a bash script

It is often useful to edit text files on a server and most servers have software installed to facilitate this.  The one I usually use is called "emacs" (http://www.gnu.org/software/emacs/); another popular text editor is "vim" (http://www.vim.org/about.php).




