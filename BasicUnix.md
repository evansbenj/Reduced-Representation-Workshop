# Basic Unix use

All of the software in this workshop is available for free and runs on a Unix operating system.  In order to use these applications, we need to know/review a few basic commands that will allow us to move through directories, copy and edit files, and see what is in a directory.

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
`ls`.
This command (short for "list") should list all of the files in the directory you are in.  

If you want more information about these files you can type this:
`ls -l`.
Now each file gets its own line with the name of the file on the right side. You should see some letters and numbers which specify who can access and manipulate each file, and also the file size and date that the file was last modified.

Sometimes we need to see what is inside a file. For text files you can use the more command to scroll through the file from the begining.  Type this:
`more ` **filename**.
where **filename** is the full name (including the suffix) of a file you can see in this directory.  You can use the spacebar to scroll down or type `q` to return to the prompt.


