# De-multiplexing Illumina Data

Most RRGS methods rely on the Illumina sequencing platform.  These machines generate data using something called a "flowcell" that is divided up into eight "lanes".  Small scale projects typically would run multiple samples (from different species or different individuals within a species) on one lane.  

A first step in our analysis pipeline is to organize data from each of our samples that were run together on an Illumina lane (De-multiplexing our data) and also to filter our data and trim off bits that have lots of errors or that have sequences from the laboratory procedures that were used to generate the data (Trimming/Quality control).

## De-Multiplexing
When samples are run on an Illumina machine, DNA is broken up into many small fragments and a small bit of DNA called an adaptor is then added on each of the fragments. This adaptor allows the sequencing process to occur, essentially by making possible high-throughput put polymerase chain reaction (ask Ben about this if you are unfamiliar). To make possible the multiplexing of samples on one Illumina lane, each sample is linked to a unique adaptor that contains a "barcode" sequence that allows us to sort out which samples each sequence came from.  For our dataset, we have nine individuals from one species (the Tonkean macaque). Each of the samples received the following barcodes (the sample name is followed by the barcode):

`PF515 CCTCTTATCA`

`PM561 TATCGTTAGT`

`PM565 TAGTGCGGTC`

`PM566 GGCCGGTAAC`

`PM567 AGGAACCTCG`

`PM582 TTATCCGTAG`

`PM584 CGCTATACGG`

`PM592 CACGCAACGA`

`PM602 ATCCGTCTAC`

We will use this information in a moment to de-multiplex our data.  Lets use emacs to make a text file that contains this information.  Please type `emacs monkey.pools` to generate a file called `monkey.pools`.  Now type the information above in the file.  When you are done type `Ctrl-X` and `Ctrl-S` to save it and then `Ctrl-X` and `Ctrl-C` to close the program.

## Fasta and Fastq format

Illumina sequence data is provided in a text file that is in a format called `fastq`.  This is a modification of another format called `fasta` in which each sequence has a header that begins with a `>` sign.  This is followed by the sequences.  Here is an example:

`>example_sequence_in_fasta_format`

`ATGCGCGCGCTAGGCTCGCGATCGGGGAGCGCGAGCTGAGCTAGCGCGATGCGCCCCGAC`

The format of `fastq` files is similar to `fasta` except that quality scores are included.  Each sequence has four lines (instead of two for `fasta` files).  The first begins with `@` followed by information about the sequence.  The second line is the nucleotide sequence. The third line is a `+` which may be followed by the same information that followed the `@` sign in the first line.  The fourth line is the quality values.  For the Illumina data we will be working with, these values range from 0–41 and are represented by single characters.  Further information is provided [here](http://en.wikipedia.org/wiki/FASTQ_format).


## Quality control and trimming

Similar to Sanger sequencing, Illumina generates sequences that have errors.  Errors typically become more common towards the end of the sequence read, and the software will sometimes (but not always) insert an "N" in positions where the base pair is difficult to call.  But sometimes it makes an incorrect call as well.  

