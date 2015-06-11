# De-multiplexing Illumina Data

(or you can go back to the Basic Unix page [here](https://github.com/evansbenj/Reduced-Representation-Workshop/blob/master/2_BasicUnix.md)).

Most RRGS methods rely on the Illumina sequencing platform.  These machines generate data using something called a "flowcell" that is divided up into eight "lanes".  Small scale projects typically would run multiple samples (from different species or different individuals within a species) on one lane.  

## Fasta and Fastq format

Illumina sequence data is provided in a text file that is in a format called `fastq`.  This is a modification of another format called `fasta` in which each sequence has a header that begins with a `>` sign.  This is followed by the sequences.  Here is an example:

```
>example_sequence_in_fasta_format`
ATGCGCGCGCTAGGCTCGCGATCGGGGAGCGCGAGCTGAGCTAGCGCGATGCGCCCCGAC
```

The format of `fastq` files is similar to `fasta` except that quality scores are included.  Each sequence has four lines (instead of two for `fasta` files).  The first begins with `@` followed by information about the sequence.  The second line is the nucleotide sequence. The third line is a `+` which may be followed by the same information that followed the `@` sign in the first line.  The fourth line is the quality values.  For the Illumina data we will be working with, these values range from 0–41 and are represented by single characters.  More details about fastq format is available [here](http://en.wikipedia.org/wiki/FASTQ_format).  Here is an example of a sequence in fastq format:

```
@HWI-ST724:202:D127MACXX:6:2114:13665:74490
TGCAGGAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAGGAAATTTTTGGGCACAAAGAACCACAGAAAAAAAATGAAAA
+HWI-ST724:202:D127MACXX:6:2114:13665:74490
AFHJJJFIJJJJIJJJJJHFDDDDDB0530&0)00&)0&05007BDD############################################
```

In this sequence the number signs indicate low quality reads at the end (right side) of the sequence.

## De-Multiplexing

A first step in our analysis pipeline is to organize data from each of our samples that were run together on an Illumina lane (De-multiplexing our data) and also to filter our data and trim off bits that have lots of errors or that have sequences from the laboratory procedures that were used to generate the data (Trimming/Quality control).  To begin please make sure you are in your home directory by typing this:

`cd ~`

When samples are run on an Illumina machine, DNA is broken up into many small fragments and a small bit of DNA called an adaptor is then added on each of the fragments. This adaptor allows the sequencing process to occur, essentially by making possible high-throughput put polymerase chain reaction (ask Ben about this if you are unfamiliar). To make possible the multiplexing of samples on one Illumina lane, each sample is linked to a unique adaptor that contains a "barcode" sequence that allows us to sort out which samples each sequence came from.  For our dataset, we have nine individuals from one species (the Tonkean macaque). Each of the samples received the following barcodes (the sample name is followed by the barcode):
```
PF515 CCTCTTATCA
PM561 TATCGTTAGT
PM565 TAGTGCGGTC
PM566 GGCCGGTAAC
PM567 AGGAACCTCG
PM582 TTATCCGTAG
PM584 CGCTATACGG
PM592 CACGCAACGA
PM602 ATCCGTCTAC
```
We will use this information in a moment to de-multiplex our data.  Lets use emacs to make a text file that contains this information.  Please type `emacs monkey.pools` to generate a file called `monkey.pools`.  Now copy and paste the information above to your emacs window.  Then type `Ctrl-X` and `Ctrl-S` to save it and then `Ctrl-X` and `Ctrl-C` to close the program.


## Quality control and trimming

A first step in analysis of Illumina data is to identify adaptor and barcode sequences in our data, sort sequneces by the barcode, and then trim off the adaptor and barcode sequences.  We can also get rid of sequences that have ambiguous barcodes due to sequencing errors.

Illumina generates sequences that have errors in base calls.  Errors typically become more common towards the end of the sequence read, and sometimes (but not always) an "N" is inserted in positions where the base pair is difficult to call.  But sometimes it makes an incorrect call as well. 

We will use software called `RADpools` to de-multiplex and trim our data.  This software is available [here](https://github.com/johnomics/RADtools/blob/master/RADpools).

To use `RADpools` we need to load some perl modules first.  Please type this:

`module load app/radtools`

The command to execute this program on our data is:

`/apps/RADtools/1.2.4/RADpools -i /home/datasets/2015_Ben_Evans/data -d ./monkey -s -f -o -t 55 -q 10`

The first part (`/apps/RADpools`) directs the computer to run the program RADpools, which is in the director called `apps`.  The `-i` flag specifies where the data are.  The `-d` flag specifies where the barcode file is that we made eariler.  The `-s` flag tells RADpools that the data have Sanger quality scores.  The `-f` flag tells RADpools to interpret the barcodes using the "fuzzy" option, which allows for errors and asigns barcodes with errors to the nearest pool. The `-o` flag directs RADpools to output the trimmed data in fastq format.  The `-t` flag tells RADpools to trim all of the sequences to a length of 55 base pairs.  The `q` flag tells RADpools to throw out any sequences that have any bases with a quality score below 10.  When the program is done sorting the data, it should generate a directory called `monkey` in your current directory. (The name is just whatever the suffix is of your `.pools` file.)  All of this information is available, of course, in the manual that comes with the program. 

## OK, now we are ready to move on to mapping reads to a reference genome.  Please click [here](https://github.com/evansbenj/Reduced-Representation-Workshop/blob/master/4_Mapping_reads_to_a_reference_genome.md).
