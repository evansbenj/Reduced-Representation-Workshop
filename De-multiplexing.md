# De-multiplexing Illumina Data

Most RRGS methods rely on the Illumina sequencing platform.  These machines generate data using something called a "flowcell" that is divided up into eight "lanes".  Small scale projects typically would run multiple samples (from different species or different individuals within a species) on one lane.  
## De-Multiplexing

## Fasta and Fastq format

## Quality control and trimming

Similar to Sanger sequencing, Illumina generates sequences that have errors.  Errors typically become more common towards the end of the sequence read, and the software will sometimes (but not always) insert an "N" in positions where the base pair is difficult to call.  But sometimes it makes an incorrect call as well.  
