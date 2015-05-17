# Automating an alignment with a bash script

Now that you have seen how to align data from one individual to a reference genome, we can automate the alignment of all individuals to the reference genome using a bash script.  We can accomplish this by defining an `array` that contains the names of all of the individuals in the analysis and then looping through this array and executing each of the commands for each individual.

Here is a bash script that should accomplish this:

```
#!/bin/bash                                                                              

path_to_bwa="/work/ben/workshop_software/bwa-0.7.12"
path_to_samtools="/work/ben/workshop_software/samtools-1.2"
path_to_data="/work/ben/2015_workshop/data/monkey"
path_to_chromosome="/work/ben/2015_workshop/rhesus_chromosomes"
chromosome="chr3.fa"

individuals="PF515                                                                       
PM561                                                                                    
PM565                                                                                    
PM566                                                                                    
PM567                                                                                    
PM582                                                                                    
PM584                                                                                    
PM592                                                                                    
PM602"

for each_individual in $individuals
do

    echo ${each_individual}
    $path_to_bwa/bwa aln $path_to_chromosome/$chromosome $path_to_data/${each_individual\
}.fastq > $path_to_data/${each_individual}.sai
    $path_to_bwa/bwa samse -r "@RG\tID:FLOWCELL1.LANE6\tSM:${each_individual}.fastq\tPL:\
illumina" $path_to_chromosome/$chromosome $path_to_data/${each_individual}.sai $path_to_\
data/${each_individual}.fastq > $path_to_data/${each_individual}.sam
    $path_to_samtools/samtools view -bt $path_to_chromosome/$chromosome -o $path_to_data\
/${each_individual}.bam $path_to_data/${each_individual}.sam
    $path_to_samtools/samtools sort $path_to_data/${each_individual}.bam $path_to_data/$\
{each_individual}_sorted

done

```

In the beginning of the script 5 variables are defined that specify, respectively, the path for the bwa and samtools programs, the path to the data, the path to the reference chromosome, and the name of the chromosome you are working on.  Please copy this section and then make a new file called `alignment_commando` using emacs by typing this:

`emacs alignment_commando`

This should open up an emacs window.  You can then paste the bash script into this window.  Now use the arrow keys to scroll up to the line that says `chromosome="chr3.fa"` and change the part that says `chr3.fa` to correspond with whatever chromosome you are working on.  For example, if youa re working on chromosome 10, please change this to instead read `chr10.fa`.

Now type `Ctrl-x` and `Ctrl-s` to save the file.

Now we need to make the file executable, so type this:

`chmod +x alignment commando`

And now we should be able to execute the file.  Type this:

`./alignment_commando`

