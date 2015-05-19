# Reduced-Representation-Workshop
This is a repository for a one day bioinformatics workshop on reduced representation genome sequencing

## Background
Many interesting organisms have big genomes, making complete genome sequencing infeasible, especially for multiple individuals.  A relatively cost-efficient solution has been developed recently called "reduced representation genome sequencing" (RRGS).  This approach enables deep sequencing of multiple genetic samples from the same (homologous) genomic regions.  It takes advantage of next generation (Illumina) sequencing technology and requires relatively simple laboratory efforts such as DNA extraction, whcih can be accomplished with a centrifuge and a heat block.  Other laboratory steps (library construction) can be outsourced or done in house, depending on the equipment and funds that are available. This approach has many applications, including phylogenomics, population genomics, linkage mapping, and analysis of gene flow.

## Goals
The goal of this workshop is to introduce students to RRGS, including how various approaches differ and advantages and disadvantages of differnet approaches. Ideally this workshop will provide a sufficient level of exposure to students that they will be able to know how to learn more and generate and analyze their own datasets.

## Specific Topics
 This workshop will be split into a morning and an afternoon section. The morning section will begin with an introduction to reduced representation genome sequencing (e.g. RADseq, GBS) and explore examples from my (Ben Evans, McMaster University, Canada) research that utilize this datatype.  The afternoon session will be interactive, and will include (a) basic Unix use, including file manipulation and text editing with Emacs, (b) quality scores, file types, and de-multiplexing of multiplexed samples, (c) aligning reduced representation data to a reference genome, (a) assembling reduced representation data without a reference genome, (d) genotype calling, and (e) applications to phylogenomics, population genomics, and linkage analysis.  

## Interactive afternoon session 
 The afternoon session will consist of a group project that is based on a RADseq dataset from the Tonkean macaque monkey (*Macaca tonkeana*).  Tonkean macaques inhabit the central Indonesian island of Sulawesi and, like other papionin monkeys, have a social system characterized by strong female philopatry and obligate male migration.  Reproductive success is thought to be more variable among males than females.  If this is true, we  expect  molecular polymorphism on the X chromosome to be elevated relative to an expectation with equal variance in reproductive success among the sexes (sounds complicated, but Ben will explain this). To test this hypothesis, in the afternoon session, each student will map RADseq data from the Tonkean macaque to one of the chromosomes of a closely related species – the rhesus macaque (*Macaca mulatta*), whose genome has been completely sequenced.  We will then evaluate diversity and polymorphism of the sex chromosomes and the autosomes of the Tonkean macaque and explore what this can tell us about the social system of these fantastic monkeys.  In doing so, my hope is that students will gain an appreciation for each of the steps of analysis of RRGS data.
 
## OK, lets begin by clicking [here](https://github.com/evansbenj/Reduced-Representation-Workshop/blob/master/2_BasicUnix.md).