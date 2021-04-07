#!/bin/bash

#$ -P bf528
#$ -cwd
#$ -pe omp 16

source /etc/bashrc

module load subread

### To get the right version of subread###
module spider subread 

module load subread/1.6.2

###Location of GTF files###
GTF=/project/bf528/project_3/reference/rn4_refGene_20180308.gtf

###Where to send countsfile to after running featureCounts###
OUT=/projectnb/bf528/users/hedgehog/project_3/project-3-hedgehog/Programmer/countsfile/allcounts.txt



###Location of Hedgehog BAM files##
BAMFILE=/projectnb/bf528/users/hedgehog/project_3/data_curator_data/bam_files/*Aligned.sortedByCoord.out.bam

featureCounts -T 16 -a $GTF  -o $OUT $BAMFILE