#!/bin/bash

#$ -P bf528
#$ -cwd

source /etc/bashrc

module load python/1.6
module load python2
module load multiqc/1.6
module load utils/4.0.2
module spider utils


multiqc #/projectnb/bf528/users/hedgehog/project_3/project-3-hedgehog/Programmer/countsfile/counts_only

echo $1

 setwd("/projectnb2/bf528/users/hedgehog/project_3/project-3-hedgehog/Programmer/countsfile/counts_only")
read.csv(SRR1177966.txt, header = TRUE, row.names = Geneid)
 -o /projectnb/bf528/users/hedgehog/project_3/project-3-hedgehog/Programmer
 
test1 == c(1:5, "6,7", "8,9,10")
tf == tempfile()
writeLines(test1, tf)