#!/bin/bash
#$ -S /bin/bash
. /etc/profile
#$ -cwd
#$ -pe threaded 10

###################################################################
#Script Name	:IQtree.sh		                                                                                              
#Description	:Construct phylogenetic tree with IQtree                                                                              
#Args:                                                                                           
#Author	:Jennifer Gerbracht                                               
#Email		:jennifer.gerbracht@gmx.de                                           
###################################################################

source activate iqtree

iqtree \
-s mafft_alignment_trimal_derep_defrag.phy \
-m Q.pfam+R5 \
-nt AUTO \
-v \
-bb 1000

conda deactivate
