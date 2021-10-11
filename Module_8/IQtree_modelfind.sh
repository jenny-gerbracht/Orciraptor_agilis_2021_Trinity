#!/bin/bash
#$ -S /bin/bash
. /etc/profile
#$ -cwd
#$ -pe threaded 10

###################################################################
#Script Name	:IQtree.sh		                                                                                              
#Description	:Find best model for IQtree                                                                              
#Args:                                                                                           
#Author	:Jennifer Gerbracht                                               
#Email		:jennifer.gerbracht@gmx.de                                           
###################################################################

source activate iqtree
iqtree \
-s mafft_alignment_trimal_derep_defrag.phy \
-m MFP \
-nt AUTO \
-v
conda deactivate
