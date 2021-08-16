#!/bin/bash

###################################################################
#Script Name	:seqkit.sh		                                                                                              
#Description	:Removes contigs from fasta that were identified by the blastn search as contamination                                                                              
#Args:                                                                                           
#Author	:Jennifer Gerbracht                                               
#Email		:jennifer.gerbracht@gmx.de                                           
###################################################################

source ../config.txt
moduledir="${mydir}/Module_3"

seqkit \
grep \
-f contaminants_orciraptor.txt \
${mydir}/Module_3/trinity_out_dir/Trinity.fasta \
-o Trinity_filtered.fasta \
-v
