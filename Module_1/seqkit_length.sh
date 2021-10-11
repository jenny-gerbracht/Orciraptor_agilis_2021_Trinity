#!/bin/bash

###################################################################
#Script Name	:seqkit_length.sh		                                                                                              
#Description	:Removes contigs from fasta that are smaller than 200 nt (for upload to TSA)                                                                       
#Args:                                                                                           
#Author	:Jennifer Gerbracht                                               
#Email		:jennifer.gerbracht@gmx.de                                           
###################################################################

source ../config.txt
moduledir="${mydir}/Module_1"

seqkit \
seq \
-m 200 \
${moduledir}/NA_transcripts.fasta > ${moduledir}/NA_transcripts_200.fasta
