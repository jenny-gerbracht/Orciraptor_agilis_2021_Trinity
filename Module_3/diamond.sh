#!/bin/sh

###################################################################
#Script Name	:diamond.sh		                                                                                              
#Description	:diamond search to remove redundant ORFs                                                                               
#Args:                                                                                           
#Author	:Jennifer Gerbracht                                               
#Email		:jennifer.gerbracht@gmx.de                                           
###################################################################

source ../config.txt
moduledir="${mydir}/Module_3"

####################################
#
# Setting up log file
#
###################################
# Get date
date=$(date "+%Y-%m-%d")

# Define log file and redirect stdout and stderr to this file
if [ ! -d "${moduledir}/Logs/" ]; then
  mytime=$(date "+%Y-%m-%d %H:%M:%S")
  echo "$mytime Make directory ${moduledir}/Logs/"
  mkdir ${moduledir}/Logs/
fi
log_file="${moduledir}/Logs/log_diamond_$date"
exec &> >(tee -a "$log_file")

####################################
#
# Diamond search
#
####################################

cd transdecoder

diamond makedb \
--in Trinity_filtered.fasta.transdecoder.pep_renamed.fas \
-d Trinity_filtered.fasta.transdecoder.pep_renamed.fas

diamond blastp \
--threads 8 \
-d Trinity_filtered.fasta.transdecoder.pep_renamed.fas.dmnd \
-q Trinity_filtered.fasta.transdecoder.pep_renamed.fas \
-a ORFsvsORFs.daa \
--max-target-seqs 2000 \
--evalue 0.00001

diamond view \
-o ORFsvsORFs.tsv \
-f 6 qseqid sseqid pident length mismatch gapopen qstart qend qlen sstart send slen evalue \
-a ORFsvsORFs.daa
