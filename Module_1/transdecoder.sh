#!/bin/bash
#$ -S /bin/bash
. /etc/profile
#$ -cwd

###################################################################
#Script Name	:transdecoder.sh		                                                                                              
#Description	:Performs ORF prediction from de novo assembly of Mougeotia sp.                                                                               
#Args:                                                                                           
#Author	:Jennifer Gerbracht                                               
#Email		:jennifer.gerbracht@gmx.de                                           
###################################################################

source ../config.txt
moduledir="${mydir}/Module_1"

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
log_file="${moduledir}/Logs/log_transdecoder_$date"
exec &> >(tee -a "$log_file")

####################################
#
# Make folders
#
####################################

if [ ! -d "${moduledir}/transdecoder/" ]; then
  mytime=$(date "+%Y-%m-%d %H:%M:%S")
  echo "$mytime Make directory ${moduledir}/transdecoder/"
  mkdir ${moduledir}/transdecoder/
fi

cd transdecoder

source activate transdecoder

TransDecoder.LongOrfs -t ${moduledir}/NA_transcripts.fasta
TransDecoder.Predict -t ${moduledir}/NA_transcripts.fasta

source deactivate
