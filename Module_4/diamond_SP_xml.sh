#!/bin/bash

###################################################################
#Script Name	:diamond_SP_xml.sh		                                                                                              
#Description	:Run blastp search vs. SwissProt database                                                       
#Args:                                                                                           
#Author	:Jennifer Gerbracht                                               
#Email		:jennifer.gerbracht@gmx.de                                           
###################################################################

# Define paths to working directory locations
source ../config.txt
moduledir="${mydir}/Module_4"

# Log file
# Get date
date=$(date "+%Y-%m-%d")

# Define log file and redirect stdout and stderr to this file
if [ ! -d "${moduledir}/Logs/" ]; then
  mytime=$(date "+%Y-%m-%d %H:%M:%S")
  echo "${mytime} Make directory ${moduledir}/Logs/"
  mkdir ${moduledir}/Logs/
fi
log_file="${moduledir}/Logs/log_diamond_$date"
exec &> >(tee -a "$log_file")

# Create necessary folders
if [ ! -d "${moduledir}/diamond/" ]; then
  mytime=$(date "+%Y-%m-%d %H:%M:%S")
  echo "${mytime} Make directory ${moduledir}/diamond/"
  mkdir ${moduledir}/diamond/
fi

echo ""
echo "###################"
echo "## Diamond SP search"
echo "###################"
echo ""

diamond \
blastp \
--query ${mydir}/Module_3/transdecoder/Orciraptor_non-redundant.faa \
--threads 10 \
--db /srv/Jenny/entap/bin/uniprot_sprot \
--out ${moduledir}/diamond/diamond_SP.xml \
--outfmt 5 \
--max-target-seqs 1
