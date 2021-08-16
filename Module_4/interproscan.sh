#!/bin/bash

###################################################################
#Script Name	:interproscan.sh		                                                                                              
#Description	:Predict functional domains with InterProScan                                                                   
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
log_file="${moduledir}/Logs/log_interproscan_$date"
exec &> >(tee -a "$log_file")

# Create necessary folders
if [ ! -d "${moduledir}/interproscan/" ]; then
  mytime=$(date "+%Y-%m-%d %H:%M:%S")
  echo "${mytime} Make directory ${moduledir}/interproscan/"
  mkdir ${moduledir}/interproscan/
fi

echo ""
echo "###################"
echo "## InterProScan"
echo "###################"
echo ""

# Remove asterisks from protein fasta

sed 's/*//g' ${mydir}/Module_3/transdecoder/Orciraptor_non-redundant.faa \
> ${mydir}/Module_3/transdecoder/Orciraptor_non-redundant.faa_nostop.fasta

${interproscan_dir}/interproscan.sh \
--input ${mydir}/Module_3/transdecoder/Orciraptor_non-redundant.faa_nostop.fasta \
--iprlookup -pa \
--output-dir ${moduledir}/interproscan \
--goterms \
--cpu 8
