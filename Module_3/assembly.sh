#!/bin/bash

###################################################################
#Script Name	:assembly.sh		                                                                                              
#Description	:Performs de novo transcriptome assembly for Orciraptor agilis                                                                               
#Args:                                                                                           
#Author	:Jennifer Gerbracht                                               
#Email		:jennifer.gerbracht@gmx.de                                           
###################################################################

source ../config.txt
moduledir="${mydir}/Module_3"
mydir="/srv/Jenny/Datasets/orciraptor/processed_github"

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
  echo "${mytime} Make directory ${moduledir}/Logs/"
  mkdir ${moduledir}/Logs/
fi
log_file="${moduledir}/Logs/log_assembly_$date"
exec &> >(tee -a "$log_file")

# Assembly with Trinity v. 2.0.6 and generate assembly summary statistics
echo ""
echo "###################"
echo "## Trinity assembly"
echo "###################"
echo ""
echo -n "Trinity version: "
Trinity --version
echo ""

Trinity \
--seqType fq \
--max_memory 50G \
--left ${mydir}/blacklist_NA_paired_unaligned_V1S1.1.fq.gz,\
${mydir}/blacklist_NA_paired_unaligned_V1S2.1.fq.gz,\
${mydir}/blacklist_NA_paired_unaligned_V1S3.1.fq.gz,\
${mydir}/blacklist_NA_paired_unaligned_V1S4.1.fq.gz,\
${mydir}/blacklist_NA_paired_unaligned_V1S5.1.fq.gz,\
${mydir}/blacklist_NA_paired_unaligned_V1S6.1.fq.gz,\
${mydir}/blacklist_NA_paired_unaligned_V1S7.1.fq.gz,\
${mydir}/blacklist_NA_paired_unaligned_V1S8.1.fq.gz,\
${mydir}/blacklist_NA_paired_unaligned_V1S9.1.fq.gz \
--right ${mydir}/blacklist_NA_paired_unaligned_V1S1.2.fq.gz,\
${mydir}/blacklist_NA_paired_unaligned_V1S2.2.fq.gz,\
${mydir}/blacklist_NA_paired_unaligned_V1S3.2.fq.gz,\
${mydir}/blacklist_NA_paired_unaligned_V1S4.2.fq.gz,\
${mydir}/blacklist_NA_paired_unaligned_V1S5.2.fq.gz,\
${mydir}/blacklist_NA_paired_unaligned_V1S6.2.fq.gz,\
${mydir}/blacklist_NA_paired_unaligned_V1S7.2.fq.gz,\
${mydir}/blacklist_NA_paired_unaligned_V1S8.2.fq.gz,\
${mydir}/blacklist_NA_paired_unaligned_V1S9.2.fq.gz \
--CPU 8 \
--normalize_reads  \
--SS_lib_type RF
