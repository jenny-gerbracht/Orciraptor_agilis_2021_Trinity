#!/bin/bash

###################################################################
#Script Name	:cazy.sh		                                                                                              
#Description	:Predict carbohydrate-active enzymes                                                                     
#Args:                                                                                           
#Author	:Jennifer Gerbracht                                               
#Email		:jennifer.gerbracht@gmx.de                                           
###################################################################

# Define paths to working directory locations
source ../../config.txt
moduledir="${mydir}/Module_4"

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
log_file="${moduledir}/Logs/log_cazy_$date"
exec &> >(tee -a "$log_file")

# Annotate carbohydrate-active enzymes with dbcan2 in HMM mode (http://bcb.unl.edu/dbCAN2/download/Tools/)
echo ""
echo "###################"
echo "## dbcan2"
echo "###################"
echo ""

python3 ${Cazy_dir}/run_dbcan.py \
${mydir}/Module_3/transdecoder/Orciraptor_non-redundant.faa \
protein \
--db_dir ${Cazy_dir}/db \
--dbCANFile dbCAN-HMMdb-V9.txt \
--hmm_eval 1e-5 \
--hmm_cpu 10 \
--out_dir cazy \
--tools hmmer 
