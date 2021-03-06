#!/bin/bash

###################################################################
#Script Name	:blastn.sh		                                                                                              
#Description	:Runs blast search against nt database, returns one hit
#		 per contig with taxonomic assignment                                                                               
#Args:                                                                                           
#Author	:Jennifer Gerbracht                                               
#Email		:jennifer.gerbracht@gmx.de                                           
###################################################################

blastn -db /srv/Jenny/NCBI/nt/nt \
       -query ./trinity_out_dir/Trinity.fasta \
       -outfmt '6 qseqid sseqid pident length evalue bitscore sgi sacc staxids sskingdoms sscinames scomnames stitle' \
       -max_target_seqs 1 \
       -max_hsps 1 \
       -evalue 1e-25 \
       -num_threads 4 \
       -out Module_3_tax_1hit.hits
