###################################################################
#Script Name	:GetORFstatus.py		                                                                                              
#Description	:Parse headers of TransDecoder output for ORF types                                                                          
#Args		:renamed TransDecoder output
#Usage		:python GetORFstatus.py <fasta file>                                                                             
#Author	:Tommy Harding                                               
#Email		:                                          
###################################################################

import sys


Fasta = sys.argv[1]


complete = 0
Fiveprime = 0
Threeprime = 0 
internal = 0

FASTA = open(Fasta)
while 1:
	line = FASTA.readline()
	if not line:
		break
	if line.startswith(">"):
		ORFstatus = line.split("type-")[1].split("_len")[0]
		if ORFstatus == "complete":
			complete += 1
		elif ORFstatus == "5prime_partial":
			Fiveprime += 1
		elif ORFstatus == "3prime_partial":
			Threeprime += 1
		elif ORFstatus == "internal":
			internal += 1
FASTA.close()

print complete, "complete ORFs (", float(complete)/(complete+Fiveprime+Threeprime+internal)*100, "%)"
print Fiveprime, "5'-partial ORFs (", float(Fiveprime)/(complete+Fiveprime+Threeprime+internal)*100, "%)"
print Threeprime, "3'-partial ORFs (", float(Threeprime)/(complete+Fiveprime+Threeprime+internal)*100, "%)"
print internal, "internal ORFs (", float(internal)/(complete+Fiveprime+Threeprime+internal)*100, "%)"
