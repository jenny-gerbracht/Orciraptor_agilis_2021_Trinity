###################################################################
#Script Name	:ParseORFsVSORFsblastp.py		                                                                                              
#Description	:This script removes redundant sequences in predicted ORFs after Trinity assembly if percent identity > 95% and ALNcov > 0.9 of smallest  
#	         sequence.                                                                            
#Args		:blastp output and fasta with ORFs
#Usage		:python ParseORFsVSORFsblastp.py <blastp> <fasta>                                                                                 
#Author	:Tommy Harding                                               
#Email		:                                          
###################################################################

import sys

blastp = sys.argv[1]
fasta = sys.argv[2]


print "Parsing BLASTP output..."
BLASTP = open(blastp)
redundant = []
qlist =[]
l = 0
while 1:
	line = BLASTP.readline()
	if not line:
		break
	q = line.split("\t")[0].strip()
	if q not in qlist:
		qlist.append(q)
		l += 1
		if l % 1000 == 0:
			print l, "queries parsed"
	s = line.split("\t")[1].strip()
	perc_ID = float(line.split("\t")[2])
	ALNlength = float(line.split("\t")[3])
	mismatch = int(line.split("\t")[4])
	query_length = int(line.split("\t")[8])
	subject_length = int(line.split("\t")[11])
	ALNcov = ALNlength/min(query_length,subject_length)
	if q != s and perc_ID > 0.95 and ALNcov > 0.9:
		if query_length == subject_length:
			if q not in redundant and s not in redundant:
				redundant.append(s)
		elif query_length > subject_length:
			if s not in redundant:
				redundant.append(s)
		elif subject_length > query_length:
			if q not in redundant:
				redundant.append(q)
BLASTP.close()


print "Extracting non-redundant sequences..."
FASTA = open(fasta)
out = open(fasta[:-4]+"_non-redundant.fas", 'w')
while 1:
	line = FASTA.readline()
	if not line:
		break
	if line.startswith(">"):
		seqname = line.split()[0].strip(">")
		if seqname not in redundant:
			out.write(line)
			line = FASTA.readline()
			while not line.startswith(">"):
				out.write(line)
				x = FASTA.tell()
				line = FASTA.readline()
				if not line:
					break				
			FASTA.seek(x)
FASTA.close()
out.close()

print len(redundant), "discarded sequences"
