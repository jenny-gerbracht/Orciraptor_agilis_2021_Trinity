###################################################################
#Script Name	:AssignRandomSeqnames.py		                                                                                              
#Description	:This script assigns random names to sequences.                                                       
#Args		: Fasta file with input sequences for alignment
#Usage		: python AssignRandomSeqnames.py <fasta> file>                                                                                          
#Author	:Tommy Harding                                               
#Email		:                                        
###################################################################

import sys, random, string

Fasta = sys.argv[1]


out = open("%s_renamed.fas" % Fasta.split(".fa")[0], 'w')
LOG = open("%s_randomNames.txt" % Fasta.split(".fa")[0], 'w')
LOG.write("Assigned random name\tOriginal seqname\n")
randomNamelist = []
FILE = open(Fasta)
while 1:
	line = FILE.readline()
	if not line:
		break
	if line.startswith(">"):
		r = ''.join(random.choice(string.ascii_lowercase + string.digits) for _ in range(8))
		while r in randomNamelist:
			r = ''.join(random.choice(string.ascii_lowercase + string.digits) for _ in range(8))
		randomNamelist.append(r)
		out.write(">%s\n" % r)
		LOG.write("%s\t%s\n" % (r, line.strip(">").strip("\n")))
	else:
		out.write(line)
out.close()
LOG.close()
