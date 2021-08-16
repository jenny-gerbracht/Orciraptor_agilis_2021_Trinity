###################################################################
#Script Name	:DereplicateALN.py		                                                                                              
#Description	:This script removes identical sequences from trimmed alignment.                                                                    
#Args		:alignment file (phy) trimmed by BGME
#Usage		:python DereplicateALN.py <alignment file> <tree>                                                                             
#Author	:Tommy Harding                                               
#Email		:                                          
###################################################################

import sys, re


ALN = sys.argv[1]

FILE = open(ALN)
seqLIST = []
ALNlist = []
FirstLine = FILE.readline()
ditchedSEQ = 0
while 1:
	line = FILE.readline()
	if not line:
		break
	name = line.split()[0]
	seq = line.split()[1]
	seqCLEAN = re.sub('-','',seq)
	if seqCLEAN not in seqLIST:
		seqLIST.append(seqCLEAN)
		info = (name,seq)
		ALNlist.append(info)
	else:
		ditchedSEQ += 1
FILE.close()

out = open("%sDerep.phy" % (ALN[:-4]), 'w')
Ntaxa = int(FirstLine.split()[0])-ditchedSEQ
out.write(" %s %s\n" % (Ntaxa,FirstLine.split()[1]))
maxLIST = []
for info in ALNlist:
	maxLIST.append(len(info[0]))
maxLEN = max(maxLIST)
for info in ALNlist:
	extraspace = (maxLEN - len(info[0]) + 1) * " "
	out.write("%s%s%s\n" % (info[0],extraspace,info[1]))
out.close()	
