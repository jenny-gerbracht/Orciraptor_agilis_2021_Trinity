###################################################################
#Script Name	:Change_Seqname_TransDecoder.py		                                                                                              
#Description	:Renames sequences from the TransDecoder output (e.g. version 2.1.0).                                                                            
#Args		:.pep output from TransDecoder
#Usage		:python Change_Seqname_TransDecoder.py <TransDecoder.pep>                                                                                 
#Author	:Tommy Harding                                               
#Email		:                                          
###################################################################

import sys

Fasta_file = sys.argv[1]


outfile=open(Fasta_file + "_renamed.fas", "w")


Fasta_infile = open(Fasta_file)

lines = Fasta_infile.readlines()
for l in range(len(lines)-1):
	line = lines[l]
	#print compline
	if line[0] == ">":
		new_seq_name = ("%s_type-%s_len%s_%s_%s" % (line.split()[0].split("|")[2].strip(), line.split()[5].split(":")[1].strip(), line.split()[6].split(":")[1].strip(), line.split()[8].split("|")[0], line.split()[8].split("|")[1].split(":")[0]))
		outfile.write(">%s\n"%(new_seq_name))
		m = 1
		while not (lines[l+m][0] == ">"):
			if (l+m) == (len(lines)-1):
				line = lines[l+m]
				outfile.write(line)
				break
			else:
				next = lines[l+m]
				outfile.write(next)
				m = m+1


Fasta_infile.close()
outfile.close()
