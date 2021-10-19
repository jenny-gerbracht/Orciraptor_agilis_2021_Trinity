###################################################################
#Script Name	:Extract_seq_using_FASTA.py		                                                                                              
#Description	:Take a fasta file and extract these sequences from another fasta file.                                                                          
#Args		:fasta file with headers of interest, fasta file with sequences to be
# 		 extracted and name of output file
#Usage		:python Extract_seq_using_FASTA.py <fasta file1> <fasta file2> <outfile name>                                                                                 
#Author	:Tommy Harding                                               
#Email		:                                          
###################################################################

import sys

LIST = sys.argv[1]
Fasta_file = sys.argv[2]
outfile_name = sys.argv[3]

print "Building seqList..."
List = open(LIST)
Seq_list = []
while 1:
	line = List.readline()
	if not line:
		break
	if line.startswith(">"):
		seqname = line.split("_")[0].strip(">")
		if seqname not in Seq_list:
			Seq_list.append(seqname)
List.close()
print len(Seq_list), 'sequences to extract'


print "Extracting sequences..."
Fasta_infile = open(Fasta_file)
lines = Fasta_infile.readlines()
Fasta_infile.close()
outfile = open(outfile_name, 'w')
n=0
for l in range(len(lines)-1):
	compline = lines[l]
	if compline[0] == ">":
		header_line = compline.split("_")[0].strip(">").strip()
		if header_line in Seq_list:
			outfile.write(compline)
			m = 1
			while not (lines[l+m][0] == ">"):
				if l+m == len(lines)-1:
					next = lines[l+m]
					outfile.write(next)
					m = m+1
					break
				else:
					next = lines[l+m]
					outfile.write(next)
					m = m+1
outfile.close()

