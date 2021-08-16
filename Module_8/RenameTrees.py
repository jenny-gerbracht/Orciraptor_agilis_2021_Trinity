###################################################################
#Script Name	:RenameTrees.py		                                                                                              
#Description	:Replace random names from tree with original ones                                                                    
#Args		:List of mapping random to original names, treefile
#Usage		:python RenameTrees.py <namelist> <tree>                                                                             
#Author	:Tommy Harding                                               
#Email		:                                          
###################################################################

from Bio import Phylo
import glob, re, sys

nameList = sys.argv[1]
tree = sys.argv[2]


def RemoveProhibitedSymbols(title):
	ProhibitedSymbols = ["\(","\)",":",";",","]
	for symbol in ProhibitedSymbols:
		newtitle = re.sub(symbol,'-',title)
		title = newtitle
	return title

##################################
dict = {}
FILE = open(nameList)
line = FILE.readline()
while 1:
	line = FILE.readline()
	if not line:
		break
	description = RemoveProhibitedSymbols(line.split("\t")[1].strip())
	dict[line.split()[0].strip()] = [description]
FILE.close()

	
TREE = open(tree)
lines = TREE.readlines()
line = lines[0]
FILE.close()
TREE = Phylo.read(tree, "newick")
tips = set([str(t) for t in TREE.get_terminals()])
for seqID in tips:
	line = re.sub(seqID,dict[seqID][0],line)

out = open("%s_Renamed.tre" % tree[:-4], 'w')
out.write(line)
out.close()

