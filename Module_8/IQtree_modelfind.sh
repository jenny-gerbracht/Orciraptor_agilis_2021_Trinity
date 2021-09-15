#!/bin/bash
#$ -S /bin/bash
. /etc/profile
#$ -cwd
#$ -pe threaded 10

source activate iqtree
iqtree \
-s mafft_alignment_trimal_derep_defrag.phy \
-m MFP \
-nt AUTO \
-v
conda deactivate
