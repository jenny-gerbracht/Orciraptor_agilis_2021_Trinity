source activate mafft
mafft-linsi mafft_input_renamed.fas > mafft_alignment
conda deactivate

trimal -in mafft_alignment -out mafft_alignment_trimal.phy -phylip -automated1
