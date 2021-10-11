# Orciraptor_agilis_2021_Trinity
Read processing and filtering, *de novo* transcriptome assembly, differential gene expression analysis and functional annotation of *Orciraptor agilis* RNA-seq data using Trinity

## Module 1: Read processing and *de novo* transcriptome assembly of prey organism *Mougeotia* sp.

1. Run readprocessing_and_assembly.sh, output is *de novo* transcriptome assembly of *Mougeotia* sp.
2. Predict ORFs to use later for decontamination: transdecoder.sh
3. Filter transcriptome for contigs smaller than 200 nt: seqkit_length.sh (only for upload)

## Module 2: Read processing of *Orciraptor agilis*

1. Run symlinks.sh
2. Run readprocessing.sh. Output are quality filtered and adapter trimmed reads.
3. Run mapping.sh. Output are reads that do not map to sequences from rRNA and/or *Mougeotia* sp.

## Module 3: *De novo* transcriptome assembly, decontamination, ORF prediction
1. Run assembly.sh to assemble the transcriptome from processed reads of all libraries. Output is *de novo* transcriptome assembly of *Orciraptor agilis* as a fasta.
2. Run blastn search with this transcriptome (nt database v5 updated on 2021-03-10): blastn.sh
  * Checked contigs with > 95% identity over a length of minimum 100 nt, saved contig identifiers of all bacterial, viral, ribosomal and algal contigs in contaminants.txt
  * Remove these sequences from transcriptome with seqkit.sh
4. ORF prediction with transdecoder.sh
5. Use Change_Seqname_TransDecoder.py to change names of the transcdecoder output
6. Perform a diamond blastp search with diamond.sh comparing all ORFs against each other
7. Run ParseORFsVSORFsblastp.py on the diamond output and the renamed transcdecoder file to
obtain a non-redundant .pep file, rename output to "Orciraptor_non-redundant.faa"
8. Use Change_Seqname_TransDecoder.py to change the names of the corresponding coding sequences
9. Call Extract_seq_using_FASTA.py to obtain the coding sequences of the non-redundant .pep file

## Module 4: Functional annotation
1. Run eggnog-mapper in diamond and hmm mode: eggnog.sh. 
2. Run InterProScan using interproscan.sh
3. Run a Diamond blastp search vs Swiss-Prot database (UniProt Release 2021_01): diamond.sh
4. Annotation of carbohydrate-active enzymes (CAZymes) with with dbcan2 in HMM mode (database dbCAN-HMMdb-V9): cazy.sh

## Module 5: Differential gene expression analysis
1. Mapping the processed reads back to the newly generated transcriptome with bowtie2 and counting with salmon in alignment-mode (bowtie2.sh).
2. Run DESeq2

## Module 6: Assembly summary statistics
1) number of Genes, isoforms, ORFs (status), Ex90 

## Module 8: GH5_5 phylogeny
1. Run AssignRandomSeqnames.py on input sequences to generate random names
2. Simplify tip labels and generate annotation file for Figtree with Renaming_and_annoration.R
3. Run alignment and trimming ?.sh
4. Remove identical sequences using DereplicateALN.py
5. Find best model with IQtree_modelfind.sh
6. Generate tree with IQtree.sh
7. Run RenameTrees.py to get tip labels
