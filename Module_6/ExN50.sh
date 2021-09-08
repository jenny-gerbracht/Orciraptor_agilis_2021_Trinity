#!/bin/bash

source ../config.txt

${Trinity_dir}/util/abundance_estimates_to_matrix.pl \
--est_method salmon \
--gene_trans_map none \
--name_sample_by_basedir \
${mydir}/Module_5/salmon/V1S1.salmon_quant/quant.sf \
${mydir}/Module_5/salmon/V1S2.salmon_quant/quant.sf \
${mydir}/Module_5/salmon/V1S3.salmon_quant/quant.sf \
${mydir}/Module_5/salmon/V1S4.salmon_quant/quant.sf \
${mydir}/Module_5/salmon/V1S5.salmon_quant/quant.sf \
${mydir}/Module_5/salmon/V1S6.salmon_quant/quant.sf \
${mydir}/Module_5/salmon/V1S7.salmon_quant/quant.sf \
${mydir}/Module_5/salmon/V1S8.salmon_quant/quant.sf \
${mydir}/Module_5/salmon/V1S9.salmon_quant/quant.sf \

${Trinity_dir}/util/misc/contig_ExN50_statistic.pl \
salmon_TPM_matrix \
${mydir}/Module_3/transdecoder/Orciraptor_non-redundant.fas \
| tee ExN50.stats

${Trinity_dir}/util/misc/plot_ExN50_statistic.Rscript \
ExN50.stats 
