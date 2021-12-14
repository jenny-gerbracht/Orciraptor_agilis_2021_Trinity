###################################################################
#Script Name	:DESeq2.R		                                                                                              
#Description	:Run Orciraptor agilis gene expression analysis and generate figures                                                       
#Author	:Jennifer Gerbracht                                               
#Email		:jennifer.gerbracht@gmx.de                                           
###################################################################

library(tximport)
library(DESeq2)
library(ggplot2)
library(ggrepel)
library(dendextend)
library(grid)
library(gridExtra)
library(tibble)
library(dplyr)
library(tidyr)
library(stringr)
library(forcats)
library(ochRe)
library(ggpubr)
library(goseq)
library(GO.db)
library(GOfuncR)
library(qvalue)
library(stringr)
library(ComplexHeatmap)
library(RColorBrewer)
library(circlize)
library(seqinr)
library(UpSetR)
library(purrr)

source(file = "../../config.txt")

###################
## Import expression data from salmon
###################

samples_j <- read.table(paste0(mydir, "/experiment.txt"), 
                        col.names = c("samples", "condition"))
samples_j <- samples_j[c(4,5,6,7,8,9,1,2,3),]
files_j <- file.path(paste0(mydir, "/Module_5"), 
                     "salmon", 
                     paste0(samples_j$samples, ".salmon_quant"), 
                     "quant.sf")
names(files_j) <- samples_j$samples
all(file.exists(files_j))
txi_i <- tximport(files_j, 
                  type = "salmon",
                  txOut = TRUE)
tpm <- txi_i$abundance
write.table(tpm,
            file = paste0(mydir, "/Module_6/salmon_TPM_matrix"),
            col.names = TRUE,
            row.names = TRUE,
            quote = FALSE,
            sep = "\t")
tpm <- rownames_to_column(as.data.frame(tpm), var = "transcript")
###################
###################
## DESeq2
###################

samples_j$condition <- factor(samples_j$condition)
dds <- DESeqDataSetFromTximport(txi_i,
                                colData = samples_j,
                                design = ~ condition)

# Pre-filtering low count genes: analogous to edgeR pre-filtering rowSums(cpm(counts)>1)>=2, gives same results as when
# using keep <- rowSums(counts(dds)) >= 10 which is described in the pre-filtering section of the DESeq2 vignette

dds <- dds[rowSums(fpm(dds) > 1) >= 2] 
dds$condition <- relevel(dds$condition, ref = "gliding")
dds <- DESeq(dds)
res_ag <- lfcShrink(dds, coef = "condition_attacking_vs_gliding")
res_dg <- lfcShrink(dds, coef = "condition_digesting_vs_gliding")
dds$condition <- relevel(dds$condition, ref = "digesting")
dds <- DESeq(dds)
resultsNames(dds)
res_ad <- lfcShrink(dds, coef = "condition_attacking_vs_digesting")
###################
###################
## Building output tables
###################

# Base table gene expression info
# TPM values
gene_expression <- as.data.frame(txi_i$abundance)
colnames(gene_expression) <- paste0(colnames(gene_expression), "_TPM")
gene_expression <- rownames_to_column(gene_expression, var = "transcript")

# Add baseMean
as.data.frame(res_ag) %>%
  rownames_to_column(var = "transcript") %>% 
  dplyr::select(transcript, baseMean) %>% 
  left_join(gene_expression, by = "transcript") -> gene_expression

# Add log2FC and padj values from each comparison
as.data.frame(res_ag) %>% 
  setNames(., paste0(colnames(res_ag), "_ag")) %>% 
  rownames_to_column(var = "transcript") %>% 
  dplyr::select(transcript, log2FoldChange_ag, padj_ag) %>% 
  left_join(gene_expression, by = "transcript") -> gene_expression

as.data.frame(res_dg) %>% 
  setNames(., paste0(colnames(res_dg), "_dg")) %>% 
  rownames_to_column(var = "transcript") %>% 
  dplyr::select(transcript, log2FoldChange_dg, padj_dg) %>% 
  left_join(gene_expression, by = "transcript") -> gene_expression

as.data.frame(res_ad) %>% 
  setNames(., paste0(colnames(res_ad), "_ad")) %>% 
  rownames_to_column(var = "transcript") %>% 
  dplyr::select(transcript, log2FoldChange_ad, padj_ad) %>% 
  left_join(gene_expression, by = "transcript") -> gene_expression

gene_expression %>% 
  dplyr::select(transcript, baseMean, V1S4_TPM, V1S5_TPM, V1S6_TPM, V1S7_TPM, V1S8_TPM, V1S9_TPM, V1S1_TPM, V1S2_TPM, V1S3_TPM, log2FoldChange_ag, padj_ag, log2FoldChange_dg, padj_dg, log2FoldChange_ad, padj_ad) -> gene_expression

write.table(gene_expression,
            file = paste0(mydir, "/Tables/gene_expression.txt"),
            row.names = FALSE,
            col.names = TRUE,
            quote = FALSE,
            sep = "\t")
###################
###################
## PCA plot
###################

# Extracting transformed count values
vsd <- vst(dds, blind = FALSE)

# Make plot
pcaData <- plotPCA(vsd, 
                   intgroup = c("condition"), 
                   returnData=TRUE)
pcaData$condition <- factor(pcaData$condition, levels = c("gliding", "attacking", "digesting"))   
percentVar <- round(100 * attr(pcaData, "percentVar"))
pdf(file = paste0(mydir, "/Figures/PCA.pdf"), width = 5, height = 7)
ggplot(pcaData, aes(PC1, PC2, color=condition)) +
  geom_point(size = 3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  coord_fixed() +
  scale_color_manual(values=c("#3333CC", "#D60093", "#1F4E79")) +
  theme_minimal()
dev.off()

###################
###################
## Clustering
###################

# Cutoff log2FC 1, padj 0.001
as.data.frame(res_ag) %>% 
  filter(padj < 0.001 & abs(log2FoldChange) >= 1) -> res_ag_diff
as.data.frame(res_ad) %>% 
  filter(padj < 0.001 & abs(log2FoldChange) >= 1) -> res_ad_diff
as.data.frame(res_dg) %>% 
  filter(padj < 0.001 & abs(log2FoldChange) >= 1) -> res_dg_diff
dge_set <- c(rownames(res_ag_diff), 
             rownames(res_ad_diff), 
             rownames(res_dg_diff))
dge_set <- unique(dge_set)

# Cutoff log2FC 2, padj 0.001
as.data.frame(res_ag) %>% 
  filter(padj < 0.001 & abs(log2FoldChange) >= 2) -> res_ag_diff_2
as.data.frame(res_ad) %>% 
  filter(padj < 0.001 & abs(log2FoldChange) >= 2) -> res_ad_diff_2
as.data.frame(res_dg) %>% 
  filter(padj < 0.001 & abs(log2FoldChange) >= 2) -> res_dg_diff_2
dge_set_2 <- c(rownames(res_ag_diff_2), 
             rownames(res_ad_diff_2), 
             rownames(res_dg_diff_2))
dge_set_2 <- unique(dge_set_2)

# Generate z-score matrix for heatmaps
vsd_matrix <- assay(vsd)
vsd_matrix_dge <- vsd_matrix[dge_set,]
heat <- t(scale(t(vsd_matrix)))
heat <- heat[dge_set,]

# Clusters rows by Pearson correlation as distance method
hc <- hclust(as.dist(1 - cor(t(as.matrix(vsd_matrix_dge)))))
my_transcript_partition_assignments <- cutree(hc, h = 80/100*max(hc$height))

# Visualise dendrogram with clusters
clust.cutree <- dendextend:::cutree(as.dendrogram(hc), h = 80/100*max(hc$height), order_clusters_as_data = FALSE)
idx <- order(names(clust.cutree))
clust.cutree <- clust.cutree[idx]
df.merge <- merge(my_transcript_partition_assignments, clust.cutree, by = "row.names")
df.merge.sorted <- df.merge[order(df.merge$y),]
lbls <- unique(df.merge.sorted$x)
dend1 <- color_branches(as.dendrogram(hc), h = 80/100*max(hc$height), groupLabels = lbls)
pdf(file = paste0(mydir, "/Figures/clustering_dendrogram.pdf"), width = 12, height = 5)
plot(dend1, leaflab = "none")
abline(h=80/100*max(hc$height), lty = 2, col="grey")
dev.off()

# Make factor_labeling.txt for goseq
data.frame(my_transcript_partition_assignments) %>% 
  rownames_to_column(var = "transcripts") -> factor_labeling
colnames(factor_labeling) <- c("transcript", "cluster")

#write.table(factor_labeling, 
#            file = paste0(mydir, "/Module_5/GO_analysis/factor_labeling.txt"),
#            row.names = FALSE,
#            col.names = FALSE,
#            quote = FALSE,
#            sep = "\t")

# Make list of clusters
clusterlist <- list()
for (i in c(1:max(my_transcript_partition_assignments))) {
  cluster <- heat[(my_transcript_partition_assignments == i),]
  clusterlist[[i]] <- cluster
}

# Plot clusters as boxplots
p <- list()
splan <- 3 - 1L

for (i in 1:max(my_transcript_partition_assignments)) {
  box <- rownames_to_column(as.data.frame(clusterlist[[i]]), var = "transcript")
  box_longer <- pivot_longer(box, cols = !transcript, names_to = "sample")
  box_longer %>% 
    filter(sample == "V1S4" | sample == "V1S5" | sample == "V1S6") %>% 
    mutate(cat = "Gliding") -> box_longer1
  box_longer %>% 
    filter(sample == "V1S7" | sample == "V1S8" | sample == "V1S9") %>% 
    mutate(cat = "Attacking") -> box_longer2
  box_longer %>% 
    filter(sample == "V1S1" | sample == "V1S2" | sample == "V1S3") %>% 
    mutate(cat = "Digesting") -> box_longer3
  comb <- rbind(box_longer1, box_longer2, box_longer3)
  comb$cat <-factor(comb$cat, levels = c("Gliding", "Attacking", "Digesting"))
  g <- ggplot(comb, aes(x = cat, y = value)) +
    geom_boxplot(outlier.size = 0,
                 outlier.shape = NA,
                 alpha = 0.5, 
                 color = "#252A52",
                 fill = "#252A52") +
    stat_smooth(aes(x = cat, y = value, group = "#252A52"), color = "#252A52",
                se = TRUE,
                method = "lm", formula = y~poly(x, splan)) +
    ggtitle(paste("cluster", i, "|", "contigs:", nrow(clusterlist[[i]]))) +
    scale_y_continuous(limits = c(-2, 2)) +
    ylab("z-score") +
    xlab(NULL) +
    theme(legend.position = "none", panel.background = element_rect(colour = "darkgrey", size=1))
  p[[i]] <- ggplotGrob(g)
}
rm(box, box_longer, box_longer1, box_longer2, box_longer3, comb)
pdf(file = paste0(mydir, "/Figures/cluster_boxplots_log2FC1.pdf"), width = 15, height = 3)
grid.arrange(grobs = p, ncol = 5)
dev.off()


###################
###################
## GO analysis
###################
## Needs blast2go.annot file 
## Needs transcript lengths obtained from Trinity scripts fasta_seq_length.pl

# Import transcript <-> GO mapping
blast2go <- read.delim(file = paste0(mydir, "/Module_5/GO_analysis/blast2go.annot"),
                       header = FALSE,
                       sep = "\t",
                       col.names = c("transcript", "GO", "text"))
blast2go %>% 
  dplyr::select(transcript, GO) %>% 
  filter(str_detect(GO, "GO:")) -> blast2go
blast2go %>% 
  group_by(transcript) %>% 
  distinct() %>% 
  summarise(all_GO = paste0(GO, collapse = ",")) -> go_annotation

#write.table(go_annotation,
#            file = paste0(mydir, "/Tables/go_annotation.txt"),
#            quote = FALSE,
#            sep = "\t")

DE_transcripts <- factor_labeling$transcript
transcript_lengths <- read.table(file = paste0(mydir, "/Module_5/GO_analysis/transcript_lengths.txt"), 
                           header = TRUE, 
                           row.names=1, 
                           com = '')
transcript_lengths <- as.matrix(transcript_lengths[,1, drop=FALSE])
background <- rownames(dds)

cluster_dotplots <- list()
cluster_dotplots_psorted <- list()
cluster_lolipop <- list()
cluster_results <- list()
clusters <- unique(factor_labeling$cluster)
sample_set_transcript_lengths <- transcript_lengths[background,]

for (i in clusters) {
  deg <- filter(factor_labeling, cluster == i)$transcript
  cat_genes_vec = as.integer(background %in% deg)
  pwf = nullp(cat_genes_vec, bias.data = sample_set_transcript_lengths)
  rownames(pwf) = background
  result <- goseq(pwf, gene2cat = blast2go)
  GOs <- result$category
  transcriptSets <- list()
  backgroundSets <- list()
  
  for (j in GOs) {
    filter(go_annotation, str_detect(all_GO, paste(j))) %>% 
      summarise(transcripts = paste0(transcript, collapse = "," )) -> tmp
    transcriptSets[[j]] <- str_split(tmp$transcripts, ",", simplify = TRUE)[str_split(tmp$transcripts, ",", simplify = TRUE) %in% filter(factor_labeling, cluster == paste(i))$transcript]
    backgroundSets[[j]] <- sum(str_detect(go_annotation$all_GO[go_annotation$transcript %in% background], j))
  }
  
  results_df <- as.data.frame(matrix(vector(), nrow(result), 8,))
  colnames(results_df) <- c("ID","Description","GeneRatio","BgRatio","pvalue","p.adjust", "transcriptID","Count")
  results_df$ID <- result$category
  results_df$Description <- result$term
  results_df$GeneRatio <- result$numDEInCat/nrow(filter(factor_labeling, cluster == i))
  results_df$BgRatio <- paste0(backgroundSets[results_df$ID], "/", length(background))
  results_df$pvalue <- result$over_represented_pvalue
  results_df$p.adjust <- p.adjust(result$over_represented_pvalue, method = "BH")
  results_df$transcriptID <- transcriptSets[results_df$ID]
  results_df$Count <- result$numDEInCat
  results_df$ontology <- Ontology(results_df$ID)
  results_df <- arrange(results_df, pvalue)
  cluster_results[[i]] <- filter(results_df, pvalue < 0.05)
}

# Check for GO numbers that have no description
cluster_results[[1]][is.na(cluster_results[[1]]$Description),]
cluster_results[[2]][is.na(cluster_results[[2]]$Description),]
cluster_results[[3]][is.na(cluster_results[[3]]$Description),]
cluster_results[[4]][is.na(cluster_results[[4]]$Description),]
cluster_results[[5]][is.na(cluster_results[[5]]$Description),]

# Fill in missing info for Description and Ontology
cluster_results[[1]][3,2] <- c("protein-L-histidine N-pros-methyltransferase activity")
cluster_results[[1]][3,9] <- c("MF")

cluster_results[[3]][45,2] <- c("opioid growth factor receptor activity")
cluster_results[[3]][45,9] <- c("MF")

#How many GO terms have a p-value < 0.05? (included in Figure)
lapply(cluster_results, nrow)

#How many GO terms have an adjusted p-value < 0.05? 
lapply(cluster_results, function(x) nrow(filter(x, p.adjust < 0.05)))

# Plots separately
for (i in clusters) {
  g <- ggplot(cluster_results[[i]][c(1:30),], aes(x = GeneRatio, y = fct_reorder(Description, GeneRatio))) + 
    geom_point(aes(size = Count, color = p.adjust)) +
    theme_bw(base_size = 14) +
    labs(size = "Count") +
    scale_colour_gradient(limits = c(0,0.05), low="red") +
    scale_size_continuous(limits = c(1, 160), range = c(2,10), breaks = c(1, 10, 50, 100, 150)) +
    ylab(NULL) +
    guides(size = guide_legend(order = 1)) +
    ggtitle(paste0("Cluster ", i))
  cluster_dotplots[[i]] <- ggplotGrob(g)
  
  h <- ggplot(cluster_results[[i]][c(1:30),], aes(x = -log10(p.adjust), y = fct_reorder(Description, -log10(p.adjust)))) + 
    geom_col(width = 0.8) +
    theme_bw(base_size = 14) +
    geom_vline(xintercept = 1.30103) +
    theme(axis.text = element_text(size = 10), axis.text.x = element_blank(), axis.ticks.x = element_blank()) +
    ylab(NULL) +
    xlab("-log10 FDR-corrected p-value") +
    geom_text(aes(x = 1.7, label = "p.adjust < 0.05", y = 5), angle = 90, size = 5)
  cluster_dotplots_psorted[[i]] <- ggplotGrob(h)
  
  l <- ggplot(cluster_results[[i]][c(1:20),], aes(x = -log10(p.adjust), y = fct_reorder(Description, -log10(p.adjust))))+
    geom_vline(xintercept = 1.30103, linetype = "dashed", color = "grey") +
    geom_segment(aes(y=fct_reorder(Description, -log10(p.adjust)), yend=fct_reorder(Description, -log10(p.adjust)), x=0, xend=-log10(p.adjust)), color = "darkgrey") +
    geom_point(aes(color = factor(ontology)), size = 5) +
    scale_color_manual(values = c("#00AFBB", "#E7B800", "#FC4E07")) +
    geom_text(aes(label = Count), size = 2.5, color = "white") +
    theme_classic() +
    #scale_x_continuous(limits = c(0, 12)) +
    theme(axis.text = element_text(size = 10)) +
    ylab(NULL) +
    xlab("-log10 FDR-corrected p-value") 
  #geom_text(aes(x = 1.7, label = "< 0.05", y = 5), angle = 90, size = 5)
  cluster_lolipop[[i]] <- ggplotGrob(l)
}

# Cluster 2 with 30 terms
l <- ggplot(cluster_results[[2]][c(1:30),], aes(x = -log10(p.adjust), y = fct_reorder(Description, -log10(p.adjust))))+
  geom_vline(xintercept = 1.30103, linetype = "dashed", color = "grey") +
  geom_segment(aes(y=fct_reorder(Description, -log10(p.adjust)), yend=fct_reorder(Description, -log10(p.adjust)), x=0, xend=-log10(p.adjust)), color = "darkgrey") +
  geom_point(aes(color = factor(ontology)), size = 6) +
  scale_color_manual(values = c("#00AFBB", "#E7B800", "#FC4E07"), guide_legend(title = "Ontology")) +
  geom_text(aes(label = Count), size = 2.5, color = "white") +
  theme_classic() +
  #scale_x_continuous(limits = c(0, 12)) +
  theme(axis.text = element_text(size = 10)) +
  ylab(NULL) +
  xlab("-log10 FDR-corrected p-value") 
cluster_lolipop[[2]] <- ggplotGrob(l)

# Cluster 4 with 40 terms
l <- ggplot(cluster_results[[4]][c(1:40),], aes(x = -log10(p.adjust), y = fct_reorder(Description, -log10(p.adjust))))+
  geom_vline(xintercept = 1.30103, linetype = "dashed", color = "grey") +
  geom_segment(aes(y=fct_reorder(Description, -log10(p.adjust)), yend=fct_reorder(Description, -log10(p.adjust)), x=0, xend=-log10(p.adjust)), color = "darkgrey") +
  geom_point(aes(color = factor(ontology)), size = 6) +
  scale_color_manual(values = c("#00AFBB", "#E7B800", "#FC4E07"), guide_legend(title = "Ontology")) +
  geom_text(aes(label = Count), size = 2.5, color = "white") +
  theme_classic() +
  #scale_x_continuous(limits = c(0, 12)) +
  theme(axis.text = element_text(size = 10)) +
  ylab(NULL) +
  xlab("-log10 FDR-corrected p-value") 
cluster_lolipop[[4]] <- ggplotGrob(l)

#pdf(file = paste0(mydir, "/Figures/cluster_dotplots_log2FC1_v2.pdf"), width = 14, height = 30)
grid::grid.newpage()
grid::grid.draw(rbind(cluster_dotplots[[1]], 
                      cluster_dotplots[[2]],
                      cluster_dotplots[[3]],
                      cluster_dotplots[[4]],
                      cluster_dotplots[[5]]))
#dev.off()


#pdf(file = paste0(mydir, "/Figures/cluster_dotplots_psorted_log2FC1_v2.pdf"), width = 14, height = 30)
grid::grid.newpage()
grid::grid.draw(rbind(cluster_dotplots_psorted[[1]], 
                      cluster_dotplots_psorted[[2]],
                      cluster_dotplots_psorted[[3]],
                      cluster_dotplots_psorted[[4]],
                      cluster_dotplots_psorted[[5]]))
dev.off()


pdf(file = paste0(mydir, "/Figures/cluster_lolipop.pdf"), width = 7.5, height = 22)
grid::grid.newpage()
grid::grid.draw(rbind(cluster_lolipop[[1]], 
                      cluster_lolipop[[2]],
                      cluster_lolipop[[3]],
                      cluster_lolipop[[4]],
                      cluster_lolipop[[5]]))
dev.off()

filter(cluster_results[[4]], p.adjust < 0.05)$Description

# Plot cluster 2 extra
pdf(file = paste0(mydir, "/Figures/cluster_lolipop_cluster2.pdf"), width = 6, height = 7)
grid::grid.newpage()
grid::grid.draw(cluster_lolipop[[2]])
dev.off()

# Plot cluster 4 extra
pdf(file = paste0(mydir, "/Figures/cluster_lolipop_cluster4.pdf"), width = 6, height = 9)
grid::grid.newpage()
grid::grid.draw(cluster_lolipop[[4]])
dev.off()

# Retrieve any GO terms significantly enriched and their ontology
sig_GOs <- rbind(filter(cluster_results[[1]], p.adjust < 0.05),
                 filter(cluster_results[[2]], p.adjust < 0.05),
                 filter(cluster_results[[3]], p.adjust < 0.05),
                 filter(cluster_results[[4]], p.adjust < 0.05),
                 filter(cluster_results[[5]], p.adjust < 0.05))
distinct(sig_GOs, ID, .keep_all = TRUE) %>% 
  dplyr::select(ID, Description) -> sig_GOs_unique

sig_GOs_unique[is.na(sig_GOs_unique$Category),]
write.table(sig_GOs_unique,
            file = paste0(mydir, "/Tables/sig_GOs_unique.txt"),
            col.names = TRUE,
            row.names = FALSE,
            quote = FALSE,
            sep = "\t")

###################
###################
## Heatmaps
###################

go_annotation %>% 
  separate_rows(all_GO, sep = ",") %>% 
  distinct() -> go_annotation_long

heat %>% 
  as.data.frame() %>% 
  rownames_to_column(var = "transcript") %>% 
  left_join(go_annotation_long, by = "transcript") %>% 
  filter(!is.na(all_GO)) -> heat_anno

# Make heatmap with GO mean z-scores, ordered
# Import manually annotated GO list
sig_GOs_annotated <- read.delim(file = paste0(mydir, "/Tables/sig_GOs_unique_anno.txt"))
sig_GOs_annotated$Category_plot <- factor(sig_GOs_annotated$Category_plot,
                                          levels = c("Signalling", "Metabolic process", "Ribosome and translation",
                                                     "Nucleic acid process and cell division", "Protein-related", "Carbohydrate-binding", 
                                                     "Cytoskeleton", "Cellular component"))
sig_GOs_annotated <- arrange(sig_GOs_annotated, Category_plot)

mylist_go_ordered <- list()
for (i in sig_GOs_annotated$ID){
  a <- filter(heat_anno, str_detect(all_GO, i))
  mylist_go_ordered[[length(mylist_go_ordered) + 1]] <- a
}
names(mylist_go_ordered) <- sig_GOs_annotated$ID
colmeans_GO <- list()
for (i in c(1:nrow(sig_GOs_annotated))) {
  colmeans_GO[[i]] <- colMeans(dplyr::select(mylist_go_ordered[[i]], c(2:10)))
}
colmeans_GO_df <- do.call(rbind, colmeans_GO)
rownames(colmeans_GO_df) <- sig_GOs_annotated$Description
colmeans_GO_df <- as.data.frame(colmeans_GO_df )

pdf(file = paste0(mydir, "/Figures/GO_mean_z_order_redcyan.pdf"), width = 10, height = 14)
draw(Heatmap(as.matrix(colmeans_GO_df), 
             name = "mean z-score", 
             cluster_rows = TRUE,
             cluster_columns = FALSE,
             cluster_row_slices = FALSE,
             show_column_names = FALSE,
             show_row_dend = FALSE,
             row_title_rot = 0,
             column_split = rep(1:3, each = 3),
             row_split = c(rep("A", 3), rep("B", 11), rep("C", 11), rep("D", 24), rep("E", 5), rep("F", 1), rep("G", 3), rep("H", 24)),
             row_title = levels(sig_GOs_annotated$Category_plot),
             bottom_annotation = HeatmapAnnotation(foo = anno_block(labels = c("G", "A", "D"))),
             col = colorRamp2(c(-1.5,-0.75,0,0.75,1.5), colorRampPalette(c("cyan", "black", "red"))(5)),
             heatmap_legend_param = list(at = c(-1.5,-0.75,0,0.75,1.5))),
     padding = unit(c(2,2,2,40), "mm"),
     heatmap_legend_side = "left")
dev.off()

###################
###################
## Heatmap cytoskeleton
###################

#Read in eggnog output
eggnog_diamond <- read.delim(paste0(mydir, "/Module_4/eggnog/orciraptor_diamond.emapper.annotations"), 
                             skip = 4,
                             quote = "")
eggnog_diamond <- eggnog_diamond[-c((nrow(eggnog_diamond)-2):(nrow(eggnog_diamond))),]
eggnog_diamond$source <- c("diamond")

# HMM
eggnog_hmm <- read.delim(paste0(mydir, "/Module_4/eggnog/orciraptor_hmm.emapper.annotations"), 
                         skip = 4,
                         quote = "")
eggnog_hmm <- eggnog_hmm[-c((nrow(eggnog_hmm)-2):(nrow(eggnog_hmm))),]
eggnog_hmm$source <- c("hmm")

# Keep all annotations from Diamond, add annotations for peptide sequences that were annotated only in HMM mode
eggnog_combined <- rbind(eggnog_diamond, (eggnog_hmm[!eggnog_hmm$X.query_name %in% eggnog_diamond$X.query_name,]))

# Write merged output
write.table(eggnog_combined,
            file = paste0(mydir, "/Module_4/eggnog/eggnog_merged.tsv"),
            col.names = TRUE,
            row.names = FALSE,
            quote = FALSE,
            sep = "\t")

eggnog_combined %>% 
  filter(narr_og_cat != "-")  %>% 
  dplyr::select(X.query_name, Preferred_name, narr_og_cat, narr_og_desc) %>% 
  rename(., transcript = X.query_name) -> eggnog_combined_select
eggnog_combined_select$transcript <- str_split(eggnog_combined_select$transcript, "_", simplify = TRUE)[,1]
eggnog_z <- filter(eggnog_combined_select, str_detect(narr_og_cat, "Z"))

heat %>% 
  as.data.frame() %>% 
  rownames_to_column(var = "transcript") -> heat_df
heat_df$transcript <- str_split(heat_df$transcript, "_", simplify = TRUE)[,1]
heat_df %>% 
  left_join(eggnog_z, by = "transcript") %>%
  filter(!is.na(narr_og_cat)) %>% 
  mutate(anno1 = paste(transcript, " | ", Preferred_name, " | ", narr_og_desc)) -> heat_eggnog_go

heat_eggnog_go %>%
  dplyr::select(14, 2:10) %>%
  column_to_rownames(var = "anno1")-> heat_eggnog_go

write.table(heat_eggnog_go,
            file = paste0(mydir, "/Tables/cyto_DE.txt"),
            sep = "\t",
            quote = FALSE)

# Automatically annotated heatmap
pdf(file = paste0(mydir, "/Figures/cyto_heatmap_order_redcyan_genes.pdf"), width = 15, height = 18)
draw(Heatmap(as.matrix(heat_eggnog_go), 
             name = "mean z-score", 
             cluster_rows = TRUE,
             cluster_columns = FALSE,
             cluster_row_slices = FALSE,
             show_column_names = FALSE,
             show_row_dend = FALSE,
             row_title_rot = 0,
             column_split = rep(1:3, each = 3),
             bottom_annotation = HeatmapAnnotation(foo = anno_block(labels = c("G", "A", "D"))),
             col = colorRamp2(c(-1.5,-0.75,0,0.75,1.5), colorRampPalette(c("cyan", "black", "red"))(5)),
             heatmap_legend_param = list(at = c(-1.5,-0.75,0,0.75,1.5))),
     padding = unit(c(2,2,2,150), "mm"),
     heatmap_legend_side = "left")
dev.off()

heat_eggnog_go_anno <- read.delim(file = paste0(mydir, "/Tables/cyto_DE_anno.txt"))
heat_eggnog_go$Category <- heat_eggnog_go_anno$Category
heat_eggnog_go$Annotation <- heat_eggnog_go_anno$Annotation
cyto_categories <- heat_eggnog_go_anno$Category
rownames(heat_eggnog_go) <- NULL

filter(heat_eggnog_go, Category != "Non-cytoskeleton") %>%
  arrange(Category) -> cyto_heat

pdf(file = paste0(mydir, "/Figures/cyto_heatmap_order_redcyan.pdf"), width = 10, height = 15)
draw(Heatmap(as.matrix(cyto_heat[,-c(10,11)]), 
             name = "mean z-score",
             row_labels = cyto_heat$Annotation,
             cluster_rows = TRUE,
             cluster_columns = FALSE,
             cluster_row_slices = FALSE,
             show_column_names = FALSE,
             show_row_dend = FALSE,
             row_title_rot = 0,
             column_split = rep(1:3, each = 3),
             row_split = c(rep("A", 6), rep("B", 6), rep("C", 2), rep("D", 1), rep("E", 28), rep("F", 3),
                           rep("G", 14), rep("H", 9), rep("I", 2)),
             row_title = c("Actin", "Actin-related", "Cell division", "Dynein",
                           "Kinesin", "Microtubule-related", "Myosin", "Regulatory (Cell cycle)", "Tubulin"),
             bottom_annotation = HeatmapAnnotation(foo = anno_block(labels = c("Gliding", "Attacking", "Digesting"))),
             col = colorRamp2(c(-1.5,-0.75,0,0.75,1.5), colorRampPalette(c("cyan", "black", "red"))(5)),
             heatmap_legend_param = list(at = c(-1.5,-0.75,0,0.75,1.5))),
     padding = unit(c(2,2,2,40), "mm"),
     heatmap_legend_side = "left")
dev.off()

###################
###################
## Upset plot
###################

# Obtain all peptide sequence identifiers 
transdecoder <- read.fasta(file = paste0(mydir, "/Module_3/transdecoder/Orciraptor_non-redundant.faa"))
peptide_ids <- attr(transdecoder, "names")
rm(transdecoder)

# Obtain sequence identifiers from diamond blast vs SP
diamond_SP <- read.delim(file = paste0(mydir, "/Module_4/diamond/diamond_SP.txt"), header = FALSE)
diamond_SP_ids <- diamond_SP$V1
rm(diamond_SP)

# Obtain sequence identifiers from ips
ips <- read.delim(file = paste0(mydir, "/Module_4/interproscan/Orciraptor_non-redundant.faa_nostop.fasta.tsv"), 
                  sep = "\t", header = FALSE, fill = TRUE,
                  col.names = c("qseqid" , "md5" , "length" , "analysis", "accession", "description", "start", "stop", "evalue", "status", "date","IPS_accession", "IPS_description", "GO", "PA"))
ips_ids <- ips$qseqid
rm(ips)

# Obtain sequence identifiers from eggnog
eggnog_ids <- eggnog_combined$X.query_name

# Obtain sequence identifiers from cazy
cazy <- read.delim(paste0(mydir, "/Module_4/cazy/cazy/overview.txt"))
cazy_ids <- cazy$Gene.ID
rm(cazy)

#Upset plot

listInput <- list("ORFs" = unique(peptide_ids),
                  "Swiss-Prot" = unique(diamond_SP_ids),
                  "InterProScan" = unique(ips_ids),
                  "EggNOG" = unique(eggnog_ids),
                  "CAZy" = unique(cazy_ids))
pdf(file = paste0(mydir, "/Figures/upset.pdf"),
    height = 8,
    width = 20)
upset(fromList(listInput), 
      nsets = 5, 
      order.by = "freq",
      mb.ratio = c(0.6, 0.4))
dev.off()

###################
###################
## CAZy plot: Figure 3A
###################
###################
## Figure ranked TPM, CBM
###################

cazy_hmmer.out <- read.delim(file = paste0(mydir, "/Module_4/cazy/cazy/hmmer.out"))
cazy_hmmer.out %>% 
  group_by(Gene.ID) %>% 
  distinct(HMM.Profile, .keep_all = TRUE) -> cazy_unique_HMM
cazy_unique_HMM$HMM.Profile <- sub(".hmm", "", cazy_unique_HMM$HMM.Profile)
cazy_unique_HMM$family <- gsub('[[:digit:]]+', '', cazy_unique_HMM$HMM.Profile)
cazy_unique_HMM$family <- str_split(cazy_unique_HMM$family, "_", simplify = TRUE)[,1]

# Add info about up / or downregulation
gene_expression %>% 
  mutate(up = ifelse(gene_expression$log2FoldChange_ag > 1 & gene_expression$padj_ag < 0.001,
                     TRUE,
                     FALSE)) %>% 
  mutate(down = ifelse(gene_expression$log2FoldChange_ag < -1 & gene_expression$padj_ag < 0.001,
                       TRUE,
                       FALSE)) %>% 
  dplyr::select(transcript, up, down) -> gene_expression_cats

cazy_unique_HMM %>% 
  dplyr::select(3, 1, 11) %>% 
  setNames(., c("transcript", "HMM", "family")) -> cazy_unique_HMM

cazy_unique_HMM %>% 
  left_join(gene_expression_cats) -> cazy_unique_HMM_expression

# Omit cazys for which there is no gene expression data 
cazy_unique_HMM_expression <- na.omit(cazy_unique_HMM_expression)

# Get number of HMMs for plotting (calculating height)
cazy_unique_HMM_expression %>% 
  group_by(family) %>% 
  summarise(number = length(unique(HMM))) -> number_family
number_family <- as.data.frame(number_family)

cazy_unique_HMM_expression %>% 
  ungroup() %>% 
  left_join(gene_expression, by = "transcript") %>% 
  group_by(transcript) %>% 
  mutate(attacking_TPM = mean(c(V1S7_TPM, V1S8_TPM, V1S9_TPM))) %>% 
  ungroup() %>% 
  dplyr::select(transcript, HMM, family, up, down, attacking_TPM) -> dotplot_data

dotplot_data %>% 
  mutate(status = case_when(dotplot_data$up == "TRUE" ~ "up",
                            dotplot_data$down == "TRUE" ~ "down",
                            dotplot_data$up == "FALSE" & dotplot_data$down == "FALSE" ~ "none")) -> dotplot_data


cols <- c("up" = "#ca0020", "down" = "#0571b0", "none" = "#595959")

dotplot_data %>% 
  filter(family == "CBM") %>% 
  group_by(transcript) %>% 
  distinct(HMM, .keep_all = TRUE) %>% 
  mutate(HMM = paste0(HMM, collapse = "+")) %>%
  distinct(HMM, .keep_all = TRUE) %>% 
  ungroup() %>% 
  arrange(-attacking_TPM) -> attacking_TPM_df

write.table(attacking_TPM_df[c(1:50),c(1:2)],
            file = paste0(mydir, "/Tables/chitin_TPM.txt"),
            row.names = FALSE,
            col.names = TRUE,
            quote = FALSE,
            sep = "\t")

# Plot the top 50 in decreasing order

#pdf(file = paste0(mydir, "/Figures/CBM_TPM.pdf"), height = 12, width = 3.5)
attacking_TPM_df[c(1:50),] %>% 
  mutate(new = fct_reorder(transcript, attacking_TPM)) %>% 
  ggplot(aes(y = new, x = attacking_TPM, fill = status)) +
  geom_col(width = 0.8) +
  scale_y_discrete(name = NULL, labels = rev(attacking_TPM_df$HMM[c(1:50)])) +
  scale_x_continuous(name = "TPM") +
  scale_fill_manual(values = cols, guide = NULL) +
  theme_bw() -> a
#dev.off()
###################
###################
## Figure ranked TPM, cleaving except AA11
###################

filter(dotplot_data, family == "GH") %>% 
  group_by(HMM) %>% 
  summarise(n = n()) %>% 
  arrange(desc(n)) -> GH_families_ranked 

filter(dotplot_data, family == "PL") %>% 
  group_by(HMM) %>% 
  summarise(n = n()) %>% 
  arrange(desc(n)) -> PL_families_ranked 

cleaving_families_ranked <- rbind(GH_families_ranked,
                                  PL_families_ranked)

cleaving_families <- cleaving_families_ranked$HMM

dotplot_data %>% 
  filter(HMM %in% cleaving_families) %>% 
  arrange(-attacking_TPM) -> attacking_TPM_df

write.table(attacking_TPM_df[c(1:50),c(1:2)],
            file = paste0(mydir, "/Tables/lytic_TPM.txt"),
            row.names = FALSE,
            col.names = TRUE,
            quote = FALSE,
            sep = "\t")

# Plot the top 50 in decreasing order

#pdf(file = paste0(mydir, "/Figures/Lytic_TPM.pdf"), height = 12, width = 3.5)
attacking_TPM_df[c(1:50),] %>% 
  mutate(new = fct_reorder(transcript, attacking_TPM)) %>% 
  ggplot(aes(y = new, x = attacking_TPM, fill = status)) +
  geom_col(width = 0.8) +
  scale_y_discrete(name = NULL, labels = rev(attacking_TPM_df$HMM[c(1:50)])) +
  scale_x_continuous(name = "TPM") +
  scale_fill_manual(values = cols, guide = NULL) +
  theme_bw() -> b
#dev.off()

###################
## Figure ranked TPM, AA11
###################

filter(dotplot_data, HMM == "AA11") %>% 
  group_by(HMM) %>% 
  summarise(n = n()) %>% 
  arrange(desc(n)) -> AA_families_ranked 

cleaving_families <- AA_families_ranked$HMM

dotplot_data %>% 
  filter(HMM %in% AA_families_ranked) %>% 
  arrange(-attacking_TPM) -> attacking_TPM_df

write.table(attacking_TPM_df[c(1:50),c(1:2)],
            file = paste0(mydir, "/Tables/AA11_TPM.txt"),
            row.names = FALSE,
            col.names = TRUE,
            quote = FALSE,
            sep = "\t")

# Plot the top 3 in decreasing order

attacking_TPM_df %>% 
  mutate(new = fct_reorder(transcript, attacking_TPM)) %>% 
  ggplot(aes(y = new, x = attacking_TPM, fill = status)) +
  geom_col(width = 0.8) +
  scale_y_discrete(name = NULL, labels = rev(attacking_TPM_df$HMM)) +
  scale_x_continuous(name = "TPM") +
  scale_fill_manual(values = cols, guide = NULL) +
  theme_bw() -> c

###################

pdf(file = paste0(mydir, "/Figures/TPM.pdf"), width = 8.2, height = 11.6)
grid::grid.newpage()
grid::grid.draw(cbind(ggplotGrob(a), 
                      ggplotGrob(b)))
dev.off()

pdf(file = paste0(mydir, "/Figures/TPM_AA11.pdf"), width = 4.1, height = 1.25)
grid::grid.newpage()
grid::grid.draw(ggplotGrob(c))
dev.off()

###################
###################
## GH5_5 expression
###################

dds$condition <- factor(c(rep("gliding", 3),
                        rep("attacking", 3),
                        rep("digesting", 3)),
                        levels = c("gliding", "attacking", "digesting"))
plotCounts(dds,
           "m.81072_type-3prime_partial_len1679_TR22360_c1_g1_i3",
           intgroup = "condition",
           normalized = TRUE,
           transform = TRUE,
           main = "Normalised counts of m.81072",
           xlab = c("gliding", "attacking", "digesting"),
           replaced = FALSE)


GH5_5_counts <- counts(dds, normalized = TRUE)["m.81072_type-3prime_partial_len1679_TR22360_c1_g1_i3",]
GH5_5_counts <- data.frame(GH5_5_counts)
GH5_5_counts$condition <- factor(c(rep("gliding", 3),
                                   rep("attacking", 3),
                                   rep("digesting", 3)),
                                 levels = c("gliding", "attacking", "digesting"))

set.seed(10)
pdf(file = paste0(mydir, "/Figures/GH5_5_counts.pdf"), height = 5)
ggplot(GH5_5_counts, aes(x = condition, y = log2(GH5_5_counts))) +
  geom_boxplot(width = 0.5) +
  stat_boxplot(geom = "errorbar",
               width = 0.3,
               linetype = "dashed") +
  geom_point(aes(color = condition),
              size = 4) +
  scale_color_manual(values=c("#3333CC", "#D60093", "#1F4E79")) +
  theme_bw()
dev.off()

###################
## CAZy plot: Figure 3?
###################

cazy_out <- read.delim(paste0(mydir, "/Module_4/cazy/cazy/overview.txt"))
cazy_out %>% 
  dplyr::select(Gene.ID, HMMER) %>% 
  rename(., transcript = Gene.ID) %>% 
  left_join(tpm, by = "transcript") -> cazy_out
cazy_out$HMMER <- gsub("\\(.*?\\)", "", cazy_out$HMMER)

gliding <- c("V1S4", "V1S5", "V1S6")
attacking <- c("V1S7", "V1S8", "V1S9")
digesting <- c("V1S1", "V1S2", "V1S3")

cazy_out %>% 
  pivot_longer(cols = starts_with("V1S"), names_to = "sample", values_to = "TPM") %>% 
  mutate(group = case_when(sample %in% gliding ~ "gliding",
                           sample %in% attacking ~ "attacking",
                           sample %in% digesting ~ "digesting")) %>% 
  group_by(transcript, group) %>% 
  mutate(mean_TPM = mean(TPM)) %>% 
  ungroup() %>% 
  mutate(family = case_when(grepl("GH", HMMER) ~ "GH",
                            grepl("GT", HMMER) ~ "GT",
                            grepl("PL", HMMER) ~ "PL",
                            grepl("CE", HMMER) ~ "CE",
                            grepl("AA", HMMER) ~ "AA",
                            grepl("CBM", HMMER) ~ "CBM")) -> cazy_expression_longer

# Attacking

cazy_expression_longer %>% 
  filter(group == "attacking") %>% 
  dplyr::select(!c(3,4)) %>% 
  distinct(transcript, .keep_all = TRUE) %>% 
  group_by(transcript) %>% 
  distinct(HMMER, .keep_all = TRUE) %>% 
  mutate(HMM = paste0(HMMER, collapse = "+")) %>%
  distinct(HMM, .keep_all = TRUE) %>% 
  ungroup()  %>% 
  arrange(-mean_TPM) -> attacking_df

# Plot the top 50 in decreasing order

attacking_df[c(1:50),] %>% 
  mutate(new = fct_reorder(transcript, mean_TPM)) %>% 
  ggplot(aes(y = new, x = mean_TPM, fill = family)) +
  geom_bar(stat = "identity") +
  scale_y_discrete(name = NULL, labels = rev(attacking_df$HMM[c(1:50)])) +
  scale_x_continuous(name = "TPM") +
  scale_fill_ochre() +
  theme_bw()

###################
## Gene expression plots cellulases
###################

as.data.frame(res_ag) %>% 
  rownames_to_column(var = "Gene.ID") -> dge_as
dge_as$cazy <- ifelse(dge_as$Gene.ID %in% cazy_out$Gene.ID, "TRUE", "FALSE")
dge_as %>% 
  left_join(cazy_out) -> dge_as

endoglucanases <- c("GH5", "GH5_1", "GH5_5", "GH5_9", "GH5_10", "GH5_11", "GH5_14", "GH5_27", "GH5_29", "GH5_31", "GH5_54", "GH6", "GH7", "GH10", "GH44")
#betaglucosidase <- c("GH1", "GH2", "GH3", "GH30")

dge_as %>% 
  mutate(endoglucanase = ifelse(dge_as$HMMER %in% c(endoglucanases), "TRUE", "FALSE")) -> dge_as

dge_as$HMMER[dge_as$endoglucanase == FALSE |
                    abs(dge_as$log2FoldChange) < 2 |
                    dge_as$padj > 0.001] <- ""

ggplot(data = dge_as, 
       aes(x = log2FoldChange, 
           y = -log10(pvalue),
           col = endoglucanase, 
           alpha = endoglucanase,
           label = HMMER)) +
  geom_point() +
  geom_point(data = filter(dge_as, endoglucanase == TRUE)) +
  geom_text_repel(color = "black",
                  max.overlaps = Inf) +
  scale_color_manual(values = c("grey", "royalblue")) +
  scale_alpha_manual(values = c(0.2, 0.8)) +
  theme_bw(base_size = 18) +
  geom_vline(xintercept=c(-2, 2), col = "black", linetype = "dashed") +
  geom_hline(yintercept=-log10(0.05), col = "black", linetype = "dashed") +
  guides(color = FALSE, alpha = FALSE)
###################

###################
## Gene expression plots chitin
###################

as.data.frame(res_ag) %>% 
  rownames_to_column(var = "transcript") -> dge_ag
dge_ag$cazy <- ifelse(dge_ag$transcript %in% cazy_out$transcript, "TRUE", "FALSE")
dge_ag %>% 
  left_join(cazy_out) -> dge_ag

chitinases <- c("GH18", "GH19")
chitin_synthase <- c("CE4")
chitin_binding <- c("CBM18", "CBM50") 
chitin_lpmo <- c("AA11")

dge_ag %>% 
  mutate(Chitinase = ifelse(dge_ag$HMMER %in% c(chitinases), "TRUE", "FALSE")) %>% 
  mutate(Chitin_Synthase = ifelse(dge_ag$transcript == c("m.77450_type-complete_len1285_TR22023_c5_g3_i1"), "TRUE", "FALSE")) %>% 
  mutate(Chitin_Binding = ifelse(dge_ag$HMMER %in% c(chitin_binding), "TRUE", "FALSE")) %>% 
  mutate(Chitin_Active_LPMO = ifelse(dge_ag$HMMER %in% c(chitin_lpmo), "TRUE", "FALSE")) -> dge_ag

dge_ag$HMMER_gene[dge_ag$Chitinase == FALSE &
                    dge_ag$Chitin_Synthase == FALSE &
                    dge_ag$Chitin_Binding == FALSE] <- ""

cairo_pdf(file = paste0(mydir, "/Figures/MAplot_chitin.pdf"), width = 8, height = 6)
ggplot(dge_ag, aes(x = (log2(baseMean)+1), y = log2FoldChange, label = HMMER_gene)) +
  #geom_point(alpha = 0.2, color = "grey", shape = 16) +
  geom_point(data = filter(dge_ag, Chitinase == TRUE), color = "#1b9e77", size = 4, shape = 16) +
  geom_point(data = filter(dge_ag, Chitin_Synthase == TRUE), color = "#d95f02", size = 4, shape = 15) +
  geom_point(data = filter(dge_ag, Chitin_Binding == TRUE), color = "#7570b3", size = 4, shape = 17) +
  geom_point(data = filter(dge_ag, Chitin_Active_LPMO == TRUE), color = "#e7298a", size = 4, shape = 8) +
  scale_x_continuous(breaks=seq(0, max((log2(dge_ag$baseMean)+1)), 5)) + 
  labs(x = "Log2 mean expression", y = "Log2 fold change") +
  geom_hline(yintercept = c(0,  1), linetype = c(1, 2),
             color = c("black", "black")) +
  geom_hline(yintercept = c(0,  2), linetype = c(1, 2),
             color = c("black",  "black")) +
  #geom_text_repel(color = "black",
  #max.overlaps = Inf,
  #force = 10,
  #segment.color = "darkgrey") +
  theme_bw(base_size = 18) 
dev.off()
