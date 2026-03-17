setwd("C:/Users/Stuart/Downloads/RNA")

library(fgsea)
library(msigdbr)
library(dplyr)
library(ggplot2)
library(biomaRt)

library(dplyr)          ## loads in dplyr package to use
library(tidyr)          ## loads in tidyr package to use
library(ggplot2)          ## loads in ggplot2 package to use
library(readr)    
library(edgeR)
library(tibble)
library(cowplot)
library(readxl)
library(limma)
library(org.Hs.eg.db)
library(clusterProfiler)
library(tibble)


data <- read_excel("raw_counts2.xlsx")

view(data)
 colnames(data)

str(data)

# Convert the data to numeric 
data_numeric <- as.data.frame(lapply(data[, -1], as.numeric))

view(data_numeric)
# Add 
# Add the gene names back to the data frame
data_numeric <- cbind(Gene = data$gene_id, data_numeric)

# Set the gene names as row names and remove the gene column
rownames(data_numeric) <- data_numeric$Gene

data_numeric <- data_numeric[, -1]

# Create a DGEList object

colnames(data_numeric)

group <- factor(rep(c("CTRL", "ERb.KO", "positiveERb", "positiveERbVH", 
"positiveErb.LY5"), each = 3))

view(group)

dge <- DGEList(counts = data_numeric, group = group)

view(dge)

# Calculate cpm values
cpm_values <- cpm(dge)

view(cpm_values)

colSums(cpm_values)

dge_cpm <- cpm(dge , log = TRUE)

view(dge_cpm)

dge_table <- as_tibble(dge_cpm, rownames = "Gene")


view(dge_table)

rownames(dge_table)

dge_table2 <- pivot_longer(dge_table, cols = -Gene,
 names_to = "Sample", values_to = "CPM")
 
 view(dge_table2)

rownames <- dge_table2$Gene

# Plot the data
plot1 <-ggplot(dge_table2, aes(x = Sample, y = CPM, fill = Sample)) +
  geom_violin(trim = FALSE, show.legend = TRUE) +
  stat_summary(fun = "median", 
  geom = "point", size = 2, color = "black") +
  theme_bw() +
  labs(title = "Violin Plot of logCPM Values by Sample",
       subtitle = "Distribution of logCPM values across samples",
       x = "Sample",
       y = "CPM (log scale)") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5)
        ) +
        coord_flip()

print(plot1)                
 view(plot1)

ggsave("violin_plot.png", 
plot = plot1, width = 10, height = 6, dpi = 300)

#Filtering lowly expressed genes

keep <- filterByExpr(dge)

view(keep)

dge_filtered <- dge[keep, , keep.lib.sizes = FALSE]

view(dge_filtered)

dge_filtered_cpm <- cpm(dge_filtered, log = TRUE)

view(dge_filtered_cpm)


dge_filtered_table <- as_tibble(dge_filtered_cpm, rownames = "Gene")
view(dge_filtered_table)
pivot_longer(dge_filtered_table, cols = -Gene,
 names_to = "Sample", values_to = "logCPM") -> PT1
view(PT1)


plot2 <- ggplot(PT1,aes(x = Sample , y = logCPM , fill = Sample)) +
 geom_violin(trim = FALSE, show.legend = TRUE) +
 stat_summary(fun = "median", geom = "point", size = 2, color = "black", shape = 95) +
 labs(title = "Violin Plot of logCPM Values by Sample (Filtered)", 
 x = "Sample", y = "logCPM", subtitle = "Distribution of logCPM values across samples (Filtered)") +
 theme_bw() +
 theme(
 plot.title = element_text(hjust = 0.5),
 plot.subtitle = element_text(hjust = 0.5),
 axis.text.x = element_text(angle = 45, hjust = 1)
 ) +
 coord_flip()
print(plot2)
ggsave("violin_plot_filtered.png",
 plot = plot2, width = 10, height = 6, dpi = 300)

#Normalization
dge_filtered_norm <- calcNormFactors(dge_filtered , method = "TMM")

view(dge_filtered_norm)

dge_filtered_norm_cpm <- cpm(dge_filtered_norm, log = TRUE)

view(dge_filtered_norm_cpm)

dge_filtered_norm_table <- as_tibble(dge_filtered_norm_cpm, 
rownames = "Gene") 
view(dge_filtered_norm_table)

pivot_longer(dge_filtered_norm_table, cols = -Gene,
 names_to = "Sample", values_to = "logCPM") -> PT2

view(PT2)

plot3 <- ggplot(PT2,aes(x = Sample , y = logCPM , fill = Sample)) +
 geom_violin(trim = FALSE, show.legend = TRUE) +
 stat_summary(fun = "median", geom = "point", size = 2, color = "black", shape = 95) +
 labs(title = "Violin Plot of logCPM Values by Sample (Normalized)", 
 x = "Sample", y = "logCPM", subtitle = "Distribution of logCPM values across samples (Normalized)") +
 theme_bw() +
 theme(
 plot.title = element_text(hjust = 0.5),
 plot.subtitle = element_text(hjust = 0.5),
 axis.text.x = element_text(angle = 45, hjust = 1)
 ) +
 coord_flip()
print(plot3)
ggsave("violin_plot_normalized.png",
 plot = plot3, width = 10, height = 6, dpi = 300)

 Overall_plot <- plot_grid(plot1, plot2, plot3, 
 labels = c("A", "B", "C"), nrow = 1 , align = "h", axis = "tb", 
 label_size = 9,label_x =c(0,0,0))

print(Overall_plot)

ggsave("Overall_violin_plots.png", 
plot = Overall_plot, width = 16, height = 10, dpi = 300)

#hierarchical clustering

dge_filtered_norm_cpm_t <- t(dge_filtered_norm_cpm)

view(dge_filtered_norm_cpm_t)

dge_dist <- dist(dge_filtered_norm_cpm_t)
dge_hclust <- hclust(dge_dist)

plot(dge_hclust, 
main = "Hierarchical Clustering of Samples", 
xlab = "", sub = "", cex = 1 , horiz = FALSE)

# Shorten long sample names

dge_hclust$labels <- abbreviate(dge_hclust$labels, 
minlength = 6)


png("Dendrogram.png", width = 10000, height = 5000, res = 500)

plot(
  as.dendrogram(dge_hclust),
  main = "Hierarchical Clustering of Samples",
  xlab = "Distance",
  ylab = "Samples",
  cex = 0.2,          # sample names
  cex.main = 2,     # title
  cex.lab = 1.5,    # axis labels
  cex.axis = 1,     # axis numbers
  horiz = TRUE
)

dev.off()


#PCA
dge_pca <- prcomp(dge_filtered_norm_cpm_t, 
scale. = FALSE , retx = TRUE)

view(dge_pca)

pca.var <- dge_pca$sdev^2

view(dge_pca$sdev)
pca.var.per <- round(pca.var/sum(pca.var)*100, 1) 

view(pca.var.per)

pca_data <- as.data.frame(dge_pca$x)

view(pca_data)

pca_data$Sample <- rownames(pca_data)

view(pca_data)

pca_plot <-ggplot(pca_data, aes(x = PC1, y = PC2, color = Sample)) +
 geom_point(size = 4) +
 xlab(paste0("PC1 (", pca.var.per[1], "%)")) +
 ylab(paste0("PC2 (", pca.var.per[2], "%)")) +
 labs(title = "PCA of Samples") +
 theme_bw() +
 theme(plot.title = element_text(hjust = 0.5))

ggsave("PCA_plot.png", plot = pca_plot, 
width = 10, height = 6, dpi = 300)
print(pca_plot)

#Differential expression analysis
view(group)

group_factor <- factor(group, levels = c("CTRL", "ERb.KO", "positiveERb", 
"positiveERbVH", "positiveErb.LY5"))

design <- model.matrix(~0 + group)

view(design)

colnames(design) <- levels(group_factor)

view(design)

#Getting mean-variance relationship
dge_filtered_norm_voom <- voom(dge_filtered_norm, design = 
design, plot = TRUE)


#Fitting linear model and empirical Bayes moderation
fit_lmfit <- lmFit(dge_filtered_norm_voom, design)
fit_eBayes <- eBayes(fit_lmfit)

#View results
view(fit_eBayes)

#Contrast matrix for Treatment vs Control
contrast_matrix <- makeContrasts(
ERb.KO_vs_CTRL = ERb.KO - CTRL,
positiveERb_vs_CTRL = positiveERb - CTRL,
positiveERbVH_vs_CTRL = positiveERbVH - CTRL,
positiveErb.LY5_vs_CTRL = positiveErb.LY5 - CTRL,
levels = design
)
view(contrast_matrix)

fit_contrasts <- contrasts.fit(fit_eBayes, contrast_matrix)

fit_contrasts <- eBayes(fit_contrasts)

#Differential expression results

dge_results <- topTable(fit_contrasts, number = Inf, 
adjust.method = "BH", sort.by = "P", coef = "positiveERbVH_vs_CTRL")  

view(dge_results)

#Diffrentially expressed genes with adjusted p-value < 0.01 and log fold change > 1 or < -1 and non-significant genes with adjusted p-value > 0.01 or log fold change between -1 and 1

diff_exp_genes <- dge_results[, c("logFC", "AveExpr", "t", "P.Value", "adj.P.Val", "B")] 

diff_exp_genes$Significant <- "Not Significant"
diff_exp_genes$Significant[diff_exp_genes$logFC > 1 & diff_exp_genes$adj.P.Val < 0.01] <- "Upregulated"
diff_exp_genes$Significant[diff_exp_genes$logFC < -1 & diff_exp_genes$adj.P.Val < 0.01] <- "Downregulated"  


view(diff_exp_genes)

#Extracting Ensembl IDs by removing the version
diff_exp_genes$Ensembl_ID <- sub("\\.\\d+$", "", rownames(diff_exp_genes))

view(diff_exp_genes)

#Converting Ensembl IDs to gene symbols using biomaRt
library(biomaRt)

mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
view(mart)

gene_symbols <- getBM(
  attributes = c("ensembl_gene_id", "hgnc_symbol"),
  filters = "ensembl_gene_id",
  values = diff_exp_genes$Ensembl_ID,
  mart = mart
)

view(gene_symbols)

#Merging gene symbols with differential expression results
diff_exp_genes2 <- merge(diff_exp_genes, gene_symbols, 
by.x = "Ensembl_ID", by.y = "ensembl_gene_id", all.x = TRUE)

view(diff_exp_genes2)

##Replace NA gene symbols with Ensembl IDs
diff_exp_genes2$hgnc_symbol[is.na(diff_exp_genes2$hgnc_symbol)] <- diff_exp_genes2$Ensembl_ID[is.na(diff_exp_genes2$hgnc_symbol)]

view(diff_exp_genes2$Ensembl_ID[is.na(diff_exp_genes2$hgnc_symbol)])

view(diff_exp_genes2$hgnc_symbol[is.na(diff_exp_genes2$hgnc_symbol)])

rownames(diff_exp_genes2)

view(diff_exp_genes2)

#Reorder columns to have gene symbol first
diff_exp_genes_order <- diff_exp_genes2[, c("hgnc_symbol", "logFC", "AveExpr", "t", "P.Value", "adj.P.Val", "B", "Significant", "Ensembl_ID")]

view(diff_exp_genes_order)

 
#Volcano plot with Siginificant and non-significant genes colored differently
diff_exp_genes_order$Significant <- "Not Significant"
diff_exp_genes_order$Significant[diff_exp_genes_order$logFC > 1 & diff_exp_genes_order$adj.P.Val < 0.01] <- "Upregulated"
diff_exp_genes_order$Significant[diff_exp_genes_order$logFC < -1 & diff_exp_genes_order$adj.P.Val < 0.01] <- "Downregulated"



view(diff_exp_genes_order)

#Volcano plot of all and merge them into one plot with significant genes colored differently

volcano_plot_positiveERb_vs_CTRL <- ggplot(diff_exp_genes_order, aes(x = logFC, y = -log10(adj.P.Val), color = Significant)) +
 geom_point(alpha = 0.6) +
 labs(title = "Volcano Plot of Differentially Expressed Genes (positiveERb_vs_CTRL)",
 x = "Log2 Fold Change",
 y = "-Log10 Adjusted P-Value") +
 labs(title = 'Volcano Plot of positiveERb_vs_CTRL', x = 'Log2 Fold Change', y = '-Log10 Adjusted P-Value') +
 theme_bw() +
 theme(plot.title = element_text(hjust = 0.5))

volcano_plot_ERb.KO_Vs_CTRL <- ggplot(diff_exp_genes_order, aes(x = logFC, y = -log10(adj.P.Val), color = Significant)) +
 geom_point(alpha = 0.6) +
 labs(title = "Volcano Plot of Differentially Expressed Genes (ERb.KO vs CTRL)",
 x = "Log2 Fold Change",
 y = "-Log10 Adjusted P-Value") +
 labs(title = 'Volcano Plot of ERb.KO vs CTRL', x = 'Log2 Fold Change', y = '-Log10 Adjusted P-Value') +
 theme_bw() +
 theme(plot.title = element_text(hjust = 0.5))

volcano_plot_positiveERbVH_vs_CTRL <- ggplot(diff_exp_genes_order, aes(x = logFC, y = -log10(adj.P.Val), color = Significant)) +
 geom_point(alpha = 0.6) +
 labs(title = "Volcano Plot of Differentially Expressed Genes (positiveERbVH vs CTRL)",
 x = "Log2 Fold Change",
 y = "-Log10 Adjusted P-Value") +
 labs(title = 'Volcano Plot of positiveERbVH vs CTRL', x = 'Log2 Fold Change', y = '-Log10 Adjusted P-Value') +
 theme_bw() +
 theme(plot.title = element_text(hjust = 0.5))


volcano_plot_positiveERb.LY5_vs_CTRL <- ggplot(diff_exp_genes_order, aes(x = logFC, y = -log10(adj.P.Val), color = Significant)) +
 geom_point(alpha = 0.6) +
 labs(title = "Volcano Plot of Differentially Expressed Genes (positiveERb.LY5 vs CTRL)",
 x = "Log2 Fold Change",
 y = "-Log10 Adjusted P-Value") +
 labs(title = 'Volcano Plot of positiveERb.LY5 vs CTRL', x = 'Log2 Fold Change', y = '-Log10 Adjusted P-Value') +
 theme_bw() +
 theme(plot.title = element_text(hjust = 0.5))


Overall_volcano_plot <- plot_grid(volcano_plot_positiveERb_vs_CTRL, volcano_plot_ERb.KO_Vs_CTRL, volcano_plot_positiveERbVH_vs_CTRL, volcano_plot_positiveERb.LY5_vs_CTRL,
labels = c("A", "B", "C", "D"), nrow = 2 , align = "h", axis = "tb", label_size = 9,label_x =c(0,0,0,0))

print(Overall_volcano_plot)

ggsave("Overall_volcano_plot.png", 
       plot = Overall_volcano_plot, 
       width = 20, height = 15, dpi = 300)


#Gene set enrichment analysis using clusterProfiler

view(diff_exp_genes_order)


library(clusterProfiler)
library(org.Hs.eg.db)


gene_conversion <- bitr(diff_exp_genes_order$Ensembl_ID, fromType = "ENSEMBL", 
toType = "ENTREZID", OrgDb = org.Hs.eg.db)  

view(gene_conversion)

table(duplicated(gene_conversion$ENSEMBL))

table(duplicated(gene_conversion$ENTREZID))

#Remove duplicated ENTREZ IDs to ensure unique mapping for GSEA

gene_conversion_unique <- gene_conversion[!duplicated(gene_conversion$ENSEMBL), ]

# Optional: also remove duplicate Entrez IDs
gene_conversion_unique <- gene_conversion_unique[!duplicated(gene_conversion_unique$ENTREZID), ]


table(duplicated(gene_conversion_unique$ENSEMBL))

table(duplicated(gene_conversion_unique$ENTREZID))




view(gene_conversion_unique)

# Make sure the vector is numeric and named with ENTREZ IDs

gene_list <- diff_exp_genes_order$logFC
names(gene_list) <- diff_exp_genes_order$Ensembl_ID

view(gene_list)
# Keep only mapped genes

gene_list <- gene_list[names(gene_list) %in% gene_conversion_unique$ENSEMBL]

# Replace names with ENTREZ IDs
names(gene_list) <- gene_conversion_unique$ENTREZID[match(names(gene_list), gene_conversion_unique$ENSEMBL)]

view(gene_list)

anyDuplicated(names(gene_list))

# 2) break ties slightly
set.seed(1)
gene_list2 <- gene_list + rnorm(length(gene_list), 0, 1e-10)
gene_list2 <- sort(gene_list2, decreasing = TRUE)

#THIS IS THE IMPORTANT STEP: sort decreasing
gene_list <- sort(gene_list, decreasing = TRUE)

#Perform GO and biological process (BP) enrichment analysis using clusterProfiler


# 3) rerun GSEA GO
gsea_GO_BP <- gseGO(
  geneList      = gene_list2,
  OrgDb         = org.Hs.eg.db,
  keyType       = "ENTREZID",
  ont           = "BP",
  pvalueCutoff  = 0.05,
  pAdjustMethod = "BH",
  nPermSimple   = 10000,
  eps           = 0,
  verbose       = FALSE
)
gsea_GO_results <- as.data.frame(gsea_GO_BP)

view(gsea_GO_results)


gsea_GO_results_final <- clusterProfiler::simplify(gsea_GO_BP)

write_delim( x = as.data.frame(gsea_GO_results_final), file = "gsea_GO_results_final.txt", delim = "\t")

#dotplot
gsea_GO_plot_BP_positiveERbVH_vs_CTRL <- dotplot(gsea_GO_BP, showCategory = 10, title = "GSEA GO Enrichment (BP) positiveERbVH_vs_CTRL") + 
  labs(x = "Gene Ratio", y = "GO Biological Process") +
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
    axis.text.y = element_text(size = 10)
  )



#Molecular function (MF) 

gsea_GO_MF <- gseGO(
  geneList      = gene_list2,
  OrgDb         = org.Hs.eg.db,
  keyType       = "ENTREZID",
  ont           = "MF",
  pvalueCutoff  = 0.05,
  pAdjustMethod = "BH",
  nPermSimple   = 10000,
  eps           = 0,
  verbose       = FALSE
)

gsea_GO_MF_results <- as.data.frame(gsea_GO_MF)

write_delim( x = gsea_GO_MF_results, file = "gsea_GO_MF_results.txt", delim = "\t")

gsea_GO_MF_plot_positiveRbVH_vs_CTRL <- dotplot(gsea_GO_MF, showCategory = 10, title = "GSEA GO Enrichment (MF) positiveRbVH_vs_CTRL") + 
  labs(x = "Gene Ratio", y = "GO Molecular Function") +
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
    axis.text.y = element_text(size = 10)
  )


#Cellular component (CC)

gsea_GO_CC <- gseGO(
  geneList      = gene_list2,
  OrgDb         = org.Hs.eg.db,
  keyType       = "ENTREZID",
  ont           = "CC",
  pvalueCutoff  = 0.05,
  pAdjustMethod = "BH",
  nPermSimple   = 10000,
  eps           = 0,
  verbose       = FALSE
)

gsea_GO_CC_results <- as.data.frame(gsea_GO_CC)

write_delim( x = gsea_GO_CC_results, file = "gsea_GO_CC_results.txt", delim = "\t")

gsea_GO_CC_plot_positiveRbVH_vs_CTRL <- dotplot(gsea_GO_CC, showCategory = 10, title = "GSEA GO Enrichment (CC) positiveERbVH_vs_CTRL") + 
  labs(x = "Gene Ratio", y = "GO Cellular Component") +
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
    axis.text.y = element_text(size = 10)
  )


Overall_GSEA_GO_positiveRbVH_vs_CTRL <- plot_grid(gsea_GO_plot_BP_positiveERbVH_vs_CTRL , gsea_GO_MF_plot_positiveRbVH_vs_CTRL, gsea_GO_CC_plot_positiveRbVH_vs_CTRL,
labels = c("A", "B", "C"), nrow = 1 , align = "h", axis = "tb", label_size = 9,label_x =c(0,0,0))


ggsave("Overall_GSEA_GO_positiveRbVH_vs_CTRL.png", 
       plot = Overall_GSEA_GO_positiveRbVH_vs_CTRL, 
       width = 20, height = 10, dpi = 300)


#KEGG pathways

gsea_Kegg <- gseKEGG(
  geneList      = gene_list2,
  organism      = "hsa",
  minGSSize     = 5,
  maxGSSize     = 500,
  pvalueCutoff  = 0.1,
  eps           = 0,
  nPermSimple   = 10000,
  verbose       = FALSE
)


view(gsea_Kegg)



Plot top pathway
# -------------------------------------------------
dotplot(gsea_Kegg, showCategory = 10)

Gsea_KEGG_positiveERbVH_vs_CTRL_dotplot <- dotplot(gsea_Kegg, showCategory = 10) + 
  labs(title = "GSEA KEGG Pathway Enrichment positiveRbVH_vs_CTRL", x = "Gene Ratio", y = "KEGG Pathway") +
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
    axis.text.y = element_text(size = 10)
  )
dev.off()

ggsave("Gsea_KEGG_positiveRbVH_vs_CTRL_dotplot.png", 
       plot = Gsea_KEGG_positiveERbVH_vs_CTRL_dotplot, 
       width = 10, height = 8, dpi = 300)


#Activated ans repressed pathways
gsea_Kegg_results <- as.data.frame(gsea_Kegg)
gsea_Kegg_results$Direction <- ifelse(gsea_Kegg_results$NES > 0, "Activated", "Repressed")
gsea_Kegg_positiveRbVH_vs_CTRL_dotplot_Activated_and_Suppresed  <- ggplot(gsea_Kegg_results, aes(x = reorder(Description, NES), y = NES, fill = Direction)) +
  geom_col() +
  coord_flip() +
  labs(title = "GSEA KEGG Pathway Enrichment positiveRbVH_vs_CTRL", x = "KEGG Pathway", y = "Normalized Enrichment Score (NES)") +
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
    axis.text.y = element_text(size = 10)
  )

ggsave("gsea_Kegg_positiveRbVH_vs_CTRL_dotplot_Activated_and_Suppresed.png", 
       plot = gsea_Kegg_positiveRbVH_vs_CTRL_dotplot_Activated_and_Suppresed, 
       width = 10, height = 8, dpi = 300)

Gsea_KEGG_positiveRbVH_vs_CTRL_overall <- plot_grid(Gsea_KEGG_positiveERbVH_vs_CTRL_dotplot, gsea_Kegg_positiveRbVH_vs_CTRL_dotplot_Activated_and_Suppresed,
labels = c("A", "B"), nrow = 1 , align = "h", axis = "tb", label_size = 9,label_x =c(0,0))

print(Gsea_KEGG_positiveRbVH_vs_CTRL_overall)

ggsave("Gsea_KEGG_positiveRbVH_vs_CTRL_overall.png", 
       plot = Gsea_KEGG_positiveRbVH_vs_CTRL_overall, 
       width = 20, height = 10, dpi = 300)








 # Reactome pathways
gene_list2 <- diff_exp_genes_significant$logFC
names(gene_list2) <- diff_exp_genes_significant$Ensembl_ID
gene_list2 <- gene_list2[names(gene_list2) %in% gene_conversion_unique$ENSEMBL]
names(gene_list2) <- gene_conversion_unique$ENTREZID[match(names(gene_list2), gene_conversion_unique$ENSEMBL)]
gene_list2 <- sort(gene_list2, decreasing = TRUE)

#Get Reactome pathways 
install.packages("ReactomePA")
library(ReactomePA)
reactome_pathways <- enrichPathway(gene = names(gene_list2), 
                                    organism = "human", 
                                    pvalueCutoff = 0.01, 
                                    qvalueCutoff = 0.05, 
                                    readable = TRUE)

reactome_ORA_ERb.KO_Vs_CTRL <- dotplot(reactome_pathways, showCategory = 10) + 
  labs(title = "Reactome Pathway Enrichment ERb.KO vs CTRL", x = "Gene Ratio", y = "Reactome Pathway") +
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
    axis.text.y = element_text(size = 10)
  )


ggsave("Reactome_ORA_ERB.KO_Vs_CTRL.png",
       plot = reactome_ORA_ERb.KO_Vs_CTRL,
       width = 10,
       height = 10,
       dpi = 300)


#Reactome pathway analysis using fgsea for a more comprehensive GSEA approach



library(msigdbr)
library(dplyr)
library(fgsea)
reactome_df <- msigdbr(species = "Homo sapiens",
                       category = "C2",
                       subcategory = "CP:REACTOME")

reactome_pathways <- reactome_df %>%
  split(x = .$entrez_gene, f = .$gs_name)


fsgea_reactome <- fgseaMultilevel(
  pathways = reactome_pathways ,
  stats = gene_list2 ,
  minSize = 5,
  maxSize = 500)

fgsea_reactome <- fgsea_reactome[order(fgsea_reactome$pval), ]
head(fgsea_reactome[, c("pathway","NES","pval","padj")])

topPathways <- fgsea_reactome %>%
  arrange(pval) %>%
  head(10)

reactome_GSEA_ERB.KO_Vs_CTRL <- ggplot(topPathways, aes(x = reorder(pathway, NES), y = NES)) +
  geom_col(aes(fill = pval)) +
  coord_flip() +
  labs(title = "Reactome GSEA – Top Pathways of ERb.KO Vs CTRL",
       x = "",
       y = "Normalized Enrichment Score (NES)") +
  theme_bw()

#How to save
ggsave("reactome_GSEA_ERB.KO_Vs_CTRL.png", 
plot = reactome_GSEA_ERB.KO_Vs_CTRL, 
width = 10, height = 13, dpi = 300)

# 