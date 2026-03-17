pwd()
setwd("C:\\Users\\Stuart\\Downloads")
setwd(dir)
dir.exists("C:/Users/Stuart/Downloads/raw counts.xlsx")
dir.exists("C:/Users/Stuart/Downloads/raw_counts")
install.packages("readxl")
library(readxl)
setwd("C://Users//Stuart//Downloads")
raw_counts <- read_excel("raw counts.xlsx")
head(raw_counts)
library("dplyr")          ## loads in dplyr package to use
library("tidyr")          ## loads in tidyr package to use
library("ggplot2")          ## loads in ggplot2 package to use
library("readr")    
library(edgeR)
library(tibble)
library(cowplot)
library(readxl)
library(limma)


if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install(c("DESeq2", "edgeR"))
install.packages(c("ggplot2", "dplyr"))



view(raw_counts)

raw_counts$Control1 <- raw_counts$BT549_CT2
raw_counts$Control2 <- raw_counts$BT540_CT3
raw_counts$Control3 <- raw_counts$BT549_CT4
raw_counts$Control4 <- raw_counts$HS578T_CT3
raw_counts$Control5 <- raw_counts$HS578T_CT4
raw_counts$Control6 <- raw_counts$HS578T_CT5
raw_counts$Treatment1 <- raw_counts$BT540_KO3
raw_counts$Treatment2 <- raw_counts$BT549_KO8
raw_counts$Treatment3 <- raw_counts$BT549_KO11
raw_counts$Treatment4 <- raw_counts$HS578T_KO3
raw_counts$Treatment5 <- raw_counts$HS578T_KO7
raw_counts$Treatment6 <- raw_counts$HS578T_KO8

raw_counts <- raw_counts %>%
  select(Geneid, Control1, Control2, Control3, Control4, Control5, Control6,
         Treatment1, Treatment2, Treatment3, Treatment4, Treatment5, Treatment6)

raw_counts_data.frame <- as.data.frame(raw_counts)

rownames(raw_counts_data.frame) <- raw_counts_data.frame$Geneid

row_names <-rownames(raw_counts_data.frame)

view(raw_counts_data.frame)

raw_counts2 <- raw_counts_data.frame[,-1]

view(raw_counts2)



group <-c("Control", "Control", "Control", "Control", "Control", "Control",
         "Treatment", "Treatment", "Treatment", "Treatment", "Treatment", "Treatment")


str(raw_counts2)

raw_counts_numeric <- as.data.frame(sapply(raw_counts2, as.numeric))

str(raw_counts_numeric)

view(raw_counts_numeric)


rownames(raw_counts_numeric) <- NULL

rownames(raw_counts_numeric) <- row_names

view(raw_counts_numeric)

dge <- DGEList(counts = raw_counts_numeric, group = group)

rownames(dge$counts)



view(dge)


dge_logCPM <- cpm(dge, log = TRUE)

view(dge_logCPM)

cpm <- cpm(dge)
colSums(cpm) 


view(dge_logCPM)


dge_logCPM_row <- as_tibble(dge_logCPM, rownames = "Geneid")


view(dge_logCPM_row)

PT1 <-pivot_longer(dge_logCPM_row, cols = -Geneid, names_to = "Sample", values_to = "logCPM")

view(PT1)

plot1 <- ggplot(PT1,aes(x = Sample , y = logCPM , fill = Sample)) +
geom_violin(trim = FALSE, show.legend = TRUE) +
stat_summary(fun = "median", geom = "point", size = 2, color = "black", shape = 95) +
labs(title = "Violin Plot of logCPM Values by Sample", x = "Sample", y = "logCPM", subtitle = "Distribution of logCPM values across samples") +
theme_bw() +
theme(
plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
plot.subtitle = element_text(hjust = 0.5, size = 12),
) +
coord_flip()

view(plot1)
print(plot1)

ggsave("plot1.png", plot = plot1, width = 10, height = 6, dpi = 300)

# Add this line to display the plot
print(plot1)


#filtering

cpm <- cpm(dge)
threshold <-rowSums(cpm > 1) >= 2
view(threshold)

DGEList_filtered <- dge[threshold,]

view(DGEList_filtered)

DGEList_filtered_2 <- cpm(DGEList_filtered, log = TRUE)

view(DGEList_filtered_2)

DGEList_filtred_row <- as_tibble(DGEList_filtered_2, rownames = "Geneid")
view(DGEList_filtred_row)

PT2 <-pivot_longer(DGEList_filtred_row, cols = -Geneid, names_to = "Sample", values_to = "logCPM")

view(PT2)
plot2 <- ggplot(PT2,aes(x = Sample , y = logCPM , fill = Sample)) +
  geom_violin(trim = FALSE, show.legend = TRUE) +
  stat_summary(fun = "median", geom = "point", size = 2, color = "black", shape = 95) +
  labs(title = "Violin Plot of logCPM Values by Sample (Filtered)", x = "Sample", y = "logCPM", subtitle = "Distribution of logCPM values across samples after filtering") +
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5, size = 12),
  ) +
  coord_flip()
view(plot2)
print(plot2)
ggsave("plot2.png", plot = plot2, width = 10, height = 6, dpi = 300)


#Normalization and filtering

Norm.filtered <- calcNormFactors(DGEList_filtered, method = "TMM")

view(Norm.filtered)

Norm.filtered_logCPM <- cpm(Norm.filtered, log = TRUE)

Norm.filtered_logCPM_row <- as_tibble(Norm.filtered_logCPM, rownames = "Geneid")

view(Norm.filtered_logCPM_row)

PT3 <-pivot_longer(Norm.filtered_logCPM_row, cols = -Geneid, names_to = "Sample", values_to = "logCPM")

view(PT3)
plot3 <- ggplot(PT3,aes(x = Sample , y = logCPM , fill = Sample)) +
  geom_violin(trim = FALSE, show.legend = TRUE) +
  stat_summary(fun = "median", geom = "point", size = 2, color = "black", shape = 95) +
  labs(title = "Violin Plot of logCPM Values by Sample (Normalized and Filtered)", x = "Sample", y = "logCPM", subtitle = "Distribution of logCPM values across samples after normalization and filtering") +
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5, size = 12),
  ) +
  coord_flip()
view(plot3)
print(plot3)
ggsave("plot3.png", plot = plot3, width = 10, height = 6, dpi = 300)


norm_plot <- plot_grid(plot1, plot2, plot3, labels = c("A", "B", "C"), ncol = 1)
print(norm_plot)
ggsave("norm_plot.png", plot = norm_plot, width = 12, height = 18, dpi = 300)

#Hierarchical clustering

distance <- dist(t(Norm.filtered_logCPM), method = "euclidean")
clusters <- hclust(distance, method = 'complete')
view(clusters)


view(distance)

plot4 <- plot(clusters, main = "Hierarchical Clustering of Samples", xlab = "Samples", ylab = "Distance", sub = "", cex = 0.8, 
labels = colnames(Norm.filtered_logCPM), hang = -1)

view(colnames(Norm.filtered_logCPM))
ggsave("plot4.png", plot = plot4, width = 10, height = 6, dpi = 300)

png("plot4_dendrogram.png", width = 1200, height = 600)
plot4_dendrogram <- plot(
  as.dendrogram(clusters),
  main = "Hierarchical Clustering of Samples",
  xlab = "Distance",
  ylab = "Samples",
  cex = 0.8,
  horiz = TRUE
)

dev.off()
ggsave("plot4_dendrogram.png", plot = plot4_dendrogram, 
width = 15, height = 6, dpi = 300)

print(plot4_dendrogram)
#PCA
pca.res <- prcomp(t(Norm.filtered_logCPM), scale. = FALSE, retx = TRUE)

view(pca.res)
pca.res$x     # PCA scores (coordinates of each sample in new PC space)
pca.res$sdev   # standard deviations of PCs (how much variation each PC explains)
pca.res$rotation  # loadings (which genes contribute to which PCs)

view(pca.res$sdev)

view(pca.res$rotation)

view(pca.res$x)

pca.var <- pca.res$sdev^2

pca.var.per <- round(pca.var/sum(pca.var)*100, 1)

pca.res.df <- as_tibble(pca.res$x)

view(pca.res$x)

print()


library(tibble)
library(ggplot2)

# ------------------------------
# 1️⃣ PCA (already calculated)
# pca.res <- prcomp(t(Norm.filtered_logCPM), scale. = FALSE, retx = TRUE)
# ------------------------------

# 2️⃣ Calculate % variance explained
pca.var <- pca.res$sdev^2
pca.var.per <- round(pca.var / sum(pca.var) * 100, 1)  # percentage

# 3️⃣ Convert PCA scores to tibble
pca.res.df <- as_tibble(pca.res$x)

# 4️⃣ Add sample labels manually
sampleLabels <- c(
  "Control1", "Control2", "Control3", "Control4", "Control5", "Control6",
  "Treatment1", "Treatment2", "Treatment3", "Treatment4", "Treatment5", "Treatment6"
)

# 5️⃣ Add group info
group <- c(
  "Control", "Control", "Control", "Control", "Control", "Control",
  "Treatment", "Treatment", "Treatment", "Treatment", "Treatment", "Treatment"
)

# 6️⃣ Add to PCA data frame
pca.res.df$sampleLabels <- sampleLabels
pca.res.df$group <- group

view(pca.res.df)

png("PCA.png", width = 800, height = 600)
ggplot(pca.res.df) +
  aes(x = PC1, y = PC2, color = group) +   # no label here
  geom_point(size = 6) +                     # only points
  xlab(paste0("PC1 (", pca.var.per[1], "%)")) +
  ylab(paste0("PC2 (", pca.var.per[2], "%)")) +
  labs(title = "PCA of Samples") +
  theme_bw()
dev.off()

#Differential expression analysis

group_DGE <- c("Control", "Control", "Control", "Control", "Control", "Control",
               "Treatment", "Treatment", "Treatment", "Treatment", "Treatment", "Treatment")

group_factor <- factor(group_DGE, levels = c("Control", "Treatment"))

view(group_factor)

design_1 <- model.matrix(~ 0 + group_factor)

view(design_1)

design_2 <-design_1

colnames(design_1) <- c("Control", "Treatment")

view(design_1)

#Getting mean variance relationship with voom

voom_dge <- voom(Norm.filtered, design_1, plot = TRUE)

voom_dge <- voom(Norm.filtered, design_1, plot = FALSE)

view(voom_dge)

view(Norm.filtered)
#how to save voom file

write.csv(voom_dge, file = "voom_dge.csv", row.names = TRUE)
#How to save voom_dge plot

png("voom_dge.png", width = 800, height = 600)
plot(voom_dge$plot)
dev.off()

ggsave("voom_dge.png", plot = voom_dge$plot, width = 10, 
height = 6, dpi = 300)


#Fitting linear model and empirical Bayes moderation
fit <- lmFit(voom_dge, design_1)

view(fit)

write.csv(fit, file = "fit.csv", row.names = TRUE)

#Contrast matrix for Treatment vs Control

contrast_matrix <- makeContrasts(Treatment - Control, levels = design_1)

view(contrast_matrix)

write.csv(contrast_matrix, file = "contrast_matrix.csv", row.names = TRUE)

#Extracting linear model fit for the contrast

fit2 <- contrasts.fit(fit, contrast_matrix)

view(fit2)

write.csv(fit2, file = "fit2.csv", row.names = TRUE)
#Empirical Bayes moderation
fit2_bayes <- eBayes(fit2)

view(fit2_bayes)

#Differential expression results
results <- topTable(fit2_bayes, adjust.method = "BH", number = Inf ,
sort.by = "P")

view(results)

write.csv(results, file = "differential_expression_results.csv", 
row.names = TRUE)

view(results)

#Diffrentially expressed genes with adjusted p-value < 0.01

# differentially expressed genes with log fold change > 1 or < -1 and adjusted p-value < 0.01
diff_exp_genes <- results[(results$adj.P.Val < 1 & results$logFC > 1) | 
                            (results$adj.P.Val < 1 & results$logFC < -1),
                            c("logFC", "AveExpr", "t", "P.Value", "adj.P.Val", "B") ]


view(diff_exp_genes)

#Extracting Ensembl IDs by removing the version
 ensembl_ids <- sub("\\.\\d+$", "", rownames(diff_exp_genes))
 
 (view(ensembl_ids))

diff_exp_genes$Ensembl_ID <- ensembl_ids

view(diff_exp_genes)

#Converting Ensembl IDs to gene symbols using biomaRt
if (!require("biomaRt", quietly = TRUE))
    install.packages("biomaRt")

library(biomaRt)
mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
gene_info <- getBM(attributes = c("ensembl_gene_id", "hgnc_symbol"),
                   filters = "ensembl_gene_id",
                   values = diff_exp_genes$Ensembl_ID,
                   mart = mart)
view(gene_info)


merge diff_exp_genes with gene_info to add gene symbols
diff_exp_genes_annotated <- merge(diff_exp_genes, gene_info, by.x = "Ensembl_ID", by.y = "ensembl_gene_id", all.x = TRUE)
view(diff_exp_genes_annotated)

#Replace NA gene symbols with Ensembl IDs
diff_exp_genes_annotated$hgnc_symbol[is.na(diff_exp_genes_annotated$hgnc_symbol)] <- diff_exp_genes_annotated$Ensembl_ID[is.na(diff_exp_genes_annotated$hgnc_symbol)]
view(diff_exp_genes_annotated)

#Reorder columns to have gene symbol first
diff_exp_genes_annotated <- diff_exp_genes_annotated[, c("hgnc_symbol", "logFC", "AveExpr", "t", "P.Value", "adj.P.Val", "B", "Ensembl_ID")]

view(diff_exp_genes_annotated)

#Extract upregulated and downregulated genes
upregulated_genes <- diff_exp_genes_annotated[diff_exp_genes_annotated$logFC > 1 & diff_exp_genes_annotated$adj.P.Val < 1, ]
downregulated_genes <- diff_exp_genes_annotated[diff_exp_genes_annotated$logFC < -1 & diff_exp_genes_annotated$adj.P.Val < 1, ]

#Save results to CSV files
write.csv(diff_exp_genes_annotated, file = "differentially_expressed_genes_annotated.csv", row.names = FALSE)
write.csv(upregulated_genes, file = "upregulated_genes.csv", row.names = FALSE)   
write.csv(downregulated_genes, file = "downregulated_genes.csv", row.names = FALSE)

#Volcano plot of differential expression results
results$Significant <- ifelse(results$adj.P.Val < 1 & results$logFC > 1, "Upregulated",
                             ifelse(results$adj.P.Val < 1 & results$logFC < -1, "Downregulated", "Not Significant"))
view(results)

volcano_plot <- ggplot(results, aes(x = logFC, y = -log10(adj.P.Val), color = Significant)) +
  geom_point(alpha = 0.6) +
  scale_color_manual(values = c("Upregulated" = "red", "Downregulated" = "blue", "Not Significant" = "grey")) +
  labs(title = "Volcano Plot of Differential Expression Results", x = "Log2 Fold Change", y = "-Log10 Adjusted P-Value") +
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
    legend.title = element_blank(),
    legend.position = "right"
  ) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "black") +
  geom_hline(yintercept = -log10(1), linetype = "dashed", color = "black")

view(volcano_plot)
print(volcano_plot)


#Save volcano plot
ggsave("volcano_plot.png", plot = volcano_plot, 
width = 10, height = 6, dpi = 300)

#Gene set enrichment analysis using clusterProfiler
if (!require("clusterProfiler", quietly = TRUE))
    install.packages("clusterProfiler")

library(clusterProfiler)
if (!require("org.Hs.eg.db", quietly = TRUE))
    install.packages("org.Hs.eg.db")

if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("org.Hs.eg.db")
library(org.Hs.eg.db)
#Prepare gene list for GSEA
gene_list <- diff_exp_genes_annotated$logFC
names(gene_list) <- diff_exp_genes_annotated$Ensembl_ID
gene_list <- sort(gene_list, decreasing = TRUE)



view(gsea_results)
head(diff_exp_genes_annotated$Ensembl_ID)

gene_conversion <- bitr(diff_exp_genes_annotated$Ensembl_ID, fromType = "ENSEMBL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)   

head(diff_exp_genes_annotated$Ensembl_ID)

diff_exp_genes_annotated$Ensembl_ID <- sub("\\..*", "",
 diff_exp_genes_annotated$Ensembl_ID)

gene_conversion <- bitr(diff_exp_genes_annotated$Ensembl_ID, 
fromType = "ENSEMBL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)

dim(gene_conversion)

# Load libraries
library(clusterProfiler)
library(org.Hs.eg.db)
library(dplyr)

# -------------------------------------------------
# 1️⃣ Remove Ensembl version numbers (IMPORTANT)
# -------------------------------------------------
diff_exp_genes_annotated$Ensembl_ID <- 
  sub("\\..*", "", diff_exp_genes_annotated$Ensembl_ID)

# -------------------------------------------------
# 2️⃣ Convert Ensembl → Entrez
# -------------------------------------------------
gene_conversion <- bitr(
  diff_exp_genes_annotated$Ensembl_ID,
  fromType = "ENSEMBL",
  toType = "ENTREZID",
  OrgDb = org.Hs.eg.db
)

# -------------------------------------------------
# 3️⃣ Merge with original data to keep logFC
# -------------------------------------------------
merged_data <- merge(
  gene_conversion,
  diff_exp_genes_annotated,
  by.x = "ENSEMBL",
  by.y = "Ensembl_ID"
)

merged_data_unique <- merged_data %>%
  group_by(ENTREZID) %>%
  summarise(logFC = max(logFC, na.rm = TRUE)) %>%  # keep max logFC per Entrez
  ungroup()
# -------------------------------------------------

# -------------------------------------------------
# 5️⃣ Create ranked gene list (VERY IMPORTANT)
# -------------------------------------------------
gene_list <- merged_data_unique$logFC
names(gene_list) <- merged_data_unique$ENTREZID

# Sort decreasing for GSEA
gene_list <- sort(gene_list, decreasing = TRUE)

# -------------------------------------------------
# 6️⃣ Run GSEA KEGG
# -------------------------------------------------
gsea_results <- gseKEGG(
  geneList = gene_list,
  organism = "hsa",
  pvalueCutoff = 1,
  verbose = FALSE
)

# -------------------------------------------------
# 7️⃣ View results
# -------------------------------------------------
head(as.data.frame(gsea_results))

# -------------------------------------------------
# 8️⃣ Plot top pathway
# -------------------------------------------------
dotplot(gsea_results, showCategory = 10)

png("gsea_dotplot.png", width = 800, height = 600)
dotplot(gsea_results, showCategory = 10) + 
  labs(title = "GSEA KEGG Pathway Enrichment", x = "Gene Ratio", y = "KEGG Pathway") +
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
    axis.text.y = element_text(size = 10)
  )
dev.off()
