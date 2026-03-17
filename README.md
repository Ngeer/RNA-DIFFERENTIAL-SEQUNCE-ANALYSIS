RNA-Seq Analysis of ERβ Signaling in Inflammatory Breast Cancer
📌 Project Overview

This project investigates the role of Estrogen Receptor β (ERβ) in aggressive breast cancer using RNA-seq data.
The goal is to identify differentially expressed genes and understand how ERβ influences gene regulation and tumor behavior.

🎯 Objectives

Perform quality control on raw RNA-seq data

Align reads to the reference genome

Quantify gene expression

Identify differentially expressed genes (DEGs)

Explore biological patterns using visualization techniques

🧪 Experimental Design

Samples were divided into the following groups:

CTRL (control)

ERβ Knockout (ERβ KO)

ERβ Positive

ERβ Positive + Treatment(s)

Each group contains biological replicates.

⚙️ Workflow Overview

The analysis pipeline includes:

Quality Control

FastQC

MultiQC

Read Trimming

Trim Galore

Alignment

STAR aligner

Quantification

featureCounts

Differential Expression Analysis

edgeR / limma-voom

📊 Data Analysis

Normalization using TMM (Trimmed Mean of M-values)

Filtering of lowly expressed genes

Statistical testing for DEGs

📈 Visualization

PCA → sample clustering

Hierarchical clustering (dendrogram) → similarity between samples

Volcano plots → significant DEGs

Violin plots → distribution of expression

🔍 Key Findings

ERβ influences gene expression linked to tumor progression

Identification of both upregulated and downregulated genes

Evidence of repressive regulatory roles of ERβ

🛠️ Tools & Technologies

Nextflow (workflow management)

R (edgeR, limma, ggplot2)

STAR

FastQC / MultiQC

featureCounts
