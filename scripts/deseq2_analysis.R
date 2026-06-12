#===============================================================================
#                       INSTALLING THE NECESSARY PACKAGES
#===============================================================================

install.packages("BiocManager")

BiocManager::install(c("GenomicFeatures"))


BiocManager::install("DESeq2")
library(DESeq2)

#===============================================================================
#                             IMPORTING DATA
#===============================================================================
# Reading counts and setting the first column ('Geneid') as row names
cnt <- read.csv("counts.csv", row.names = 1)

# Reading the metadata 
met <- read.csv("metadata.csv")

# Assigning the SampleID column directly to the row names of met
rownames(met) <- met$SampleID

# Round the counts to the nearest whole integer to handle any fractional 
# mapping values generated during the HISAT2 alignment pipeline
cnt_rounded <- round(cnt)

# Converting the experimental 'Condition' column to a factor, and setting 
# 'wildType' as the baseline
met$Condition <- factor(met$Condition, levels = c("WildType", "Mutant", "Complemented"))

# Making sure that the colnames in cnt matches the rownames in met
all(colnames(cnt) %in% rownames(met))

# Making sure that both colnames in cnt and rownames in met align 
all(colnames(cnt) == rownames(met))

#===============================================================================
#                   DIFFERENTIAL EXPRESSION OF GENE ANALYSIS
#===============================================================================
# Building DESeq2 Dataset
dds <- DESeqDataSetFromMatrix(countData = cnt_rounded,
                              colData = met,
                              design = ~Condition)

# Checking whether dataset has been generated
dds

# Removal of low count reads
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]
dds

# Setting reference for DEG Analysis to be doubly sure
dds$Condition <- relevel(dds$Condition, ref = "WildType")

# Performing the differential analysis of the gene
deg <- DESeq(dds)

# Extracting differential expression results for Mutant vs. WildType baseline
# Syntax: contrast = c("Factor_Name", "Test_Condition", "Reference_Baseline")
res <- results(deg, contrast = c("Condition", "Mutant", "WildType"))

# Saving results in csv format
write.csv(res, "DESeq2_results_Mutant_vs_WildType.csv")

# Summary statistics of results
summary(res)

# Extracting Mutant vs WildType results optimized for an adjusted p-value (FDR) < 0.05
res0.05 <- results(deg, contrast = c("Condition", "Mutant", "WildType"), alpha = 0.05)

# View the summary 
summary(res0.05)

# Filtering for upregulated genes with p adjusted value of 0.05 and LFC > 1
res_up <- subset(res,padj < 0.05 & log2FoldChange > 1)
summary(res_up)

# Filtering for upregulated genes with p adjusted value of 0.05 and LFC > 1
res_down <- subset(res,padj < 0.05 & log2FoldChange < -1)
summary(res_down)

# Converting res0.05 to dataframe
res0.05.df <- as.data.frame(res0.05)
str(res0.05.df)

#===============================================================================
#                     GENE FILTERING AND ANNOTATION
#===============================================================================
# Converting ensembl gene Ids of S aureus to gene name
# Loading the package and reading the local GFF file
library(rtracklayer)
gff_table <- as.data.frame(import("saureus.gff"))

# Extracting column 13 (Gene Symbol) and column 10 (Locus ID) by name
mapping_key <- unique(gff_table[gff_table$type == "gene", c("Name", "ID")])

# Merging to append the systematic ID column cleanly to the results
final_data <- merge(as.data.frame(res0.05), mapping_key, by.x = "row.names", by.y = "Name", all.x = TRUE)

# Adding and renaming columns 
colnames(final_data)[colnames(final_data) == "Row.names"] <- "Gene_Symbol"
colnames(final_data)[colnames(final_data) == "ID"] <- "Locus_ID"

# Viewing the final table
head(final_data)

# Saving final table dataframe into CSV format
write.csv(final_data, "final_DESeq2_results_Mutant_vs_WildType.csv", row.names = FALSE)

#===============================================================================
#                                 VISUALIZATION
#===============================================================================
install.packages("dplyr")
install.packages("ggplot2")

library(dplyr)
library(ggplot2)
library(DESeq2)

# Generating the PCA plot
vsd <- vst(deg, blind = FALSE)
plotPCA(vsd, intgroup = "Condition") + labs(color = "Condition")

# Estimating size factors
sizeFactors(deg)

# Calculating the Dispersion Estimates
plotDispEsts(deg)

# Building MA plot
plotMA(res0.05)

# Filtering for the top 10 most significant genes
top_genes <- final_data%>%
  arrange(padj)%>%
  head(10)
top_genes
write.csv(top_genes, "top_10_genes.csv", row.names = FALSE)

# Generating Volcano Plot
vol <- final_data%>%
  filter(!is.na(padj))
ggplot(vol, aes(x = log2FoldChange, y = -log10(padj), color = padj < 0.05 & abs(log2FoldChange) > 1)) + 
  geom_point() + 
  labs(
    title = "Differential Gene Expression: Mutant vs Wild-Type",
    x = "Log2 Fold Change",
    y = "-Log10 Adjusted P-value",
    color = "Significance Status"
  )

# Visualizing top Differential Gene Expression via annotated Volcano Plot
install.packages("ggrepel")
library(ggrepel)

ggplot(vol, aes(x = log2FoldChange, y = -log10(padj), color = padj < 0.05 & abs(log2FoldChange) > 1)) + 
  geom_point(alpha = 0.6) +
  
  geom_text_repel(data = top_genes, aes(label = Gene_Symbol), size = 3, max.overlaps = 10) +
  
  xlim(-7.5, 5) + 
 
  labs(title = "Mutant vs Wild-Type", x = "Log2 Fold Change", y = "-Log10 Padj", color = "Significant")


# Building the Heat Map visualization for the  top 5 up-regulated genes and top 5 down-regulated genes
BiocManager::install("ComplexHeatmap")
library(ComplexHeatmap)
library(circlize) # Used for color mapping

# Filtering for the top 5 Upregulated and top 5 Downregulated genes 
top5_up <- final_data %>% filter(padj < 0.05 & log2FoldChange > 1) %>% arrange(padj) %>% head(5)
top5_down <- final_data %>% filter(padj < 0.05 & log2FoldChange < -1) %>% arrange(padj) %>% head(5)

# Combining the top genes
top10_genes <- rbind(top5_up, top5_down)

# Extracting the normalized counts using the Gene_Symbol mapping
mat <- counts(deg, normalized = TRUE)[top10_genes$Gene_Symbol, ]

# # Getting the z score
mat_scaled <- t(scale(t(mat)))

# Generating the heatmap
Heatmap(
  mat_scaled, 
  name = "Z-Score", 
  column_title = "S. aureus Transcriptomic Expression Profiles",
  row_title = "Top Up/Down-Regulated Genes",
  col = colorRamp2(c(-2, 0, 2), c("#313695", "#ffffbf", "#a50026")), # Classic Blue-Yellow-Red palette
  show_row_dend = TRUE, 
  show_column_dend = TRUE
)