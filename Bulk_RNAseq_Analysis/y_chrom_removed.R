#y-linked genes removed pipeline

y_chrom <- read.csv("y_chrom_genes.txt", header = TRUE)

sBAT_filtered_1 <- sBAT_filtered
sBAT_filtered_1$ID <- rownames(sBAT_filtered_1)
rownames(sBAT_filtered_1) <- NULL

#nrow = 22552; took out 8 rows
sBAT_filtered_y <- sBAT_filtered_1[!sBAT_filtered_1$ID %in% y_chrom$Gene.stable.ID, ]
rownames(sBAT_filtered_y) <- sBAT_filtered_y$ID
sBAT_filtered_y$ID <- NULL

#deseq2
sBAT_dds_y <- DESeqDataSetFromMatrix(countData = sBAT_filtered_y, colData = sBAT_column_data, design = ~condition)
sBAT_dds_y$condition <- relevel(sBAT_dds_y$condition, ref = "WT")
sBAT_dds_analysis_y <- DESeq(sBAT_dds_y)
sBAT_dds_result_y <- results(sBAT_dds_analysis_y)
sBAT_deseq_y <- as.data.frame(sBAT_dds_result_y)

#edgeR
sBAT_er_data_y <- DGEList(counts=sBAT_filtered_y, genes=rownames(sBAT_filtered_y), group=sBAT_condition_factor)
sBAT_er_data_y <- calcNormFactors(sBAT_er_data_y)
sBAT_er_data_y <- estimateDisp(sBAT_er_data_y)
sBAT_er_y <- exactTest(sBAT_er_data_y)
sBAT_er_y <- as.data.frame(topTags(sBAT_er_y, n=nrow(sBAT_filtered_y)))
#number of genes with FDR less than 0.05
table(sBAT_er_y$FDR < 0.05)

#top genes list < 0.05
sBAT_er_y












































