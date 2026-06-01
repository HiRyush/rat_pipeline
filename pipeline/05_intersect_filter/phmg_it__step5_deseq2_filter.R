#!/usr/bin/env Rscript
# Step 5: DESeq2 Expression-aware Filter
# Identify genes newly expressed in Treatment → flag variants in those genes

suppressMessages(library(methods))
if (!requireNamespace("DESeq2", quietly = TRUE)) {
    if (!requireNamespace("BiocManager", quietly = TRUE))
        install.packages("BiocManager", repos="https://cloud.r-project.org")
    BiocManager::install("DESeq2", ask = FALSE)
}

library(DESeq2)

cat("================================================================\n")
cat("Step 5: DESeq2 Expression-aware Filter\n")
cat("================================================================\n\n")

# 1. Read featureCounts output
counts_file <- "/home/yusanghyeon/RAT_project/pipeline_candidates/results/deseq2/gene_counts.txt"
counts_raw <- read.table(counts_file, header = TRUE, row.names = 1, skip = 1, sep = "\t")

# Extract count columns (columns 6 onward = sample BAM counts)
count_matrix <- counts_raw[, 6:ncol(counts_raw)]
# Clean column names
colnames(count_matrix) <- gsub(".*markdup\\.", "", colnames(count_matrix))
colnames(count_matrix) <- gsub("\\.dedup\\.bam", "", colnames(count_matrix))

cat("Samples:", paste(colnames(count_matrix), collapse = ", "), "\n")
cat("Genes:", nrow(count_matrix), "\n\n")

# 2. Sample metadata
condition <- factor(c(rep("Control", 5), rep("Treatment", 10)),
                    levels = c("Control", "Treatment"))
coldata <- data.frame(condition = condition, row.names = colnames(count_matrix))

# 3. Run DESeq2
dds <- DESeqDataSetFromMatrix(countData = round(count_matrix),
                               colData = coldata,
                               design = ~ condition)

# Pre-filter: remove genes with total count < 10
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep, ]
cat("Genes after filtering (total count >= 10):", nrow(dds), "\n")

dds <- DESeq(dds)
res <- results(dds, contrast = c("condition", "Treatment", "Control"))

# 4. Identify "newly expressed" genes in Treatment
# Definition: genes where Control mean count < 5 AND Treatment mean count >= 10
control_mean <- rowMeans(counts(dds, normalized = TRUE)[, 1:5])
treatment_mean <- rowMeans(counts(dds, normalized = TRUE)[, 6:15])

newly_expressed <- names(which(control_mean < 5 & treatment_mean >= 10))
cat("\nNewly expressed genes in Treatment (Control<5, Treatment>=10):",
    length(newly_expressed), "\n")

# 5. Save results
out_dir <- "/home/yusanghyeon/RAT_project/pipeline_candidates/results/deseq2"

# Full DESeq2 results
res_df <- as.data.frame(res)
res_df$gene <- rownames(res_df)
write.csv(res_df, file.path(out_dir, "deseq2_results.csv"), row.names = FALSE)

# Newly expressed gene list
writeLines(newly_expressed, file.path(out_dir, "newly_expressed_genes.txt"))

# Also save significantly DE genes for reference
sig_genes <- rownames(res_df[!is.na(res_df$padj) & res_df$padj < 0.05, ])
writeLines(sig_genes, file.path(out_dir, "significant_DE_genes.txt"))

cat("Significantly DE genes (padj < 0.05):", length(sig_genes), "\n")
cat("\nOutput saved to:", out_dir, "\n")
cat("Done.\n")
