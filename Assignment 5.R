# ============================================================
# Microarray Data Analysis - GSE10072
# AI & Omics Research Internship - Assignment
# ============================================================

# --- Install & Load Packages ---
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")

BiocManager::install(c("affy", "limma", "Biobase", "GEOquery", "arrayQualityMetrics", "hgu133a.db"),
                     ask = FALSE, update = FALSE)
install.packages("pheatmap")
#-------------------------------
# Load Required Libraries
#-------------------------------
library(GEOquery)
library(AnnotationDbi)
library(hgu133a.db)
library(limma)
library(tibble)  
library(dplyr)
library(ggplot2)
library(pheatmap)  

# ============================================================
# Step 1: Load GEO dataset
# ============================================================

# Download Series Matrix file

options(download.file.method = "libcurl")
gset <- getGEO("GSE10072", GSEMatrix = TRUE)


# If more than one platform, pick GPL96 (HG-U133A)
if (length(gset) > 1) idx <- grep("GPL96", attr(gset, "names")) else idx <- 1

# Extract actual Expression Set
eset <- gset[[idx]]

# Now use exprs()
expr_data <- exprs(eset)

# Phenotype/sample info
pheno <- pData(eset)

# Feature data (probe annotation)
feature_data <- fData(eset)


#------------------------------------------
# Step 2: Define groups and create design matrix
#------------------------------------------

# View phenotype info to see how groups are labeled
head(pheno[, grep("characteristics", colnames(pheno))])

# In GSE10072, the sample info is stored like:
# "tissue: normal" or "tissue: tumor"
# We'll extract those values

group_list <- ifelse(grepl("normal", pheno$source_name_ch1, ignore.case = TRUE),
                     "Normal", "Cancer")

# Check distribution
table(group_list)

# Convert to factor (required by limma)
group_list <- factor(group_list, levels = c("Normal", "Cancer"))


# Create design matrix
design <- model.matrix(~ 0 + group_list)
colnames(design) <- levels(group_list)
design

# Define contrast (Cancer vs Normal)
contrast_matrix <- makeContrasts(Cancer_vs_Normal = Cancer - Normal, levels = design)
contrast_matrix


#------------------------------------------
# Step 3: Differential Gene Expression Analysis using limma
#------------------------------------------

# Apply log2 transform if data isn’t already log2
if (max(expr_data) > 100) {
  expr_data <- log2(expr_data + 1)
}

# Fit linear model
fit <- lmFit(expr_data, design)
colnames(design)

# Apply contrast (Cancer vs Normal)
fit2 <- contrasts.fit(fit, contrast_matrix)

# Compute statistics with empirical Bayes
fit2 <- eBayes(fit2)

# View top differentially expressed genes
topTable(fit2, adjust = "fdr", number = 10)


#------------------------------------------
# Step 4 : Map Probe IDs to Gene Symbols
#------------------------------------------

# Map probe IDs to gene symbols using hgu133a.db
probe2gene <- AnnotationDbi::select(hgu133a.db,
                                    keys = rownames(expr_data),
                                    columns = c("SYMBOL", "GENENAME"),
                                    keytype = "PROBEID")

# Merge gene symbols with differential expression results
deg_results <- topTable(fit2, adjust = "fdr", number = Inf)
deg_results$PROBEID <- rownames(deg_results)
deg_results <- merge(deg_results, probe2gene, by = "PROBEID", all.x = TRUE)

#------------------------------------------
# Step 5: Handle duplicate probes
#------------------------------------------

# Some genes have multiple probes — keep the one with highest |logFC|
deg_results <- deg_results %>%
  group_by(SYMBOL) %>%
  slice_max(order_by = abs(logFC), n = 1) %>%
  ungroup()

# Remove rows without a gene symbol
deg_results <- na.omit(deg_results)

#------------------------------------------
# Step 6: Save DEG results as CSV files
#------------------------------------------

# Create results folder if it doesn’t exist
if (!dir.exists("Results")) dir.create("Results")

# Save full DEG list
write.csv(deg_results, "Results/All_DEGs.csv", row.names = FALSE)

# Filter up- and downregulated genes (adjusted p < 0.05)
upregulated <- deg_results %>% filter(adj.P.Val < 0.05 & logFC > 0)
downregulated <- deg_results %>% filter(adj.P.Val < 0.05 & logFC < 0)

write.csv(upregulated, "Results/Upregulated.csv", row.names = FALSE)
write.csv(downregulated, "Results/Downregulated.csv", row.names = FALSE)

# Display summary
cat("Total upregulated genes:", nrow(upregulated), "\n")
cat("Total downregulated genes:", nrow(downregulated), "\n")

#------------------------------------------
# Step 6: Save DEG results as CSV files
#------------------------------------------


# Volcano plot
volcano_file <- "Results/Volcano_Plot.png"

ggplot(deg_results, aes(x = logFC, y = -log10(adj.P.Val))) +
  geom_point(aes(color = ifelse(logFC > 0 & adj.P.Val < 0.05, "Up",
                                ifelse(logFC < 0 & adj.P.Val < 0.05, "Down", "NS"))),
             alpha = 0.6) +
  scale_color_manual(values = c("Up" = "red", "Down" = "blue", "NS" = "grey")) +
  theme_minimal() +
  xlab("log2 Fold Change") +
  ylab("-log10 Adjusted P-value") +
  ggtitle("Volcano Plot: Cancer vs Normal") +
  theme(legend.title = element_blank())

# Save as PNG
ggsave(volcano_file, width = 8, height = 6)

# ============================================================
# Step 8: Heatmap of top 25 DEGs
# ============================================================

# Select top 25 genes by absolute logFC
top25_genes <- deg_results %>% 
  arrange(desc(abs(logFC))) %>% 
  slice(1:25)

# Extract expression data for these genes
top25_expr <- expr_data[top25_genes$PROBEID, ]

# Row names as gene symbols
rownames(top25_expr) <- top25_genes$SYMBOL

# Prepare annotation for columns
annotation_col <- data.frame(Group = group_list)
rownames(annotation_col) <- colnames(top25_expr)  


# Save as PNG
heatmap_file <- "Results/Top25_DEGs_Heatmap.png"
png(heatmap_file, width = 1000, height = 800)

pheatmap::pheatmap(top25_expr,
                   scale = "row",
                   annotation_col = annotation_col,  
                   show_rownames = TRUE,
                   show_colnames = FALSE,
                   main = "Top 25 DEGs Heatmap")
dev.off()

