# ============================================================
# Microarray Preprocessing - GSE10072
# AI & Omics Research Internship - Assignment
# ============================================================

# --- Install & Load Packages ---
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")

BiocManager::install(c("affy", "limma", "Biobase", "GEOquery", "arrayQualityMetrics", "hgu133a.db"),
                     ask = FALSE, update = FALSE)

library(affy)
library(limma)
library(Biobase)
library(GEOquery)
library(arrayQualityMetrics)
library(hgu133a.db)
library(ggplot2)

# Download Series Matrix file
gset <- getGEO("GSE10072", GSEMatrix = TRUE)

# If more than one platform, pick GPL96 (HG-U133A)
if (length(gset) > 1) idx <- grep("GPL96", attr(gset, "names")) else idx <- 1

# Extract actual ExpressionSet
eset <- gset[[idx]]

# Now you can use exprs()
expr_data <- exprs(eset)

# Phenotype/sample info
pheno <- pData(eset)

# Feature data (probe annotation)
feature_data <- fData(eset)


# ============================================================
# STEP 1: Load raw CEL files
# ============================================================
setwd("C:/Users/Dell/Documents/AI_Omics_Internship_2025/GSE10072/raw")
rawData <- ReadAffy()
rawData


# ============================================================
# STEP 2: QC Before Normalization
# ============================================================
# Full QC report (HTML)
arrayQualityMetrics(expressionset = rawData,
                    outdir = "../QC_BeforeNorm",
                    force = TRUE)


# ============================================================
# STEP 3: Normalize with RMA
# ============================================================
normData <- rma(rawData)

# ============================================================
# STEP 4: QC After Normalization
# ============================================================
arrayQualityMetrics(expressionset = normData,
                    outdir = "../QC_AfterNorm",
                    force = TRUE)


# Extract normalized expression
exprs_df <- as.data.frame(exprs(normData))
dim(exprs_df)


# ============================================================
# STEP 5: Filter Low Variance Transcripts
# ============================================================

# Calculate median intensity for each probe
row_median <- apply(exprs_df, 1, median)

# Plot the Median Intensities Probe Distribution
hist(row_median, breaks = 100, freq = FALSE, main= "Median Intensity Distribution")

Threshold = 5.5 

# Add a vertical red line at threshold = 5.5
abline(v = 5.5, col = "red", lwd = 2)

index = row_median > Threshold
exprs_filtered <- exprs_df[index, ]

summary(row_median)

# ============================================================
# STEP 6: Define Phenotype Groups (Normal vs Cancer)
# ============================================================

# Extract phenotype/sample information
pheno <- pData(gset[[1]])

# View the first few rows to understand columns
head(pheno)

# Check which column contains the group info (e.g., source_name_ch1 or characteristics_ch1)
unique(pheno$source_name_ch1)

# Relabel samples into groups (modify the pattern if needed)
pheno$Group <- ifelse(grepl("normal", pheno$source_name_ch1, ignore.case = TRUE),
                      "Normal",
                      "Cancer")

# Check the group distribution
table(pheno$Group)

# Confirm everything looks correct
head(pheno[, c("source_name_ch1", "Group")])

# ============================================================
# STEP 7A: Boxplot After Normalization
# ============================================================

# Create boxplot of normalized data
boxplot(exprs_df,
        main = "Boxplot of Normalized Expression Data",
        las = 2,              # rotate sample names
        col = "lightblue",
        outline = FALSE,
        ylab = "Log2 Expression Intensity")

# save the plot
dev.copy(png, "../Boxplot_After_Normalization.png", width = 1200, height = 600)
dev.off()


# ============================================================
# STEP 7B: PCA Plot After Normalization
# ============================================================

# Perform PCA
pca <- prcomp(t(exprs_df), scale. = TRUE)

# Basic PCA scatter plot
plot(pca$x[,1], pca$x[,2],
     main = "PCA Plot (After Normalization)",
     xlab = "PC1",
     ylab = "PC2",
     col = "steelblue",
     pch = 19)

# Optionally save the plot
dev.copy(png, "../PCA_After_Normalization.png", width = 800, height = 600)
dev.off()

cat("After filtering:", nrow(exprs_filtered))

# Save all objects in your current R session
save.image("C:/Users/Dell/Documents/AI_Omics_Internship_2025/SaniaKhurshid_Class_4_Assignment.RData")

