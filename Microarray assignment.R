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
exprs_df <- as.data.frame(exprs_norm)

exprs_out <- cbind(PROBEID = rownames(exprs_norm), exprs_norm)
