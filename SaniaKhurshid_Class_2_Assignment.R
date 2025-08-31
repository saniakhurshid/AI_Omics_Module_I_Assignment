# ---------------------------
# Differential Expression Analysis Assignment
# ---------------------------

# 1. Define paths
input_path <- "C:/Users/Dell/Downloads"   
output_path <- "C:/Users/Dell/Documents/AI_Omics_Internship_2025/Results" 

# 2. Function to classify genes
classify_gene <- function(logFC, padj) {
  if (is.na(padj)) padj <- 1
  if (logFC > 1 & padj < 0.05) {
    return("Upregulated")
  } else if (logFC < -1 & padj < 0.05) {
    return("Downregulated")
  } else {
    return("Not_Significant")
  }
}

# 3. Function to process a file
process_file <- function(filename, outdir) {
  df <- read.csv(file.path(input_path, filename))
  
  # Handle missing padj values
  df$padj[is.na(df$padj)] <- 1
  
  # Classify genes
  df$status <- mapply(classify_gene, df$logFC, df$padj)
  
  # Create Results folder if it doesn't exist
  if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)
  
  # Save processed file
  out_file <- file.path(outdir, paste0("processed_", filename))
  write.csv(df, out_file, row.names = FALSE)
  
  # Print summary
  cat("\nSummary for", filename, ":\n")
  print(table(df$status))
  
  return(df)
}

# 4. Run analysis for both datasets
df1 <- process_file("DEGs_Data_1.csv", output_path)
df2 <- process_file("DEGs_Data_2.csv", output_path)

# 5. Visualization - Barplots
barplot(table(df1$status), 
        col = c("blue", "grey", "red"), 
        main = "DEGs_Data_1 Gene Classification", 
        ylab = "Number of Genes")

barplot(table(df2$status), 
        col = c("blue", "grey", "red"), 
        main = "DEGs_Data_2 Gene Classification", 
        ylab = "Number of Genes")

# 6. Extract top 10 Upregulated and Downregulated genes (example for df2)
top_up <- head(df2[order(-df2$logFC), c("Gene_Id", "logFC", "padj")], 10)
top_down <- head(df2[order(df2$logFC), c("Gene_Id", "logFC", "padj")], 10)

cat("\nTop 10 Upregulated Genes (DEGs_Data_2):\n")
print(top_up)

cat("\nTop 10 Downregulated Genes (DEGs_Data_2):\n")
print(top_down)

# Save all objects in current R session
save.image("C:/Users/Dell/Documents/AI_Omics_Internship_2025/SaniaKhurshid_Class_2_Assignment.RData")

