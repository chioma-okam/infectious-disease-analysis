if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("clusterProfiler")
BiocManager::install("DOSE")
BiocManager::install("enrichplot")
BiocManager::install("clusterProfiler")
install.packages("readxl")
library("readxl")
install.packages("pheatmap")
library(pheatmap)
# Load your data
malaria_data <- read_excel("GSE265864_salmon_abundances.xlsx")
# Check the first few rows to make sure it looks right
head(malaria_data)
#summary of how samples are grouped
colnames(malaria_data)
# 1. Load data and convert to a standard data frame
counts_matrix <- as.data.frame(malaria_data)
# 2. Fix duplicate GENE symbols (Rows)
unique_genes <- make.unique(as.character(counts_matrix[, 1]))
rownames(counts_matrix) <- unique_genes
counts_matrix <- counts_matrix[, -1] # Drop the original text column
# 3. Fix duplicate MONKEY names (Columns)
unique_monkeys <- make.unique(colnames(counts_matrix))
colnames(counts_matrix) <- unique_monkeys

# 4. Grab all Brainstem columns from your matrix
brainstem_columns <- grep("Brainstem", colnames(counts_matrix), value = TRUE)
counts_brainstem_all <- counts_matrix[, brainstem_columns]
# 5. Exclude the first 4 columns (The 4 Healthy Controls shown in the paper)
# This leaves exactly the 14 infected columns for your analysis!
counts_brainstem <- counts_brainstem_all[, -c(1:4)]

# 6. Double-check to ensure R sees exactly 14 columns now
print(ncol(counts_brainstem))

# 7. Build your metadata table matching the paper's exact groups
sample_metadata <- data.frame(
  row.names = colnames(counts_brainstem),
  condition = c(rep("Untreated", 5), 
                rep("MB_6mg", 3), 
                rep("MB_10mg", 3), 
                rep("MB_17mg", 3))) # Changed to 3 to perfectly match the paper!

# 8. Set 'Untreated' as your statistical baseline reference
sample_metadata$condition <- factor(sample_metadata$condition, 
                                    levels = c("Untreated", "MB_6mg", "MB_10mg", "MB_17mg"))
# 9. Load DESeq2 and build the dataset object
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install("DESeq2")
library(DESeq2)
dds <- DESeqDataSetFromMatrix(countData = round(counts_brainstem),
                              colData = sample_metadata,
                              design = ~ condition)

# 10. Run the master statistical calculations
dds <- DESeq(dds)
# 11. Extract results for the 17mg dose vs Untreated baseline
res_17mg <- results(dds, contrast=c("condition", "MB_17mg", "Untreated"))

# 12. Print the statistical breakdown
summary(res_17mg)

BiocManager::install("EnhancedVolcano")

library(EnhancedVolcano)

# 1. Identify your top 50 most statistically significant genes, here i gave instruction for 50 genes
top50_genes <- head(rownames(res_17mg[order(res_17mg$padj, na.last = TRUE), ]), 50)

# 2. Create a custom color vector for every single gene
keyvals <- rep('grey30', nrow(res_17mg))
names(keyvals) <- rep('Not Significant', nrow(res_17mg))

# 3. Assign Royal Blue to Downregulated targets (padj < 0.05 and LFC <= -1)
keyvals[which(res_17mg$padj < 0.05 & res_17mg$log2FoldChange <= -1)] <- 'royalblue'
names(keyvals)[which(res_17mg$padj < 0.05 & res_17mg$log2FoldChange <= -1)] <- 'Downregulated by Drug'

# 4. Assign Red3 to Upregulated targets (padj < 0.05 and LFC >= 1)
keyvals[which(res_17mg$padj < 0.05 & res_17mg$log2FoldChange >= 1)] <- 'firebrick3'
names(keyvals)[which(res_17mg$padj < 0.05 & res_17mg$log2FoldChange >= 1)] <- 'Upregulated by Drug'

# 5. Generate the plot forcing exactly those 50 labels
EnhancedVolcano(res_17mg,
                lab = rownames(res_17mg),
                x = 'log2FoldChange',
                y = 'padj',
                selectLab = top50_genes,     # This forces R to label ONLY your top 50 genes
                title = 'Infected Baseline vs. 17mg Methylene Blue Treatment',
                subtitle = 'Showcasing the Top 50 Brainstem Differentially Expressed Genes',
                pCutoff = 0.05,            
                FCcutoff = 1.0,             
                pointSize = 2.0,            # size of the dots
                labSize = 3.5,               # Slightly smaller text so all 50 names fit neatly
                drawConnectors = TRUE,       # Draws elegant lines connecting the text to the dot
                widthConnectors = 0.5,
                colConnectors = 'grey50',
                max.overlaps = 30,           # <-- THIS LINE fixes the text stacking over DEPP1 or any gene!
                # tells R to "Spend up to 30 iterative attempts trying to shift, nudge, and rearrange these words(gene names) until they no longer touch each other."
                colCustom = keyvals,          # <-- This activates your new custom color rules!
                legendPosition = 'right', 
                legendLabSize = 10,
                legendIconSize = 4.0,
                caption = 'Right side = Upregulated by Drug | Left side = Downregulated by Drug')

# 1. Convert your results to a clean data frame
res_17mg_df <- as.data.frame(res_17mg)
res_17mg_df$Gene_Symbol <- rownames(res_17mg_df)
res_17mg_df <- res_17mg_df[, c("Gene_Symbol", "baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj")]

# 2. Export ALL 25,979 genes to a master CSV
write.csv(res_17mg_df, file = "Brainstem_Malaria_17mg_vs_Untreated_All_Genes.csv", row.names = FALSE)

# 3. Filter and export ONLY your statistically significant genes (p < 0.05, absolute LFC >= 1)
significant_genes <- subset(res_17mg_df, padj < 0.05 & abs(log2FoldChange) >= 1)
significant_genes <- significant_genes[order(-abs(significant_genes$log2FoldChange)), ]
write.csv(significant_genes, file = "Brainstem_Malaria_17mg_Significant_Targets.csv", row.names = FALSE)

# 4. Extract just the top 50 genes based on statistical significance
top50_df <- subset(res_17mg_df, Gene_Symbol %in% top50_genes)

# Sort them so they match the order of significance
top50_df <- top50_df[order(top50_df$padj), ]

# Export to a separate, dedicated Excel-friendly file
write.csv(top50_df, file = "Brainstem_Malaria_17mg_Top_50_Lobe_Genes.csv", row.names = FALSE)



library(pheatmap)
# 1. Transform the raw counts to a normalized log2 scale (stabilizes variance)
vsd <- vst(dds, blind = FALSE)
mat <- assay(vsd)

# 2. Filter rows for your top 50 genes
mat_top50 <- mat[top50_genes, ]

# 3. Create a clean design table to label the top columns (Samples)
# Make sure "condition" matches the column name in your colData(dds)
df_annotation <- as.data.frame(colData(dds)[, c("condition"), drop = FALSE])

# 4. Define custom publication colors: Blue (Downregulated) to Red (Upregulated)
heatmap_colors <- colorRampPalette(c("royalblue3", "white", "firebrick3"))(50)

# 5. Render the Heatmap
pheatmap(mat_top50, 
         scale = "row",                      # Centers and scales data gene-by-gene
         clustering_distance_cols = "euclidean", 
         clustering_method = "complete",     # Groups similar samples together
         annotation_col = df_annotation,     # Adds the color bar at the top for groups
         color = heatmap_colors,             # Uses your custom blue-white-red palette
         show_rownames = TRUE,               # Shows your 50 gene symbols
         show_colnames = FALSE,              # Hides messy raw sample IDs
         main = "Expression Profile of Top 50 Brainstem DEGs",
         border_color = NA                  # Removes harsh borders for a smooth look
         )

# Preparing for GO and KEGG in R  

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager") 
BiocManager::install("fgsea")
# Install the core pathway and mapping packages
BiocManager::install(c("clusterProfiler", "org.Mmu.eg.db", "enrichplot"))

# Load the libraries to make sure they open cleanly
library(clusterProfiler)
library(org.Mmu.eg.db)
library(enrichplot)

# A. Force convert your DESeq results into a clean, queryable data frame
res_17mg_df <- as.data.frame(res_17mg)
# B. Add the Gene Symbols as an actual text column instead of just rownames
res_17mg_df$Gene_Symbol <- rownames(res_17mg_df)
# 1. Now extract all detected genes to use as your statistical background
background_genes <- res_17mg_df$Gene_Symbol
# 2. Extract only your statistically significant targets (p < 0.05 and absolute LFC >= 1)
sig_genes_df <- subset(res_17mg_df, padj < 0.05 & abs(log2FoldChange) >= 1)
significant_genes <- sig_genes_df$Gene_Symbol       
# 3. Convert your text Gene Symbols into official Entrez ID numbers
bg_entrez <- bitr(background_genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Mmu.eg.db)
sig_entrez <- bitr(significant_genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Mmu.eg.db)       
# Look at your mapped significant IDs to confirm it worked!
head(sig_entrez)

# 1. Run the GO algorithm with relaxed, discovery-phase thresholds
go_results <- enrichGO(
  gene          = sig_entrez$ENTREZID,
  universe      = bg_entrez$ENTREZID,
  OrgDb         = org.Mmu.eg.db,
  ont           = "BP",
  pAdjustMethod = "none",              # Look at raw p-values first to capture trends
  pvalueCutoff  = 0.05, # Only pathways with an unadjusted p-value less than 0.05 are kept.
  qvalueCutoff  = 0.2, # 
  readable      = TRUE)
# 2. Convert to a data frame again
go_results_df <- as.data.frame(go_results)
# Check your Environment window now—it should show a bunch of observations!
print(paste("Pathways found:", nrow(go_results_df)))

# Save the master pathway table to your project folder
write.csv(go_results_df, "Brainstem_Malaria_GO_Biological_Processes.csv", row.names = FALSE)

# Let's peek at the top 5 most enriched pathways!
head(go_results_df[, c("ID", "Description", "GeneRatio", "p.adjust")], 5)
# This will render your Dotplot
dotplot(go_results, showCategory = 5, title = "Top Enriched Biological Processes") # determines the no. of pathways to appear on your plot

# This will render your Barplot
barplot(go_results, showCategory = 5, title = "GO Significance Levels") # determines the no. of enriched pathways you want to appear

# 1. Run the KEGG Pathway Enrichment Analysis
kegg_results <- enrichKEGG(
  gene         = sig_entrez$ENTREZID, # the function looks at your entire pool of significant genes mixed together, irrespective of if MB turned it up or down will show on graph
  universe     = bg_entrez$ENTREZID, #Establishes your "background universe. This represents the complete list of all genes detected across your entire brainstem RNA-seq run. 
  organism     = "mcc",                # 'mcc' is the official KEGG code for Rhesus Macaque
  pvalueCutoff = 0.05, #Drops any calculated pathway from your final dataset if its unadjusted enrichment probability score is worse than 0.05.
  qvalueCutoff = 0.2   # Applies a false discovery rate control layer. It keeps pathways only if they maintain a 20% or lower chance of being a false positive.
)

# 2. Convert the KEGG results to a readable data frame
kegg_results_df <- as.data.frame(kegg_results)

# Print the number of pathways found to your console
print(paste("KEGG Pathways found:", nrow(kegg_results_df)))

# 3. If pathways are found, save the table and generate the dotplot
if(nrow(kegg_results_df) > 0) {
  # Save the master KEGG spreadsheet to your project folder
  write.csv(kegg_results_df, "Brainstem_Malaria_KEGG_Pathways.csv", row.names = FALSE)
  
  # Render and print the KEGG dotplot
  print(dotplot(kegg_results, showCategory = 10, title = "Top Enriched KEGG Pathways"))
}
