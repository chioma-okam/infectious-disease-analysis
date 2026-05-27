setwd()
install.packages("writexl")
library(writexl)
library(readxl)
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
require(DOSE)
# Load data
dataframe <- read_excel("VL vs Asymptomatic GSE77528.top.table.xlsx")
dataframe <-readr::read_tsv("VL vs Asymptomatic GSE77528.top.table.tsv")
# Make sure log2FC is numeric: when from excel convert to numeric ist
dataframe$logFC <- as.numeric(dataframe$logFC)

# Ensure Gene.symbol is character
dataframe$Gene.symbol <- as.character(dataframe$Gene.symbol)

# Filter out entries with missing gene symbols
dataframe <- dataframe[!is.na(dataframe$Gene.symbol) & dataframe$Gene.symbol != "", ]

# Remove duplicates (some genes appear more than once)
dataframe <- dataframe[!duplicated(dataframe$Gene.symbol), ]

# Now create a named numeric vector
gene_list <- dataframe$logFC
names(gene_list) <- dataframe$Gene.symbol

# Remove NAs and infinite values
gene_list <- gene_list[!is.na(gene_list) & is.finite(gene_list)]

# Sort in decreasing order
gene_list <- sort(gene_list, decreasing = TRUE)

# Check names
head(gene_list, 5)
# Check top few
head(gene_list)

sum(is.na(names(gene_list)) | names(gene_list) == "") # If this returns anything > 0, 
#that's the problem. You have unnamed values (or blank gene symbols).


library(clusterProfiler)
library(org.Hs.eg.db)

##Run GSEA using gseGO()
gse <- gseGO(
  geneList = gene_list,
  ont = "MF",
  keyType = "SYMBOL",
  minGSSize = 20, # to increase or decrease no of mappings to a gene ID
  maxGSSize = 800,
  pvalueCutoff = 0.05,
  verbose = TRUE,
  OrgDb = org.Hs.eg.db,
  pAdjustMethod = "fdr")## BH can also be used


#Dot plot (split by enrichment direction)
library(enrichplot)
require(DOSE)

dotplot(gse, showCategory = 20, split = ".sign") + 
  facet_grid(. ~ .sign) +
  ggtitle("Dot Plot – GSEA GO Terms (VL vs Asymptomatic)") #show category 20 means showing 20 significant genes, 
# the no can be increased or decreased depending on what you want

## Recommended easier step for saving result in excel
install.packages("writexl")
library(writexl)
write_xlsx(gse@result, "GSEA_results_VL_vs_Asymp_BP.xlsx")



# How to prepare data for Ridge plot
# Ridge plot
ridgeplot(gse, showCategory = 20) +
  ggtitle("Ridge Plot – GSEA GO Terms (VL vs Asymptomatic_ML_Top20)")

print(ridgeplot(gse, showCategory = 20)) ## when you print this way it won't have any caption
## it is better to export using PNG format 

## KEGG PATHWAY
library(readxl)
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
## Step 1: Load and prepare data
dataframe <- read_excel("VL vs Asymptomatic GSE77528.top.table.xlsx")

# Make sure logFC is numeric
dataframe$logFC <- as.numeric(dataframe$logFC)
# Ensure Gene.symbol is character
dataframe$Gene.symbol <- as.character(dataframe$Gene.symbol)

# Remove missing or empty gene symbols
dataframe <- dataframe[!is.na(dataframe$Gene.symbol) & dataframe$Gene.symbol != "", ]

# Remove duplicates
dataframe <- dataframe[!duplicated(dataframe$Gene.symbol), ]
## Step 2: Convert SYMBOL → ENTREZID for KEGG
gene_df <- bitr(dataframe$Gene.symbol,
                fromType = "SYMBOL",
                toType = "ENTREZID",
                OrgDb = org.Hs.eg.db)
# Merge Entrez IDs back into the dataframe to keep logFC order
dataframe <- merge(dataframe, gene_df, by.x = "Gene.symbol", by.y = "SYMBOL")
# Create a named vector for GSEA
gene_list <- dataframe$logFC
names(gene_list) <- dataframe$ENTREZID
gene_list <- sort(gene_list, decreasing = TRUE)

## Step 3: Run GSEA KEGG
gse_kegg <- gseKEGG(
  geneList = gene_list,
  organism = "hsa",      # human
  minGSSize = 20,
  maxGSSize = 800,
  pvalueCutoff = 0.05,
  pAdjustMethod = "fdr",
  verbose = TRUE)
## Step 4: Visualizations
# Dot plot
dotplot(gse_kegg, showCategory = 20, split = ".sign") + 
  facet_grid(. ~ .sign) +
  ggtitle("Dot Plot – GSEA KEGG Pathways (VL vs Asymptomatic)")
# Ridge plot
ridgeplot(gse_kegg, showCategory = 20) +
  ggtitle("Ridge Plot – GSEA KEGG Pathways (VL vs Asymptomatic)")
# Save as Excel file
kegg_kf <- as.data.frame(gse_kegg)

write.xlsx(kegg_kf, file = "KEGG_GSEA_VL_Asymptomatic_results.xlsx", rowNames = FALSE)


