# DE analysis - this script is for plotting DE data
# import ggplot package
library(tidyverse)
library(ggpubr)

# import repel text package
library(ggrepel)

dirPath <- '/Users/timyeh/Library/CloudStorage/GoogleDrive-timyeh007@gmail.com/My Drive/Postdoc Zhang Lab/RNAseq/Yoshimi et al Nature 2019/quantas analysis/human'

DEdir <- file.path(dirPath, 'hs.g1vsg2.expr.diff.txt')

# get gene sets with significant expression changes in each AS type
DEdata <- read.table(DEdir, header = 1, fill = TRUE)

# filter undefined values
DEdata_NoNA <- na.omit(DEdata)
DEdata_filtered <- DEdata_NoNA[DEdata_NoNA$gene.names != "No",] 

# convert factor to numeric for calculation
DEdata_filtered$logFC <- as.numeric(as.character(DEdata_filtered$logFC))
DEdata_filtered$logCPM <- as.numeric(as.character(DEdata_filtered$logCPM))

# create a new column for coloring differentially-expressed genes
DEdata_filtered$diffExp <- "No"
DEdata_filtered$diffExp[(DEdata_filtered$p.adj < 0.01) & (DEdata_filtered$logFC > 0.6)] <- "Upregulated"
DEdata_filtered$diffExp[(DEdata_filtered$p.adj < 0.01) & (DEdata_filtered$logFC < -0.6)] <- "Downregulated"

uR <- DEdata_filtered[DEdata_filtered$diffExp == "Upregulated",]
dR <- DEdata_filtered[DEdata_filtered$diffExp == "Downregulated",]

# sort the rows based on diffExp column to rearrange layers for later plotting
DEdata_filtered <- DEdata_filtered %>% arrange(diffExp)

# create a column that get gene name labels for DE genes
DEdata_filtered$DElabel <- NA
# need to covert from factor to string
DEdata_filtered$gene.names <- as.character(DEdata_filtered$gene.names)
#DEdata_filtered$DElabel[DEdata_filtered$diffExp != "No"] <- DEdata_filtered$gene.names[DEdata_filtered$diffExp != "No"]

# label the top 20 genes with smallest p-values
DEdata_sorted <- DEdata_filtered %>% arrange(p.adj)
DEdata_sorted$DElabel[1:20] <- DEdata_sorted$gene.names[1:20]

# plot data
ggplot(data = DEdata_sorted, aes(x = logFC, y = -log10(p.adj), col = diffExp, label = DElabel)) + 
    geom_point() + 
    theme_classic() +
    geom_text_repel() +
    scale_colour_manual(values = c("blue", "snow3", "red")) +
    geom_vline(xintercept = c(-0.6, 0.6), linetype = "dotted") +
    geom_hline(yintercept = -log10(0.01), linetype = "dotted") +
    rremove("legend") +
    theme(axis.text.x = element_text(size=15, face="bold", color = "black"),
        axis.text.y = element_text(size=15, face="bold", color = "black"))

# import packages for GO analysis
library(clusterProfiler)
library(org.Mm.eg.db)
library(org.Hs.eg.db)
#library(AnnotationDbi)

# GO analysis
# unregulated genes
upGOresults <- enrichGO(gene = uR$gene.names, OrgDb = "org.Hs.eg.db", keyType = "SYMBOL", ont = "BP")
dfupGO <- as.data.frame(upGOresults@result)
rownames(dfupGO) <- 1:nrow(dfupGO)
dfupGO_sorted <- dfupGO %>% arrange(p.adjust)

ggplot(data = dfupGO_sorted[1:20,], aes(x = Count, y = reorder(Description, +Count), fill = p.adjust)) + 
  geom_bar(stat = "identity") +
  theme_classic() +
  scale_x_continuous(expand = c(0, 0)) +
  scale_fill_distiller(palette = "YlOrRd", 
                       limits = c(min(dfupGO_sorted[1:20,]$p.adjust), max(dfupGO_sorted[1:20,]$p.adjust)), 
                       breaks = c(min(dfupGO_sorted[1:20,]$p.adjust), max(dfupGO_sorted[1:20,]$p.adjust))) +
  guides(fill = guide_colourbar(barwidth = 1, barheight = 5, ticks = FALSE)) +
  theme(axis.text.x = element_text(size=15, face="bold", color = "black"),
        axis.text.y = element_text(size=15, face="bold", color = "black"))

# downregulated genes
dnGOresults <- enrichGO(gene = dR$gene.names, OrgDb = "org.Hs.eg.db", keyType = "SYMBOL", ont = "BP")
dfdnGO <- as.data.frame(dnGOresults@result)
rownames(dfdnGO) <- 1:nrow(dfdnGO)
dfdnGO_sorted <- dfdnGO %>% arrange(p.adjust)

ggplot(data = dfdnGO_sorted[1:20,], aes(x = Count, y = reorder(Description, +Count), fill = p.adjust)) + 
  geom_bar(stat = "identity") +
  theme_classic() +
  scale_x_continuous(expand = c(0, 0)) +
  scale_fill_distiller(palette = "Blues", 
                       limits = c(min(dfdnGO_sorted[1:20,]$p.adjust), max(dfdnGO_sorted[1:20,]$p.adjust)), 
                       breaks = c(min(dfdnGO_sorted[1:20,]$p.adjust), max(dfdnGO_sorted[1:20,]$p.adjust))) +
  guides(fill = guide_colourbar(barwidth = 1, barheight = 5, ticks = FALSE)) +
  theme(axis.text.x = element_text(size=15, face="bold", color = "black"),
        axis.text.y = element_text(size=15, face="bold", color = "black"))


