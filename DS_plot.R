# import ggplot package
library(tidyverse)
library(ggpubr)

# import repel text package
library(ggrepel)

dirPath <- '/Users/timyeh/Library/CloudStorage/GoogleDrive-timyeh007@gmail.com/My Drive/Postdoc Zhang Lab/RNAseq/Yoshimi et al Nature 2019/quantas analysis/human'

cassDir <- file.path(dirPath, 'hs.cass_dataset.diff.txt')
mutxDir <- file.path(dirPath, 'hs.mutx_dataset.diff.txt')
iretDir <- file.path(dirPath, 'hs.iret_dataset.diff.txt')
alt5Dir <- file.path(dirPath, 'hs.alt5_dataset.diff.txt')
alt3Dir <- file.path(dirPath, 'hs.alt3_dataset.diff.txt')
tacaDir <- file.path(dirPath, 'hs.taca_dataset.diff.txt')

AStype <- list(cassDir, mutxDir, iretDir, alt5Dir, alt3Dir, tacaDir)
typeNames <- c("Cass", "Mutx", "Iret", "Alt5", "Alt3", "Taca")

# create empty lists
AS_list <- list()
inList <- list()
exList <- list()

# counter
i = 1


# get gene sets with significant splicing changes in each AS type
for (x in AStype){
  ASdata <- read.table(x, header = 1)
  
  # store filtered datasets inside a list
  AS_list[[i]] <- na.omit(ASdata)
  i = i + 1
}

for (t in 1:6){
# filter the dataset
AS_list[[t]] <- AS_list[[t]][AS_list[[t]]$FDR != 0,]

# create a column that get gene name labels for DE genes
AS_list[[t]]$DSlabel <- NA
# need to covert from factor to string
AS_list[[t]]$gene <- as.character(AS_list[[t]]$gene)

# create a new column for coloring differentially-spliced genes
AS_list[[t]]$diffSpl <- "No"
AS_list[[t]]$diffSpl[(AS_list[[t]]$FDR < 0.01) & (AS_list[[t]]$dI_g1_vs_g2 > 0.1)] <- "Inclusion"
AS_list[[t]]$diffSpl[(AS_list[[t]]$FDR < 0.01) & (AS_list[[t]]$dI_g1_vs_g2 < -0.1)] <- "Exclusion"

inGroup <- AS_list[[t]][AS_list[[t]]$diffSpl == "Inclusion",]
exGroup <- AS_list[[t]][AS_list[[t]]$diffSpl == "Exclusion",]

# label the top 10 genes with smallest p-values
inGroup <- inGroup %>% arrange(FDR)
exGroup <- exGroup %>% arrange(FDR)

if(nrow(inGroup) == 0){}

else if (nrow(inGroup) < 10){
  inGroup$DSlabel[1:nrow(inGroup)] <- inGroup$gene[1:nrow(inGroup)]
}

else {inGroup$DSlabel[1:10] <- inGroup$gene[1:10]}


if(nrow(exGroup) == 0){}

else if (nrow(exGroup) < 10){
  exGroup$DSlabel[1:nrow(exGroup)] <- inGroup$gene[1:nrow(exGroup)]
}

else {exGroup$DSlabel[1:10] <- exGroup$gene[1:10]}


# store differentially-spliced data
inList[[t]] <- inGroup
exList[[t]] <- exGroup

# create a column that get gene name labels for DE genes
AS_list[[t]]$DSlabel <- NA
# need to covert from factor to string
AS_list[[t]]$gene <- as.character(AS_list[[t]]$gene)
AS_list[[t]]$DSlabel[AS_list[[t]]$diffSpl != "No"] <- AS_list[[t]]$gene[AS_list[[t]]$diffSpl != "No"]

# sort the rows based on diffSpl column to rearrange layers for later plotting
AS_list[[t]] <- AS_list[[t]] %>% arrange(desc(AS_list[[t]]$diffSpl))

# plot data
AS_list[[6+t]] <- ggplot(data = AS_list[[t]], aes(x = dI_g1_vs_g2, y = -log10(FDR), col = diffSpl, label = DSlabel)) + 
                  geom_point() + 
                  theme_classic() +
                  scale_colour_manual(values = c("blue", "red", "snow3")) +
                  geom_vline(xintercept = c(-0.1, 0.1), linetype = "dotted") +
                  geom_hline(yintercept = -log10(0.01), linetype = "dotted") +
                  geom_text_repel() +
                  rremove("xlab") +
                  rremove("ylab") +
                  xlim(c(-0.5, 0.5))
}

# name elements in each list
names(inList) <- typeNames
names(exList) <- typeNames

ggarrange(AS_list[[7]], AS_list[[8]], AS_list[[9]], AS_list[[10]], AS_list[[11]], AS_list[[12]],
          legend = "none",
          labels = typeNames,
          nrow = 2, ncol = 3)


# import packages for GO analysis
library(clusterProfiler)
#library(org.Mm.eg.db)
library(org.Hs.eg.db)
library(AnnotationDbi)

# GO analysis
# unregulated genes

#for (tp in typeNames){
  # Exon inclusion
  in_strsplit <- str_split(inList[["Cass"]]$gene, "//")
  in_dfsplit <- data.frame(in_strsplit)
  in_geneGroup <- in_dfsplit[2,]
  
  inGOresults <- enrichGO(gene = in_geneGroup, OrgDb = "org.Hs.eg.db", keyType = "SYMBOL", ont = "BP")
  dfinGO <- as.data.frame(inGOresults@result)
  rownames(dfinGO) <- 1:nrow(dfinGO)
  dfinGO_sorted <- dfinGO %>% arrange(p.adjust)
  
  ggplot(data = dfinGO_sorted[1:20,], aes(x = Count, y = reorder(Description, +Count), fill = p.adjust)) + 
    geom_bar(stat = "identity") +
    theme_classic() +
    scale_x_continuous(expand = c(0, 0)) +
    scale_fill_distiller(palette = "YlOrRd", 
                         limits = c(min(dfinGO_sorted[1:20,]$p.adjust), max(dfinGO_sorted[1:20,]$p.adjust)), 
                         breaks = c(min(dfinGO_sorted[1:20,]$p.adjust), max(dfinGO_sorted[1:20,]$p.adjust))) +
    guides(fill = guide_colourbar(barwidth = 1, barheight = 5, ticks = FALSE)) +
    theme(axis.text.x = element_text(size=15, face="bold", color = "black"),
          axis.text.y = element_text(size=15, face="bold", color = "black"))
  
  # Exon skipping
  ex_strsplit <- str_split(exList[["Cass"]]$gene, "//")
  ex_dfsplit <- data.frame(ex_strsplit)
  ex_geneGroup <- ex_dfsplit[2,]
  
  exGOresults <- enrichGO(gene = ex_geneGroup, OrgDb = "org.Hs.eg.db", keyType = "SYMBOL", ont = "BP")
  dfexGO <- as.data.frame(exGOresults@result)
  rownames(dfexGO) <- 1:nrow(dfexGO)
  dfexGO_sorted <- dfexGO %>% arrange(p.adjust)
  
  ggplot(data = dfexGO_sorted[1:20,], aes(x = Count, y = reorder(Description, +Count), fill = p.adjust)) + 
    geom_bar(stat = "identity") +
    theme_classic() +
    scale_x_continuous(expand = c(0, 0)) +
    scale_fill_distiller(palette = "Blues", 
                         limits = c(min(dfexGO_sorted[1:20,]$p.adjust), max(dfexGO_sorted[1:20,]$p.adjust)), 
                         breaks = c(min(dfexGO_sorted[1:20,]$p.adjust), max(dfexGO_sorted[1:20,]$p.adjust))) +
    guides(fill = guide_colourbar(barwidth = 1, barheight = 5, ticks = FALSE)) +
    theme(axis.text.x = element_text(size=15, face="bold", color = "black"),
          axis.text.y = element_text(size=15, face="bold", color = "black"))
  
#}


