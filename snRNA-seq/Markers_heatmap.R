library(Matrix)
library(DoubletFinder)
library(Seurat)
library(dplyr)
library(stringr)
library(cowplot)
library(ggplot2)
library(RColorBrewer)
library(tidyverse)
library(ggforce) 
library(ggpubr)
library(pheatmap)
library(reshape2)

Lung.combine.CT <- readRDS( "/data1/mashuai/data/COVID-19/Seurat_v5/RDS/COVID-19.celltype.rds")

Idents(Lung.combine.CT) <- Lung.combine.CT@meta.data$stage1
Lung.combine.CT.downsample <- subset(Lung.combine.CT,downsample = 20000)
Idents(Lung.combine.CT.downsample) <- Lung.combine.CT.downsample@meta.data$celltype

CT_markers <- FindAllMarkers(Lung.combine.CT.downsample, only.pos = TRUE, min.pct = 0.25,logfc.threshold=0.25)
write.csv(CT_markers, '/data1/mashuai/data/COVID-19/Seurat_v5/Lung_FindAllMarkers_celltype.csv')

Lung.combine.CT.downsample@meta.data$celltype <- factor(Lung.combine.CT.downsample@meta.data$celltype,levels=c("AT1","AT2","AT2.trans","Basal","Club","Goblet","Ciliated","Art.EC","Vei.EC","Cap.EC.g","Cap.EC.a","Lym.EC","Alv.Fib","Adv.Fib","Myofib","Air.SMC","Vas.SMC","Peri","CD4.T","CD8.T","Pro.T","NK","BC","Plasmo","AM","Mono","DC","Mast"))

Lung.meta <- Lung.combine.CT.downsample@meta.data
celltypes <- c("AT1","AT2","AT2.trans","Basal","Club","Goblet","Ciliated","Art.EC","Vei.EC","Cap.EC.g","Cap.EC.a","Lym.EC","Alv.Fib","Adv.Fib","Myofib","Air.SMC","Vas.SMC","Peri","CD4.T","CD8.T","Pro.T","NK","BC","Plasmo","AM","Mono","DC","Mast")

ct.barc <- list()
for (i in 1:length(celltypes)){
  ct.barc[[i]] <- rownames(subset(Lung.meta, celltype==celltypes[i]))
}
names(ct.barc) <- celltypes

exp_mat <- as.matrix(GetAssayData(Lung.combine.CT.downsample, slot='data'))
cellmarker <- read.csv('/data1/mashuai/data/COVID-19/Seurat_v5/Lung_FindAllMarkers_celltype.csv')
top30 <- cellmarker %>% group_by(cluster) %>% top_n(n = 30, wt = avg_logFC)
top30$cluster <- factor(top30$cluster,levels=c("AT1","AT2","AT2.trans","Basal","Club","Goblet","Ciliated","Art.EC","Vei.EC","Cap.EC.g","Cap.EC.a","Lym.EC","Alv.Fib","Adv.Fib","Myofib","Air.SMC","Vas.SMC","Peri","CD4.T","CD8.T","Pro.T","NK","BC","Plasmo","AM","Mono","DC","Mast"))

ct.exp <- data.frame(matrix(nrow=28366))

for (cell in celltypes){
  barc <- ct.barc[[cell]]
  tmp.mat <- exp_mat[,barc]
  tmp <- data.frame(Mean=rowMeans(tmp.mat))
  colnames(tmp) <- cell
  ct.exp <- cbind(ct.exp, tmp)
}

ct.exp.tp30 <- ct.exp[as.character(top30$gene),]
ct.exp.tp30[is.na(ct.exp.tp30)] <- 0
ct.exp.tp30 <- subset(ct.exp.tp30, select = -matrix.nrow...28366. )

my_palette <- colorRampPalette(c('#6495ED', 'white', '#EE3B3B'))(50)
p1 <- pheatmap(ct.exp.tp30, show_rownames = F, show_colnames = T, scale='row', cluster_cols = F, cluster_rows = F,
         fontsize_row = 10, fontsize_col = 10, 
         treeheight_row = 0, treeheight_col = 20, family = 'Arial', color = my_palette)
		 
pdf("/data1/mashuai/data/COVID-19/Seurat_v5/Celltype_top30_heatmap.pdf ",height = 12, width = 7)
print(p1)
dev.off()	


