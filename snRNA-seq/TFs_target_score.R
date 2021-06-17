library(DoubletFinder)
library(Seurat)
library(dplyr)
library(stringr)
library(cowplot)
library(ggplot2)
library(ggpubr)
library(RColorBrewer)
library(tidyverse)
library(ggforce) 
library(ggsci)
library(pheatmap)

Lung.combine.CT <- readRDS("/data1/mashuai/data/COVID-19/Seurat_v5/RDS/COVID-19.celltype.rds")
Lung.combine.CT@meta.data$celltype.stage1 <- paste0(Lung.combine.CT@meta.data$celltype, "_", Lung.combine.CT@meta.data$stage1)
Idents(Lung.combine.CT) <- Lung.combine.CT@meta.data$stage1
Lung.combine.CT.trim <- subset(Lung.combine.CT,idents=c('O-Control','O-COVID'))
Idents(Lung.combine.CT) <- Lung.combine.CT@meta.data$celltype
Idents(Lung.combine.CT.trim) <- Lung.combine.CT.trim@meta.data$celltype.stage1


NFKB1.list <- read.table('/data1/mashuai/data/GeneSet/NFKB1_target_total.list')
NFKB1.list  <- as.character(NFKB1.list[,1])
HIF1A.list <- read.table('/data1/mashuai/data/GeneSet/HIF1A_target_total.list')
HIF1A.list  <- as.character(HIF1A.list[,1])
FOXO3.list <- read.table('/data1/mashuai/data/GeneSet/FOXO3_target_total.list')
FOXO3.list  <- as.character(FOXO3.list[,1])

Lung.combine.CT.trim <- AddModuleScore( object = Lung.combine.CT.trim, features=list(NFKB1.list),name = 'NFKB1_target')
Lung.combine.CT.trim <- AddModuleScore( object = Lung.combine.CT.trim, features=list(HIF1A.list),name = 'HIF1A_target')
Lung.combine.CT.trim <- AddModuleScore( object = Lung.combine.CT.trim, features=list(FOXO3.list),name = 'FOXO3_target')
								 
Lung.combine.CT.trim@meta.data$celltype.stage1 <- factor(Lung.combine.CT.trim@meta.data$celltype.stage1,levels=rev(c("AT1_O-Control","AT2_O-Control","AT2.trans_O-Control","Basal_O-Control","Club_O-Control","Goblet_O-Control","Ciliated_O-Control","Art.EC_O-Control","Vei.EC_O-Control","Cap.EC.g_O-Control","Cap.EC.a_O-Control","Lym.EC_O-Control","Alv.Fib_O-Control","Adv.Fib_O-Control","Myofib_O-Control","Air.SMC_O-Control","Vas.SMC_O-Control","Peri_O-Control","CD4.T_O-Control","CD8.T_O-Control","Pro.T_O-Control","NK_O-Control","BC_O-Control","Plasmo_O-Control","AM_O-Control","Mono_O-Control","DC_O-Control","Mast_O-Control",
"AT1_O-COVID","AT2_O-COVID","AT2.trans_O-COVID","Basal_O-COVID","Club_O-COVID","Goblet_O-COVID","Ciliated_O-COVID","Art.EC_O-COVID","Vei.EC_O-COVID","Cap.EC.g_O-COVID","Cap.EC.a_O-COVID","Lym.EC_O-COVID","Alv.Fib_O-COVID","Adv.Fib_O-COVID","Myofib_O-COVID","Air.SMC_O-COVID","Vas.SMC_O-COVID","Peri_O-COVID","CD4.T_O-COVID","CD8.T_O-COVID","Pro.T_O-COVID","NK_O-COVID","BC_O-COVID","Plasmo_O-COVID","AM_O-COVID","Mono_O-COVID","DC_O-COVID","Mast_O-COVID")
))


colourCount = length(unique(Lung.combine.CT@meta.data$celltype.stage1))
getPalette = colorRampPalette(brewer.pal(12, "Paired"))

pdf("/data1/mashuai/data/COVID-19/Seurat_v5/RidgePlot.target.score.pdf",height=10,width=10)
RidgePlot(Lung.combine.CT.trim,features=c('NFKB1_target1'),group.by='celltype.stage1',cols=getPalette(colourCount))	
RidgePlot(Lung.combine.CT.trim,features=c('HIF1A_target1'),group.by='celltype.stage1',cols=getPalette(colourCount))	
RidgePlot(Lung.combine.CT.trim,features=c('FOXO3_target1'),group.by='celltype.stage1',cols=getPalette(colourCount))	
dev.off()



