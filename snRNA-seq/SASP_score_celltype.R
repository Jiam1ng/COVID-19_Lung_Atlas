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

Lung.combine.CT <- readRDS("/data1/mashuai/data/COVID-19/Seurat_v4/RDS/Lung.combine.CT.rds")
Lung.combine.CT@meta.data$celltype.stage1 <- paste0(Lung.combine.CT@meta.data$celltype, "_", Lung.combine.CT@meta.data$stage1)
Idents(Lung.combine.CT) <- Lung.combine.CT@meta.data$stage1
Lung.combine.CT.trim <- subset(Lung.combine.CT,idents=c('O-Control','O-COVID'))
Idents(Lung.combine.CT) <- Lung.combine.CT@meta.data$celltype
Idents(Lung.combine.CT.trim) <- Lung.combine.CT.trim@meta.data$celltype.stage1
	
SASP.list <- read.table('/data1/mashuai/data/GeneSet/SASP_V2.list')
SASP.list  <- SASP.list$V1

#ANG, CXCL8, MMP12,

SASP.list <- SASP.list[-1] 
SASP.list <- SASP.list[-26] 
SASP.list <- SASP.list[-59]
SASP.list <- SASP.list[-5]
SASP.list <- as.character(SASP.list)

Lung.combine.CT.trim <- AddModuleScore( object = Lung.combine.CT.trim, features=list(SASP.list),name = 'SASP')

SASP.expression <- AverageExpression(Lung.combine.CT.trim,features=SASP.list,assays='RNA')
			


my_palette <- colorRampPalette(c('#6495ED', 'white', '#EE3B3B'))(50)


Celltype = data.frame(
Group = factor(rep(c("Control", "COVID"), c(28, 28))),
Type = c("Epihelial_cell","Epihelial_cell","Epihelial_cell","Epihelial_cell","Epihelial_cell","Epihelial_cell","Epihelial_cell","Endothelial_cell","Endothelial_cell","Endothelial_cell","Endothelial_cell","Endothelial_cell","Stromal_cell","Stromal_cell","Stromal_cell","Stromal_cell","Stromal_cell","Stromal_cell","Immune_cell","Immune_cell","Immune_cell","Immune_cell","Immune_cell","Immune_cell","Immune_cell","Immune_cell","Immune_cell","Immune_cell","Epihelial_cell","Epihelial_cell","Epihelial_cell","Epihelial_cell","Epihelial_cell","Epihelial_cell","Epihelial_cell","Endothelial_cell","Endothelial_cell","Endothelial_cell","Endothelial_cell","Endothelial_cell","Stromal_cell","Stromal_cell","Stromal_cell","Stromal_cell","Stromal_cell","Stromal_cell","Immune_cell","Immune_cell","Immune_cell","Immune_cell","Immune_cell","Immune_cell","Immune_cell","Immune_cell","Immune_cell","Immune_cell"),
Celltype = c("AT1","AT2","AT2.trans","Basal","Club","Goblet","Ciliated","Art.EC","Vei.EC","Cap.EC.g","Cap.EC.a","Lym.EC","Alv.Fib","Adv.Fib","Myofib","Air.SMC","Vas.SMC","Peri","CD4.T","CD8.T","Pro.T","NK","BC","Plasmo","AM","Mono","DC","Mast","AT1","AT2","AT2.trans","Basal","Club","Goblet","Ciliated","Art.EC","Vei.EC","Cap.EC.g","Cap.EC.a","Lym.EC","Alv.Fib","Adv.Fib","Myofib","Air.SMC","Vas.SMC","Peri","CD4.T","CD8.T","Pro.T","NK","BC","Plasmo","AM","Mono","DC","Mast")
)

rownames(Celltype) = paste(c("AT1_O-Control","AT2_O-Control","AT2.trans_O-Control","Basal_O-Control","Club_O-Control","Goblet_O-Control","Ciliated_O-Control","Art.EC_O-Control","Vei.EC_O-Control","Cap.EC.g_O-Control","Cap.EC.a_O-Control","Lym.EC_O-Control","Alv.Fib_O-Control","Adv.Fib_O-Control","Myofib_O-Control","Air.SMC_O-Control","Vas.SMC_O-Control","Peri_O-Control","CD4.T_O-Control","CD8.T_O-Control","Pro.T_O-Control","NK_O-Control","BC_O-Control","Plasmo_O-Control","AM_O-Control","Mono_O-Control","DC_O-Control","Mast_O-Control",
"AT1_O-COVID","AT2_O-COVID","AT2.trans_O-COVID","Basal_O-COVID","Club_O-COVID","Goblet_O-COVID","Ciliated_O-COVID","Art.EC_O-COVID","Vei.EC_O-COVID","Cap.EC.g_O-COVID","Cap.EC.a_O-COVID","Lym.EC_O-COVID","Alv.Fib_O-COVID","Adv.Fib_O-COVID","Myofib_O-COVID","Air.SMC_O-COVID","Vas.SMC_O-COVID","Peri_O-COVID","CD4.T_O-COVID","CD8.T_O-COVID","Pro.T_O-COVID","NK_O-COVID","BC_O-COVID","Plasmo_O-COVID","AM_O-COVID","Mono_O-COVID","DC_O-COVID","Mast_O-COVID"), sep = "")


ann_colors = list(Group = c(Control = "grey85", COVID = "#4B1C00"),
                  Type = c(Epihelial_cell = "#21637F", Endothelial_cell = "#E71F16", Stromal_cell = '#FAD1AC',Immune_cell = '#D396C0'),
				  Celltype = c(
                  AT1 = "#A6CEE3",AT2 = "#6FAACF", AT2.trans = "#3887BC", Basal = "#3F8EAA",Club = "#7BB899",Goblet = "#ADDC86",Ciliated = "#79C360",
                  Art.EC = "#45A939",Vei.EC = "#669E48",Cap.EC.g = "#B89B74",Cap.EC.a =  "#F9908F",Lym.EC = "#EF5C5C",
                  Alv.Fib = "#E52829",Adv.Fib = "#EA4A34", Myofib = "#F58E56",Air.SMC = "#FDB762",Vas.SMC = "#FE9D35", Peri = "#FE8308",
                  CD4.T = "#ED8F47",CD8.T = "#D7A49E", Pro.T = "#BBA0CD", NK = "#9471B4", BC = "#6D419C",Plasmo = "#A18499", AM = "#DDD399", Mono = "#F0E084", DC = "#D09C56",Mast = "#B15928")

				)


p_SASP <- pheatmap(SASP.expression$RNA, show_rownames = T, show_colnames = T, scale='row', cluster_cols = T, cluster_rows = T,
         fontsize_row = 10, fontsize_col = 10, border_color='NA',kmeans_k= 8,cutree_cols=2,
		 annotation_col = Celltype,annotation_colors = ann_colors,
         treeheight_row = 20, treeheight_col = 20, family = 'Arial', color = my_palette,clustering_method = "ward.D")

	
	
	
	