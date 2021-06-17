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

DEGs_snRNA.list <- read.table('/data1/mashuai/data/COVID-19/Seurat_v5/DEGs_0.5_CT/gene_freq_all.trim.list')
DEGs_ID  <- as.character(DEGs_snRNA.list$V1)


DEGs.expression <- AverageExpression(Lung.combine.CT.trim,features=DEGs_ID,assays='RNA')

my_palette <- colorRampPalette(c('#6495ED', 'white', '#EE3B3B'))(50)


Celltype = data.frame(
Group = factor(rep(c("Control", "COVID"), c(21, 21))),
Type = c("Epihelial_cell","Epihelial_cell","Epihelial_cell","Epihelial_cell","Epihelial_cell","Epihelial_cell","Epihelial_cell","Endothelial_cell","Endothelial_cell","Endothelial_cell","Endothelial_cell","Stromal_cell","Stromal_cell","Stromal_cell","Stromal_cell","Immune_cell","Immune_cell","Immune_cell","Immune_cell","Immune_cell","Proliferative_cell","Epihelial_cell","Epihelial_cell","Epihelial_cell","Epihelial_cell","Epihelial_cell","Epihelial_cell","Epihelial_cell","Endothelial_cell","Endothelial_cell","Endothelial_cell","Endothelial_cell","Stromal_cell","Stromal_cell","Stromal_cell","Stromal_cell","Immune_cell","Immune_cell","Immune_cell","Immune_cell","Immune_cell","Proliferative_cell"),
Celltype = c("AT1","AT2","AT1-AT2","Basal","Club","Gob","Ciliated","AEC","VEC","CEC","LEC","Fib","Myofib","SMC","Pericyte","TC","BC","Mac","DC","Mast","Pro","AT1","AT2","AT2.s","Basal","Club","Gob","Ciliated","AEC","VEC","CEC","LEC","Fib","Myofib","SMC","Per","TC","BC","Mac","DC","Mast","Pro")
)

rownames(Celltype) = paste(c("AT1_O-Control","AT2_O-Control","AT2-s_O-Control","Basal_O-Control","Club_O-Control","Gob_O-Control","Ciliated_O-Control","AEC_O-Control","VEC_O-Control","CEC_O-Control","LEC_O-Control","Fib_O-Control","Myofib_O-Control","SMC_O-Control","Per_O-Control","TC_O-Control","BC_O-Control","Mac_O-Control","DC_O-Control","Mast_O-Control","Pro_O-Control",
"AT1_O-COVID","AT2_O-COVID","AT2-s_O-COVID","Basal_O-COVID","Club_O-COVID","Gob_O-COVID","Ciliated_O-COVID","AEC_O-COVID","VEC_O-COVID","CEC_O-COVID","LEC_O-COVID","Fib_O-COVID","Myofib_O-COVID","SMC_O-COVID","Per_O-COVID","TC_O-COVID","BC_O-COVID","Mac_O-COVID","DC_O-COVID","Mast_O-COVID","Pro_O-COVID"), sep = "")


ann_colors = list(Group = c(Control = "grey85", COVID = "#4B1C00"),
                  Type = c(Epihelial_cell = "#21637F", Endothelial_cell = "#E71F16", Stromal_cell = '#FAD1AC',Immune_cell = '#D396C0', Proliferative_cell = "#B15928"),
                  log2FoldChange = c('grey95', '#EE3B3B'),
				  Celltype = c(AT1 = "#A6CEE3",AT2 = "#5B9EC9", AT2.s = "#2D82AF", Basal = "#7EBA98",Club = "#98D277",Gob = "#52AF43",Ciliated = "#6F9E4C",AEC = "#DD9A88",VEC = "#F16667",CEC = "#E42022",LEC =  "#F06C45",Fib = "#FDBB69",Myofib = "#FE982C",SMC = "#F78620",Per = "#D9A295",TC = "#B294C7",BC = "#7D54A5",Mac = "#9E8099",DC = "#F0EB99",Mast = "#DBB466",Pro="#B15928")
				)


p5 <- pheatmap(DEGs.expression$RNA, show_rownames = T, show_colnames = T, scale='row', cluster_cols = T, cluster_rows = T,
         fontsize_row = 10, fontsize_col = 10, border_color='NA',kmeans_k= 10,cutree_cols=2,
		 annotation_col = Celltype,annotation_colors = ann_colors,
         treeheight_row = 20, treeheight_col = 20, family = 'Arial', color = my_palette,clustering_method = "ward.D")


pdf("/data1/mashuai/data/COVID-19/Seurat_v4/DEGs_snRNA_heatmap.pdf",height = 4, width = 10)
print(p5)
dev.off()

cluster.id <- p5$kmeans$cluster

write.csv(cluster.id,"/data1/mashuai/data/COVID-19/Seurat_v4/DEGs_snRNA_heatmap_clusterid.csv")


