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

EC.CT <- readRDS("/data1/mashuai/data/COVID-19/Seurat_v5/RDS/COVID-19_EC.CT.rds")
Cell_adhesion.list <- read.table('/data1/mashuai/data/GeneSet/EC.cell_adhesion.list')
Cell_adhesion.list  <- as.character(Cell_adhesion.list$V1)


Cell_adhesion.expression <- AverageExpression(EC.CT,features=Cell_adhesion.list,assays='RNA')

my_palette <- colorRampPalette(c('#6495ED', 'white', '#EE3B3B'))(50)


Celltype = data.frame(
Group = factor(rep(c("Control", "COVID"), c(6, 6))),
Celltype = c("Art.EC","Vei.EC","Cap.EC.g.1","Cap.EC.g.2","Cap.EC.a","Lym.EC","Art.EC","Vei.EC","Cap.EC.g.1","Cap.EC.g.2","Cap.EC.a","Lym.EC")
)

rownames(Celltype) = paste(c("Art.EC_O-Control","Vei.EC_O-Control","Cap.EC.g.1_O-Control","Cap.EC.g.2_O-Control","Cap.EC.a_O-Control","Lym.EC_O-Control",
"Art.EC_O-COVID","Vei.EC_O-COVID","Cap.EC.g.1_O-COVID","Cap.EC.g.2_O-COVID","Cap.EC.a_O-COVID","Lym.EC_O-COVID"), sep = "")


ann_colors = list(Group = c(Control = "grey85", COVID = "#4B1C00"),
				  Celltype = c(
                  Art.EC = "#EC9A91",Vei.EC = "#E93E3F",Cap.EC.g.1 = "#F06C45",Cap.EC.g.2 = "#FDAC4F",Cap.EC.a =  "#FB820F",Lym.EC = "#D1AAB7")

				)


p_Cell_adhesion <- pheatmap(Cell_adhesion.expression$RNA, show_rownames = T, show_colnames = T, scale='row', cluster_cols = T, cluster_rows = T,
         fontsize_row = 10, fontsize_col = 10, border_color='NA',cutree_cols=2,
		 annotation_col = Celltype,annotation_colors = ann_colors,
         treeheight_row = 20, treeheight_col = 20, family = 'Arial', color = my_palette,clustering_method = "ward.D")


pdf("/data1/mashuai/data/COVID-19/Seurat_v5/EC/DEGs_Cell_adhesion_heatmap_up.pdf",height = 9.2, width = 4.5)
print(p_Cell_adhesion)
dev.off()
