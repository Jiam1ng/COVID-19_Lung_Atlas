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
library(tidyverse)
library(ggforce) 
library(ggsci)

Lung.combine.CT<-readRDS( "/data1/mashuai/data/COVID-19/Seurat_v5/RDS/COVID-19.celltype.rds")
DefaultAssay(Lung.combine.CT) <- "RNA"
Lung.combine.CT@meta.data$celltype.stage1 <- paste0(Lung.combine.CT@meta.data$celltype, "_", Lung.combine.CT@meta.data$stage1)
Idents(Lung.combine.CT) <- Lung.combine.CT@meta.data$stage1
Lung.combine.CT.trim <- subset(Lung.combine.CT,idents=c('O-Control','O-COVID'))

Apoptosis.list <- read.table('/data1/mashuai/data/GeneSet/Cell_death.list')
Apoptosis.list  <- as.character(Apoptosis.list[,1])
Pyroptosis.list <- read.table('/data1/mashuai/data/GeneSet/Pyroptosis.list')
Pyroptosis.list  <- as.character(Pyroptosis.list[,1])
Ferroptosis.list <- read.table('/data1/mashuai/data/GeneSet/Ferroptosis.list')
Ferroptosis.list  <- as.character(Ferroptosis.list[,1])
Necrocytosis.list <- read.table('/data1/mashuai/data/GeneSet/Necrocytosis.list')
Necrocytosis.list  <- as.character(Necrocytosis.list[,1])

Death.list <- c(Apoptosis.list,Pyroptosis.list,Ferroptosis.list,Necrocytosis.list)

Death.list <- Death.list[!duplicated(Death.list)]

Lung.combine.CT.trim <- AddModuleScore( object = Lung.combine.CT.trim, features=list(Death.list),name = 'Death_score')



colourCount = length(unique(Lung.combine.CT.trim@meta.data$celltype.stage1))
getPalette = colorRampPalette(brewer.pal(12, "Paired"))

RidgePlot(Lung.combine.CT.trim,features=c('Death_score1'),group.by='celltype.stage1',cols=getPalette(colourCount))	
pdf("/data1/mashuai/data/COVID-19/Seurat_v5/RidgePlot.Death.score.pdf",height=3,width=10)
RidgePlot(Lung.combine.CT.trim,features=c('Death_score1'),group.by='stage1',cols=c("grey80", "#4B1C00")	)
dev.off()


my_comparisons <-  list(c("O-COVID", "O-Control"))
	
	
Lung_Death <- data.frame(Lung.combine.CT.trim@meta.data$celltype,Lung.combine.CT.trim@meta.data$stage1,Lung.combine.CT.trim@meta.data$Death_score1)	

colnames(Lung_Death) <- c("Celltype","stage1","Death_score")	
Lung_Death$stage1 <- factor(Lung_Death$stage1,levels=c('O-Control','O-COVID'))
	
Lung_Death$Celltype <- factor(Lung_Death$Celltype,levels= c("AT1","AT2","AT2.trans","Basal","Club","Goblet","Ciliated","Art.EC","Vei.EC","Cap.EC.g","Cap.EC.a","Lym.EC","Alv.Fib","Adv.Fib","Myofib","Air.SMC","Vas.SMC","Peri","CD4.T","CD8.T","Pro.T","NK","BC","Plasmo","AM","Mono","DC","Mast"))

  

	
p4 <- ggplot(Lung_Death,aes(stage1,Death_score,fill=stage1))+
  facet_wrap(~Celltype, scales='free') +
  geom_boxplot(outlier.colour = NA,notch = T,size = 0.4)+
  theme_classic()+
  theme(legend.position = 'none',
        axis.title.y = element_text(size=12),
        axis.text = element_text(size=12),
        axis.title.x = element_blank())+
  scale_fill_npg()+
  scale_fill_manual(values=c("grey80", "#4B1C00"))+
  stat_compare_means(label='p.signif', label.x=1.5, method='wilcox.test',comparisons=my_comparisons)+
  stat_compare_means(label.y = max(Lung_Death$Death)+0.1)	



pdf("/data1/mashuai/data/COVID-19/Seurat_v5/BoxlinPlot.Death_score.celltype.pdf",height=20,width=13)
print(p4)
dev.off()

Epi.score <- c('CLIC5', 'CAV1', 'CAV2', 'SFTPC', 'PGC', 'WIF1','SFTA3')

Idents(Lung.combine.CT.trim) <- Lung.combine.CT.trim@meta.data$celltype
Epi.sub <- subset(Lung.combine.CT.trim,idents=c('AT1','AT2','AT2.trans'))

DotPlot(Epi.sub,features=Epi.score,group.by='celltype',split.by='stage1')

Epi.sub <- AddModuleScore( object = Epi.sub, features=list(Epi.score),name = 'Epi_score')
