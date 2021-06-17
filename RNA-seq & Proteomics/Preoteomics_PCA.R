library(missMDA)
library(ggplot2)
library(RColorBrewer)
library(pheatmap)

matr <- read.csv("D:/master/project/COVID-19/proteome/HA007DA初分析结果-总/COVID_MS_final.csv")
p.val <- matr$O.Y.P.value
p.adj <- p.adjust(p.val, method = 'BH')

matr$OvsY_P_adj <- p.adj

write.csv(matr, "D:/master/project/COVID-19/proteome/HA007DA初分析结果-总/COVID_MS_final_padj.csv", row.names=F)
dep <- subset(matr, abs(RM_LD_RD.O.Log2FC) > 1.5 & P_adj < 0.05)
dep$DE <- ifelse(dep$RM_LD_RD.O.Log2FC>0, "up", "down")
#down   up 
#796 1503

# PCA
protein <- read.csv("D:/master/project/01_COVID-19/proteome/HA007DA初分析结果-总/COVID_MS_final_padj.csv")
protein <- protein[!is.na(protein$PvsO_P),]
protein <- protein[,c(10:66)]

pca.dat <- imputePCA(protein)
pca.dat.rec <- pca.dat$fittedX
colnames(pca.dat.rec) <- colnames(protein)
protein.pca <- prcomp(pca.dat.rec,scale=FALSE)
pca.res <- protein.pca$rotation


pca.res <- as.data.frame(pca.res)
pca.res$sample <- rownames(pca.res)
pca.res$condition <- factor(c(rep("Patient", 45), rep("Old", 12)), levels = c("Old", "Patient"))

write.csv(pca.res, 'D:/master/project/COVID-19/proteome/HA007DA初分析结果-总/COVID19_protein_PCA_PvsO.csv', row.names=F)

# 4*3
ggplot(pca.res, aes(PC1, PC2, color=condition, shape=condition)) + 
  geom_point(size=3) +
  theme_bw() +
  scale_color_manual(values = c('#D96846', '#4B2415'))
#  geom_text_repel(aes(label=sample), point.padding = unit(0.8, "lines"),size=4, color='black') +

eu.dist <-as.matrix(dist(t(pca.dat.rec), method="euclidean"))
sampleDist_colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(eu.dist, cellheight = 10, cellwidth = 10, fontsize = 10, border_color = 'white',
         col=sampleDist_colors, treeheight_row = 10, treeheight_col = 10)
