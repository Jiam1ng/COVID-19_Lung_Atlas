library(DESeq2)
library(ggplot2)
library(dplyr)
library(stringr)
library(biomaRt)
library(RColorBrewer)
library(pheatmap)
library(ggthemes)

build_mat <- function(Count_dir, Out_dir, Sample, range.seq){
  count_files <- dir(Count_dir)[range.seq]
  for (i in 1:length(count_files)){
    if (i==1){
      count_mat <- read.table(paste0(Count_dir, '/', count_files[i], '/', count_files[i], '_count.txt'), sep='\t', head=F)
      colnames(count_mat) <- c('Gene_ID', Sample[i])
    }else{
      tmp <- read.table(paste0(Count_dir, '/', count_files[i], '/', count_files[i], '_count.txt'), sep='\t', head=F)
      colnames(tmp) <- c('Gene_ID', Sample[i])
      count_mat <- cbind(count_mat,tmp[,-1])
    }
  }
  rownames(count_mat) <- count_mat$Gene_ID
  count_mat <- count_mat[1:(dim(count_mat)[1]-6),-1]
  return(as.matrix(count_mat))
}

tmp2 <- read.table('D:/3.projects/1.single-cell/4.projects/9.COVID-19/02_bulkseq/01_count/1-1CGRLD_BKDL202583814-1a/1-1CGRLD_BKDL202583814-1a_count.txt', header=F)

#5:(dim(count_mat)[1])
#samples <- samples[c(1:18, 33:44,19:32)]
patient <- samples[c(1:18,33:44)]
#young <- samples[35:38]
young <- samples[c(21,22,25:28,31,32)]
#old <- samples[31:34]
old <- samples[c(19,20,23,24,29,30,45,46)]

out_dir <- 'D:/3.projects/1.single-cell/4.projects/9.COVID-19/02_bulkseq/03_DESeq2/final'

sample <- c(young, old)
condition <- factor(c(rep('Young', length(young)), rep('Old', length(old))), levels=c('Young', 'Old'))
group <- 'OvsY'
ran.seq <- c(c(21,22,25:28,31,32), c(19,20,23,24,29,30,45,46))

sample <- c(old, patient)
condition <- factor(c(rep('Old', length(old)), rep('Patient', length(patient))), levels=c('Old', 'Patient'))
group <- 'PvsO'
ran.seq <- c(c(19,20,23,24,29,30,45,46), c(1:18,33:44))

sample <- c(young, old, patient)
condition <- c(rep('A', 23), rep('B', 23))
group <- 'PvsO'
ran.seq <- c(c(21,22,25:28,31,32), c(19,20,23,24,29,30,45,46), c(1:18,33:44))


count_mat <- build_mat(count.dir, out_dir, sample, ran.seq)
colnames(count_mat) <- sample
head(count_mat)
write.csv(count_mat, paste0(out_dir, paste0('/COVID-19_', group, '_RNAseq_count_matrix.csv')))

colData <- data.frame(sample, condition)


dds <- DESeqDataSetFromMatrix(count_mat, colData, design= ~ condition)

#Pre-filtering low count genes
# While it is not necessary to pre-filter low count genes before running the DESeq2 functions,there are two reasons which
# make pre-filtering useful: by removing rows in which there are no reads or nearly no reads, we reduce the memory size of
# the dds data object and we increase the speed of the transformation and testing functions within DESeq2.
# Here we perform a minimal pre-filtering to remove rows that have only 0 or 1 read.
# Note that more strict filtering to increase power is automatically applied via independent filtering on the mean of
# normalized counts within the results function.
dds <- dds[rowSums(counts(dds)) > 1,]

#normalize
dds2 <- DESeq(dds)

# acquire the results using function results(), and assign to res
res <- results(dds2)

# transforms the count data to the log2 scale in a way which minimizes 
# differences between samples for rows with small counts, and which 
# normalizes with respect to library size
rld <- rlog(dds2)
write.csv(assay(rld), paste0(out_dir,paste0('/COVID-19_', group, '_RNAseq_rlog_data.csv')))

# We can order our results table by the smallest BH-adjusted p values:
res <- res[order(res$padj),]
resdata <- merge(as.data.frame(res),as.data.frame(counts(dds2,normalize=TRUE)),by="row.names",sort=FALSE)
names(resdata)[1]<-"GeneID"


#Exporting results to CSV files ####

geneID_symbol <- read.table('D:/5.documents/data/genome/human/hg19_virusAdded/hg19_gene_info_df.txt', sep='\t', header = F)
colnames(geneID_symbol) <- c('GeneID', 'Symbol', "GeneType")

#ensembl <- useEnsembl("ENSEMBL_MART_ENSEMBL",dataset="hsapiens_gene_ensembl", mirror="www")
#hs_id <- getBM(attributes = c('hgnc_symbol','ensembl_gene_id'), mart = ensembl)
#colnames(hs_id) <- c('Symbol', 'GeneID')

resdata_anno<-merge(geneID_symbol,resdata,by="GeneID")
#resdata_anno[resdata_anno$Symbol=='ensembl',]$Symbol <- as.character(resdata_anno[resdata_anno$Symbol=='ensembl',]$GeneID)
write.csv(resdata_anno, paste0(out_dir,'/COVID-19_', group, '_all_gene_anno.csv'),row.names = F)

#resdata_anno <- read.csv(paste0(out_dir,'/COVID-19_', 'PvsO', '_all_gene_anno.csv'))
#sub.gene <- read.table('D:/5.documents/data/genome/human/hg19_virusAdded/IDtoGENENAME.txt', sep='\t', header = F)
#resdata_anno <- subset(resdata_anno, GeneID %in% sub.gene$V1)
diff_gene <-subset(resdata_anno, padj < 0.05 & abs(log2FoldChange) >= 1.5)
write.csv(diff_gene, paste0(out_dir,'/COVID-19_', group, '_diff_gene_anno_sub.csv'),row.names = F)

p.diff <- diff_gene

up.share <- data.frame()
shared <- intersect(subset(diff_gene, log2FoldChange>0)$GeneID, subset(test, log2FoldChange>0)$GeneID)
p.diff$type <- 'PvsO'
p.diff[p.diff$GeneID %in% shared,]$type <- 'shared'
shared <- intersect(subset(diff_gene, log2FoldChange<0)$GeneID, subset(p.diff, log2FoldChange<0)$GeneID)
p.diff[p.diff$GeneID %in% shared,]$type <- 'shared'

write.csv(p.diff, paste0(out_dir,'/COVID-19_', group, '_diff_gene_anno.csv'),row.names = F)
# plot ####
plot_dir <- 'D:/3.projects/1.single-cell/4.projects/1.SC/bulkseq/03_plot/final'

## Euclidean distance plot
## Another useful way to assess overall similarity between samples is sample clustering. To this end,
## We can apply the dist function to the transpose of the transformed count matrix to get sample-to-sample distances.
## Here, We use the R function dist to calculate the Euclidean distance between samples.
## To avoid that the distance measure is dominated by a few highly variable genes, and have a roughly equal contribution
## from all genes, we use the regularized log transformed( rlog-transformed) data by DESeq2.
## Heatmap showing the Euclidean distances between the samples as calculated from the regularized log transformation

test <- read.csv("D:/3.projects/1.single-cell/4.projects/9.COVID-19/02_bulkseq/03_DESeq2/final/COVID-19_PvsO_RNAseq_rlog_data.csv")
test <- test[,-1]
sampleDists <- dist(t(test))
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- colnames(sampleDists)
colnames(sampleDistMatrix) <- NULL
sampleDist_colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)



#6*5
p <- pheatmap(sampleDistMatrix,
              clustering_distance_rows=sampleDists,
              clustering_distance_cols=sampleDists,
              col=sampleDist_colors,
              cellwidth = 10, cellheight = 10, fontsize = 8,
              border_color = 'white',
              treeheight_row = 20, treeheight_col = 20)

## principal components analysis and plot
## PCA plot shows the 2D plane spanned by their first two principal components.
## This type of plot is useful for visualizing the overall effect of experimental covariates and batch effects.

pca_data <- plotPCA(rld, intgroup=c("condition"), returnData=TRUE)
pca_data$group <- factor(c(rep('Young', 12), rep('Old', 12)), levels=c('Young', 'Old'))
write.csv(pca_data,paste(out_dir,"COVID-19_PvsO_PCA_data.csv",sep=""),row.names = F)
percentVar <- round(100 * attr(pca_data, "percentVar"))

#4*3
# young #9FC7DA
# old #E39996
# patient #BFAFCC
ggplot(pca_data, aes(PC1, PC2, color=group, shape=group)) +
  geom_point(size=3) +
  geom_text_repel(data=pca_data, aes(label=name),
                  color='black') +
  scale_color_manual(values=c('#A2C5D7', '#DA6843')) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance"))+
  theme_bw() +
  theme(axis.text = element_text(color='black'), axis.ticks = element_line(color='black'))

# Boxplot of the expression value
#5*5
expression_value_data<-assay(rld)
expression_value_color<-rainbow(7)
boxplot(expression_value_data,col=expression_value_color,main="Expression value of rlog data")

########################## MA Plot ##########################
## The MA-plot represents each gene with a dot. The x axis is the average expression over all samples, and the y axis is the
## log2 fold change of normalized counts between treatment and control.(M means log fold change, A means the count mean.)
## Genes with an adjusted p value below a threshold (here 0.1, the default) are shown in red.Points which fall out of the
## window are plotted as open triangles pointing either up or down.
plotMA(res, main="Spinal_cord_MA_Plot.", ylim=c(-5,5))



Volcano_data <- resdata_anno
Volcano_data$threshold <- factor(ifelse(Volcano_data$padj < 0.05 & abs(Volcano_data$log2FoldChange) > 1,ifelse(Volcano_data$log2FoldChange > 1 ,'Upregulated','Downregulated'),'Unchanged'),
                                 levels = c('Upregulated', 'Unchanged', 'Downregulated'))

volcano_sum<-as.data.frame(t(summary(Volcano_data$threshold)))

Volcano_tile <- paste0('\nCutoff for pvalue is 0.05','\nThe number of upregulated gene is ',volcano_sum$`Up-regulated`,'\nThe number of downregulated gene is ',volcano_sum$`Down-regulated`,'\nThe number of unchanged gene is ',volcano_sum$Unchanged )


Volcano_data$lg10<-log10(Volcano_data$padj)

#b21f1f
#1a2a6c
pdf(paste(out_dir,"COVID19_PvsO_volcano_plot.pdf",sep="/"), width = 5, height = 6)
ggplot( data=Volcano_data,
        aes(x=log2FoldChange, y =lg10,colour=threshold,fill=threshold)) +
  scale_color_manual(values=c('#c00000', 'grey90', '#4372c3'))+
  geom_point(alpha=0.8, size=1.8) +
#  xlim(c(-6.5,6.5)) +ylim(c(0,20))+
#  ggtitle(Volcano_tile) +
  theme_few() +
  geom_vline(xintercept=c(-1,1),lty=4,col="grey",lwd=0.6)+
  geom_hline(yintercept = -log10(0.05),lty=4,col="grey",lwd=0.6)+
#  theme(plot.title = element_text(hjust = 0.5,size=12),legend.title = element_blank())+
  labs(x="log2 (fold change)",y="-log10 (Pvalue)")
dev.off()

heatmap_data<-diff_gene
row.names(heatmap_data)<-heatmap_data$GeneID
heatmap_data<-heatmap_data[,seq(10,dim(diff_gene)[2])]
heatmap_data<-data.matrix(heatmap_data)


#3.5*6
phm <- pheatmap(heatmap_data,
         cluster_rows = T,
         cluster_cols =T,
         clustering_distance_rows = "correlation",#euclidean
         clustering_distance_cols = "correlation",#euclidean
         clustering_method = "complete",
         scale = "row",#"none", "row", "column"
         color = colorRampPalette(c('#0089CE', 'white', '#C63C3C'))(100),
         #color = colorRampPalette(c("blue", "white", "red"))(100),alpha=0.8,
         #c("skyblue","yellow","tomato")
         #color=colorRampPalette(c('green','yellow','red'), bias=1)(50),
         border_color = NA, #"grey60",
         fontsize = 14,
         fontsize_row = 10, #fontsize,
         fontsize_col = 10, #fontsize,
         fontsize_number = 0.8* fontsize,
         #margins=c(5,10)
         #clustering_callback = identity2,
         kmeans_k = NA,
         breaks = NA,
         cutree_rows = NA,
         cutree_cols = NA,
         treeheight_row = 0,
         treeheight_col = 20,
         show_rownames = FALSE,
         show_colnames = TRUE,
         legend = TRUE,
         legend_breaks = NA,
         legend_labels = NA,
         annotation = NA,
         annotation_col= NA,
         #annotation_row = annotation_row,#NA,
         #annotation_colors = ann_colors,
         annotation_legend = TRUE,
         drop_levels = TRUE,
         #display_numbers = F,
         #display_numbers = TRUE, number_format = "%.2f", number_color="purple")
         #number_format = "%.2f",
         #number_color = "grey30",
         gaps_row = NULL,
         gaps_col = NULL,
         labels_row = NULL,
         labels_col = NULL,
         cellwidth = NA,
         cellheight = NA,
         #revC=TRUE,
         width =8,
         height = 10,
)

gene.ord <- phm[['tree_row']]$labels
hm.mat <- tmp
rownames(hm.mat) <- hm.mat$GeneID
hm.mat <- hm.mat[gene.ord, seq(10,dim(hm.mat)[2])]
hm.mat <- hm.mat[!is.na(rowSums(hm.mat)),]

pheatmap(hm.mat,
         cluster_rows = T,
         cluster_cols =F,
         clustering_distance_rows = "correlation",#euclidean
         clustering_distance_cols = "correlation",#euclidean
         clustering_method = "complete",
         scale = "row",#"none", "row", "column"
         color = colorRampPalette(c('#0089CE', 'white', '#C63C3C'))(100),
         #color = colorRampPalette(c("blue", "white", "red"))(100),alpha=0.8,
         #c("skyblue","yellow","tomato")
         #color=colorRampPalette(c('green','yellow','red'), bias=1)(50),
         border_color = NA, #"grey60",
         fontsize = 14,
         fontsize_row = 10, #fontsize,
         fontsize_col = 10, #fontsize,
         fontsize_number = 0.8* fontsize,
         #margins=c(5,10)
         #clustering_callback = identity2,
         kmeans_k = NA,
         breaks = NA,
         cutree_rows = NA,
         cutree_cols = NA,
         treeheight_row = 0,
         treeheight_col = 20,
         show_rownames = FALSE,
         show_colnames = TRUE,
         legend = TRUE,
         legend_breaks = NA,
         legend_labels = NA,
         annotation = NA,
         annotation_col= NA,
         #annotation_row = annotation_row,#NA,
         #annotation_colors = ann_colors,
         annotation_legend = TRUE,
         drop_levels = TRUE,
         #display_numbers = F,
         #display_numbers = TRUE, number_format = "%.2f", number_color="purple")
         #number_format = "%.2f",
         #number_color = "grey30",
         gaps_row = NULL,
         gaps_col = NULL,
         labels_row = NULL,
         labels_col = NULL,
         cellwidth = NA,
         cellheight = NA,
         #revC=TRUE,
         width =8,
         height = 10,
)

# GO plot ####
term <- read.csv('D:/3.projects/1.single-cell/4.projects/9.COVID-19/02_bulkseq/04_GO/final/PvsO_select.csv')
colnames(term)[1] <- 'group'
term$Description <- factor(term$Description, levels = rev(as.character(term$Description)))
term$LogP <- term$LogP*(-1)
#up F8CAAA
#down AFC8E1
#7*3
ggplot(subset(term, group=='down'), aes(Description, LogP)) +
  geom_bar(stat='identity', width=0.8, fill='#AFC8E1')+
  #  scale_fill_gradient(low='#f5af19', high='#C23C32')+
  coord_flip()+
  scale_y_continuous(expand = c(0,0))+
  theme(panel.background = element_blank(), panel.grid = element_blank(),
        axis.line = element_line(colour = "black"), axis.text = element_text(color='black'), axis.ticks = element_line(color='black'))

term <- read.csv('D:/3.projects/1.single-cell/4.projects/9.COVID-19/02_bulkseq/04_GO/200922/COVID-19_bulkseq_PvsO_termSelect.csv')
colnames(term)[1] <- 'group'
term$Description <- factor(term$Description, levels = rev(as.character(term$Description)))
term$LogP <- term$LogP*(-1)
#up C5591C
#down 2F5498
ggplot(subset(term, group=='up'), aes(Description, LogP)) +
  geom_bar(stat='identity', width=0.8, fill='#C5591C')+
  #  scale_fill_gradient(low='#f5af19', high='#C23C32')+
  coord_flip()+
  scale_y_continuous(expand = c(0,0))+
  theme(panel.background = element_blank(), panel.grid = element_blank(),
        axis.line = element_line(colour = "black"), axis.text = element_text(color='black'), axis.ticks = element_line(color='black'))

ace2 <- read.csv('D:/3.projects/1.single-cell/4.projects/9.COVID-19/02_bulkseq/03_DESeq2/200922/ACE2_OvsY.csv')
ace2$age <- factor(ace2$age, levels=c('Young', 'Old'))
# young #9FC7DA
# old #E39996
# patient #BFAFCC

#3*3.5
ggviolin(ace2, x='age', y='ACE2', add = 'boxplot', color='age', add.params = list(fill = "white")) +
  scale_color_manual(values=c('#E39996', '#BFAFCC'))+
  geom_jitter(aes(color=age)) +
  theme(panel.grid = element_blank(), panel.background = element_blank(),
        axis.line = element_line(color='black'), axis.ticks = element_line(color='black'), legend.position = 'none',
        axis.text.y = element_text(color='black'), axis.text.x = element_text(color='black'))

ggplot(ace2, aes(age, ACE2, color=age)) +
  geom_boxplot(outlier.colour = NA) +
  geom_jitter() +
  scale_color_manual(values=c('#9FC7DA', '#E39996'))+
  theme(panel.grid = element_blank(), panel.background = element_blank(),
        axis.line = element_line(color='black'), axis.ticks = element_line(color='black'), legend.position = 'none',
        axis.text.y = element_text(color='black'), axis.text.x = element_text(color='black'))

# gene violin ####
covid.all <- read.csv('D:/3.projects/1.single-cell/4.projects/9.COVID-19/02_bulkseq/03_DESeq2/200922/COVID-19_PvsO_all_gene_anno.csv')
tmp <- subset(covid.all, Symbol %in% c('CDKN1A', 'CDKN2A'))
tmp <- tmp[,c(2, 10:45)]
tmp <- melt(tmp, id.vars='Symbol')
colnames(tmp) <- c('gene', 'sample', 'exp')
tmp$group <- factor(c(rep('Old', 16), rep('Patient', 56)), levels=c('Old', 'Patient'))

#3*4
ggviolin(subset(tmp, gene=='CDKN1A'), x='group', y='exp', add = 'boxplot', color='group', add.params = list(fill = "white")) +
  scale_color_manual(values=c('#D75D36', '#491D0D'))+
#  geom_jitter() +
  stat_compare_means(label='p.format', label.x=1.5, method='wilcox.test') +
  theme(panel.grid = element_blank(), panel.background = element_blank(),
        axis.line = element_line(color='black'), axis.ticks = element_line(color='black'), legend.position = 'none',
        axis.text.y = element_text(color='black'), axis.text.x = element_text(color='black'))

#3*3.5
ggplot(subset(tmp, gene=='CDKN2A'), aes(group, exp, color=group)) +
  geom_boxplot(outlier.colour = NA, width=0.5) +
  scale_color_manual(values=c('#D75D36', '#491D0D'))+
  geom_jitter() +
  stat_compare_means(label='p.format', label.x=1.5, method='wilcox.test') +
  theme(panel.grid = element_blank(), panel.background = element_blank(),
        axis.line = element_line(color='black'), axis.ticks = element_line(color='black'), legend.position = 'none',
        axis.text.y = element_text(color='black'), axis.text.x = element_text(color='black'))

# COVID enrichment ####
term <- read.csv('D:/3.projects/1.single-cell/4.projects/9.COVID-19/02_bulkseq/04_GO/final/PvsO_COVID_select.csv')
colnames(term)[1] <- 'DE'
term <- subset(term, DE=='down')
tmp  <- term$Description
term <- term[,4]
term <- as.matrix(term)
rownames(term) <- tmp
term <- term*(-1)

pheatmap(term, show_rownames = T, show_colnames = T, scale='none', cluster_cols = F, cluster_rows = T,
         border_color='white', fontsize_row = 10, fontsize_col = 10, cellwidth = 15, cellheight = 15, angle_col = 90,
         treeheight_row = 0, family = 'Arial', treeheight_col = 0, color = colorRampPalette(c('white', '#0089CE'))(50)[20:50])

# snRNA-seq ####
deg <- read.csv('D:\3.projects\1.single-cell\4.projects\9.COVID-19\02_bulkseq\03_DESeq2\final/COVID-19_PvsO_diff_gene_anno_sub.csv')

# OvsY ####
oy.all = read.csv('D:/3.projects/1.single-cell/4.projects/9.COVID-19/02_bulkseq/03_DESeq2/final/COVID-19_OvsY_all_gene_anno.csv')

oy.sub <- subset(oy.all, Symbol %in% c('CDKN1A', 'CDKN1B'))
oy.sub <- oy.sub[,c(2,10:25)]
oy.sub <- melt(oy.sub, id.vars='Symbol')
colnames(oy.sub) <- c('gene', 'sample', 'exp')
oy.sub$age <- factor(c(rep('Young', 16),rep('Old', 16)),, levels=c('Young', 'Old'))
#3*4
ggviolin(oy.sub, x='age', y='exp', add = 'boxplot', fill='age', add.params = list(fill = "white")) +
  scale_fill_manual(values=c('#585f86', '#f5a069'))+
  geom_jitter() +
  stat_compare_means(label='p.format', label.x=1.5, method='wilcox.test') +
  facet_wrap(~gene) +
  theme(panel.grid = element_blank(), panel.background = element_blank(),
        axis.line = element_line(color='black'), axis.ticks = element_line(color='black'), legend.position = 'none',
        axis.text.y = element_text(color='black'), axis.text.x = element_text(color='black'))

library(Hmisc)

all.mat = read.csv('D:/3.projects/1.single-cell/4.projects/9.COVID-19/02_bulkseq/03_DESeq2/final/COVID-19_PvsO_all_gene_anno_virusAdded.csv')
all.mat = all.mat[!is.na(all.mat$padj),]
rownames(all.mat) = all.mat$GeneID
all.mat = all.mat[,-c(1:9)]
all.mat = t(all.mat)

viral.cor = rcorr(all.mat, type = 'spearman')

corRes <- function(cor,p, n) {
  ut <- upper.tri(cor) 
  data.frame(row = rownames(cor)[row(cor)[ut]] ,
             column = rownames(cor)[col(cor)[ut]], 
             cor =(cor)[ut], 
             p = p[ut],
             n = n[ut])
}


cor.res <- corRes(viral.cor$r, viral.cor$P, viral.cor$n)

y.out <- subset(y.out, p<0.01)


# UPR pathway ###
all.mat = read.csv('D:/3.projects/1.single-cell/4.projects/9.COVID-19/02_bulkseq/03_DESeq2/final/COVID-19_PvsO_all_gene_anno.csv')
#upr = read.table('D:/3.projects/1.single-cell/4.projects/9.COVID-19/UPR/apoptosis.txt', header=T, sep = '\t')
upr = read.table('D:/3.projects/1.single-cell/4.projects/9.COVID-19/UPR/GO_INFLAMMATORY_RESPONSE.txt', header=T)
upr.mat = subset(all.mat, Symbol %in% as.character(upr$gene))
rownames(upr.mat) = upr.mat$GeneID
upr.sum = upr.mat[,-c(1:9)]
upr.sum = apply(upr.sum, 2, sum)
upr.sum = data.frame(sample = names(upr.sum),
                     group = factor(c(rep('Control', 8), rep('COVID-19', 30)), levels=c('Control', 'COVID-19')),
                     sum = upr.sum)
upr.sum$log10 = log10(upr.sum$sum)
ggviolin(upr.sum, x='group', y='log10', add = 'boxplot', fill='group', add.params = list(fill = "white")) +
  scale_fill_manual(values=c('#d8d8d8', '#4a1c00'))+
  geom_jitter() +
  stat_compare_means(label = "p.format", label.x=1.5, method='wilcox.test') +
  #  facet_wrap(~gene, scales = 'free', ncol=5) +
  theme(panel.grid = element_blank(), panel.background = element_blank(),
        axis.line = element_line(color='black'), axis.ticks = element_line(color='black'), legend.position = 'right',
        axis.text.y = element_text(color='black'), axis.text.x = element_blank())

all.protein = read.csv('D:/3.projects/1.single-cell/4.projects/9.COVID-19/04_Proteome/COVID_MS_final_padj.csv')
all.protein = all.protein[!is.na(all.protein$PvsO_P_adj),]
upr.mat = subset(all.protein, Protein.Gene %in% as.character(upr$gene))
rownames(upr.mat) = upr.mat$Protein.Gene
upr.sum = upr.mat[,-c(1:9)]
upr.sum = apply(upr.sum, 2, sum)
upr.sum = upr.sum[1:57]
upr.sum = data.frame(sample = names(upr.sum),
                     group = factor(c(rep('COVID-19', 45), rep('Control', 12)), levels=c('Control', 'COVID-19')),
                     sum = upr.sum)
upr.sum$log10 = log10(upr.sum$sum)
ggviolin(upr.sum, x='group', y='log10', add = 'boxplot', fill='group', add.params = list(fill = "white")) +
  scale_fill_manual(values=c('#d8d8d8', '#4a1c00'))+
  geom_jitter() +
  stat_compare_means(label = "p.format", label.x=1.5, method='wilcox.test') +
  #  facet_wrap(~gene, scales = 'free', ncol=5) +
  theme(panel.grid = element_blank(), panel.background = element_blank(),
        axis.line = element_line(color='black'), axis.ticks = element_line(color='black'), legend.position = 'right',
        axis.text.y = element_text(color='black'), axis.text.x = element_blank())


upr <- read.table('/data5/lijiaming/projects/01_single-cell/10.COVID-19/03_result/UPR/GO_INFLAMMATORY_RESPONSE.txt', header=T)
Idents(lung) = lung$stage1
genes <- as.character(upr$gene)
genes = genes[genes!=""]
genes = list(genes)
lung <- AddModuleScore(lung, features=genes, seed=15, name="inflammation")
meta = lung[[]]
pdf('inflammation_overall_density.pdf', height=3, width=7)
RidgePlot(lung, cols=c('#d8d8d8', '#4a1c00'), idents=c('O-Control', 'O-COVID'), features='inflammation1') +
  geom_vline(xintercept = median(subset(meta, stage1=='O-Control')$inflammation1), color='black', linetype='dashed') +
  geom_vline(xintercept = median(subset(meta, stage1=='O-COVID')$inflammation1), color='black', linetype='dashed')
dev.off()

pdf('apoptosis_vln_celltype.pdf', height=3, width=10)
VlnPlot(lung, cols=c('#d8d8d8', '#4a1c00'), idents=c('O-Control', 'O-COVID'), features='apoptosis1', 
        split.by='stage1', pt.size=0, group.by='celltype')
dev.off()

p = ggviolin(subset(meta, stage1!='Y-Control'), x='stage1', y='IRE1', add = 'boxplot', color='stage1', add.params = list(fill = "white")) +
  scale_color_manual(values=c('#d8d8d8', '#4a1c00'))+
#  geom_jitter() +
  stat_compare_means(label = "p.format", label.x=1.5, method='wilcox.test') +
  facet_wrap(~celltype, scales = 'free', ncol=5) +
  theme(panel.grid = element_blank(), panel.background = element_blank(),
        axis.line = element_line(color='black'), axis.ticks = element_line(color='black'), legend.position = 'right',
        axis.text.y = element_text(color='black'), axis.text.x = element_blank())
pdf('vlnplot_IRE1.pdf', heigh=10, width=10)
print(p)
dev.off()

