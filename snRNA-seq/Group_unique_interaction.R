library(dplyr)
library(reshape)

read_in <- function(group){
  Mean <- read.table(paste0(outpath, group, '/significant_means.txt'), sep='\t', header = T)
  Mean <- Mean[, -c(3:4,10:12)]
  Mean <- melt(Mean, id=c('id_cp_interaction', 'interacting_pair', 'gene_a', 'gene_b','secreted', 'receptor_a', 'receptor_b'))
  Mean <- plyr::rename(Mean, c(variable='cell_pair', value='mean'))
  Mean <- Mean[!is.na(Mean$mean),]
  
  pvalue <- read.table(paste0(outpath, group, '/pvalues.txt'), sep='\t', header = T)
  pvalue <- pvalue[, -c(3:11)]
  pvalue <- melt(pvalue, id=c('id_cp_interaction', 'interacting_pair'))
  pvalue <- plyr::rename(pvalue, c(variable='cell_pair', value='pvalue'))
  
  df <- left_join(Mean, pvalue[], by=c('id_cp_interaction', 'interacting_pair', 'cell_pair')) %>% mutate(Group=group)
  df[df$pvalue==0,]$pvalue <- 0.001
  
  return(df)
}

old = read_in('O')
patient = read_in('P')


get_group_unique <- function(df1, df2){
  df1$ident <- paste0(df1$id_cp_interaction, '_', df1$cell_pair)
  df2$ident <- paste0(df2$id_cp_interaction, '_', df2$cell_pair)
  overlap_id <- intersect(df1$ident, df2$ident)
  shared <- subset(df1, ident %in% overlap_id)
  unique1 <- subset(df1, ident %in%
                      setdiff(as.character(df1$ident), overlap_id))
  unique1$logP <- log10(unique1$pvalue)*(-1)
  unique1$logmean <- log2(unique1$mean)*(-1)
  
  unique2 <- subset(df2, ident %in%
                      setdiff(as.character(df2$ident), overlap_id))
  unique2$logP <- log10(unique2$pvalue)*(-1)
  unique2$logmean <- log2(unique2$mean)*(-1)
  
  share <- inner_join(df1, df2, by='ident')
  
  return(list(unique1, unique2, share))
}

tmp <- get_group_unique(old, patient)

#shared <- tmp[[3]]
o_unique <- tmp[[1]]
p_unique <- tmp[[2]]
write.csv(o_unique, paste0(result, 'COVID_EC_interaction_Ounique_celltype.csv'), row.names = F)
write.csv(p_unique, paste0(result, 'COVID_EC_interaction_Punique_celltype.csv'), row.names = F)
