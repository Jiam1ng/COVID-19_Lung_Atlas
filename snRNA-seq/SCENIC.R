library(SCENIC)
library(AUCell)
library(RcisTarget)
library(Seurat)

covid <- readRDS("/data1/mashuai/data/COVID-19/Seurat_v4/RDS/Lung.combine.CT.rds")
celltypes <- levels(covid)
deg.up <- c()
for (cell in celltypes){
  tmp <- read.table(paste0("/data1/mashuai/data/COVID-19/Seurat_v4/EC/DEGs_0.5_CT/up/", cell, "_PvsC.up.list"))
  tmp$gene <- rownames(tmp)
  tmp$celltype <- cell
  deg.up <- rbind(deg.up, tmp)
}

deg.down <- c()
for (cell in celltypes){
  tmp <- read.table(paste0("/data1/mashuai/data/COVID-19/Seurat_v4/EC/DEGs_0.5_CT/down/", cell, "_PvsC.down.list"))
  tmp$gene <- rownames(tmp)
  tmp$celltype <- cell
  deg.down <- rbind(deg.down, tmp)
}

used <- subset(covid, stage1 != "Y-Control")

for (cell in celltypes){
  tmp <- subset(used, EC.celltype==cell)
  exp.mat <- as.matrix(GetAssayData(tmp, slot='data'))
  genes <- c(as.character(subset(deg.up, celltype==cell)$gene), as.character(subset(deg.down, celltype==cell)$gene))
  exp.mat <- exp.mat[genes,]
  
  dir.create(paste0('/data5/lijiaming/projects/01_single-cell/10.COVID-19/03_result/SCENIC/EC/PvsO/', cell))
  setwd(paste0('/data5/lijiaming/projects/01_single-cell/10.COVID-19/03_result/SCENIC/EC/PvsO/', cell))
  
  org="hgnc"
  dbDir="/data2/zhengyandong/SCENIC_cisTarget_databases/hg19_cisTarget_databases"
  myDatasetTitle="SCENIC analysis of COVID"
  data(defaultDbNames)
  dbs <- defaultDbNames[[org]]
  scenicOptions <- initializeScenic(org=org, dbDir=dbDir, dbs=dbs, datasetTitle=myDatasetTitle, nCores=10) 
  
  saveRDS(scenicOptions, file="int/scenicOptions.Rds") 
  
  # Gene filter
  genesKept <- geneFiltering(exp.mat, scenicOptions=scenicOptions,
                             minCountsPerGene=3*.01*ncol(exp.mat),
                             minSamples=ncol(exp.mat)*.01)
  exprMat_filtered <- exp.mat[genesKept, ]
  
  # Run Genie3
  runCorrelation(exprMat_filtered, scenicOptions)
  runGenie3(exprMat_filtered, scenicOptions)
  
  # Run the remaining
  scenicOptions@settings$verbose <- TRUE
  scenicOptions@settings$nCores <- 10
  scenicOptions@settings$seed <- 15
  runSCENIC_1_coexNetwork2modules(scenicOptions)
  runSCENIC_2_createRegulons(scenicOptions)
  runSCENIC_3_scoreCells(scenicOptions, exprMat_filtered)
  
  print(paste0(cell, " finished"))
}
