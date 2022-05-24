chicken_data <- list(read.delim('/mnt/cephfs3_rw/exchange/cs02/raid0/tanyagrig/mg_data/chicken_GSM3963889_AB3013.txt', header = TRUE, sep = "\t"),
                     read.delim('/mnt/cephfs3_rw/exchange/cs02/raid0/tanyagrig/mg_data/chicken_GSM3963890_AB5658.txt', header = TRUE, sep = "\t"))

chicken_data <- lapply(chicken_data, FUN = function(x){
  x[,mapply(x, 1, FUN = sum) > 500]})
names(chicken_data) <- paste0('chicken', as.character(1:length(chicken_data)))
chicken_table <- do.call(cbind, chicken_data)

chicken_table <- chicken_table[rownames(chicken_table) %like% 'MT' == FALSE,]
chicken_table <- chicken_table[rownames(chicken_table) %like% 'ERCC' == FALSE,]
chicken_table <- chicken_table[rownames(chicken_table) %like% 'RP' == FALSE,]
chicken_table <- chicken_table[rownames(chicken_table) %like% 'RNA' == FALSE,]
chicken_table <- chicken_table[rownames(chicken_table) %like% 'COX' == FALSE,]

chicken_seurat <- CreateSeuratObject(chicken_table)
chicken_seurat$dataset <- c(rep('chicken1', ncol(chicken_data$chicken1)),
                            rep('chicken2', ncol(chicken_data$chicken2)))

chicken_seurat <- NormalizeData(chicken_seurat)
chicken_seurat <- FindVariableFeatures(chicken_seurat, selection.method = "dispersion", nfeatures = 2000)
chicken_seurat <- ScaleData(chicken_seurat)
variable.genes <- VariableFeatures(object = chicken_seurat)
chicken_seurat <- RunPCA(chicken_seurat, features = variable.genes)
chicken_seurat <- FindTopFeatures(chicken_seurat, min.cutoff = 'q0')
chicken_seurat <- RunSVD(chicken_seurat)
chicken_seurat <- RunUMAP(object = chicken_seurat, reduction = 'lsi', dims = 2:30)
chicken_seurat <- FindNeighbors(object = chicken_seurat, reduction = 'lsi', dims = 2:30)
chicken_seurat <- FindClusters(object = chicken_seurat, algorithm = 3, resolution = 0.5)


chicken_harmony <- RunHarmony(
  object = chicken_seurat,
  group.by.vars = 'dataset',
  reduction = 'pca',
  assay.use = 'RNA'
)

chicken_harmony <- chicken_harmony %>% 
  RunUMAP(reduction = "harmony", dims = 1:20) %>% 
  FindNeighbors(reduction = "harmony", dims = 1:20) %>% 
  FindClusters(resolution = 0.5)
DimPlot(chicken_harmony, reduction = 'umap', group.by = 'dataset', pt.size = 1) +
  DimPlot(chicken_harmony, reduction = 'umap')

chicken_markers <- FindAllMarkers(chicken_harmony)

moudules_chicken <- list.files("/home/PAK-CSPMZ/vakimov/hotspot_tables", pattern="chicken*", full.names = T)
names(moudules_chicken) <- list.files("/home/PAK-CSPMZ/vakimov/hotspot_tables", pattern="chicken*", full.names = F)
modules_chicken <- lapply(moudules_chicken, read.csv, header=T)

features1 <- list(moudules_chicken$chicken1.csv$Gene)
features2 <- list(moudules_chicken$chicken1.csv$Gene)
features3 <- list(moudules_chicken$chicken1.csv$Gene)


chicken_moduled <- AddModuleScore(object = chicken_harmony, features = features1, ctrl = 5, name = 'module1')
chicken_moduled <- AddModuleScore(object = chicken_moduled, features = features2, ctrl = 5, name = 'module2')
chicken_moduled <- AddModuleScore(object = chicken_moduled, features = features3, ctrl = 5, name = 'module3')

FeaturePlot(object = chicken_moduled, features = 'module11', cols = c('green', 'red'), pt.size = 4) + ggtitle(" chicken module 1") 
FeaturePlot(object = chicken_moduled, features = 'module21', cols = c('green', 'red'), pt.size = 4) + ggtitle(" chicken module 2") 
FeaturePlot(object = chicken_moduled, features = 'module31', cols = c('green', 'red'), pt.size = 4) + ggtitle(" chicken module 3")
