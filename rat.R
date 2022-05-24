rat_data <- list(read.delim('/mnt/cephfs3_rw/exchange/cs02/raid0/tanyagrig/mg_data/rat_GSM3963900_AB6890.txt', header = TRUE, sep = "\t"),
                 read.delim('/mnt/cephfs3_rw/exchange/cs02/raid0/tanyagrig/mg_data/rat_GSM3963901_AB6891.txt', header = TRUE, sep = "\t"),
                 read.delim('/mnt/cephfs3_rw/exchange/cs02/raid0/tanyagrig/mg_data/rat_GSM3963902_AB6892.txt', header = TRUE, sep = "\t"),
                 read.delim('/mnt/cephfs3_rw/exchange/cs02/raid0/tanyagrig/mg_data/rat_GSM3963903_AB6893.txt', header = TRUE, sep = "\t"),
                 read.delim('/mnt/cephfs3_rw/exchange/cs02/raid0/tanyagrig/mg_data/rat_GSM3963904_AB6894.txt', header = TRUE, sep = "\t"),
                 read.delim('/mnt/cephfs3_rw/exchange/cs02/raid0/tanyagrig/mg_data/rat_GSM3963905_AB6895.txt', header = TRUE, sep = "\t"))

rat_data <- lapply(rat_data, FUN = function(x){
  x[,mapply(x, 1, FUN = sum) > 500]})
names(rat_data) <- paste0('rat', as.character(1:length(rat_data)))
rat_table <- do.call(cbind, rat_data)

rat_table <- rat_table[rownames(rat_table) %like% 'Mt' == FALSE,]
rat_table <- rat_table[rownames(rat_table) %like% 'mt-' == FALSE,]
rat_table <- rat_table[rownames(rat_table) %like% 'MT-' == FALSE,]
rat_table <- rat_table[rownames(rat_table) %like% 'Ercc' == FALSE,]
rat_table <- rat_table[rownames(rat_table) %like% 'ERCC' == FALSE,]
rat_table <- rat_table[rownames(rat_table) %like% 'ercc' == FALSE,]
rat_table <- rat_table[rownames(rat_table) %like% 'Rp' == FALSE,]
rat_table <- rat_table[rownames(rat_table) %like% 'RP' == FALSE,]
rat_table <- rat_table[rownames(rat_table) %like% 'rp' == FALSE,]
rat_table <- rat_table[rownames(rat_table) %like% 'Rna' == FALSE,]
rat_table <- rat_table[rownames(rat_table) %like% 'rna' == FALSE,]
rat_table <- rat_table[rownames(rat_table) %like% 'RNA' == FALSE,]
rat_table <- rat_table[rownames(rat_table) %like% 'Cox' == FALSE,]
rat_table <- rat_table[rownames(rat_table) %like% 'COX' == FALSE,]

rat_seurat <- CreateSeuratObject(rat_table)
rat_seurat$dataset <- c(rep('rat1', ncol(rat_data$rat1)),
                        rep('rat2', ncol(rat_data$rat2)),
                        rep('rat3', ncol(rat_data$rat3)),
                        rep('rat4', ncol(rat_data$rat4)),
                        rep('rat5', ncol(rat_data$rat5)),
                        rep('rat6', ncol(rat_data$rat6)))

rat_seurat <- NormalizeData(rat_seurat)
rat_seurat <- FindVariableFeatures(rat_seurat, selection.method = "dispersion", nfeatures = 2000)
rat_seurat <- ScaleData(rat_seurat)
variable.genes <- VariableFeatures(object = rat_seurat)
rat_seurat <- RunPCA(rat_seurat, features = variable.genes)
rat_seurat <- FindTopFeatures(rat_seurat, min.cutoff = 'q0')
rat_seurat <- RunSVD(rat_seurat)
rat_seurat <- RunUMAP(object = rat_seurat, reduction = 'lsi', dims = 2:30)
rat_seurat <- FindNeighbors(object = rat_seurat, reduction = 'lsi', dims = 2:30)
rat_seurat <- FindClusters(object = rat_seurat, algorithm = 3, resolution = 0.5)


rat_harmony <- RunHarmony(
  object = rat_seurat,
  group.by.vars = 'dataset',
  reduction = 'pca',
  assay.use = 'RNA'
)

rat_harmony <- rat_harmony %>% 
  RunUMAP(reduction = "harmony", dims = 1:20) %>% 
  FindNeighbors(reduction = "harmony", dims = 1:20) %>% 
  FindClusters(resolution = 0.5)
DimPlot(rat_harmony, reduction = 'umap', group.by = 'dataset', pt.size = 1) +
  DimPlot(rat_harmony, reduction = 'umap')

rat_markers <- FindAllMarkers(rat_harmony)

# rat_counts <- rat_harmony@assays$RNA@counts
# write.csv(mouse_counts, '/home/PAK-CSPMZ/vakimov/mg_counts/rat_counts.csv')



moudules_rat <- list.files("/home/PAK-CSPMZ/vakimov/hotspot_tables", pattern="rat1*", full.names = T)
names(moudules_rat) <- list.files("/home/PAK-CSPMZ/vakimov/hotspot_tables", pattern="rat1*", full.names = F)
moudules_rat <- lapply(moudules_rat, read.csv, header=T)

#modules_human <- do.call(rbind, modules_human)
features1 <- list(moudules_rat$rat1.csv$Gene)
features2 <- list(moudules_rat$rat2.csv$Gene)
features3 <- list(moudules_rat$rat3.csv$Gene)
features4 <- list(moudules_rat$rat4.csv$Gene)
features5 <- list(moudules_rat$rat5.csv$Gene)
features6 <- list(moudules_rat$rat6.csv$Gene)
features7 <- list(moudules_rat$rat7.csv$Gene)

moudules_rat <- AddModuleScore(object = rat_harmony, features = features1, ctrl = 5, name = 'module1')
moudules_rat <- AddModuleScore(object = moudules_rat, features = features2, ctrl = 5, name = 'module2')
moudules_rat <- AddModuleScore(object = moudules_rat, features = features3, ctrl = 5, name = 'module3')
moudules_rat <- AddModuleScore(object = moudules_rat, features = features4, ctrl = 5, name = 'module4')
moudules_rat <- AddModuleScore(object = moudules_rat, features = features5, ctrl = 5, name = 'module5')
moudules_rat <- AddModuleScore(object = moudules_rat, features = features6, ctrl = 5, name = 'module6')
moudules_rat <- AddModuleScore(object = moudules_rat, features = features7, ctrl = 5, name = 'module7')

FeaturePlot(object = moudules_rat, features = 'module11', cols = c('green', 'red'), pt.size = 4) + ggtitle(" rat module 1")
FeaturePlot(object = moudules_rat, features = 'module21', cols = c('green', 'red'), pt.size = 4) + ggtitle(" rat module 2")
FeaturePlot(object = moudules_rat, features = 'module31', cols = c('green', 'red'), pt.size = 4) + ggtitle(" rat module 3")
FeaturePlot(object = moudules_rat, features = 'module41', cols = c('green', 'red'), pt.size = 4) + ggtitle(" rat module 4")
FeaturePlot(object = moudules_rat, features = 'module51', cols = c('green', 'red'), pt.size = 4) + ggtitle(" rat module 5")
FeaturePlot(object = moudules_rat, features = 'module61', cols = c('green', 'red'), pt.size = 4) + ggtitle(" rat module 6")
FeaturePlot(object = moudules_rat, features = 'module71', cols = c('green', 'red'), pt.size = 4) + ggtitle(" rat module 7")
