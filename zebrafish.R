zebrafish_data <- list(read.delim('/mnt/cephfs3_rw/exchange/cs02/raid0/tanyagrig/mg_data/zebra_fish_GSM3963891_AB3014.txt', header = TRUE, sep = "\t"),
                       read.delim('/mnt/cephfs3_rw/exchange/cs02/raid0/tanyagrig/mg_data/zebra_fish_GSM3963892_AB5657.txt', header = TRUE, sep = "\t"),
                       read.delim('/mnt/cephfs3_rw/exchange/cs02/raid0/tanyagrig/mg_data/zebra_fish_GSM3963893_AB7153.txt', header = TRUE, sep = "\t"),
                       read.delim('/mnt/cephfs3_rw/exchange/cs02/raid0/tanyagrig/mg_data/zebra_fish_GSM3963894_AB7154.txt', header = TRUE, sep = "\t"),
                       read.delim('/mnt/cephfs3_rw/exchange/cs02/raid0/tanyagrig/mg_data/zebra_fish_GSM3963895_AB7155.txt', header = TRUE, sep = "\t"))

zebrafish_data <- lapply(zebrafish_data, FUN = function(x){
  x[,mapply(x, 1, FUN = sum) > 500]})
names(zebrafish_data) <- paste0('zebrafish', as.character(1:length(zebrafish_data)))
zebrafish_table <- do.call(cbind, zebrafish_data)


zebrafish_table <- zebrafish_table[rownames(zebrafish_table) %like% 'mt-' == FALSE,]
zebrafish_table <- zebrafish_table[rownames(zebrafish_table) %like% 'ercc' == FALSE,]
zebrafish_table <- zebrafish_table[rownames(zebrafish_table) %like% 'rps' == FALSE,]
zebrafish_table <- zebrafish_table[rownames(zebrafish_table) %like% 'RPS' == FALSE,]
zebrafish_table <- zebrafish_table[rownames(zebrafish_table) %like% 'RPL' == FALSE,]
zebrafish_table <- zebrafish_table[rownames(zebrafish_table) %like% 'rpl' == FALSE,]
zebrafish_table <- zebrafish_table[rownames(zebrafish_table) %like% 'RNA' == FALSE,]
zebrafish_table <- zebrafish_table[rownames(zebrafish_table) %like% 'rna' == FALSE,]
zebrafish_table <- zebrafish_table[rownames(zebrafish_table) %like% 'ERCC' == FALSE,]
zebrafish_table <- zebrafish_table[rownames(zebrafish_table) %like% 'cox' == FALSE,]

zebrafish_seurat <- CreateSeuratObject(zebrafish_table)
zebrafish_seurat$dataset <- c(rep('zebrafish1', ncol(zebrafish_data$zebrafish1)),
                              rep('zebrafish2', ncol(zebrafish_data$zebrafish2)),
                              rep('zebrafish3', ncol(zebrafish_data$zebrafish3)),
                              rep('zebrafish4', ncol(zebrafish_data$zebrafish4)),
                              rep('zebrafish5', ncol(zebrafish_data$zebrafish5)))

zebrafish_seurat <- NormalizeData(zebrafish_seurat)
zebrafish_seurat <- FindVariableFeatures(zebrafish_seurat, selection.method = "dispersion", nfeatures = 2000)
zebrafish_seurat <- ScaleData(zebrafish_seurat)
variable.genes <- VariableFeatures(object = zebrafish_seurat)
zebrafish_seurat <- RunPCA(zebrafish_seurat, features = variable.genes)
zebrafish_seurat <- FindTopFeatures(zebrafish_seurat, min.cutoff = 'q0')
zebrafish_seurat <- RunSVD(zebrafish_seurat)
zebrafish_seurat <- RunUMAP(object = zebrafish_seurat, reduction = 'lsi', dims = 2:30)
zebrafish_seurat <- FindNeighbors(object = zebrafish_seurat, reduction = 'lsi', dims = 2:30)
zebrafish_seurat <- FindClusters(object = zebrafish_seurat, algorithm = 3, resolution = 0.5)


zebrafish_harmony <- RunHarmony(
  object = zebrafish_seurat,
  group.by.vars = 'dataset',
  reduction = 'pca',
  assay.use = 'RNA'
)

zebrafish_harmony <- zebrafish_harmony %>% 
  RunUMAP(reduction = "harmony", dims = 1:20) %>% 
  FindNeighbors(reduction = "harmony", dims = 1:20) %>% 
  FindClusters(resolution = 0.5)
DimPlot(zebrafish_harmony, reduction = 'umap', group.by = 'dataset', pt.size = 1) +
  DimPlot(zebrafish_harmony, reduction = 'umap')

zebrafish_markers <- FindAllMarkers(zebrafish_harmony)
# 
# zebrafish_counts <- zebrafish_harmony@assays$RNA@counts
# write.csv(zebrafish_counts, '/home/PAK-CSPMZ/vakimov/mg_counts/zebrafish_counts.csv')

modules_zebrafish <- list.files("/home/PAK-CSPMZ/vakimov/hotspot_tables", pattern="zebrafish*", full.names = T)
names(modules_zebrafish) <- list.files("/home/PAK-CSPMZ/vakimov/hotspot_tables", pattern="zebrafish*", full.names = F)
modules_zebrafish <- lapply(modules_zebrafish, read.csv, header=T)


features1 <- list(modules_zebrafish$zebrafish1.csv$Gene)
features2 <- list(modules_zebrafish$zebrafish2.csv$Gene)
features3 <- list(modules_zebrafish$zebrafish3.csv$Gene)
features4 <- list(modules_zebrafish$zebrafish4.csv$Gene)
features5 <- list(modules_zebrafish$zebrafish5.csv$Gene)
features6 <- list(modules_zebrafish$zebrafish6.csv$Gene)
features7 <- list(modules_zebrafish$zebrafish7.csv$Gene)
features8 <- list(modules_zebrafish$zebrafish8.csv$Gene)
features9 <- list(modules_zebrafish$zebrafish9.csv$Gene)
features10 <- list(modules_zebrafish$zebrafish10.csv$Gene)
features11 <- list(modules_zebrafish$zebrafish11.csv$Gene)
features12 <- list(modules_zebrafish$zebrafish12.csv$Gene)
features13 <- list(modules_zebrafish$zebrafish13.csv$Gene)
features14 <- list(modules_zebrafish$zebrafish14.csv$Gene)
features15 <- list(modules_zebrafish$zebrafish15.csv$Gene)
features16 <- list(modules_zebrafish$zebrafish16.csv$Gene)
features17 <- list(modules_zebrafish$zebrafish17.csv$Gene)

modules_zebrafish <- AddModuleScore(object = zebrafish_harmony, features = features1, ctrl = 5, name = 'module1')
modules_zebrafish <- AddModuleScore(object = modules_zebrafish, features = features2, ctrl = 5, name = 'module2')
modules_zebrafish <- AddModuleScore(object = modules_zebrafish, features = features3, ctrl = 5, name = 'module3')
modules_zebrafish <- AddModuleScore(object = modules_zebrafish, features = features4, ctrl = 5, name = 'module4')
modules_zebrafish <- AddModuleScore(object = modules_zebrafish, features = features5, ctrl = 5, name = 'module5')
modules_zebrafish <- AddModuleScore(object = modules_zebrafish, features = features6, ctrl = 5, name = 'module6')
modules_zebrafish <- AddModuleScore(object = modules_zebrafish, features = features7, ctrl = 5, name = 'module7')
modules_zebrafish <- AddModuleScore(object = modules_zebrafish, features = features8, ctrl = 5, name = 'module8')
modules_zebrafish <- AddModuleScore(object = modules_zebrafish, features = features9, ctrl = 5, name = 'module9')
modules_zebrafish <- AddModuleScore(object = modules_zebrafish, features = features10, ctrl = 5, name = 'module10')
modules_zebrafish <- AddModuleScore(object = modules_zebrafish, features = features11, ctrl = 5, name = 'module11')
modules_zebrafish <- AddModuleScore(object = modules_zebrafish, features = features12, ctrl = 5, name = 'module12')
modules_zebrafish <- AddModuleScore(object = modules_zebrafish, features = features13, ctrl = 5, name = 'module13')
modules_zebrafish <- AddModuleScore(object = modules_zebrafish, features = features14, ctrl = 5, name = 'module14')
modules_zebrafish <- AddModuleScore(object = modules_zebrafish, features = features15, ctrl = 5, name = 'module15')
modules_zebrafish <- AddModuleScore(object = modules_zebrafish, features = features16, ctrl = 5, name = 'module16')
modules_zebrafish <- AddModuleScore(object = modules_zebrafish, features = features17, ctrl = 5, name = 'module17')

FeaturePlot(object = modules_zebrafish, features = 'module11', cols = c('green', 'red'), pt.size = 4) + ggtitle(" zebrafish module 1")
FeaturePlot(object = modules_zebrafish, features = 'module21', cols = c('green', 'red'), pt.size = 4) + ggtitle(" zebrafish module 2")
FeaturePlot(object = modules_zebrafish, features = 'module31', cols = c('green', 'red'), pt.size = 4) + ggtitle(" zebrafish module 3")
FeaturePlot(object = modules_zebrafish, features = 'module41', cols = c('green', 'red'), pt.size = 4) + ggtitle(" zebrafish module 4")
FeaturePlot(object = modules_zebrafish, features = 'module51', cols = c('green', 'red'), pt.size = 4) + ggtitle(" zebrafish module 5")
FeaturePlot(object = modules_zebrafish, features = 'module61', cols = c('green', 'red'), pt.size = 4) + ggtitle(" zebrafish module 6")
FeaturePlot(object = modules_zebrafish, features = 'module71', cols = c('green', 'red'), pt.size = 4) + ggtitle(" zebrafish module 7")
FeaturePlot(object = modules_zebrafish, features = 'module81', cols = c('green', 'red'), pt.size = 4) + ggtitle(" zebrafish module 8")
FeaturePlot(object = modules_zebrafish, features = 'module91', cols = c('green', 'red'), pt.size = 4) + ggtitle(" zebrafish module 9")
FeaturePlot(object = modules_zebrafish, features = 'module101', cols = c('green', 'red'), pt.size = 4) + ggtitle(" zebrafish module 10")
FeaturePlot(object = modules_zebrafish, features = 'module111', cols = c('green', 'red'), pt.size = 4) + ggtitle(" zebrafish module 11")
FeaturePlot(object = modules_zebrafish, features = 'module121', cols = c('green', 'red'), pt.size = 4) + ggtitle(" zebrafish module 12")
FeaturePlot(object = modules_zebrafish, features = 'module131', cols = c('green', 'red'), pt.size = 4) + ggtitle(" zebrafish module 13")
FeaturePlot(object = modules_zebrafish, features = 'module141', cols = c('green', 'red'), pt.size = 4) + ggtitle(" zebrafish module 14")
FeaturePlot(object = modules_zebrafish, features = 'module151', cols = c('green', 'red'), pt.size = 4) + ggtitle(" zebrafish module 15")
FeaturePlot(object = modules_zebrafish, features = 'module161', cols = c('green', 'red'), pt.size = 4) + ggtitle(" zebrafish module 16")
FeaturePlot(object = modules_zebrafish, features = 'module171', cols = c('green', 'red'), pt.size = 4) + ggtitle(" zebrafish module 17")
