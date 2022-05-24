human_data <- list(read.delim('/mnt/cephfs3_rw/exchange/cs02/raid0/tanyagrig/mg_data/human_GSM3963869_AB3501.txt', header = TRUE, sep = "\t"),
                   read.delim('/mnt/cephfs3_rw/exchange/cs02/raid0/tanyagrig/mg_data/human_GSM3963870_AB3018.txt', header = TRUE, sep = "\t"),
                   read.delim('/mnt/cephfs3_rw/exchange/cs02/raid0/tanyagrig/mg_data/human_GSM3963871_AB5662.txt', header = TRUE, sep = "\t"),
                   read.delim('/mnt/cephfs3_rw/exchange/cs02/raid0/tanyagrig/mg_data/human_GSM3963872_AB5663.txt', header = TRUE, sep = "\t"),
                   read.delim('/mnt/cephfs3_rw/exchange/cs02/raid0/tanyagrig/mg_data/human_GSM3963873_AB5664.txt', header = TRUE, sep = "\t"),
                   read.delim('/mnt/cephfs3_rw/exchange/cs02/raid0/tanyagrig/mg_data/human_GSM3963874_AB5665.txt', header = TRUE, sep = "\t"))

human_data <- lapply(human_data, FUN = function(x){
  x[,mapply(x, 1, FUN = sum) > 500]})
names(human_data) <- paste0('human', as.character(1:length(human_data)))
human_table <- do.call(cbind, human_data)

human_table <- human_table[rownames(human_table) %like% 'MT' == FALSE,]
human_table <- human_table[rownames(human_table) %like% 'ERCC' == FALSE,]
human_table <- human_table[rownames(human_table) %like% 'RP' == FALSE,]
human_table <- human_table[rownames(human_table) %like% 'RNA' == FALSE,]
human_table <- human_table[rownames(human_table) %like% 'COX' == FALSE,]

human_seurat <- CreateSeuratObject(human_table)
human_seurat$dataset <- c(rep('human1', ncol(human_data$human1)),
                          rep('human2', ncol(human_data$human2)),
                          rep('human3', ncol(human_data$human3)),
                          rep('human4', ncol(human_data$human4)),
                          rep('human5', ncol(human_data$human5)),
                          rep('human6', ncol(human_data$human6)))

human_seurat <- NormalizeData(human_seurat)
human_seurat <- FindVariableFeatures(human_seurat, selection.method = "dispersion", nfeatures = 2000)
human_seurat <- ScaleData(human_seurat)
variable.genes <- VariableFeatures(object = human_seurat)
human_seurat <- RunPCA(human_seurat, features = variable.genes)
human_seurat <- FindTopFeatures(human_seurat, min.cutoff = 'q0')
human_seurat <- RunSVD(human_seurat)
human_seurat <- RunUMAP(object = human_seurat, reduction = 'lsi', dims = 2:30)
human_seurat <- FindNeighbors(object = human_seurat, reduction = 'lsi', dims = 2:30)
human_seurat <- FindClusters(object = human_seurat, algorithm = 3, resolution = 0.5)

human_harmony <- RunHarmony(
  object = human_seurat,
  group.by.vars = 'dataset',
  reduction = 'pca',
  assay.use = 'RNA'
)

human_harmony <- human_harmony %>% 
  RunUMAP(reduction = "harmony", dims = 1:20) %>% 
  FindNeighbors(reduction = "harmony", dims = 1:20) %>% 
  FindClusters(resolution = 0.5)

human_markers <- FindAllMarkers(human_harmony)

modules_human <- list.files("/home/PAK-CSPMZ/vakimov/hotspot_tables", pattern="human*", full.names = T)
names(modules_human) <- list.files("/home/PAK-CSPMZ/vakimov/hotspot_tables", pattern="human*", full.names = F)
modules_human <- lapply(modules_human, read.csv, header=T)

#modules_human <- do.call(rbind, modules_human)
features1 <- list(modules_human$human1.csv$Gene)
features2 <- list(modules_human$human2.csv$Gene)
features3 <- list(modules_human$human3.csv$Gene)
features4 <- list(modules_human$human4.csv$Gene)
features5 <- list(modules_human$human5.csv$Gene)
features6 <- list(modules_human$human6.csv$Gene)
features7 <- list(modules_human$human7.csv$Gene)
features8 <- list(modules_human$human8.csv$Gene)
features9 <- list(modules_human$human9.csv$Gene)


human_moduled <- AddModuleScore(object = human_harmony, features = features1, ctrl = 5, name = 'module1')
human_moduled <- AddModuleScore(object = human_moduled, features = features2, ctrl = 5, name = 'module2')
human_moduled <- AddModuleScore(object = human_moduled, features = features3, ctrl = 5, name = 'module3')
human_moduled <- AddModuleScore(object = human_moduled, features = features4, ctrl = 5, name = 'module4')
human_moduled <- AddModuleScore(object = human_moduled, features = features5, ctrl = 5, name = 'module5')
human_moduled <- AddModuleScore(object = human_moduled, features = features6, ctrl = 5, name = 'module6')
human_moduled <- AddModuleScore(object = human_moduled, features = features7, ctrl = 5, name = 'module7')
human_moduled <- AddModuleScore(object = human_moduled, features = features8, ctrl = 5, name = 'module8')
human_moduled <- AddModuleScore(object = human_moduled, features = features9, ctrl = 5, name = 'module9')

FeaturePlot(object = human_moduled, features = 'module11', cols = c('green', 'red'), pt.size = 4) + ggtitle(" human module 1")
FeaturePlot(object = human_moduled, features = 'module21', cols = c('green', 'red'), pt.size = 4) + ggtitle(" human module 2")
FeaturePlot(object = human_moduled, features = 'module31', cols = c('green', 'red'), pt.size = 4) + ggtitle(" human module 3")
FeaturePlot(object = human_moduled, features = 'module41', cols = c('green', 'red'), pt.size = 4) + ggtitle(" human module 4")
FeaturePlot(object = human_moduled, features = 'module51', cols = c('green', 'red'), pt.size = 4) + ggtitle(" human module 5")
FeaturePlot(object = human_moduled, features = 'module61', cols = c('green', 'red'), pt.size = 4) + ggtitle(" human module 6")
FeaturePlot(object = human_moduled, features = 'module71', cols = c('green', 'red'), pt.size = 4) + ggtitle(" human module 7")
FeaturePlot(object = human_moduled, features = 'module81', cols = c('green', 'red'), pt.size = 4) + ggtitle(" human module 8")
FeaturePlot(object = human_moduled, features = 'module91', cols = c('green', 'red'), pt.size = 4) + ggtitle(" human module 9")
