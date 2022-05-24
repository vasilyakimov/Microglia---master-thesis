sheep_data <- list(read.delim('/mnt/cephfs3_rw/exchange/cs02/raid0/tanyagrig/mg_data/sheep_GSM3963879_AB3498.txt', header = TRUE, sep = "\t"),
                   read.delim('/mnt/cephfs3_rw/exchange/cs02/raid0/tanyagrig/mg_data/sheep_GSM3963880_AB5661.txt', header = TRUE, sep = "\t"))

sheep_data <- lapply(sheep_data, FUN = function(x){
  x[,mapply(x, 1, FUN = sum) > 500]})
names(sheep_data) <- paste0('sheep', as.character(1:length(sheep_data)))
sheep_table <- do.call(cbind, sheep_data)

sheep_table <- sheep_table[rownames(sheep_table) %like% 'MT-' == FALSE,]
sheep_table <- sheep_table[rownames(sheep_table) %like% 'ERCC' == FALSE,]
sheep_table <- sheep_table[rownames(sheep_table) %like% 'RP' == FALSE,]
sheep_table <- sheep_table[rownames(sheep_table) %like% 'RNA' == FALSE,]
sheep_table <- sheep_table[rownames(sheep_table) %like% 'COX' == FALSE,]

sheep_seurat <- CreateSeuratObject(sheep_table)
sheep_seurat$dataset <- c(rep('sheep1', ncol(sheep_data$sheep1)),
                          rep('sheep2', ncol(sheep_data$sheep2)))


sheep_seurat <- NormalizeData(sheep_seurat)
sheep_seurat <- FindVariableFeatures(sheep_seurat, selection.method = "dispersion", nfeatures = 2000)
sheep_seurat <- ScaleData(sheep_seurat)
variable.genes <- VariableFeatures(object = sheep_seurat)
sheep_seurat <- RunPCA(sheep_seurat, features = variable.genes)
sheep_seurat <- FindTopFeatures(sheep_seurat, min.cutoff = 'q0')
sheep_seurat <- RunSVD(sheep_seurat)
sheep_seurat <- RunUMAP(object = sheep_seurat, reduction = 'lsi', dims = 2:30)
sheep_seurat <- FindNeighbors(object = sheep_seurat, reduction = 'lsi', dims = 2:30)
sheep_seurat <- FindClusters(object = sheep_seurat, algorithm = 3, resolution = 0.5)


sheep_harmony <- RunHarmony(
  object = sheep_seurat,
  group.by.vars = 'dataset',
  reduction = 'pca',
  assay.use = 'RNA'
)

sheep_harmony <- sheep_harmony %>% 
  RunUMAP(reduction = "harmony", dims = 1:20) %>% 
  FindNeighbors(reduction = "harmony", dims = 1:20) %>% 
  FindClusters(resolution = 0.5)
DimPlot(sheep_harmony, reduction = 'umap', group.by = 'dataset', pt.size = 1) +
  DimPlot(sheep_harmony, reduction = 'umap')

sheep_markers <- FindAllMarkers(sheep_harmony)

sheep_counts <- sheep_harmony@assays$RNA@counts
write.csv(sheep_counts, '/home/PAK-CSPMZ/vakimov/mg_counts/sheep_counts.csv')