macaque_data <- list(read.delim('/mnt/cephfs3_rw/exchange/cs02/raid0/tanyagrig/mg_data/macaca_GSM3963875_AB3499.txt', header = TRUE, sep = "\t"),
                     read.delim('/mnt/cephfs3_rw/exchange/cs02/raid0/tanyagrig/mg_data/macaca_GSM3963876_AB5660.txt', header = TRUE, sep = "\t"))

macaque_data <- lapply(macaque_data, FUN = function(x){
  x[,mapply(x, 1, FUN = sum) > 500]})
names(macaque_data) <- paste0('macaque', as.character(1:length(macaque_data)))
macaque_table <- do.call(cbind, macaque_data)

macaque_table <- macaque_table[rownames(macaque_table) %like% 'MT-' == FALSE,]
macaque_table <- macaque_table[rownames(macaque_table) %like% 'ERCC' == FALSE,]
macaque_table <- macaque_table[rownames(macaque_table) %like% 'RP' == FALSE,]
macaque_table <- macaque_table[rownames(macaque_table) %like% 'RNA' == FALSE,]
macaque_table <- macaque_table[rownames(macaque_table) %like% 'COX' == FALSE,]

macaque_seurat <- CreateSeuratObject(macaque_table)
macaque_seurat$dataset <- c(rep('macaque1', ncol(macaque_data$macaque1)),
                            rep('macaque2', ncol(macaque_data$macaque2)))

macaque_seurat <- NormalizeData(macaque_seurat)
macaque_seurat <- FindVariableFeatures(macaque_seurat, selection.method = "dispersion", nfeatures = 2000)
macaque_seurat <- ScaleData(macaque_seurat)
variable.genes <- VariableFeatures(object = macaque_seurat)
macaque_seurat <- RunPCA(macaque_seurat, features = variable.genes)
macaque_seurat <- FindTopFeatures(macaque_seurat, min.cutoff = 'q0')
macaque_seurat <- RunSVD(macaque_seurat)
macaque_seurat <- RunUMAP(object = macaque_seurat, reduction = 'lsi', dims = 2:30)
macaque_seurat <- FindNeighbors(object = macaque_seurat, reduction = 'lsi', dims = 2:30)
macaque_seurat <- FindClusters(object = macaque_seurat, algorithm = 3, resolution = 0.5)

macaque_harmony <- RunHarmony(
  object = macaque_seurat,
  group.by.vars = 'dataset',
  reduction = 'pca',
  assay.use = 'RNA'
)

macaque_harmony <- macaque_harmony %>% 
  RunUMAP(reduction = "harmony", dims = 1:20) %>% 
  FindNeighbors(reduction = "harmony", dims = 1:20) %>% 
  FindClusters(resolution = 0.5)
DimPlot(macaque_harmony, reduction = 'umap', group.by = 'dataset', pt.size = 1) +
  DimPlot(macaque_harmony, reduction = 'umap')

macaque_markers <- FindAllMarkers(macaque_harmony)

macaque_counts <- macaque_harmony@assays$RNA@counts
write.csv(macaque_counts, '/home/PAK-CSPMZ/vakimov/mg_counts/macaque_counts.csv')
