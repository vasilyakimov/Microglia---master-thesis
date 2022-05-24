hamster_data <- list(read.delim('/mnt/cephfs3_rw/exchange/cs02/raid0/tanyagrig/mg_data/hamster_GSM3963887_AB3016.txt', header = TRUE, sep = "\t"),
                     read.delim('/mnt/cephfs3_rw/exchange/cs02/raid0/tanyagrig/mg_data/hamster_GSM3963888_AB5656.txt', header = TRUE, sep = "\t"))

hamster_data <- lapply(hamster_data, FUN = function(x){
  x[,mapply(x, 1, FUN = sum) > 500]})
names(hamster_data) <- paste0('hamster', as.character(1:length(hamster_data)))
hamster_table <- do.call(cbind, hamster_data)

hamster_table <- hamster_table[rownames(hamster_table) %like% 'Mt' == FALSE,]
hamster_table <- hamster_table[rownames(hamster_table) %like% 'Ercc' == FALSE,]
hamster_table <- hamster_table[rownames(hamster_table) %like% 'ercc' == FALSE,]
hamster_table <- hamster_table[rownames(hamster_table) %like% 'Rp' == FALSE,]
hamster_table <- hamster_table[rownames(hamster_table) %like% 'rp' == FALSE,]
hamster_table <- hamster_table[rownames(hamster_table) %like% 'Rna' == FALSE,]
hamster_table <- hamster_table[rownames(hamster_table) %like% 'rna' == FALSE,]
hamster_table <- hamster_table[rownames(hamster_table) %like% 'ERCC' == FALSE,]
hamster_table <- hamster_table[rownames(hamster_table) %like% 'Cox' == FALSE,]
hamster_table <- hamster_table[rownames(hamster_table) %like% 'COX' == FALSE,]

hamster_seurat <- CreateSeuratObject(hamster_table)
hamster_seurat$dataset <- c(rep('hamster1', ncol(hamster_data$hamster1)),
                            rep('hamster2', ncol(hamster_data$hamster2)))

hamster_seurat <- NormalizeData(hamster_seurat)
hamster_seurat <- FindVariableFeatures(hamster_seurat, selection.method = "dispersion", nfeatures = 2000)
hamster_seurat <- ScaleData(hamster_seurat)
variable.genes <- VariableFeatures(object = hamster_seurat)
hamster_seurat <- RunPCA(hamster_seurat, features = variable.genes)
hamster_seurat <- FindTopFeatures(hamster_seurat, min.cutoff = 'q0')
hamster_seurat <- RunSVD(hamster_seurat)
hamster_seurat <- RunUMAP(object = hamster_seurat, reduction = 'lsi', dims = 2:30)
hamster_seurat <- FindNeighbors(object = hamster_seurat, reduction = 'lsi', dims = 2:30)
hamster_seurat <- FindClusters(object = hamster_seurat, algorithm = 3, resolution = 0.5)


hamster_harmony <- RunHarmony(
  object = hamster_seurat,
  group.by.vars = 'dataset',
  reduction = 'pca',
  assay.use = 'RNA'
)

hamster_harmony <- hamster_harmony %>% 
  RunUMAP(reduction = "harmony", dims = 1:20) %>% 
  FindNeighbors(reduction = "harmony", dims = 1:20) %>% 
  FindClusters(resolution = 0.5)
DimPlot(hamster_harmony, reduction = 'umap', group.by = 'dataset', pt.size = 1) +
  DimPlot(hamster_harmony, reduction = 'umap')

hamster_markers <- FindAllMarkers(hamster_harmony)

hamster_counts <- hamster_harmony@assays$RNA@counts
write.csv(hamster_counts, '/home/PAK-CSPMZ/vakimov/mg_counts/hamster_counts.csv')
