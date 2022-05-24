bmr_data <- list(read.delim('/mnt/cephfs3_rw/exchange/cs02/raid0/tanyagrig/mg_data/bmr_GSM3963896_AB6896.txt', header = TRUE, sep = "\t"),
                 read.delim('/mnt/cephfs3_rw/exchange/cs02/raid0/tanyagrig/mg_data/bmr_GSM3963897_AB6897.txt', header = TRUE, sep = "\t"),
                 read.delim('/mnt/cephfs3_rw/exchange/cs02/raid0/tanyagrig/mg_data/bmr_GSM3963898_AB6898.txt', header = TRUE, sep = "\t"),
                 read.delim('/mnt/cephfs3_rw/exchange/cs02/raid0/tanyagrig/mg_data/bmr_GSM3963899_AB6899.txt', header = TRUE, sep = "\t"))

bmr_data <- lapply(bmr_data, FUN = function(x){
  x[,mapply(x, 1, FUN = sum) > 500]})
names(bmr_data) <- paste0('bmr', as.character(1:length(bmr_data)))
bmr_table <- do.call(cbind, bmr_data)


bmr_table <- bmr_table[rownames(bmr_table) %like% 'Mt' == FALSE,]
bmr_table <- bmr_table[rownames(bmr_table) %like% 'Ercc' == FALSE,]
bmr_table <- bmr_table[rownames(bmr_table) %like% 'ercc' == FALSE,]
bmr_table <- bmr_table[rownames(bmr_table) %like% 'Rp' == FALSE,]
bmr_table <- bmr_table[rownames(bmr_table) %like% 'rp' == FALSE,]
bmr_table <- bmr_table[rownames(bmr_table) %like% 'Rna' == FALSE,]
bmr_table <- bmr_table[rownames(bmr_table) %like% 'rna' == FALSE,]
bmr_table <- bmr_table[rownames(bmr_table) %like% 'ERCC' == FALSE,]
bmr_table <- bmr_table[rownames(bmr_table) %like% 'Cox' == FALSE,]
bmr_table <- bmr_table[rownames(bmr_table) %like% 'COX' == FALSE,]

bmr_seurat <- CreateSeuratObject(bmr_table)
bmr_seurat$dataset <- c(rep('bmr1', ncol(bmr_data$bmr1)),
                        rep('bmr2', ncol(bmr_data$bmr2)),
                        rep('bmr3', ncol(bmr_data$bmr3)),
                        rep('bmr4', ncol(bmr_data$bmr4)))

bmr_seurat <- NormalizeData(bmr_seurat)
bmr_seurat <- FindVariableFeatures(bmr_seurat, selection.method = "dispersion", nfeatures = 2000)
bmr_seurat <- ScaleData(bmr_seurat)
variable.genes <- VariableFeatures(object = bmr_seurat)
bmr_seurat <- RunPCA(bmr_seurat, features = variable.genes)
bmr_seurat <- FindTopFeatures(bmr_seurat, min.cutoff = 'q0')
bmr_seurat <- RunSVD(bmr_seurat)
bmr_seurat <- RunUMAP(object = bmr_seurat, reduction = 'lsi', dims = 2:30)
bmr_seurat <- FindNeighbors(object = bmr_seurat, reduction = 'lsi', dims = 2:30)
bmr_seurat <- FindClusters(object = bmr_seurat, algorithm = 3, resolution = 0.5)

bmr_harmony <- RunHarmony(
  object = bmr_seurat,
  group.by.vars = 'dataset',
  reduction = 'pca',
  assay.use = 'RNA'
)

bmr_harmony <- bmr_harmony %>% 
  RunUMAP(reduction = "harmony", dims = 1:20) %>% 
  FindNeighbors(reduction = "harmony", dims = 1:20) %>% 
  FindClusters(resolution = 0.5)
DimPlot(bmr_harmony, reduction = 'umap', group.by = 'dataset', pt.size = 1) +
  DimPlot(bmr_harmony, reduction = 'umap')

bmr_markers <- FindAllMarkers(bmr_harmony)
bmr_counts <- bmr_harmony@assays$RNA@counts
write.csv(bmr_counts, '/home/PAK-CSPMZ/vakimov/mg_counts/bmr_counts.csv')
