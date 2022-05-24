marmoset_data <- list(read.delim('/mnt/cephfs3_rw/exchange/cs02/raid0/tanyagrig/mg_data/marmoset_GSM3963877_AB3015.txt', header = TRUE, sep = "\t"),
                      read.delim('/mnt/cephfs3_rw/exchange/cs02/raid0/tanyagrig/mg_data/marmoset_GSM3963878_AB5659.txt', header = TRUE, sep = "\t"))

marmoset_data <- lapply(marmoset_data, FUN = function(x){
  x[,mapply(x, 1, FUN = sum) > 500]})
names(marmoset_data) <- paste0('marmoset', as.character(1:length(marmoset_data)))
marmoset_table <- do.call(cbind, marmoset_data)

# marmoset_table <- marmoset_table[rownames(marmoset_table) %like% 'MT' == FALSE,]
# marmoset_table <- marmoset_table[rownames(marmoset_table) %like% 'ERCC' == FALSE,]
# marmoset_table <- marmoset_table[rownames(marmoset_table) %like% 'RP' == FALSE,]
# marmoset_table <- marmoset_table[rownames(marmoset_table) %like% 'RNA' == FALSE,]
# marmoset_table <- marmoset_table[rownames(marmoset_table) %like% 'COX' == FALSE,]

marmoset_seurat[["percent.mt"]] <- PercentageFeatureSet(marmoset_seurat, pattern = "^MT-")
# Visualize QC metrics as a violin plot
VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)


marmoset_ensemble90 = useMart("ENSEMBL_MART_ENSEMBL", 'cjacchus_gene_ensembl', host="https://aug2017.archive.ensembl.org")
marmoset_symbol = 'hgnc_symbol'
marmoset_to_human <- convert_genes(marmoset_table, marmoset_ensemble90, marmoset_symbol, human_ensemble90)

marmoset_seurat <- CreateSeuratObject(marmoset_table)
marmoset_seurat$dataset <- c(rep('marmoset1', ncol(marmoset_data$marmoset1)),
                             rep('marmoset2', ncol(marmoset_data$marmoset2)))

marmoset_seurat <- NormalizeData(marmoset_seurat)
marmoset_seurat <- FindVariableFeatures(marmoset_seurat, selection.method = "dispersion", nfeatures = 2000)
marmoset_seurat <- ScaleData(marmoset_seurat)
variable.genes <- VariableFeatures(object = marmoset_seurat)
marmoset_seurat <- RunPCA(marmoset_seurat, features = variable.genes)
marmoset_seurat <- FindTopFeatures(marmoset_seurat, min.cutoff = 'q0')
marmoset_seurat <- RunSVD(marmoset_seurat)
marmoset_seurat <- RunUMAP(object = marmoset_seurat, reduction = 'lsi', dims = 2:30)
marmoset_seurat <- FindNeighbors(object = marmoset_seurat, reduction = 'lsi', dims = 2:30)
marmoset_seurat <- FindClusters(object = marmoset_seurat, algorithm = 3, resolution = 0.5)

marmoset_harmony <- RunHarmony(
  object = marmoset_seurat,
  group.by.vars = 'dataset',
  reduction = 'pca',
  assay.use = 'RNA'
)

marmoset_harmony <- marmoset_harmony %>% 
  RunUMAP(reduction = "harmony", dims = 1:20) %>% 
  FindNeighbors(reduction = "harmony", dims = 1:20) %>% 
  FindClusters(resolution = 0.5)
DimPlot(marmoset_harmony, reduction = 'umap', group.by = 'dataset', pt.size = 1) +
  DimPlot(marmoset_harmony, reduction = 'umap')

marmoset_markers <- FindAllMarkers(marmoset_harmony)

marmoset_counts <- marmoset_harmony@assays$RNA@counts
write.csv(marmoset_counts, '/home/PAK-CSPMZ/vakimov/mg_counts/marmoset_counts.csv')

