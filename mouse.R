library(Signac)
library(Seurat)
library(GenomeInfoDb)
library(EnsDb.Hsapiens.v86)
library(BSgenome.Hsapiens.UCSC.hg38)
library(ggplot2)
library(future)
library(pbapply)
library(parallel)
library(dplyr)
library(rliger)
library(harmony)
library(clusterProfiler)
library(proxy)
library(ComplexHeatmap)
library(data.table)

mouse_data <- list(read.delim('/mnt/cephfs3_rw/exchange/cs02/raid0/tanyagrig/mg_data/mouse_GSM3963881_AB3010.txt', header = TRUE, sep = "\t"),
                   read.delim('/mnt/cephfs3_rw/exchange/cs02/raid0/tanyagrig/mg_data/mouse_GSM3963882_AB3011.txt', header = TRUE, sep = "\t"),
                   read.delim('/mnt/cephfs3_rw/exchange/cs02/raid0/tanyagrig/mg_data/mouse_GSM3963883_AB3012.txt', header = TRUE, sep = "\t"),
                   read.delim('/mnt/cephfs3_rw/exchange/cs02/raid0/tanyagrig/mg_data/mouse_GSM3963884_AB4708.txt', header = TRUE, sep = "\t"),
                   read.delim('/mnt/cephfs3_rw/exchange/cs02/raid0/tanyagrig/mg_data/mouse_GSM3963885_AB4712.txt', header = TRUE, sep = "\t"),
                   read.delim('/mnt/cephfs3_rw/exchange/cs02/raid0/tanyagrig/mg_data/mouse_GSM3963886_AB4713.txt', header = TRUE, sep = "\t"))

mouse_data <- lapply(mouse_data, FUN = function(x){
  x[,mapply(x, 1, FUN = sum) > 500]})
names(mouse_data) <- paste0('mouse', as.character(1:length(mouse_data)))
mouse_table <- do.call(cbind, mouse_data)

mouse_table <- mouse_table[rownames(mouse_table) %like% 'Mt' == FALSE,]
mouse_table <- mouse_table[rownames(mouse_table) %like% 'mt-' == FALSE,]
mouse_table <- mouse_table[rownames(mouse_table) %like% 'Ercc' == FALSE,]
mouse_table <- mouse_table[rownames(mouse_table) %like% 'ERCC' == FALSE,]
mouse_table <- mouse_table[rownames(mouse_table) %like% 'ercc' == FALSE,]
mouse_table <- mouse_table[rownames(mouse_table) %like% 'Rp' == FALSE,]
mouse_table <- mouse_table[rownames(mouse_table) %like% 'RP' == FALSE,]
mouse_table <- mouse_table[rownames(mouse_table) %like% 'rp' == FALSE,]
mouse_table <- mouse_table[rownames(mouse_table) %like% 'Rna' == FALSE,]
mouse_table <- mouse_table[rownames(mouse_table) %like% 'rna' == FALSE,]
mouse_table <- mouse_table[rownames(mouse_table) %like% 'RNA' == FALSE,]
mouse_table <- mouse_table[rownames(mouse_table) %like% 'Cox' == FALSE,]
mouse_table <- mouse_table[rownames(mouse_table) %like% 'COX' == FALSE,]


mouse_seurat <- CreateSeuratObject(mouse_table)
mouse_seurat$dataset <- c(rep('mouse1', ncol(mouse_data$mouse1)),
                          rep('mouse2', ncol(mouse_data$mouse2)),
                          rep('mouse3', ncol(mouse_data$mouse3)),
                          rep('mouse4', ncol(mouse_data$mouse4)),
                          rep('mouse5', ncol(mouse_data$mouse5)),
                          rep('mouse6', ncol(mouse_data$mouse6)))

mouse_seurat <- NormalizeData(mouse_seurat)
mouse_seurat <- FindVariableFeatures(mouse_seurat, selection.method = "dispersion", nfeatures = 2000)
mouse_seurat <- ScaleData(mouse_seurat)
variable.genes <- VariableFeatures(object = mouse_seurat)
mouse_seurat <- RunPCA(mouse_seurat, features = variable.genes)
mouse_seurat <- FindTopFeatures(mouse_seurat, min.cutoff = 'q0')
mouse_seurat <- RunSVD(mouse_seurat)
mouse_seurat <- RunUMAP(object = mouse_seurat, reduction = 'lsi', dims = 2:30)
mouse_seurat <- FindNeighbors(object = mouse_seurat, reduction = 'lsi', dims = 2:30)
mouse_seurat <- FindClusters(object = mouse_seurat, algorithm = 3, resolution = 0.5)

mouse_harmony <- RunHarmony(
  object = mouse_seurat,
  group.by.vars = 'dataset',
  reduction = 'pca',
  assay.use = 'RNA'
)

mouse_harmony <- mouse_harmony %>% 
  RunUMAP(reduction = "harmony", dims = 1:20) %>% 
  FindNeighbors(reduction = "harmony", dims = 1:20) %>% 
  FindClusters(resolution = 0.5)
DimPlot(mouse_harmony, reduction = 'umap', group.by = 'dataset', pt.size = 1) +
  DimPlot(mouse_harmony, reduction = 'umap')

mouse_markers <- FindAllMarkers(mouse_harmony)

mouse_counts <- mouse_harmony@assays$RNA@counts
write.csv(mouse_counts, '/home/PAK-CSPMZ/vakimov/mg_counts/mouse_counts.csv')

par()
plot(mtcars)
VlnPlot(human_harmony, features = c("nFeature_RNA", "nCount_RNA"), group.by = 'orig.ident', y.max = 4000, pt.size = 0)
VlnPlot(mouse_harmony, features = c("nFeature_RNA", "nCount_RNA"), group.by = 'orig.ident', y.max = 4000, pt.size = 0)
VlnPlot(rat_harmony, features = c("nFeature_RNA", "nCount_RNA"), group.by = 'orig.ident', y.max = 4000, pt.size = 0)
VlnPlot(macaque_harmony, features = c("nFeature_RNA", "nCount_RNA"), group.by = 'orig.ident', y.max = 4000, pt.size = 0)
VlnPlot(marmoset_harmony, features = c("nFeature_RNA", "nCount_RNA"), group.by = 'orig.ident', y.max = 4000, pt.size = 0)
VlnPlot(hamster_harmony, features = c("nFeature_RNA", "nCount_RNA"), group.by = 'orig.ident', y.max = 4000, pt.size = 0)
VlnPlot(bmr_harmony, features = c("nFeature_RNA", "nCount_RNA"), group.by = 'orig.ident', y.max = 4000, pt.size = 0)
VlnPlot(sheep_harmony, features = c("nFeature_RNA", "nCount_RNA"), group.by = 'orig.ident', y.max = 4000, pt.size = 0)
VlnPlot(zebrafish_harmony, features = c("nFeature_RNA", "nCount_RNA"), group.by = 'orig.ident', y.max = 4000, pt.size = 0)
VlnPlot(chicken_harmony, features = c("nFeature_RNA", "nCount_RNA"), group.by = 'orig.ident', y.max = 4000, pt.size = 0)

moudules_mouse <- list.files("/home/PAK-CSPMZ/vakimov/hotspot_tables", pattern="mouse*", full.names = T)
names(moudules_mouse) <- list.files("/home/PAK-CSPMZ/vakimov/hotspot_tables", pattern="mouse*", full.names = F)
moudules_mouse <- lapply(moudules_mouse, read.csv, header=T)

#modules_human <- do.call(rbind, modules_human)
features1 <- list(moudules_mouse$mouse1.csv$Gene)
features2 <- list(moudules_mouse$mouse2.csv$Gene)
features3 <- list(moudules_mouse$mouse3.csv$Gene)
features4 <- list(moudules_mouse$mouse4.csv$Gene)
features5 <- list(moudules_mouse$mouse5.csv$Gene)
features6 <- list(moudules_mouse$mouse6.csv$Gene)
features7 <- list(moudules_mouse$mouse7.csv$Gene)

moudules_mouse <- AddModuleScore(object = mouse_harmony, features = features1, ctrl = 5, name = 'module1')
moudules_mouse <- AddModuleScore(object = moudules_mouse, features = features2, ctrl = 5, name = 'module2')
moudules_mouse <- AddModuleScore(object = moudules_mouse, features = features3, ctrl = 5, name = 'module3')
moudules_mouse <- AddModuleScore(object = moudules_mouse, features = features4, ctrl = 5, name = 'module4')
moudules_mouse <- AddModuleScore(object = moudules_mouse, features = features5, ctrl = 5, name = 'module5')
moudules_mouse <- AddModuleScore(object = moudules_mouse, features = features6, ctrl = 5, name = 'module6')
moudules_mouse <- AddModuleScore(object = moudules_mouse, features = features7, ctrl = 5, name = 'module7')

FeaturePlot(object = moudules_mouse, features = 'module11', cols = c('green', 'red'), pt.size = 4) + ggtitle(" mouse module 1")
FeaturePlot(object = moudules_mouse, features = 'module21', cols = c('green', 'red'), pt.size = 4) + ggtitle(" mouse module 2")
FeaturePlot(object = moudules_mouse, features = 'module31', cols = c('green', 'red'), pt.size = 4) + ggtitle(" mouse module 3")
FeaturePlot(object = moudules_mouse, features = 'module41', cols = c('green', 'red'), pt.size = 4) + ggtitle(" mouse module 4")
FeaturePlot(object = moudules_mouse, features = 'module51', cols = c('green', 'red'), pt.size = 4) + ggtitle(" mouse module 5")
FeaturePlot(object = moudules_mouse, features = 'module61', cols = c('green', 'red'), pt.size = 4) + ggtitle(" mouse module 6")
FeaturePlot(object = moudules_mouse, features = 'module71', cols = c('green', 'red'), pt.size = 4) + ggtitle(" mouse module 7")

