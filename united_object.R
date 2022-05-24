
library("data.table")                              # Load data.table
library("Seurat") 

mouse_data <- cbind(read.delim('/mnt/cephfs3_rw/exchange/cs02/raid0/tanyagrig/mg_data/mouse_GSM3963881_AB3010.txt', header = TRUE, sep = "\t"),
                   read.delim('/mnt/cephfs3_rw/exchange/cs02/raid0/tanyagrig/mg_data/mouse_GSM3963882_AB3011.txt', header = TRUE, sep = "\t"),
                   read.delim('/mnt/cephfs3_rw/exchange/cs02/raid0/tanyagrig/mg_data/mouse_GSM3963883_AB3012.txt', header = TRUE, sep = "\t"),
                   read.delim('/mnt/cephfs3_rw/exchange/cs02/raid0/tanyagrig/mg_data/mouse_GSM3963884_AB4708.txt', header = TRUE, sep = "\t"),
                   read.delim('/mnt/cephfs3_rw/exchange/cs02/raid0/tanyagrig/mg_data/mouse_GSM3963885_AB4712.txt', header = TRUE, sep = "\t"),
                   read.delim('/mnt/cephfs3_rw/exchange/cs02/raid0/tanyagrig/mg_data/mouse_GSM3963886_AB4713.txt', header = TRUE, sep = "\t"))

human_data <- cbind(read.delim('/mnt/raid0/tanyagrig/mg_data/human_GSM3963869_AB3501.txt', header = TRUE, sep = "\t"),
                   read.delim('/mnt/raid0/tanyagrig/mg_data/human_GSM3963870_AB3018.txt', header = TRUE, sep = "\t"),
                   read.delim('/mnt/raid0/tanyagrig/mg_data/human_GSM3963871_AB5662.txt', header = TRUE, sep = "\t"),
                   read.delim('/mnt/raid0/tanyagrig/mg_data/human_GSM3963872_AB5663.txt', header = TRUE, sep = "\t"),
                   read.delim('/mnt/raid0/tanyagrig/mg_data/human_GSM3963873_AB5664.txt', header = TRUE, sep = "\t"),
                   read.delim('/mnt/raid0/tanyagrig/mg_data/human_GSM3963874_AB5665.txt', header = TRUE, sep = "\t"))

bmr_data <- cbind(read.delim('/mnt/raid0/tanyagrig/mg_data/bmr_GSM3963896_AB6896.txt', header = TRUE, sep = "\t"),
                 read.delim('/mnt/raid0/tanyagrig/mg_data/bmr_GSM3963897_AB6897.txt', header = TRUE, sep = "\t"),
                 read.delim('/mnt/raid0/tanyagrig/mg_data/bmr_GSM3963898_AB6898.txt', header = TRUE, sep = "\t"),
                 read.delim('/mnt/raid0/tanyagrig/mg_data/bmr_GSM3963899_AB6899.txt', header = TRUE, sep = "\t"))

chicken_data <- cbind(read.delim('/mnt/raid0/tanyagrig/mg_data/chicken_GSM3963889_AB3013.txt', header = TRUE, sep = "\t"),
                     read.delim('/mnt/raid0/tanyagrig/mg_data/chicken_GSM3963890_AB5658.txt', header = TRUE, sep = "\t"))

hamster_data <- cbind(read.delim('/mnt/raid0/tanyagrig/mg_data/hamster_GSM3963887_AB3016.txt', header = TRUE, sep = "\t"),
                     read.delim('/mnt/raid0/tanyagrig/mg_data/hamster_GSM3963888_AB5656.txt', header = TRUE, sep = "\t"))

macaque_data <- cbind(read.delim('/mnt/raid0/tanyagrig/mg_data/macaca_GSM3963875_AB3499.txt', header = TRUE, sep = "\t"),
                     read.delim('/mnt/raid0/tanyagrig/mg_data/macaca_GSM3963876_AB5660.txt', header = TRUE, sep = "\t"))

marmoset_data <- cbind(read.delim('/mnt/raid0/tanyagrig/mg_data/marmoset_GSM3963877_AB3015.txt', header = TRUE, sep = "\t"),
                      read.delim('/mnt/raid0/tanyagrig/mg_data/marmoset_GSM3963878_AB5659.txt', header = TRUE, sep = "\t"))

sheep_data <- cbind(read.delim('/mnt/raid0/tanyagrig/mg_data/sheep_GSM3963879_AB3498.txt', header = TRUE, sep = "\t"),
                   read.delim('/mnt/raid0/tanyagrig/mg_data/sheep_GSM3963880_AB5661.txt', header = TRUE, sep = "\t"))

zebrafish_data <- cbind(read.delim('/mnt/raid0/tanyagrig/mg_data/zebra_fish_GSM3963891_AB3014.txt', header = TRUE, sep = "\t"),
                       read.delim('/mnt/raid0/tanyagrig/mg_data/zebra_fish_GSM3963892_AB5657.txt', header = TRUE, sep = "\t"),
                       read.delim('/mnt/raid0/tanyagrig/mg_data/zebra_fish_GSM3963893_AB7153.txt', header = TRUE, sep = "\t"),
                       read.delim('/mnt/raid0/tanyagrig/mg_data/zebra_fish_GSM3963894_AB7154.txt', header = TRUE, sep = "\t"),
                       read.delim('/mnt/raid0/tanyagrig/mg_data/zebra_fish_GSM3963895_AB7155.txt', header = TRUE, sep = "\t"))

rat_data <- cbind(read.delim('/mnt/raid0/tanyagrig/mg_data/rat_GSM3963900_AB6890.txt', header = TRUE, sep = "\t"),
                 read.delim('/mnt/raid0/tanyagrig/mg_data/rat_GSM3963901_AB6891.txt', header = TRUE, sep = "\t"),
                 read.delim('/mnt/raid0/tanyagrig/mg_data/rat_GSM3963902_AB6892.txt', header = TRUE, sep = "\t"),
                 read.delim('/mnt/raid0/tanyagrig/mg_data/rat_GSM3963903_AB6893.txt', header = TRUE, sep = "\t"),
                 read.delim('/mnt/raid0/tanyagrig/mg_data/rat_GSM3963904_AB6894.txt', header = TRUE, sep = "\t"),
                 read.delim('/mnt/raid0/tanyagrig/mg_data/rat_GSM3963905_AB6895.txt', header = TRUE, sep = "\t"))

# filtering cells with less than 500 UMIs

human_data <- human_data[,mapply(human_data, 1, FUN = sum) > 500]
mouse_data <- mouse_data[,mapply(mouse_data, 1, FUN = sum) > 500]
chicken_data <- chicken_data[,mapply(chicken_data, 1, FUN = sum) > 500]
marmoset_data <- marmoset_data[,mapply(marmoset_data, 1, FUN = sum) > 500]
macaque_data <- macaque_data[,mapply(macaque_data, 1, FUN = sum) > 500]
sheep_data <- sheep_data[,mapply(sheep_data, 1, FUN = sum) > 500]
zebrafish_data <- zebrafish_data[,mapply(zebrafish_data, 1, FUN = sum) > 500]
hamster_data <- hamster_data[,mapply(hamster_data, 1, FUN = sum) > 500]
rat_data <- rat_data[,mapply(rat_data, 1, FUN = sum) > 500]


# Basic function to convert mouse to human gene names
convert_genes <- function(x, y, z, h){
  require("biomaRt")
  genesV2 = getLDS(attributes = z, filters = z, values = rownames(x) , mart = y, attributesL = c("hgnc_symbol"), martL = h, uniqueRows=T)
  x1 <- genesV2[,1]
  d1 <- unique(x1[-pmatch(unique(genesV2[,1]),genesV2[,1])])
  genes <- genesV2[genesV2[,1] %in% d1 == FALSE,]
  x2 <- genes[,2]
  d2 <- unique(x2[-pmatch(unique(genes[,2]),genes[,2])])
  genes <- genes[genes[,2] %in% d2 == FALSE,]
  to_human <- x[genes[,1],]
  rownames(to_human) <- genes[,2]
  to_human
}

human_ensemble95 = useMart("ENSEMBL_MART_ENSEMBL", 'hsapiens_gene_ensembl', host="https://jan2019.archive.ensembl.org")
human_ensemble90 = useMart("ENSEMBL_MART_ENSEMBL", 'hsapiens_gene_ensembl', host="https://aug2017.archive.ensembl.org")
human_ensemble96 = useMart("ENSEMBL_MART_ENSEMBL", 'hsapiens_gene_ensembl', host="https://apr2019.archive.ensembl.org")
mouse_ensemble95 = useMart("ENSEMBL_MART_ENSEMBL", 'mmusculus_gene_ensembl', host="https://jan2019.archive.ensembl.org")
mouse_ensemble90 = useMart("ENSEMBL_MART_ENSEMBL", 'mmusculus_gene_ensembl', host="https://aug2017.archive.ensembl.org")
chicken_ensemble90 = useMart("ENSEMBL_MART_ENSEMBL", 'ggallus_gene_ensembl', host="https://aug2017.archive.ensembl.org")
marmoset_ensemble90 = useMart("ENSEMBL_MART_ENSEMBL", 'cjacchus_gene_ensembl', host="https://aug2017.archive.ensembl.org")
macaque_ensemble90 = useMart("ENSEMBL_MART_ENSEMBL",'mmulatta_gene_ensembl', host="https://aug2017.archive.ensembl.org")
sheep_ensemble90 = useMart("ENSEMBL_MART_ENSEMBL", 'oaries_gene_ensembl', host="https://aug2017.archive.ensembl.org")
zebrafish_ensemble90 = useMart("ENSEMBL_MART_ENSEMBL", 'drerio_gene_ensembl', host="https://aug2017.archive.ensembl.org")
hamster_ensemble90 =useMart("ENSEMBL_MART_ENSEMBL", 'mauratus_gene_ensembl', host="https://aug2017.archive.ensembl.org")
rat_ensemble95 = useMart("ENSEMBL_MART_ENSEMBL", 'rnorvegicus_gene_ensembl', host="https://jan2019.archive.ensembl.org")
bmr_ensemble96 = useMart("ENSEMBL_MART_ENSEMBL", 'ngalili_gene_ensembl', host="https://apr2019.archive.ensembl.org")
bmr_ensemble95 = useMart("ENSEMBL_MART_ENSEMBL", 'ngalili_gene_ensembl', host="https://jan2019.archive.ensembl.org")

human_symbol = 'hgnc_symbol'
mouse_symbol <- 'mgi_symbol'
rat_symbol <- 'rgd_symbol'
chicken_symbol = 'hgnc_symbol'
marmoset_symbol = 'hgnc_symbol'
macaque_symbol = 'hgnc_symbol'
sheep_symbol = 'hgnc_symbol'
zebrafish_symbol = 'zfin_id_symbol'
hamster_symbol <- 'mgi_symbol'
bmr_symbol <- 'mgi_symbol'

mouse_to_human <- convert_genes(mouse_data, mouse_ensemble95, mouse_symbol, human_ensemble95)
chicken_to_human <- convert_genes(chicken_data, chicken_ensemble90, chicken_symbol, human_ensemble90)
marmoset_to_human <- convert_genes(marmoset_data, marmoset_ensemble90, marmoset_symbol, human_ensemble90)
macaque_to_human <- convert_genes(macaque_data, macaque_ensemble90, macaque_symbol, human_ensemble90)
sheep_to_human <- convert_genes(sheep_data, sheep_ensemble90, sheep_symbol, human_ensemble90)
zebrafish_to_human <- convert_genes(zebrafish_data, zebrafish_ensemble90, zebrafish_symbol, human_ensemble90)
zebrafish_to_human <- zebrafish_to_human[complete.cases(zebrafish_to_human),]
hamster_to_human <- convert_genes(hamster_data, hamster_ensemble90, hamster_symbol, human_ensemble90)
rat_to_human <- convert_genes(rat_data, rat_ensemble95, rat_symbol, human_ensemble95)

mouse_seurat <- CreateSeuratObject(mouse_to_human)
human_seurat <- CreateSeuratObject(human_data)
chicken_seurat <- CreateSeuratObject(chicken_to_human)
marmoset_seurat <- CreateSeuratObject(marmoset_to_human)
macaque_seurat <- CreateSeuratObject(macaque_to_human)
sheep_seurat <- CreateSeuratObject(sheep_to_human)
zebrafish_seurat <- CreateSeuratObject(zebrafish_to_human)
hamster_seurat <- CreateSeuratObject(hamster_to_human)
rat_seurat <- CreateSeuratObject(rat_to_human)

mouse_seurat$dataset <- 'mouse'
human_seurat$dataset <- 'human'
chicken_seurat$dataset <- 'chicken'
marmoset_seurat$dataset <- 'marmoset'
macaque_seurat$dataset <- 'macaque'
sheep_seurat$dataset <- 'sheep'
zebrafish_seurat$dataset <- 'zebrafish'
hamster_seurat$dataset <- 'hamster'
rat_seurat$dataset <- 'rat'

seurat <- merge(x = human_seurat,
                y = c(mouse_seurat, chicken_seurat, marmoset_seurat, macaque_seurat, sheep_seurat, zebrafish_seurat, hamster_seurat, rat_seurat))

seurat <- NormalizeData(seurat)
seurat <- FindVariableFeatures(seurat, selection.method = "dispersion", nfeatures = 2000)
seurat <- ScaleData(seurat)
variable.genes <- VariableFeatures(object = seurat)
seurat <- RunPCA(seurat, features = variable.genes)
seurat <- FindTopFeatures(seurat, min.cutoff = 'q0')
seurat <- RunSVD(seurat)
seurat <- RunUMAP(object = seurat, reduction = 'lsi', dims = 2:30)
seurat <- FindNeighbors(object = seurat, reduction = 'lsi', dims = 2:30)
seurat <- FindClusters(object = seurat, algorithm = 3, resolution = 0.5)

harmony <- RunHarmony(
  object = seurat,
  group.by.vars = 'dataset',
  reduction = 'pca',
  assay.use = 'RNA'
)

harmony <- harmony %>% 
  RunUMAP(reduction = "harmony", dims = 1:20) %>% 
  FindNeighbors(reduction = "harmony", dims = 1:20) %>% 
  FindClusters(resolution = 0.5)
DimPlot(harmony, reduction = 'umap', group.by = 'dataset') +
  DimPlot(harmony, reduction = 'umap')

gene_list <- c('CD74', 'CX3CR1', 'CSF1R', 'CSTB', 'C3', 'B2M', 'A2M', 'P2RY12', 'P2RY13', 'APOE', 'SPP1', 'CD81',
               "CSF1R", 'CD68', 'C1QA', 'C1QB', 'C1QC', 'FCER1G', 'CTSS', 'P2RY12', 'SELENOP', 'FCER1G', 'LAPTM5',
               'PLD4', 'PSAP', 'TPT1', 'SAMSN1', 'BTGI', 'ITM2B', 'HEXB', 'SPARC', 'GRP34', 'IFNGR1')

DoHeatmap(
  harmony,
  features = gene_list,
  group.by = "dataset",
  slot = 'data'
)
