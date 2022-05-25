# Download module tables

modules_chicken <- list.files("/home/PAK-CSPMZ/vakimov/hotspot_tables", pattern="chicken*", full.names = T)
names(modules_chicken) <- list.files("/home/PAK-CSPMZ/vakimov/hotspot_tables", pattern="*chicken*", full.names = F)
modules_chicken <- lapply(modules_chicken, read.csv, header=T)

modules_mouse <- list.files("/home/PAK-CSPMZ/vakimov/hotspot_tables", pattern="mouse*", full.names = T)
names(modules_mouse) <- list.files("/home/PAK-CSPMZ/vakimov/hotspot_tables", pattern="mouse*", full.names = F)
modules_mouse <- lapply(modules_mouse, read.csv, header=T)

modules_zebrafish <- list.files("/home/PAK-CSPMZ/vakimov/hotspot_tables", pattern="zebrafish*", full.names = T)
names(modules_zebrafish) <- list.files("/home/PAK-CSPMZ/vakimov/hotspot_tables", pattern="zebrafish*", full.names = F)
modules_zebrafish <- lapply(modules_zebrafish, read.csv, header=T)

modules_rat <- list.files("/home/PAK-CSPMZ/vakimov/hotspot_tables", pattern="rat *", full.names = T)
names(modules_rat) <- list.files("/home/PAK-CSPMZ/vakimov/hotspot_tables", pattern="rat1*", full.names = F)
modules_rat <- lapply(modules_rat, read.csv, header=T)

modules_human <- list.files("/home/PAK-CSPMZ/vakimov/hotspot_tables", pattern="human*", full.names = T)
names(modules_human) <- list.files("/home/PAK-CSPMZ/vakimov/hotspot_tables", pattern="human*", full.names = F)
modules_human <- lapply(modules_human, read.csv, header=T)

# Basic function to convert gene names
convert_genes <- function(x, y, z, h){
  require("biomaRt")
  genesV2 = getLDS(attributes = z, filters = z, values = x$Gene , mart = y, attributesL = c("hgnc_symbol"), martL = h, uniqueRows=T)
  human_genes <- genesV2[,2]
  human_genes
}

# Loading ensembl data

human_ensemble95 = useMart("ENSEMBL_MART_ENSEMBL", 'hsapiens_gene_ensembl', host="https://jan2019.archive.ensembl.org")
human_ensemble90 = useMart("ENSEMBL_MART_ENSEMBL", 'hsapiens_gene_ensembl', host="https://aug2017.archive.ensembl.org")
mouse_ensemble95 = useMart("ENSEMBL_MART_ENSEMBL", 'mmusculus_gene_ensembl', host="https://jan2019.archive.ensembl.org")
mouse_ensemble90 = useMart("ENSEMBL_MART_ENSEMBL", 'mmusculus_gene_ensembl', host="https://aug2017.archive.ensembl.org")
chicken_ensemble90 = useMart("ENSEMBL_MART_ENSEMBL", 'ggallus_gene_ensembl', host="https://aug2017.archive.ensembl.org")
zebrafish_ensemble90 = useMart("ENSEMBL_MART_ENSEMBL", 'drerio_gene_ensembl', host="https://aug2017.archive.ensembl.org")
rat_ensemble95 = useMart("ENSEMBL_MART_ENSEMBL", 'rnorvegicus_gene_ensembl', host="https://jan2019.archive.ensembl.org")

human_symbol = 'hgnc_symbol'
mouse_symbol <- 'mgi_symbol'
rat_symbol <- 'rgd_symbol'
chicken_symbol = 'hgnc_symbol'
zebrafish_symbol = 'zfin_id_symbol'
bmr_symbol <- 'mgi_symbol'

# converting genes

human_converted <- lapply(modules_human, function(x){ x$Gene })
mouse_converted <- lapply(modules_mouse, function(x){ convert_genes(x, mouse_ensemble95, mouse_symbol, human_ensemble95) })
rat_converted <- lapply(modules_rat, function(x){ convert_genes(x, rat_ensemble95, rat_symbol, human_ensemble95) })
zebrafish_converted <- lapply(modules_zebrafish, function(x){ convert_genes(x, zebrafish_ensemble90, zebrafish_symbol, human_ensemble90) })
chicken_converted <- lapply(modules_chicken, function(x){ convert_genes(x, chicken_ensemble90, chicken_symbol, human_ensemble90) })

converted <- c(human_converted, mouse_converted, rat_converted, zebrafish_converted, chicken_converted)
converted <- lapply(converted, FUN = unique)
# creating jaccard similarity matrix

jaccard <- function(a, b) {
  intersection = length(intersect(a, b))
  union = length(a) + length(b) - intersection
  return (intersection/union)
}

jaccard_matrix <- matrix(ncol = length(converted), nrow = length(converted))
i = 1
while (i <= length(converted)){
  jaccard_matrix[,i] = as.numeric(lapply(converted, function(x){ jaccard(converted[[i]], x) }))
  i = i + 1
}

rownames(jaccard_matrix) <- names(converted)
colnames(jaccard_matrix) <- names(converted)

library(corrplot)

plot_color <- colorRampPalette(c('#F7F7F7', '#F4A582', '#D6604D', '#D6604D','#D6604D','#D6604D','#D6604D','#D6604D','#D6604D','#D6604D','#D6604D','#D6604D','#D6604D','#D6604D','#B2182B', '#67001F'))

corrplot(jaccard_matrix, method = 'color',
         col = plot_color(200), col.lim = c(0,1), is.corr = F, order = 'hclust')