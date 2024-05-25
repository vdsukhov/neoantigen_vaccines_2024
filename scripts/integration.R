# Code written by Maria Firulyova

suppressMessages(library(data.table))
suppressMessages(library(dplyr))
suppressMessages(library(functools))
suppressMessages(library(ggplot2))
suppressMessages(library(gridExtra))
suppressMessages(library(Matrix))
suppressMessages(library(magrittr))
suppressMessages(library(sctransform))
suppressMessages(library(Seurat))
suppressMessages(library(argparse))
suppressMessages(library(glmGamPoi))
suppressMessages(library(SeuratDisk))
suppressMessages(library(SeuratWrappers))
suppressMessages(library(miQC))
suppressMessages(library(flexmix))
suppressMessages(library(Matrix))

set.seed(1)


parser <-
  ArgumentParser(description = 'Get scRNA-seq related figures from the paper')
parser$add_argument('--data',
                    type = "character",
                    help = 'Path to directory with 10x output subdirs')
parser$add_argument('--tsne',
                    type = "character",
                    help = 'Path to binary fast-tsne')
parser$add_argument('--out_dir',
                    type = "character",
                    help = 'Path to output directory')
## SET VARIABLES

args <- parser$parse_args()

## FUNCTIONS


add_metadata <- function(data) {
  mito.genes <-
    grep(pattern = "^Mt\\.|^MT\\.|^mt\\.|^Mt-|^MT-|^mt-",
         x = rownames(x = GetAssayData(object = data)),
         value = TRUE)
  percent.mito <-
    Matrix::colSums(GetAssayData(object = data, slot = "counts")[mito.genes, ]) /
    Matrix::colSums(GetAssayData(object = data, slot = "counts"))
  data[['percent.mito']] <- percent.mito
  data[['percent.mt']] <- percent.mito
  data[['percent.mito_log10']] <- log10(data[['percent.mito']] + 1)
  data[['nCount_RNA_log10']] <- log10(data[['nCount_RNA']] + 1)
  data[['nFeature_RNA_log10']] <- log10(data[['nFeature_RNA']] + 1)
  data[['nCount_RNA_log2']] <- log2(data[['nCount_RNA']] + 1)
  data[['nFeature_RNA_log2']] <- log2(data[['nFeature_RNA']] + 1)
  data[['scaled_mito']] <- scale(percent.mito)
  data[['scaled_nCount_RNA']] <- scale(data[['nCount_RNA_log10']])
  data[['nCount_ADT_log10']] <- log10(data[['nCount_ADT']] + 1)
  data[['nFeature_ADT_log10']] <- log10(data[['nFeature_ADT']] + 1)
  data[['nCount_ADT_log2']] <- log2(data[['nCount_ADT']] + 1)
  data[['nFeature_ADT_log2']] <- log2(data[['nFeature_ADT']] + 1)
  data[['scaled_nCount_ADT']] <- scale(data[['nCount_ADT_log10']])
  attr(data$scaled_nCount_ADT, "scaled:center") <- NULL
  attr(data$scaled_nCount_ADT, "scaled:scale") <- NULL
  attr(data$scaled_nCount_RNA, "scaled:center") <- NULL
  attr(data$scaled_nCount_RNA, "scaled:scale") <- NULL
  attr(data$scaled_mito, "scaled:center") <- NULL
  attr(data$scaled_mito, "scaled:scale") <- NULL
  data
}

get_conf_interval <- function(dataset, parameter) {
  left <- mean(dataset[[parameter]][[1]]) - qnorm(0.975)
  right <- mean(dataset[[parameter]][[1]]) + qnorm(0.975)
  return(c(left, right))
}

filter_mito <- function(dataset){
  expr <- FetchData(object = dataset, vars = 'scaled_mito')
  dataset <- dataset[, which(x = expr < get_conf_interval(dataset, 'scaled_mito')[2])]
  dataset
}


## GATHERING DATA TOGETHER
    

options(future.globals.maxSize = 10000 * 1024^2)

get_df <- function(path) {
  data <- Read10X(path)
  rownames(x = data[["Antibody Capture"]]) <- gsub(pattern = "_.*", replacement = "",
                                                   x = rownames(x = data[["Antibody Capture"]]))
  obj <- CreateSeuratObject(counts = data[["Gene Expression"]], min.cells = 3, min.features = 200)
  obj[["ADT"]] <- CreateAssayObject(data[["Antibody Capture"]][, colnames(x = obj)])
  #obj <- subset(obj, cells = sample(colnames(obj), 5e2)) ## TODO: TEST OPTION
  obj$sample <- gsub('.*/', '', path)
  obj$dosage <- gsub('.*-', '', path)
  obj <- RenameCells(object = obj, add.cell.id = unique(obj$sample))
  obj <- add_metadata(obj)
  obj
}

get_whole_obj <- function(pathes) {
  objects <- lapply(pathes, get_df)
  names(objects) <- sapply(objects, function(x) unique(x$sample))
  objects
}

pathes <- list.files(args$data, full.names = T)

whole <- get_whole_obj(pathes)

setwd(args$out_dir)

## Number of cells before

cells.before <- sapply(whole, function(x) dim(GetAssayData(object = x, slot = "counts"))[2])

## MIQC

whole <- sapply(whole, function(x) RunMiQC(x, nFeature_RNA = "nFeature_RNA",
                                           posterior.cutoff = 0.75,
                                           model.slot = "flexmix_model"))

## FILTER MT CONTENT

whole <- sapply(whole, function(x) filter_mito(x))


## NORMALIZATION

whole <- sapply(whole, function(x) SCTransform(
  x,
  ncells=min(100000, ncol(x)),
  vars.to.regress = c("percent.mito"),
  verbose = T,
  conserve.memory = T
))


whole <- sapply(whole, function(x) {
  x <- NormalizeData(x, normalization.method = 'CLR', margin = 2, assay = "ADT")
  x <- ScaleData(x, assay = 'ADT')
  x <- RunPCA(x, features = rownames(x[["ADT"]]), reduction.name = 'apca', assay = 'ADT', approx = F)
  x
})

## CELL CYCLE SCORING

s.genes <- stringr::str_to_title(cc.genes$s.genes)
g2m.genes <- stringr::str_to_title(cc.genes$g2m.genes)


whole <- sapply(whole, function(x) CellCycleScoring(x, s.features = s.genes,
                                                    g2m.features = g2m.genes, set.ident = F,
                                                    assay = 'SCT'))
## INTEGRATION: RNA


whole.features <- SelectIntegrationFeatures(object.list = whole, nfeatures = 2000)

whole <- lapply(X = whole, FUN = function(x) {
  x <- RunPCA(x, features = whole.features)
})

whole <- PrepSCTIntegration(object.list = whole, anchor.features = whole.features,
                            verbose = FALSE)
whole.anchors <- FindIntegrationAnchors(object.list = whole, normalization.method = "SCT",
                                        anchor.features = whole.features, verbose = FALSE, reduction = 'rpca')
whole.integrated <- IntegrateData(anchorset = whole.anchors, normalization.method = "SCT",
                                  verbose = FALSE)

## PCA
gc()

whole.integrated <- RunPCA(whole.integrated, verbose = FALSE)

scale_data <- cbind(cbind(whole$`MGI2258_M_KC-1-P-plus-PBS`@assays$ADT@scale.data,
      whole$`MGI2258_M_KC-2-H-plus-HD`@assays$ADT@scale.data), whole$`MGI2258_M_KC-3-L-plus-LD`@assays$ADT@scale.data)

whole.integrated@assays$ADT@scale.data <- scale_data

whole.integrated <- RunPCA(whole.integrated, features = rownames(whole.integrated[["ADT"]]),
                           reduction.name = 'apca', assay = 'ADT', approx = F)

## FindMultiModalNeighbors

whole.integrated <- FindMultiModalNeighbors(
  whole.integrated, reduction.list = list("pca", "apca"), 
  dims.list = list(1:30, 1:2), modality.weight.name = "RNA.weight"
)

## TSNE

whole.integrated <- RunTSNE(whole.integrated, dims = 1:20, tsne.method = "FIt-SNE",
                            fast_tsne_path = args$tsne, nthreads = 4, max_iter = 2000)
## UMAP

whole.integrated <- RunUMAP(whole.integrated, dims = 1:20)

whole.integrated <- RunUMAP(whole.integrated, nn.name = "weighted.nn", reduction.name = "wnn.umap", reduction.key = "wnnUMAP_")

## CLUSTERING

whole.integrated <- FindNeighbors(object = whole.integrated, dims = 1:20)
whole.integrated <- FindClusters(whole.integrated, graph.name = "wsnn", resolution = seq(0.2, 0.6, 0.2),
                                 algorithm = 3, verbose = FALSE)
whole.integrated <- FindClusters(object = whole.integrated, resolution = c(0.2, 0.4, 0.6, 0.8, 1))

## SAVING: DATASET


save(list = c('whole.integrated', 'whole.features', 'whole.anchors'), file = "object.RData")

## FINDING ANS SAVING MARKERS

analyze_object <- function(object, ident) {
  Idents(object) <- object[[ident]]
  if (length(levels(object)) == 1) {
    return(message(sprintf('%s: since only one cluster was identified, markers can not be found', ident)))
  }
  out_dir <- paste0('markers/', ident)
  dir.create(out_dir, recursive = T)
  whole.markers <- FindAllMarkers(object = object,
                                  assay='SCT',
                                  only.pos = TRUE,
                                  min.pct = 0.10,
                                  test.use = 'wilcox',
                                  max.cells.per.ident = 3e3,
                                  random.seed = 42)
  write.table(whole.markers, paste(out_dir, "markers.tsv", sep = '/'), sep="\t", quote=F, row.names=F)
}

idents <- grep('.*snn.*', colnames(whole.integrated@meta.data), value = T)
sapply(idents, function(ident) analyze_object(object = whole.integrated, ident = ident))


## Number of cells after

cells.after <- sapply(whole, function(x) length(colnames(x = x)))
cells.diff <- cells.before-cells.after
rbind(cells.before, cells.after, cells.diff)
