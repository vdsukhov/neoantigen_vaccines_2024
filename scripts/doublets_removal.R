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

set.seed(1)


parser <-
  ArgumentParser(description = 'Get scRNA-seq related figures from the paper')
parser$add_argument('--data',
                    type = "character",
                    help = 'Path to seurat rda')
parser$add_argument('--tsne',
                    type = "character",
                    help = 'Path to binary fast-tsne')
parser$add_argument('--out_dir',
                    type = "character",
                    help = 'Path to output directory')
parser$add_argument('--clusters',
                    type = "integer", nargs='+',
                    help = 'target clusters')
parser$add_argument('--ident',
                    type = "character",
                    help = 'resolution name')
parser$add_argument('--invert_clustering',
                    type = "logical")
parser$add_argument('--adt',
                    type = "character",
                    help = 'ADT name for cells filtering')
## SET VARIABLES

args <- parser$parse_args()

if (is.null(args$invert_clustering)) {
  args$invert_clustering <- F
}

if (is.null(args$adt)) {
  args$adt <- F
}

if (is.null(args$clusters)) {
  args$clusters <- F
  args$ident <- F
}

print(args)

## GATHERING DATA TOGETHER


options(future.globals.maxSize = 15000 * 1024^2)

get_data <- function(data, res, clusters, invert_clustering, adt) {
  load(data)
  whole.integrated@meta.data <- whole.integrated@meta.data %>% 
    mutate(`ADT_TCR-beta` = whole.integrated@assays$ADT@data['TCR-beta',]) %>% 
    mutate(`ADT_CD90-2` = whole.integrated@assays$ADT@data['CD90-2',]) %>% 
    mutate(`scaled_ADT_CD90-2` = scale(`ADT_CD90-2`)) %>% 
    mutate(`scaled_ADT_TCR-beta` = scale(`ADT_TCR-beta`))
  whole.integrated <- subset(whole.integrated,
                             cells = rownames(whole.integrated@meta.data %>% filter(nCount_RNA_log2 > 10)))
  if (!is.null(res)) {
    Idents(whole.integrated) <- res
    whole.integrated <- subset(whole.integrated, idents = clusters, invert = invert_clustering)
  }
  if (!is.null(adt)) {
    adt <- sprintf('scaled_ADT_%s', adt)
    target_cells <- rownames(whole.integrated@meta.data %>%
                               dplyr::filter(!!as.symbol(adt) > 0))
    whole.integrated <- subset(whole.integrated, cells = target_cells)
  }
  DefaultAssay(whole.integrated) <- 'RNA'
  whole.integrated[['integrated']] <- NULL
  whole.integrated[['SCT']] <- NULL
  SplitObject(whole.integrated, split.by = 'sample')
}

whole <- get_data(args$data, args$ident, args$clusters, args$invert_clustering, args$adt)
setwd(args$out_dir)

## Number of cells before

cells.before <- sapply(whole, function(x) dim(GetAssayData(object = x, slot = "counts"))[2])

## MIQC

whole <- sapply(whole, function(x) RunMiQC(x, nFeature_RNA = "nFeature_RNA",
                                           posterior.cutoff = 0.75,
                                           model.slot = "flexmix_model"))

## NORMALIZATION

whole <- sapply(whole, function(x) SCTransform(
  x,
  ncells=min(100000, ncol(x)),
  vars.to.regress = c("percent.mito"),
  method = "glmGamPoi",
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

## TSNE

whole.integrated <- RunTSNE(whole.integrated, dims = 1:20, tsne.method = "FIt-SNE",
                            fast_tsne_path = args$tsne, nthreads = 4, max_iter = 2000)
## UMAP

whole.integrated <- RunUMAP(whole.integrated, dims = 1:20)


## CLUSTERING

whole.integrated <- FindNeighbors(object = whole.integrated, dims = 1:20)
whole.integrated <- FindClusters(object = whole.integrated, resolution = c(0.2, 0.4, 0.6, 0.8, 1))

## SAVING: DATASET


save(list = c('whole.integrated', 'whole.features', 'whole.anchors'), file = "object.RData")

## AVERAGING

cluster.averages <- AverageExpression(object = whole.integrated, assays = 'SCT', slot = 'data')
sapply(names(cluster.averages), 
       function(x) write.table(cluster.averages[[x]], file=paste0(x, "_clusters.tsv")))

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

idents <- paste0('integrated_snn_res.', c(0.2, 0.4, 0.6, 0.8, 1))
sapply(idents, function(ident) analyze_object(object = whole.integrated, ident = ident))


## Number of cells after

cells.after <- sapply(whole, function(x) length(colnames(x = x)))
cells.diff <- cells.before-cells.after
rbind(cells.before, cells.after, cells.diff)
