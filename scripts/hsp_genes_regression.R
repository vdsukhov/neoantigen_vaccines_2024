# Code written by Maria Firulyova

suppressMessages(library(data.table))
suppressMessages(library(dplyr))
suppressMessages(library(ggplot2))
suppressMessages(library(Matrix))
suppressMessages(library(magrittr))
suppressMessages(library(sctransform))
suppressMessages(library(Seurat))
suppressMessages(library(argparse))
suppressMessages(library(glmGamPoi))
suppressMessages(library(SeuratWrappers))
suppressMessages(library(miQC))
suppressMessages(library(flexmix))
suppressMessages(library(SCNPrep))
suppressMessages(library(DoubletFinder))
suppressMessages(library(jsonlite))


set.seed(1)


parser <-
  ArgumentParser(description = 'Get scRNA-seq related figures from the paper')
parser$add_argument('--data',
                    type = "character",
                    help = 'Path to seurat rda')
parser$add_argument('--out_dir',
                    type = "character",
                    help = 'Path to output directory')
parser$add_argument('--json_file',
                    type = "character",
                    help = 'path to json with given filters for data')
parser$add_argument('--sample_id',
                    type = "character",
                    help = 'Sample ID')
parser$add_argument('--tsne',
                    type = "character",
                    help = 'path to fast tsne')
## SET VARIABLES

args <- parser$parse_args()

print(args)

## GATHERING DATA TOGETHER

options(future.globals.maxSize = 20000 * 1024^2)

get_data <- function(data, out_dir, js_f=NULL) {
  load(data)
  meta <- whole.integrated@meta.data
  if ("thresholds" %in% names(js_f)) {
    for (thr in names(js_f[['thresholds']])) {
      info <- strsplit(js_f[['thresholds']][thr][[1]], " ")[[1]]
      if (info[[1]] == '<') {
        meta <- meta %>%
          filter(!!as.symbol(thr) < as.numeric(info[2]))
      } else{
        meta <- meta %>% filter(!!as.symbol(thr) > as.numeric(info[2]))
      }
    }
  }
  
  if ("doublets" %in% names(js_f)) {
    for (thr in names(js_f[['doublets']])) {
      meta <-meta %>% filter(!!as.symbol(thr) != js_f[['doublets']][[thr]])
    }
  }
  
  if ("target_clusters" %in% names(js_f)) {
    meta <- meta %>% filter(!!as.symbol(names(js_f[['target_clusters']])) %in% strsplit(js_f[['target_clusters']][1][[1]], " ")[[1]])
  }
  if ("gene" %in% names(js_f)) {
    for (thr in names(js_f[['gene']])) {
      info <- js_f[['gene']][thr][[1]]
      cells_to_save <- colnames(whole.integrated)[log10(whole.integrated@assays$SCT@counts[thr,] + 1)
                                                  > js_f[['gene']][thr][[1]]]
      meta <- meta[cells_to_save,]
    }
  }
  whole.integrated <- subset(whole.integrated, cells = rownames(meta))
  if ("genes_to_exclude" %in% names(js_f)) {
    genes <- fread(js_f$genes_to_exclude)$gene
    if (sum(genes %in% rownames(whole.integrated@assays$RNA@counts)) > 0) {
      excluded_genes <- subset(whole.integrated, features = genes)
      target_genes <- rownames(whole.integrated@assays$RNA@counts)[!rownames(whole.integrated@assays$RNA@counts) %in% genes]
      whole.integrated <- subset(whole.integrated, features = target_genes)
      save(list=c('excluded_genes'), file = sprintf('%s/obj_with_excluded_genes.RData', out_dir))
    }
  }
  hsp.genes <-
    grep(pattern = "^Hsp",
         x = rownames(whole.integrated@assays$RNA@counts),
         value = TRUE)
  percent.hsp <- colSums(whole.integrated@assays$RNA@counts[hsp.genes,]) / colSums(whole.integrated@assays$RNA@counts)
  whole.integrated[['percent.hsp']] <- percent.hsp
  whole.integrated[['percent.hsp_log10']] <- log10(whole.integrated[['percent.hsp']] + 1)
  write.table(whole.integrated@meta.data, file = sprintf('%s/meta_of_remained_cells.csv', out_dir), quote = F)
  DefaultAssay(whole.integrated) <- 'RNA'
  whole.integrated[['integrated']] <- NULL
  whole.integrated[['SCT']] <- NULL
  SplitObject(whole.integrated, split.by = 'sample')
}

add_excluded_genes <- function(object, path) {
  obj_with_excluded_genes <- get(load(path))
  obj_with_excluded_genes <- subset(obj_with_excluded_genes, cells = colnames(object))
  object@assays$SCT@counts <- rbind(object@assays$SCT@counts, obj_with_excluded_genes@assays$SCT@counts)
  object
}

if (!is.null(args$json_file)) {
  js_f <- jsonlite::read_json(args$json_file, simplifyVector = T)
  
  whole <- get_data(data = args$data, js_f = js_f,
                    out_dir = args$out_dir)
} else{
  whole <- get_data(data = args$data, out_dir = args$out_dir)
}

setwd(args$out_dir)

## Number of cells before

cells.before <- sapply(whole, function(x)
  dim(GetAssayData(object = x, slot = "counts"))[2])

## MIQC

whole <- sapply(whole, function(x) RunMiQC(x, percent.mt = "percent.mito",
                                           nFeature_RNA = "nFeature_RNA",
                                           posterior.cutoff = 0.75,
                                           model.slot = "flexmix_model"))

## NORMALIZATION

whole <- sapply(whole, function(x) SCTransform(
  x,
  ncells=min(100000, ncol(x)),
  vars.to.regress = c("percent.mito", "percent.hsp"),
  method = "glmGamPoi",
  verbose = T,
  conserve.memory = T
))


## CELL CYCLE SCORING

s.genes <- stringr::str_to_title(cc.genes$s.genes)
g2m.genes <- stringr::str_to_title(cc.genes$g2m.genes)

whole <- sapply(whole, function(x) CellCycleScoring(x, s.features = s.genes,
                                                    g2m.features = g2m.genes,
                                                    set.ident = F,
                                                    assay = 'SCT'))

## INTEGRATION

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


## TSNE

# whole.integrated <- RunTSNE(whole.integrated, dims = 1:20, tsne.method = "FIt-SNE",
#                             nthreads = 4, max_iter = 2000)
## UMAP

whole.integrated <- RunUMAP(whole.integrated, dims = 1:20)


## CLUSTERING

whole.integrated <- FindNeighbors(object = whole.integrated, dims = 1:20)
whole.integrated <- FindClusters(object = whole.integrated, resolution = c(0.2, 0.4, 0.6, 0.8, 1))

## SAVING: DATASET

if (!is.null(args$json_file)) {
  if ('genes_to_exclude' %in% names(js_f) & !('obj_with_excluded_genes' %in% names(js_f))) {
    whole.integrated <- add_excluded_genes(whole.integrated, 'obj_with_excluded_genes.RData')
  }
  
  if (('genes_to_exclude' %in% names(js_f)) & ('obj_with_excluded_genes' %in% names(js_f))) {
    whole.integrated <- add_excluded_genes(whole.integrated, js_f$obj_with_excluded_genes)
  }
}


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

idents <- paste0('integrated_snn_res.', c(0.2, 0.4, 0.6, 0.8, 1))
sapply(idents, function(ident) analyze_object(object = whole.integrated, ident = ident))

## Number of cells after

cells.after <- sapply(whole, function(x) length(colnames(x = x)))
cells.diff <- cells.before-cells.after
rbind(cells.before, cells.after, cells.diff)

## Conversion

markers_files <- list.files(pattern = 'markers.tsv', recursive = T, full.names = T)
markers <- lapply(markers_files, function(x) fread(x) %>% dplyr::rename(avg_logFC = 'avg_log2FC'))
names(markers) <- gsub('/|markers', '', stringr::str_extract(markers_files, 'markers/*.*/'))

whole.integrated@meta.data <- whole.integrated@meta.data %>%
  mutate_if(is.character, as.factor) %>% 
  dplyr::select(!contains('wsnn'))

migrateSeuratObject(whole.integrated,
                    assay="SCT",
                    species='mm',
                    outdir = args$sample_id,
                    public = F,
                    curated = T,
                    markers = markers,
                    generateMarkers = F,
                    generateGMTS = F,
                    name=args$sample_id,
                    token=args$sample_id)
