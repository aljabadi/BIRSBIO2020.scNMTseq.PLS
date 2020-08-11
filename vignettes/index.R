params <-
list(on_my_mac = c(user = TRUE), save_output = FALSE, local_data = FALSE, 
    mini_run = FALSE, matching_rna_for_umap = FALSE, drop_lineages = "c('Primitive_endoderm','Visceral_endoderm', 'ExE_ectoderm')", 
    umap_params = c(run.seed = 42, n_neighbors = 15, n_components = 2, 
    min_dist = 0.55))

## ---- include = FALSE---------------------------------------------------------
library(knitr)
knitr::opts_chunk$set(collapse = TRUE, 
                      comment = "#>", 
                      fig.width = 8, 
                      cache = FALSE)

## ---- eval=params$on_my_mac, include=FALSE------------------------------------
sapply(list.files('../R', full.names = TRUE), source)
local.lib <- '../lib'
dir.create(local.lib)
.libPaths(local.lib)
## uncomment this if you are not building the vignette along with the package
# remotes::install_github('mixOmicsTeam/mixOmics@MultiAssayExperiment', upgrade = 'never')

## ----load packages, warning=FALSE, message=FALSE------------------------------
library(BIRSBIO2020.scNMTseq.PLS)
library(MultiAssayExperiment)
library(scater)
library(scran)
library(mixOmics)
library(ggplot2)
library(magrittr)
library(reshape2)

## ---- eval=params$save_output, include=params$save_output---------------------
#  ## create directories to save the figures and output data:
#  if (!dir.exists('figures')) {
#    cat('creating "figures" folder ...\n')
#    dir.create('figures')
#  }
#  
#  if (!dir.exists('savedata')) {
#     cat('creating "savedata" folder ...\n')
#    dir.create('savedata')
#  }

## ---- eval=!params$mini_run---------------------------------------------------
cat('loading data from Cloudstor ...\n')
gastru.mae_path <- url('https://cloudstor.aarnet.edu.au/plus/s/jsW7nh4YwThw8Q5/download')

## ---- eval=params$on_my_mac & !params$mini_run, include=FALSE-----------------
cat('loading data from local folder ...\n')
gastru.mae_path <- 'savedata/scnmtseq_gastrulation_mae-sce.rds'

## ---- eval=params$on_my_mac & params$mini_run, include=FALSE------------------
#  cat('loading mini data from local folder ...\n')
#  gastru.mae_path <- 'savedata/mini.gastru.mae.rds'

## ----load data----------------------------------------------------------------
gastru.mae <- readRDS(gastru.mae_path)

## -----------------------------------------------------------------------------
# Subset RNA expression and DNA methylation modalities
keep_assays <- grep("rna|met",names(assays(gastru.mae)))
# Remove putative extraembryonic cells
keep_cells <- !(gastru.mae$lineage %in% params$drop_lineages)
# Keep cells that pass QC for DNA methylation and RNA expression
keep_cells <- keep_cells & (gastru.mae$pass_metQC==TRUE & gastru.mae$pass_rnaQC==TRUE)
## keep full rna SCE for UMAP
rna.sce <- gastru.mae@ExperimentList$rna
## subset MAE
gastru.mae <- gastru.mae[,keep_cells, keep_assays]
## rna SCE for matching cells
rna.sce.matching <- gastru.mae@ExperimentList$rna

## -----------------------------------------------------------------------------
gastru.mae@ExperimentList
gastru.mae@ExperimentList$rna <- logcounts(gastru.mae@ExperimentList$rna)

## -----------------------------------------------------------------------------
get_pct_missing <- function(arr) {
  round(100*sum(is.na(arr))/prod(dim(arr)))
}

df <- lapply(experiments(gastru.mae), function(w)
     {
     c(
       'N (# cells)' = dim(w)[2],
       'P (# features)' = dim(w)[1],
       '% data missing' = get_pct_missing(w)
     )
     })
df <- data.frame(df)
kable(t(df), align = 'l')

## -----------------------------------------------------------------------------
table(gastru.mae$stage) %>% as.data.frame() %>% t() %>% set_rownames(c('stage', '# of cells')) %>% kable()
## create a data.frame from cell metadata
coldata <- data.frame(colData(gastru.mae))

## ---- message=FALSE-----------------------------------------------------------
# get the methylation assays
met_assays <- grep(names(gastru.mae), pattern = '^met', value = TRUE)
# add dimensions to labels for ggplot
dims <- lapply(experiments(gastru.mae[,,met_assays]), dim)
dims <- sapply(dims, function(x) sprintf(' (%s, %s)', x[2], x[1]))
names(met_assays) <- paste0(met_assays, dims) %>% as.list()
# calculate the feature detection in a data.frame for methylation assays
coverages <- lapply(met_assays, function(assay_name) {
  mat <- assay(gastru.mae, assay_name)
  NAs <- rowSums(!is.na(mat))/dim(mat)[2]*100
  data.frame(pct_NAs=NAs)
})
# create a long data.frame containing the assay name for plot
coverages <- rbindListWithNames(coverages)
coverages$dataset <- factor(coverages$dataset, levels = unique(coverages$dataset), ordered = TRUE)

## ---- fig.width=8, fig.asp=0.4------------------------------------------------
cov_plot <- ggplot(coverages, aes(x = pct_NAs)) + geom_density(fill = 'lightblue', show.legend = FALSE) +
  geom_vline(aes(xintercept=mean(pct_NAs)),
             color="blue", linetype="dashed", size=0.5) +
  labs(x = '% of cells detecting the feature') + facet_wrap(.~dataset, nrow = 2) +
  theme_bw() + theme(strip.text.x = element_text(size = 10, face = 'bold', color = 'purple'))
  
cov_plot

## ---- eval=params$save_output, include=params$save_output---------------------
#  ggsave(cov_plot, filename = 'figures/covplots.pdf', width = 8, height = 4)

## -----------------------------------------------------------------------------
# Create colour palettes for stages:
# col_pallete <-dput( viridisLite::viridis(n = 4) )
col_pallete <- c("#440154FF", "#31688EFF", "#35B779FF", "#FDE725FF")
names(col_pallete) <- c("E4.5", "E5.5", "E6.5", "E7.5")

## helper function to create UMAP plots coloured by stage
plot_umap_by_stage <- function(df, dims = c(1,2)) {
  axes <- paste0('UMAP_', dims)
  ggplot(df, aes_string(axes[1], axes[2])) + geom_point(aes(col=stage)) +
    theme_classic()+ scale_color_manual(values = col_pallete) + labs(col = 'Stage')
}

## -----------------------------------------------------------------------------
decomp <- modelGeneVar(rna.sce)
## filter by mean expression and significance of biological variation signal
hvgs <- rownames(decomp)[decomp$p.value<0.01 & decomp$mean > 0.01]
length(hvgs)

## -----------------------------------------------------------------------------
## PCA first: retrieve 25 PCs
rna.sce <- runPCA(rna.sce,  ncomponents = 25, subset_row=hvgs)

RNGversion('4.0.0')
set.seed(params$umap_params['run.seed'])
rna.sce <- runUMAP(rna.sce, dimred="PCA", 
                   ncomponents = params$umap_params['n_components'], 
                   n_neighbors = params$umap_params['n_neighbors'], 
                   min_dist = params$umap_params['min_dist'])

## -----------------------------------------------------------------------------
# data.frame of embeddings with cell metadata
df <-  data.frame(reducedDim(rna.sce,"UMAP"))
colnames(df) <- paste0('UMAP_', seq_along(df))
df <- cbind(df, colData(rna.sce))

## ---- fig.asp=0.7-------------------------------------------------------------
plot_umap_by_stage(df = df, dims = c(1,2))

## -----------------------------------------------------------------------------
decomp <- modelGeneVar(rna.sce.matching)
## filter by mean expression and significance of biological variation signal
hvgs <- rownames(decomp)[decomp$p.value<0.01 & decomp$mean > 0.01]
length(hvgs)

## -----------------------------------------------------------------------------
## PCA first: retrieve 25 PCs
rna.sce.matching <- runPCA(rna.sce.matching,  ncomponents = 25, subset_row=hvgs)

RNGversion('4.0.0')
set.seed(params$umap_params['run.seed'])
rna.sce.matching <- runUMAP(rna.sce.matching, dimred="PCA", 
                   ncomponents = params$umap_params['n_components'], 
                   n_neighbors = params$umap_params['n_neighbors'], 
                   min_dist = params$umap_params['min_dist'])

## -----------------------------------------------------------------------------
# data.frame of embeddings with cell metadata
df <-  data.frame(reducedDim(rna.sce.matching,"UMAP"))
colnames(df) <- paste0('UMAP_', seq_along(df))
df <- cbind(df, colData(rna.sce.matching))

## ---- fig.asp=0.7-------------------------------------------------------------
plot_umap_by_stage(df = df, dims = c(1,2))

## -----------------------------------------------------------------------------
# data dimensions
lapply(experiments(gastru.mae), dim)

## -----------------------------------------------------------------------------
## number of components
ncomp <- 2
## feature scaling
scale <- TRUE

## ----run pls------------------------------------------------------------------
cat(sprintf('Running PLS with %s components performing variable selection\n', ncomp))
st <- system.time({
  mmspls <-
    multimodal_sPLS_wrapper(mae = gastru.mae, study_assays = NULL,
                            ncomp = ncomp, scale = TRUE, design = 'null', lineages = NULL,
                            stages = NULL, DA = NULL, keepX = NULL, save = FALSE)
})['elapsed']
st <- round(st/60, 1)

cat('\nPLS run finished. Runtime: ', st, 'min\n')

## ---- eval=params$save_output, echo=params$save_output------------------------
#  saveRDS(mmspls, file = 'savedata/MultiModalSparsePLS-All.rds')
#  # mmspls <- readRDS('savedata/MultiModalSparsePLS-All.rds')

## -----------------------------------------------------------------------------
plotIndiv_wrapper <- function(pls_obj, coldata, col_cell = 'stage', legend.title = 'Stage', comp = c(1,2), ...) {
  ## order coldata cells based on matching cells in pls modalities
  coldata <- coldata[rownames(pls_obj$X[[1]]),]
  cell_class <- coldata[,col_cell]  ## col factor for cells
  cell_class <- factor(cell_class)
  suppressWarnings({
    plotIndiv(pls_obj, pch = 16, comp = comp,group = cell_class, legend = TRUE, legend.title = legend.title, size.subtitle = 10, cex=0.8, ...)
  })
}

## ---- fig.width=8, fig.height=8-----------------------------------------------
plotIndiv_wrapper(pls_obj = mmspls, coldata = coldata, col_cell = 'stage', legend.title = 'Stage', comp = c(1,2), subtitle = names(mmspls$X), col.per.group = col_pallete)

## ---- eval=params$save_output, echo=FALSE-------------------------------------
#  ggsave(filename = 'figures/plotIndiv-MultiModalSparsePLS-stage.pdf', width = 8, height = 6)

## -----------------------------------------------------------------------------
concordance <- function(pls_obj, comp = 1) {
  
  ## correlation of components with RNA components
  cor_with_y <- sapply(pls_obj$variates[-1], function(x_variates) {
    y_variates <- pls_obj$variates[[1]]
    y_variates <- y_variates[, comp, drop = FALSE]
    cor <- cor(y_variates, x_variates)
    cor <- diag(cor)
    names(cor) <- paste0('Component ', comp)
    cor['mean'] <- mean(cor)
    cor <- round(cor, 2)
    cor
  })
  ## make a long data.frame
  cor_with_y <- cor_with_y %>% as.data.frame() %>% t() %>% 
    reshape2::melt(id.vars = 1:3) %>% 
    set_colnames(c('Modality', 'Component', 'Concordance'))
  ## print concordance measures for each componet and the average:
  # sapply(unique(as.character(cor_with_y$Component)), function(z){
  #   cat('#> Component:', z )
  #   dplyr::filter(cor_with_y, Component == z)[,c('Modality', 'Concordance')] %>% 
  #     kable() 
  # })

  ggplot(cor_with_y, aes(x = Modality, y = Concordance)) + geom_col(aes(fill = Component), position = 'dodge2') + theme_classic() +
    labs(y = sprintf('Correlation with RNA (comp %s) ', paste(comp, collapse = ', ')), x = '', fill = '') + ylim(c(0, 1)) + theme(axis.text.x = element_text(angle = 45, vjust = 0.8, hjust=1))
}

## -----------------------------------------------------------------------------
concordance(pls_obj = mmspls, comp = c(1, 2))

