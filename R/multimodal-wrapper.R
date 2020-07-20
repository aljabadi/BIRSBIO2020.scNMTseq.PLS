#' MuliModalSparsePLS(DA)
#'
#' @param mae A MultiAssayExperiment object
#' @param study_assays study assays. first must be 'rna' for unsupervised analysis
#' @param lineages Character, lineages to keep, from unique(colData(mae)$lineage10x_2)
#' @param stages Character, stages to keep, from unique(colData(mae)$stage)
#' @param ncomp Integer, number of components
#' @param scale Logical, whether to scale the features
#' @param DA NULL for unsupervised, or from names(colData(mae)) for supervised (PLSDA)
#'
multimodal_analysis_wrapper <- function(mae, study_assays, lineages=NULL, stages = NULL, ncomp = 2, scale = TRUE, design = 'null', DA=NULL, keepX = NULL, save = TRUE) {
  #
  # if (is.null(study_assays)) {
  #   study_assays <- names(experiments(mae))
  #   study_assays <- study_assays[!grepl(pattern = '^wt_*', x = study_assays)]
  # }
  # if (is.integer(study_assays)) {
  #   study_assays <- names(experiments(mae))[study_assays]
  # }
  # mae <- mae[,,study_assays]
  # mae <- mixOmics::make_unique_feature_names(mae)
  # if (!is.null(lineages)) {
  #   mae <- mae[,mae$lineage10x_2 %in% lineages,]
  #   lin <- 'AllLineages'
  # } else{
  #   lineages <- unique(mae$lineage10x_2 )
  #   lin <- paste0(lineages, collapse = '-')
  # }
  #
  # if (is.null(stages)) {
  #   stages <- unique(mae$stage)
  # }
  #
  # stag <- paste(stages, collapse = '-')
  # mae <- mae[,mae$stage %in% stages,]
  #
  # if (is.null(keepX)) {
  #   keepX <- lapply(named_list(study_assays), function(z){
  #     p <- dim(assay(mae, z))[1]
  #     keepx <- ifelse(p > 500, 50, max(25, ceiling(p/10)))
  #     rep(keepx, ncomp)
  #   })
  # }
  # if (!is.null(DA)) {
  #   pls_mode = 'PLSDA'
  #   formu <- as.formula(sprintf('%s ~ %s', DA, paste(study_assays, collapse = ' + ')))
  #   out <- MultiModalSparsePLSDA(data = mae, formula = formu, keep_features = keepX, ncomp = ncomp, design = design, scale = scale)
  # } else {
  #   if (study_assays[1] != 'rna') {
  #     stop('First assay should be "rna"')
  #   }
  #   pls_mode = 'PLS'
  #   formu <- as.formula(sprintf('rna ~ %s', paste(study_assays[-1], collapse = ' + ')))
  #   out <- MultiModalSparsePLS(data = mae, formula = formu, keep_features = keepX, ncomp = ncomp, design = design, scale = scale)
  # }
  # NAME <- sprintf('savedata/MultiModalSparse%s_%s__%s__%s__%s.rds', pls_mode, ifelse(is.null(DA), 'rna', DA), paste0(study_assays[-1], collapse = '-'), lin, stag)
  #
  #
  # if (isTRUE(save)) {
  #   saveRDS(out, file = NAME)
  # }
  # out
}
