## helper function to create UMAP plots coloured by stage
#' @export
plot_reducedDim <- function(sce, reducedDim = 'UMAP',
                            dims = c(1,2), col_palette,
                            col = 'stage') {
  comps.prefix <- ifelse(reducedDim == 'UMAP', 'UMAP_', 'PC_')
  df <- data.frame(reducedDim(sce, reducedDim))[,dims]
  axes <- colnames(df) <- paste0(comps.prefix, dims)
  df[,col] <- colData(sce)[,col]
  p <- ggplot(df, aes_string(axes[1], axes[2])) + geom_point(aes_string(col=col), alpha = 0.75) +
    theme_classic()+
    labs(col = c(stage = 'Stage', lineage = 'Lineage')[col]) +
    theme(
      # legend.position=c(1,1),
      # legend.justification=c(1,1),
      legend.box.background = element_rect(colour = "black"),
      legend.text=element_text(size=7, face='bold'),
      legend.title=element_text(size=8, face='bold'),
      strip.text.x = element_text(size = 10, face = 'bold'))

  if (!is.null(col_palette)) {
    p <- p + scale_color_manual(values = col_palette)
  }
  p
}
