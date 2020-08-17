#' @export
plot_pls <- function(pls_obj, stage, legend.title = 'Stage', comp = c(1,2), col_palette) {
  variates <- lapply(pls_obj$variates, function(arr){
    df <- as.data.frame(arr[,comp])
    df$stage <- stage
    df
  })
  variates <- rbindListWithNames(variates, new_col = 'Modality')
  axes <- colnames(variates)[1:2]
  axes_labels <- paste0('component ', comp)
  ggplot(variates) +
    theme_classic() +
    geom_point(aes_string(x = axes[1], y = axes[2], col = 'stage'), alpha = 0.75) +
    facet_wrap(.~Modality, ncol = 3, scales = "free") +
    scale_color_manual(values = col_palette) +
    labs(x = axes_labels[1], y = axes_labels[2], col = legend.title) +
    theme(legend.position=c(0.9,0.15), legend.justification=c(0.6,.15),
          legend.text=element_text(size=9),
          legend.box.background = element_rect(colour = "black"),
          legend.title=element_text(size=9),
          strip.text.x = element_text(size = 10, face = 'bold'))

}
