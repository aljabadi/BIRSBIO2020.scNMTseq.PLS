## -------- helper functions to create nipals plots coloured by stage
plot_nipals_redcuedDims <- function(nipals_comps_df, red_dim = 'PC', col_palette, dims=c(1,2),
                                    facet_ncol = 2) {
  axes <- paste0(red_dim, '_', dims)
  ggplot(nipals_comps_df, aes_string(axes[1], axes[2])) +
    geom_point(aes(col = stage)) +
    scale_color_manual(values = col_palette) +
    facet_wrap(.~Modality, scales = 'free', ncol = facet_ncol) +
    theme_classic() +
    theme(legend.box.background = element_rect(colour = "black"),
          legend.title = element_blank())
}
## ------------------------------------------------------------------------ ##

get_nipals_scores <- function(nipals_comps, cell_order, stage) {
  ## get scores
  nipals_comps_df <- lapply(named_list(names(nipals_comps)), function(x){
    y <- nipals_comps[[x]]
    df <- data.frame(y$scores)
    ## ensure cells order
    df <- df[cell_order,]
    colnames(df) <- paste0('PC_', seq_along(df))
    df$stage <- stage
    df
  })
  nipals_comps_df <- rbindListWithNames(nipals_comps_df, new_col = 'Modality')
  nipals_comps_df
}
## ------------------------------------------------------------------------ ##

#' @export
plot_nipals <-
  function(nipals_comps,
           cell_order,
           stage,
           red_dim = 'PC',
           col_palette,
           dims = c(1, 2),
           facet_ncol = 2) {

    nipals_comps_df <- get_nipals_scores(nipals_comps, cell_order, stage)
    plot_nipals_redcuedDims(
      nipals_comps_df,
      red_dim = red_dim,
      col_palette = col_palette,
      dims = dims,
      facet_ncol = facet_ncol
    )
  }
