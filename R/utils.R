## ----------- all_identical -----------
## check if elements in a list are identical
#' @export
all_identical <- function(lst) {
  for (j in seq_along(lst[-1])) {
    if (!identical(lst[[j]], lst[[j+1]]))
      stop(sprintf("not identical elements: %s and %s",j , j+1 ), call. = FALSE)
  }
  TRUE
}

## ----------- rbindListWithNames -----------
## base::rbind a named list of data.frames adding a new column
## indicating the name of the dataset in the list
#' rbind a list of arrays
#'
#' @param lst A named list of arrays
#' @param new_col Character, the name of the new column which adds the list names
#' to the resulting array
#'
#' @name utils
#' @export
rbindListWithNames <- function(lst, new_col = "dataset") {
  lst_with_newcol <- mapply(x=names(lst), y=lst, FUN = function(x, y){
    y[,new_col] <- x
    y
  }, SIMPLIFY = FALSE)
  Reduce(rbind, lst_with_newcol)
}

## ----------- named_list -----------
## create a named list from a character vector
## that can be used with apply family and
## return a named list
#' @export
named_list <- function(char) {
  out <- as.list(char)
  names(out) <- char
  out
}

## ----------- ggplot color hue -----------
#' @export
gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  grDevices::hcl(h = hues, l = 65, c = 100)[1:n]
}

## ----------- get_pct_missing -----------
## for array 'arr', output the percentage of NAs (0-100)
#' @export
get_pct_missing <- function(arr) {
  100*sum(is.na(arr))/prod(dim(arr))
}

