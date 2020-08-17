#' @export
kmeans_accuracy <- function(variates, labels, nstart=50, iter.max=2000) {
  names(labels) <- rownames(variates)
  labels <- as.data.frame(labels)
  labels$cell <- rownames(labels)
  kmeans_res <- kmeans(x = variates, centers = length(unique(labels$labels)), iter.max=iter.max, nstart = nstart)
  labels$cluster <- kmeans_res$cluster
  tbl <- table(labels$labels, labels$cluster)
  stage_cluster <- rep(NA, length(unique(labels$cluster)))
  names(stage_cluster) <- rownames(tbl)
  for (stage in rownames(tbl)) {
    stage_cluster[stage] <- which.max(tbl[stage,])
  }
  ## order columns
  tbl <- tbl[,stage_cluster[rownames(tbl)]]
  total_cells <- rowSums(tbl)
  class_weight <- total_cells / sum(total_cells)
  misclassified <- rowSums(tbl) - diag(tbl)
  err_class <- misclassified / total_cells
  balanced_err <- mean(err_class)
  1 - balanced_err
}

# kmeans_perf <- function(variates, labels, times=10, ...) {
#   balanced_accuracies <- sapply(seq_len(times), function(f){
#     kmeans_accuracy(variates, labels, ...)
#   })
# }
