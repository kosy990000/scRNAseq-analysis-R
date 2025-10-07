

# input x (matrix, nmad, na_to_false, verbose=출력 시 나올 것)
# return
#   : res$outlier -> boolean vector
#   : res$stats -> the list of varname, med, lower, upper
is_outlier <- function(x, nmads = 5, na_to_false = TRUE) {
  #들어오는 변수 이름 column$varname
  varname <- sub(".*\\$", "", deparse(substitute(x)))
  med <- median(x, na.rm = TRUE)
  dev <- mad(x, center = med, constant = 1, na.rm = TRUE)
  
  if (!is.finite(med) || !is.finite(dev) || dev == 0) {
    res <- rep(FALSE, length(x))
    stats <- list(
      varname = varname,
      median  = med,
      mad     = dev,
      lower   = NA,
      upper   = NA
    )
    if (verbose) print(stats)
    return(list(outlier = res, stats = stats))
  }
  
  lower <- med - nmads * dev
  upper <- med + nmads * dev
  res <- (x < lower) | (x > upper)
  
  if (na_to_false) res[is.na(res)] <- FALSE
  
  
  stats <- setNames(
    list(med, lower, upper),
    paste0(varname, c("_median", "_lower", "_upper"))
  )
  
  return(list(outlier = res, stats = stats))
}

# top_n_genes의 비율을 계산
# input :
#   counts_mat = gene x cell
#   n = the number of top_n genes
# output: 
#   pct = proportion of the top_n_genes
pct_counts_in_top_n_genes <- function(counts_mat, n = 20) {
  # counts_mat: dgCMatrix (gene × cell)
  total_counts <- Matrix::colSums(counts_mat)
  
  top_counts <- vapply(
    seq_len(ncol(counts_mat)),
    function(j) {
      col <- counts_mat[, j]
      if (sum(col) == 0) return(0)
      sum(sort(col, decreasing = TRUE)[seq_len(min(n, length(col)) )])
    },
    numeric(1)
  )
  
  pct <- ifelse(total_counts > 0, 100 * top_counts / total_counts, NA_real_)
  return(pct)
}