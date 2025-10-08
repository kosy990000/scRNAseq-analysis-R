suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(dplyr)
  library(scDblFinder)           # 더블릿 탐지
  library(SingleCellExperiment)  # as.SingleCellExperiment()
})

# input x (matrix, nmad, na_to_false, verbose=출력 시 나올 것)
# return
#   : res$outlier -> boolean vector
#   : res$stats -> the list of varname, med, lower, upper
is_outlier <- function(x, nmads = 5, na_to_false = TRUE, varname = NULL, verbose = TRUE) {
  # --- 변수명 자동 처리: substitute()는 Fallback용 ---
  if (is.null(varname)) {
    varname <- tryCatch(
      {
        raw_name <- deparse(substitute(x))
        # md[[v]] 같은 입력 제거
        gsub(".*\\[\\[|\\]\\].*", "", raw_name)
      },
      error = function(e) "var"
    )
  }
  
  # --- 중앙값 & MAD 계산 ---
  med <- median(x, na.rm = TRUE)
  dev <- mad(x, center = med, constant = 1, na.rm = TRUE)
  
  # --- 이상치 계산 불가한 경우 ---
  if (!is.finite(med) || !is.finite(dev) || dev == 0) {
    res <- rep(FALSE, length(x))
    stats <- list(
      var    = varname,
      median = med,
      mad    = dev,
      lower  = NA,
      upper  = NA
    )
    if (verbose) print(stats)
    return(list(outlier = res, stats = stats))
  }
  
  # --- MAD 기반 cutoff ---
  lower <- med - nmads * dev
  upper <- med + nmads * dev
  res   <- (x < lower) | (x > upper)
  if (na_to_false) res[is.na(res)] <- FALSE
  
  stats <- data.frame(
    var    = varname,
    median = med,
    mad    = dev,
    lower  = lower,
    upper  = upper,
    stringsAsFactors = FALSE
  )
  
  if (verbose) {
    cat(sprintf("[QC] %-15s median=%.4f  MAD=%.4f  lower=%.4f  upper=%.4f\n",
                varname, med, dev, lower, upper))
  }
  
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



# 1) 샘플 라벨 & 셀ID 접두사 부여
add_sample_and_prefix <- function(seu, sample_id) {
  stopifnot(inherits(seu, "Seurat"))
  seu$sample <- sample_id
  colnames(seu) <- paste0(sample_id, "_", colnames(seu))
  seu
}

# 2) QC 지표 계산(Seurat 내부). species_pattern: human "^MT-", mouse "^mt-"
attach_qc_metrics <- function(seu, topn = 20, species_pattern = "^MT-") {
  stopifnot(inherits(seu, "Seurat"))
  
  # (1) 미토콘드리아 비율 계산 — SCE 스케일(0~1)로 변환
  seu[["percent.mt"]] <- PercentageFeatureSet(seu, pattern = species_pattern)
  seu$percent.mt <- seu$percent.mt / 100     # ← fraction으로 변환 (SCE와 동일)
  
  # (2) 상위 20개 유전자의 비율 계산 — 0~100 스케일 유지
  cmat <- GetAssayData(seu, assay = DefaultAssay(seu), slot = "counts")
  seu$top20_pct <- pct_counts_in_top_n_genes(cmat, n = topn)  # 이미 % 단위
  
  # (3) log1p 지표 추가 — SCE 기준에 맞춤
  seu$log1p_nCount_RNA   <- log1p(seu$nCount_RNA)
  seu$log1p_nFeature_RNA <- log1p(seu$nFeature_RNA)
  
  return(seu)
}

# 3) scDblFinder 실행: Seurat -> SCE 변환 후 라벨만 meta에 붙여서 반환(제거 X)
attach_doublet_labels <- function(seu, samples_col = "sample", dbr = 0.06, seed = 9999) {
  set.seed(seed)
  sce <- as.SingleCellExperiment(seu)          # Seurat -> SCE
  # 샘플 컬럼 없으면 NULL
  samples_vec <- if (samples_col %in% colnames(colData(sce))) sce[[samples_col]] else NULL
  sce_dbl <- scDblFinder(sce, dbr = dbr, samples = samples_vec)
  lab <- as.character(colData(sce_dbl)$scDblFinder.class)
  seu$scDblFinder.class <- lab
  seu$doublet_flag <- lab %in% c("doublet", "ambiguous")
  seu
}

# 4) 데이터 주도형 QC 뷰 생성기(원본 보존, 조건만 바꿔 subset)
make_qc_view <- function(
    seu,
    min_features = 500,
    max_mt = 15,        # percent 기준
    min_counts = 0,
    drop_doublet = TRUE,
    keep_ambiguous = FALSE
) {
  md <- seu@meta.data
  stopifnot(all(c("nFeature_RNA","nCount_RNA","percent.mt") %in% names(md)))
  keep <- md$nFeature_RNA >= min_features &
    md$nCount_RNA   >= min_counts   &
    md$percent.mt   <= max_mt
  if ("doublet_flag" %in% names(md) && drop_doublet) {
    if (keep_ambiguous && "scDblFinder.class" %in% names(md)) {
      keep <- keep & !(md$scDblFinder.class %in% c("doublet"))
    } else {
      keep <- keep & !md$doublet_flag
    }
  }
  subset(seu, cells = rownames(md)[keep])
}