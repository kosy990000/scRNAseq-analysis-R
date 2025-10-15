# BiocManager 설치 (없다면)

suppressPackageStartupMessages({
  library(scran)
  library(ggplot2)
  library(Seurat)
  library(Matrix)
  library(patchwork)
  library(stringr)
  library(scDblFinder)
  library(dplyr)
  # library(readr)  # readr 미사용; base R write.csv 사용
})


is_outlier <- function(x, nmads = 5, na_to_false = TRUE, verbose = TRUE) {
  varname <- deparse(substitute(x))  # 인자로 들어온 변수명
  
  med <- median(x, na.rm = TRUE)
  dev <- mad(x, center = med, constant = 1, na.rm = TRUE)
  
  if (!is.finite(med) || !is.finite(dev) || dev == 0) {
    res <- rep(FALSE, length(x))
    if (verbose) cat("[", varname, "] MAD=0 or NA → 모두 FALSE\n")
    return(res)
  }
  
  lower <- med - nmads * dev
  upper <- med + nmads * dev
  res <- (x < lower) | (x > upper)
  
  if (na_to_false) res[is.na(res)] <- FALSE
  
  if (verbose) {
    cat("[", varname, "] median:", med,
        " MAD:", dev,
        " Lower:", lower,
        " Upper:", upper, "\n")
  }
  return(res)
}

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



# ------------------------- folder path ---------------------- #

setwd("/home/ko/scRNA-analysis")
data_path <- "/home/ko/scRNA-analysis/data/GSE226225_RAW"

file_prefixes <- unique(gsub("_(barcodes|features|matrix).*", "", basename(list.files(data_path))))

out_dir <- file.path(data_path, "QC_outputs")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)


all_cells_qc <- list()   # 파일별 per-cell QC 테이블 누적
summary_rows <- list()   # 파일별 요약 누적

for (prefix in file_prefixes) {
  message("\n[Processing] ", prefix)
  fname = prefix
  mtx_file  <- file.path(data_path, paste0(prefix, "_matrix.mtx.gz"))
  feat_file <- file.path(data_path, paste0(prefix, "_features.tsv.gz"))
  barc_file <- file.path(data_path, paste0(prefix, "_barcodes.tsv.gz"))
  
  # Seurat ReadMtx
  mat <- ReadMtx(
    mtx = mtx_file,
    cells = barc_file,
    features = feat_file,
    feature.column = 2   # 보통 2번째 column이 gene symbol
  )
  
  # Seurat object 생성
  
  seurat_obj <- CreateSeuratObject(counts = mat, assay = "RNA", project = "scRNA")
  sce <- as.SingleCellExperiment(seurat_obj)
  
  # Remove the original Seurat object and free memory
  rm(seurat_obj, mtx, mtx_file, feat_file, barc_file)
  gc()
  
  # QC 지표 (pct는 0–1 스케일)
  sce$total_counts             <- Matrix::colSums(counts(sce))
  sce$log1p_total_counts       <- log(sce$total_counts + 1)
  sce$n_genes_by_counts        <- Matrix::colSums(counts(sce) > 0)
  sce$log1p_n_genes_by_counts  <- log(sce$n_genes_by_counts + 1)
  sce$pct_counts_in_top_20_genes <- pct_counts_in_top_n_genes(counts(sce), n = 20)
  
  mt_mask <- grepl("^mt-", rownames(sce), ignore.case = TRUE)
  # ← 0–1 스케일 유지 (퍼센트 ×100 안 함)
  sce$pct_counts_mt <- Matrix::colSums(counts(sce)[mt_mask, , drop = FALSE]) / sce$total_counts
  rm(mt_mask)
  
  # per-cell QC (Before)
  qc_before <- as.data.frame(colData(sce))
  qc_before$barcode    <- rownames(qc_before)
  qc_before$file       <- fname
  qc_before$condition  <- "Before"
  
  # outlier 플래그
  mt_outlier        <- is_outlier(qc_before$pct_counts_mt, 3)
  total_outlier     <- is_outlier(qc_before$log1p_total_counts, 5)
  genes_outlier     <- is_outlier(qc_before$log1p_n_genes_by_counts, 5)
  toptwenty_outlier <- is_outlier(qc_before$pct_counts_in_top_20_genes, 5)
  
  qc_before$mt_outlier        <- mt_outlier
  qc_before$total_outlier     <- total_outlier
  qc_before$genes_outlier     <- genes_outlier
  qc_before$toptwenty_outlier <- toptwenty_outlier
  
  qc_before$outlier_not_mt       <- total_outlier | genes_outlier | toptwenty_outlier
  qc_before$outlier_any       <- mt_outlier | total_outlier | genes_outlier | toptwenty_outlier
  
  # 집계(베이직 QC 전후)
  n_cells_before <- nrow(qc_before)
  n_not_mt_outliers  <- sum(qc_before$outlier_not_mt, na.rm = TRUE)
  n_any_outliers <- sum(qc_before$outlier_any, na.rm = TRUE)
  n_pass_basic   <- n_cells_before - n_any_outliers
  
  # basic QC 적용
  keep_cells <- qc_before$barcode[!qc_before$outlier_any]
  sce_after  <- sce[, keep_cells, drop = FALSE]
  rm(sce); gc()
  
  # Doublet 제거
  set.seed(9999)
  sce_dbl <- scDblFinder(
    sce_after,
    dbr = 0.06,
    samples = if ("sample" %in% colnames(colData(sce_after))) sce_after$sample else NULL
  )
  sce_clean   <- sce_dbl[, sce_dbl$scDblFinder.class == "singlet"]
  n_after_all <- ncol(sce_clean)
  message(sprintf("  Cells(after all QC incl. doublet removal): %d", n_after_all))
  
  # 더블릿 정보 병합 → 상태 라벨
  doublet_info <- data.frame(
    barcode = colnames(sce_dbl),
    scDblFinder.class = as.character(colData(sce_dbl)$scDblFinder.class),
    stringsAsFactors = FALSE
  )
  
  qc_merged <- qc_before %>%
    dplyr::left_join(doublet_info, by = "barcode") %>%
    dplyr::mutate(
      doublet_flag = ifelse(is.na(scDblFinder.class), FALSE,
                            scDblFinder.class %in% c("doublet","ambiguous")),
      outlier_with_doublet = outlier_any | doublet_flag,
      qc_status = dplyr::case_when(
        outlier_any ~ "fail_metric",
        !outlier_any & doublet_flag ~ "fail_doublet",
        TRUE ~ "pass_all"
      )
    )
  
  # 요약 업데이트
  n_cells_before          <- nrow(qc_merged)
  n_mt_outliers           <- sum(qc_merged$mt_outlier,   na.rm = TRUE)
  n_any_outliers_basic    <- sum(qc_merged$outlier_any,  na.rm = TRUE)
  n_any_outliers_w_double <- sum(qc_merged$outlier_with_doublet, na.rm = TRUE)
  n_pass_basic            <- n_cells_before - n_any_outliers_basic
  n_pass_all              <- sum(qc_merged$qc_status == "pass_all", na.rm = TRUE)
  
  
  # per-cell QC 행렬(모든 지표/라벨 포함)
  ## [수정 2] select는 '컬럼명'만. 지금은 정상적으로 컬럼이 생겼으니 그대로 사용 가능
  qc_matrix_file <- qc_merged %>%
    dplyr::mutate(file = fname) %>%
    dplyr::select(
      barcode, file,
      total_counts, log1p_total_counts,
      n_genes_by_counts, log1p_n_genes_by_counts,
      pct_counts_mt,
      outlier_not_mt,      # total 기반 outlier
      mt_outlier,          # mt 기반 outlier
      outlier_any,         # metric 기반 종합
      outlier_with_doublet,# metric + doublet
      doublet_flag,
      qc_status
    )
  
  # 누적
  all_cells_qc[[fname]] <- qc_matrix_file
  summary_rows[[fname]] <- tibble::tibble(
    file                  = fname,
    cells_before          = n_cells_before,
    mt_outliers           = n_mt_outliers,
    not_mt_outliers       = n_not_mt_outliers,
    any_outliers_basic    = n_any_outliers_basic,
    cells_after_basic     = n_pass_basic,
    any_outliers_w_double = n_any_outliers_w_double,
    cells_after_allQC     = n_pass_all
  )
  
  # 메모리 정리: 정의 안 된 변수는 rm 대상에서 제거
  rm(sce_after, sce_dbl, sce_clean, qc_before, qc_merged, qc_matrix_file,
     keep_cells)  # <- 존재하지 않는 genes_outlier, toptwenty_outlier 제거
  gc()
  
}


# ---------- 통합 저장 ----------

summary_df <- dplyr::bind_rows(summary_rows) |> dplyr::arrange(file)
write.csv(summary_df, file = file.path(out_dir, "dataDriven_QC_summary_by_file.csv"), row.names = FALSE)

qc_all_cells <- dplyr::bind_rows(all_cells_qc)
write.csv(qc_all_cells, file = file.path(out_dir, "dataDriven_ALL_files_QC_matrix_with_doublet.csv"), row.names = FALSE)

message("\n[Done] Outputs written to: ", out_dir)
