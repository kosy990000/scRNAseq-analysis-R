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

# ---------- 유틸 함수 ----------  #

compute_qc_limits <- function(df, use_quantile = FALSE, q = c(0.01, 0.99)) {
  rng <- function(x) {
    if (use_quantile) range(quantile(x, q, na.rm = TRUE)) else range(x, na.rm = TRUE)
  }
  list(
    total_counts      = rng(df$total_counts),        # for x/y log scales
    n_genes_by_counts = rng(df$n_genes_by_counts),   # for y log scale
    pct_counts_mt     = rng(df$pct_counts_mt)        # for color & y (linear)
  )
}

# 입력 데이터를 받아 qc_before와 qc_After를 보는 함수
qc_plots <- function(df) {
  # 로그 스케일 안전: 0/NA 제거용
  df_log <- subset(
    df,
    is.finite(total_counts) & is.finite(n_genes_by_counts) &
      total_counts > 0 & n_genes_by_counts > 0
  )
  
  # 1) total_counts 바이올린 (log Y)
  p1 <- ggplot(df_log, aes(x = condition, y = total_counts, fill = condition)) +
    geom_violin(trim = TRUE) +
    geom_jitter(width = 0.2, size = 0.3, alpha = 0.4) +
    scale_y_log10(labels = scales::comma) +
    labs(x = NULL, y = "n_total_counts") +
    theme_bw() +
    theme(legend.position = "none")
  
  # 2) pct_counts_mt 바이올린 (linear Y)
  p2 <- ggplot(df, aes(x = condition, y = pct_counts_mt, fill = condition)) +
    geom_violin(trim = TRUE) +
    geom_jitter(width = 0.2, size = 0.3, alpha = 0.4) +
    labs(x = NULL, y = "pct_counts_mt") +
    theme_bw() +
    theme(legend.position = "none")
  
  # 3) total_counts 히스토그램 (facet, 축 공유)
  p3 <- ggplot(df, aes(x = total_counts)) +
    geom_histogram(bins = 100, fill = "steelblue", color = "black") +
    scale_x_log10(labels = scales::comma) +
    labs(x = "total_counts", y = "Frequency") +
    theme_bw() +
    facet_wrap(~condition, nrow = 1, scales = "fixed")
  
  # 4) total_counts vs n_genes_by_counts (facet, log-log, 축 공유)
  df_log2 <- subset(
    df,
    is.finite(total_counts) & is.finite(n_genes_by_counts) &
      total_counts > 0 & n_genes_by_counts > 0
  )
  p4 <- ggplot(df_log2, aes(x = total_counts, y = n_genes_by_counts)) +
    geom_point(alpha = 0.5, size = 0.6, color = "steelblue") +
    scale_x_log10(labels = scales::comma) +
    scale_y_log10(labels = scales::comma) +
    labs(x = "total_counts", y = "n_genes_by_counts") +
    theme_bw() +
    facet_wrap(~condition, nrow = 1, scales = "fixed")
  
  # 5) total_counts vs mt_counts (facet, log-log, 축 공유)
  df$mt_counts <- df$pct_counts_mt * df$total_counts
  df_log_mt <- subset(
    df,
    is.finite(total_counts) & is.finite(mt_counts) &
      total_counts > 0 & mt_counts > 0
  )
  p5 <- ggplot(df_log_mt, aes(x = total_counts, y = mt_counts)) +
    geom_point(alpha = 0.5, size = 0.6, color = "steelblue") +
    scale_x_log10(labels = scales::comma) +
    scale_y_log10(labels = scales::comma) +
    labs(x = "total_counts", y = "mt_counts") +
    theme_bw() +
    facet_wrap(~condition, nrow = 1, scales = "fixed") +   # ✅ 여기 괄호 닫기
    geom_abline(slope = 1, intercept = log10(0.1), 
                color = "red", linetype = "dashed") +
    geom_abline(slope = 1, intercept = log10(0.2), 
                color = "blue", linetype = "dashed") +
    geom_abline(slope = 1, intercept = log10(0.3), 
                color = "green", linetype = "dashed")
  
  list(p1 = p1, p2 = p2, p3 = p3, p4 = p4, p5 = p5)
}


is_outlier <- function(x, nmads = 5,  na_to_false = TRUE) {
  varname <- deparse(substitute(x))  # 인자로 들어온 변수명 추출
  
  if (length(x) == 0 || all(is.na(x))) {
    return(rep(FALSE, length(x)))
  }
  med <- median(x, na.rm = TRUE)
  dev <- mad(x, center = med, constant = 1, na.rm = TRUE)
  
  if (!is.finite(med) || !is.finite(dev) || dev == 0) {
    res <- rep(FALSE, length(x))
    if (verbose) cat("[", varname, "] Lower/Upper 미정 (dev=0 또는 NA). 모두 FALSE 처리.\n")
    return(res)
  }
  
  lower <- med - nmads * dev
  upper <- med + nmads * dev
  res <- (x < lower) | (x > upper)
  
  if (na_to_false) res[is.na(res)] <- FALSE
  
  cat("[", varname, "] Lower:", lower, "Upper:", upper, "\n")
  return(res)
}


pct_counts_in_top_n_genes <- function(counts_mat, n = 20) {
  total_counts <- Matrix::colSums(counts_mat)
  top_counts <- vapply(
    seq_len(ncol(counts_mat)),
    function(j) {
      col <- counts_mat[, j]
      if (sum(col) == 0) return(0)
      sum(sort(col, decreasing = TRUE)[seq_len(min(n, length(col)))])
    },
    numeric(1)
  )
  ifelse(total_counts > 0, top_counts / total_counts, NA_real_)  # ← 0–1 스케일 유지
}

# ---------- 입력/출력 ----------
folder_path <- "scRNA-analysis"
file_list <- list.files(path = folder_path, pattern = "countmatrix.txt$", full.names = TRUE) |> sort()

out_dir <- file.path(folder_path, "QC_outputs")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)


all_cells_qc <- list()   # 파일별 per-cell QC 테이블 누적
summary_rows <- list()   # 파일별 요약 누적

# ---------- 메인 루프 ----------
for (f in file_list) {
  
  message("\n[Processing] ", f)
  fname <- basename(f)
  # count matrix 불러오기
  counts <- read.table(f, 
                       header = TRUE, row.names = 1, sep="\t")
  
  # annotation 불러오기 임시적인 연결
  anno <- read.table(paste0(folder_path,"/YoungOld_Liver_tissue_annotation.txt"), 
                     header = TRUE, sep="\t")
  
  # Seurat object 생성
  seurat_obj <- CreateSeuratObject(counts = counts)
  
  rownames(anno) <- anno$cell_id
  seurat_obj <- AddMetaData(seurat_obj, metadata = anno)
  
  # 읽기 → SCE 변환
  sce <- as.SingleCellExperiment(seurat_obj)
  rm(serat_obj); gc()
  
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
  set.seed(123)
  sce_dbl <- scDblFinder(
    sce_after,
    dbr = 0.06,
    samples = if ("sample" %in% colnames(colData(sce_after))) sce_after$sample else NULL
  )
  sce_clean <- sce_dbl[, sce_dbl$scDblFinder.class == "singlet"]
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
  
  message(sprintf(
    paste0("  Cells(before): %d | not_mt_outliers %d |mt-outliers: %d | any-outliers(basic): %d | ",
           "pass-basic: %d | any-outliers(+doublet): %d | pass-all: %d"),
    n_cells_before, n_not_mt_outliers ,n_mt_outliers, n_any_outliers_basic, n_pass_basic,
    n_any_outliers_w_double, n_pass_all
  ))
  
  # per-cell QC 행렬(모든 지표/라벨 포함)
  qc_matrix_file <- qc_merged %>%
    dplyr::mutate(file = fname) %>%
    dplyr::select(
      barcode, file,
      total_counts, log1p_total_counts,
      n_genes_by_counts, log1p_n_genes_by_counts,
      pct_counts_in_top_20_genes, pct_counts_mt,
      # 원하는 QC flag들
      outlier_not_mt,
      mt_outlier,               # 2번
      outlier_any,              # 3번
      outlier_with_doublet,     # 4번
      doublet_flag,             # 5번
      qc_status                 # 6번 (최종)
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
  
  # 메모리 정리 (전역과 이름 충돌 방지)
  rm(sce_after, sce_dbl, sce_clean, qc_before, qc_merged, qc_matrix_file,
     keep_cells, mt_outlier, total_outlier, genes_outlier, toptwenty_outlier)
  gc()
}

# ---------- 통합 저장 ----------
summary_df <- dplyr::bind_rows(summary_rows) |> dplyr::arrange(file)
write.csv(summary_df, file = file.path(out_dir, "QC_summary_by_file.csv"), row.names = FALSE)

qc_all_cells <- dplyr::bind_rows(all_cells_qc)
write.csv(qc_all_cells, file = file.path(out_dir, "ALL_files_QC_matrix_with_doublet.csv"), row.names = FALSE)

message("\n[Done] Outputs written to: ", out_dir)

