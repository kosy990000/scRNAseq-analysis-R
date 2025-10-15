# BiocManager 설치 (없다면)


library(scran)
library(ggplot2)
library(zellkonverter)     # readH5AD
library(SingleCellExperiment)
library(Seurat)
library(Matrix)
library(scater)
library(patchwork)


# before 데이터로 축 범위를 뽑는 헬퍼
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




qc_plots <- function(df) {
  # 로그 스케일 안전: 0/NA 제거
  df_log <- subset(df,
                   is.finite(total_counts) & is.finite(n_genes_by_counts) &
                     total_counts > 0 & n_genes_by_counts > 0)
  
  # 1) total_counts 바이올린
  p1 <- ggplot(df_log, aes(x = condition, y = total_counts, fill = condition)) +
    geom_violin(trim = TRUE) +
    geom_jitter(width = 0.2, size = 0.3, alpha = 0.4) +
    scale_y_log10(labels = scales::comma) +
    labs(x = NULL, y = "n_total_counts") +
    theme_bw() +
    theme(legend.position = "none")
  
  # 2) pct_counts_mt 바이올린 (선형축)
  p2 <- ggplot(df, aes(x = condition, y = pct_counts_mt, fill = condition)) +
    geom_violin(trim = TRUE) +
    geom_jitter(width = 0.2, size = 0.3, alpha = 0.4) +
    labs(x = NULL, y = "pct_counts_mt") +
    theme_bw() +
    theme(legend.position = "none")
  
  # 히스토그램: facet으로 before/after 나눔
  p3 <- ggplot(df, aes(x = total_counts)) +
    geom_histogram(bins = 100, fill = "steelblue", color = "black") +
    scale_x_log10(labels = scales::comma) +
    labs(x = "total_counts", y = "Frequency") +
    theme_bw() +
    facet_wrap(~condition, nrow = 1, scales = "fixed")  # VS 비교
  
  # 산점도: facet으로 before/after 나눔
  df_log <- subset(df, total_counts > 0 & n_genes_by_counts > 0 & 
                     is.finite(total_counts) & is.finite(n_genes_by_counts))
  
  p4 <- ggplot(df_log, aes(x = total_counts, y = n_genes_by_counts)) +
    geom_point(alpha = 0.5, size = 0.6, color = "steelblue") +
    scale_x_log10(labels = scales::comma) +
    scale_y_log10(labels = scales::comma) +
    labs(x = "total_counts", y = "n_genes_by_counts") +
    theme_bw() +
    facet_wrap(~condition, nrow = 1, scales = "fixed")
  
  list(p1 = p1, p2 = p2, p3 = p3, p4 = p4)
}

is_outlier <- function(x, nmads = 5) {
  med <- median(x, na.rm = TRUE)
  dev <- mad(x, center = med, constant = 1, na.rm = TRUE)
  (x < (med - nmads * dev)) | (x > (med + nmads * dev))
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

sce <- readH5AD("/GSE182416_old_2.h5ad")  # 필요시 절대경로 권장

# SCE 내부 확인
assayNames(sce)   # 예: "counts", "logcounts" 유무 확인

assays(sce)$counts <- assays(sce)$X


assays(sce)$X <- NULL

# 4 quality control matrix
sce$total_counts <- Matrix::colSums(counts(sce)) # n_cells_count
sce$log1p_total_counts <- log(sce$total_counts + 1)   
sce$n_genes_by_counts <- Matrix::colSums(counts(sce) > 0)
sce$log1p_n_genes_by_counts <- log(sce$n_genes_by_counts + 1)
sce$pct_counts_in_top_20_genes <- pct_counts_in_top_n_genes(counts(sce), n = 20)

mt_mask <- grepl("^mt-", rownames(sce), ignore.case = TRUE)
sce$pct_counts_mt <- Matrix::colSums(counts(sce)[mt_mask, , drop=FALSE]) / sce$total_counts * 100

rm(mt_mask)


qc_before <- as.data.frame(colData(sce))


qc_before$outlier <- (
  is_outlier(qc_before$log1p_total_counts, 5) |
    is_outlier(qc_before$log1p_n_genes_by_counts, 5) |
    is_outlier(qc_before$pct_counts_in_top_20_genes, 5)
)

table(qc_before$outlier)

qc_after <- subset(qc_before, outlier == FALSE)

qc_before$condition <- "Before"
qc_after$condition  <- "After"

qc_all <- bind_rows(qc_before, qc_after) %>%
  mutate(condition = factor(condition, levels = c("Before", "After")))

# QC 패널 저장
plots  <- qc_plots(qc_all)
figure <- (plots$p1 | plots$p2) / (plots$p3 | plots$p4)
ggsave("qc_panel_before_after.pdf", figure, width = 12, height = 10, dpi = 300)

# ---- 여기서 '먼저' keep_cells 계산하고 정리하자 ----
keep_cells <- intersect(colnames(sce), rownames(qc_after))
sce_after  <- sce[, keep_cells, drop = FALSE]

# Seurat 변환
cnt <- counts(sce_after)
if (!inherits(cnt, "dgCMatrix")) {
  cnt <- as(Matrix::Matrix(cnt, sparse = TRUE), "dgCMatrix")
}
meta <- as.data.frame(colData(sce_after))
stopifnot(identical(rownames(meta), colnames(cnt)))  # 또 한 번 가드

s_after <- CreateSeuratObject(counts = cnt, meta.data = meta, assay = "RNA", project = "after")
saveRDS(s_after, file = "after_seurat.rds")

# 이제 정리
rm(qc_before, qc_after)