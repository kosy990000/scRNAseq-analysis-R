

source("utills/utillsFileIO.R")  # file_reader_10x_Genom_cell_ranger(): Seurat 반환
source("utills/utillsQC.R")    # 필요시




data_name <- "GSE239591_RAW"
dir_list  <- get_folder_list(data_name)
file_list <- dir_list$file_list
out_dir   <- dir_list$out_dir

seurat_merged <- NULL
qc_tables     <- list()
qc_thr_list   <- list()



for (f in file_list) {
  message("\n[Processing] ", f)
  fname      <- basename(f)
  sample_id  <- tools::file_path_sans_ext(basename(tools::file_path_sans_ext(fname)))
  
  # 1) 로드: 이미 Seurat 객체를 반환(네 구현)
  seu <- file_reader_10x_Genom_cell_ranger(f, project = "scRNA")
  
  # 2) 라벨/접두사
  seu <- add_sample_and_prefix(seu, sample_id = sample_id)
  
  
  
  # 3) QC 지표(Seurat 내부)
  #  - human이면 "^MT-", mouse면 "^mt-"
  species_pat <- if (any(grepl("^MT-", rownames(seu)))) "^MT-" else "^mt-"
  seu <- attach_qc_metrics(seu, species_pattern = species_pat)

  md  <- seu@meta.data
  
  vars <- list(
    log1p_nFeature_RNA = list(nmads = 5),
    log1p_nCount_RNA   = list(nmads = 5),
    percent.mt   = list(nmads = 3),
    top20_pct    = list(nmads = 5)
  )
  
  # MAD 기반 threshold 계산
  thr_df <- do.call(rbind, list(
    is_outlier(md$log1p_nFeature_RNA, 5, varname="log1p_nFeature_RNA")$stats,
    is_outlier(md$log1p_nCount_RNA, 5, varname="log1p_nCount_RNA")$stats,
    is_outlier(md$percent.mt, 3, varname="percent.mt")$stats,
    is_outlier(md$top20_pct, 5, varname="top20_pct")$stats
  ))
  thr_df
  
  keep <- with(md,
               (log1p_nFeature_RNA >= thr_df$lower[thr_df$var=="log1p_nFeature_RNA"] &
                  log1p_nFeature_RNA <= thr_df$upper[thr_df$var=="log1p_nFeature_RNA"]) &
                 
                 (log1p_nCount_RNA >= thr_df$lower[thr_df$var=="log1p_nCount_RNA"] &
                    log1p_nCount_RNA <= thr_df$upper[thr_df$var=="log1p_nCount_RNA"]) &
                 
                 (percent.mt >= thr_df$lower[thr_df$var=="percent.mt"] &
                    percent.mt <= thr_df$upper[thr_df$var=="percent.mt"]) &
                 
                 (top20_pct >= thr_df$lower[thr_df$var=="top20_pct"] &
                    top20_pct <= thr_df$upper[thr_df$var=="top20_pct"])
  )
  # QC 테이블용 메타데이터에 pass flag 추가
  qc_md <- seu@meta.data %>%
    tibble::rownames_to_column("barcode_prefixed") %>%
    mutate(
      sample = seu$sample,
      qc_pass_basic = keep  # 1차 QC 통과 여부 (TRUE/FALSE)
    )
  

  seu_basic <- subset(seu, cells = rownames(md)[keep])
  rm(seu)
  
  message(sprintf("[QC filter] %d / %d cells retained (%.1f%%)",
                  sum(keep), nrow(md), 100 * mean(keep)))
  
  
  seu_qc <- attach_doublet_labels(seu_basic, dbr = 0.06, seed = 9999)
  rm(seu_basic)
  
  n_before <- ncol(seu_qc)
  seu_qc <- subset(seu_qc, subset = scDblFinder.class == "singlet")
  n_after <- ncol(seu_qc)
  message(sprintf("[Doublet removal] %d → %d cells (%.1f%% kept)",
                  n_before, n_after, 100 * n_after / n_before))
  
  # Doublet pass flag 생성
  doublet_pass <- colnames(seu_qc)
  
  qc_md <- qc_md %>%
    mutate(
      qc_pass_doublet = barcode_prefixed %in% doublet_pass,
      qc_pass_all = qc_pass_basic & qc_pass_doublet
    )
  
  # QC 결과 테이블 저장
  qc_tables[[fname]] <- qc_md %>%
    select(barcode_prefixed, sample, nCount_RNA, nFeature_RNA, percent.mt,
           qc_pass_basic, qc_pass_doublet, qc_pass_all)
  
  # 병합
  seurat_merged <- if (is.null(seurat_merged)) seu_qc else merge(seurat_merged, seu_qc)
  
  rm(seu_qc, md, thr_df, keep, n_before, n_after, fname, sample_id)
  gc()
}

# (1) 통합 threshold 파일

thr_all <- dplyr::bind_rows(qc_thr_list, .id = "sample")
thr_path <- file.path(out_dir, "all_samples_QC_thresholds.csv")
write.csv(thr_all, thr_path, row.names = FALSE)

# (2) 통합 QC metrics 파일
qc_all <- dplyr::bind_rows(qc_tables, .id = "sample")
qc_path <- file.path(out_dir, "all_samples_QC_metrics.csv")
write.csv(qc_all, qc_path, row.names = FALSE)

# (3) QC 이후 Seurat 객체 저장 (필요 시)
seu_rds_path <- file.path(out_dir, sprintf("passBasicQC.rds"))
saveRDS(seurat_merged , seu_rds_path, compress = "xz")



