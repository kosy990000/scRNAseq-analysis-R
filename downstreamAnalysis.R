source("utills/utillsFileIO.R")

library(Seurat)
library(SingleCellExperiment)
library(zellkonverter)  # as.Seurat()


#------------------------------------------------------------------------------#
data_name <- "GSE161340_RAW"

file_list <- get_QCpass_data_rds(data_name)

# 1) SCE 읽기
sce_list <- lapply(file_list, readRDS)
names(sce_list) <- tools::file_path_sans_ext(basename(file_list))

# 2) 샘플 라벨/셀 이름 접두사
for (nm in names(sce_list)) {
  colData(sce_list[[nm]])$sample <- nm
  colnames(sce_list[[nm]]) <- paste0(nm, "_", colnames(sce_list[[nm]]))
}

# 3) 각 SCE를 Seurat로 변환
seurat_list <- lapply(sce_list, function(sce) {
  assn <- assayNames(sce)
  counts_layer    <- if ("counts" %in% assn) "counts" else assn[1]
  logcounts_layer <- if ("logcounts" %in% assn) "logcounts" else counts_layer
  
  zellkonverter::as.Seurat(sce, counts = counts_layer, data = logcounts_layer)
})

# 4) Seurat 병합 (합집합 feature, 0 padding; 셀ID에 파일명 접두사 유지)
if (length(seurat_list) == 1) {
  seu_merged <- seurat_list[[1]]
} else {
  seu_merged <- merge(
    x = seurat_list[[1]],
    y = seurat_list[-1],
    add.cell.ids = names(seurat_list)
  )
}

# 5) 후처리(선택)
seu_merged <- NormalizeData(seu_merged)
seu_merged <- FindVariableFeatures(seu_merged, nfeatures = 3000)
seu_merged <- ScaleData(seu_merged, verbose = FALSE)

# 6) 확인
print(seu_merged)
table(seu_merged$sample)


for (file in file_list) {
  
  fname <- basename(f)
  
  
}