# 10x Genomics(Cell Ranger)
# 입출력 즉 tar.gz 내부가 matrix, barcode, feature or genes 일때
# 입력은 단일 파일로 받음
# 즉 tar 구조가 tar 풀기 -> 단일 폴더 -> 내부 matrix, barcode, feature or genes 일때
file_reader_10x_Genom_cell_ranger <- function(f, project = "scRNA") {
  # Seruat 설치 확인 여부
  
  # f = tar.gz 파일 경로
  extract_dir <- file.path(tempdir(), tools::file_path_sans_ext(basename(f)))
  
  
  # 압축해제
  utils::untar(f, exdir = extract_dir)
  
  # 자동으로 matrix.mtx 파일 위치 찾기
  mtx_file <- list.files(extract_dir, pattern = "^matrix\\.mtx$", 
                         full.names = TRUE, recursive = TRUE)
  
  # 그 폴더가 data_dir
  data_dir <- dirname(mtx_file[1])
  
  # 같은 폴더 안에서 barcodes / genes 파일도 자동으로 탐색
  feat_file <- list.files(data_dir, pattern = "^(features|genes)\\.tsv$", full.names = TRUE)
  barc_file <- list.files(data_dir, pattern = "^barcodes\\.tsv$", full.names = TRUE)
  
  # 이제 ReadMtx에 path만 넘기면 됨
  mat <- ReadMtx(
    mtx      = mtx_file[1],
    cells    = barc_file[1],
    features = feat_file[1],
    feature.column = 2   # 보통 2열(gene symbol) 사용
  )
  
  seurat_obj <- CreateSeuratObject(counts = mat, assay = "RNA", project = project)
  
  return (seurat_obj)
}


# set the project folder path
# I want setting the path with
# set folder
# excuse dir -> data -> GSE_data_RAW -> QC_outputs -> figure
# if select the regular data like normalization
# set the file_path with strage 
set_dir_with_10x_cellRanger <- function(data_name){
  folder_path <- getwd()
  
  
  data_path <- paste0(folder_path, "/data/" ,data_name)
  file_list <- list.files(path = data_path, pattern = "tar.gz$", full.names = TRUE) |> sort()
  
  out_dir <- file.path(data_path, "QC_outputs")
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  
  return(list(file_list = file_list, out_dir = out_dir))
}


get_QCpass_data_rds <- function(data_name){
  folder_path <- getwd()
  
  data_path <- paste0(folder_path, "/data/", data_name, "/QC_outputs")
  file_list <- list.files(path= data_path, pattern = "sce.rds$", full.name = TRUE) |> sort()
  
  return(file_list)
}