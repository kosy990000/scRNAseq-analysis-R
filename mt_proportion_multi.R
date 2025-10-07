# 필요한 패키지
library(dplyr)
library(ggplot2)
library(ggpubr)
library(patchwork)
library(stringr)

#-------------------------------#
# 1) 방향(큰 그룹) 계산 헬퍼
#-------------------------------#
.dir_label_df <- function(df, group_col = "group", y_col = "pct_counts_mt",
                          metric = c("median","mean"),
                          y_pad_ratio = 0.06) {
  metric <- match.arg(metric)
  # 그룹별 요약
  dir_tbl <- df %>%
    group_by(.data[[group_col]]) %>%
    summarise(
      med = median(.data[[y_col]], na.rm = TRUE),
      mn  = mean(.data[[y_col]],   na.rm = TRUE),
      .groups = "drop"
    )
  
  # 비교 기준 선택
  score_col <- if (metric == "median") "med" else "mn"
  g_hi <- dir_tbl[[group_col]][which.max(dir_tbl[[score_col]])]
  g_lo <- dir_tbl[[group_col]][which.min(dir_tbl[[score_col]])]
  
  # 라벨 위치
  y_max   <- max(df[[y_col]], na.rm = TRUE)
  y_range <- diff(range(df[[y_col]], na.rm = TRUE))
  y_pos   <- y_max + y_pad_ratio * ifelse(is.finite(y_range) && y_range > 0, y_range, abs(y_max) + 1)
  
  data.frame(
    group1 = g_lo, group2 = g_hi,
    y.position = y_pos,
    label = paste0(g_hi, " > ", g_lo, " (by ", metric, ")"),
    stringsAsFactors = FALSE
  )
}

#---------------------------------------------#
# 2) 메인 플로팅 함수 (한 판에 모든 주석/검정 포함)
#---------------------------------------------#
plot_violin_with_stats <- function(
    df,
    group_col = "group",
    y_col = "pct_counts_mt",
    comparisons = list(c("Quercetin","Vector")),  # Wilcoxon 비교쌍
    metric = c("median","mean"),                  # 방향 판단 기준
    violin_fill = "red", violin_alpha = 0.55,
    box_fill    = "skyblue", box_alpha = 0.8,
    mean_pt = list(shape = 16, size = 2.5, color="black"),
    median_pt = list(shape = 15, size = 3, color="black"),
    ylim = NULL,                 # c(ymin, ymax) 주면 coord_cartesian 적용
    y_pad_ratio = 0.06,          # 방향 라벨 여백
    title = NULL,
    order_levels = NULL
) {
  metric <- match.arg(metric)
  
  # 방향 라벨 데이터
  lab_df <- .dir_label_df(df, group_col, y_col, metric = metric, y_pad_ratio = y_pad_ratio)
  
  df[[group_col]] <- factor(df[[group_col]], levels = order_levels)
  
  
  p <- ggplot(df, aes(x = .data[[group_col]], y = .data[[y_col]])) +
    geom_violin(trim = TRUE, width = 1, fill = violin_fill, alpha = violin_alpha) +
    geom_boxplot(width = 0.18, outlier.shape = NA, fill = box_fill, color = "black", alpha = box_alpha) +
    do.call(stat_summary, c(list(fun = mean,   geom = "point"), mean_pt)) +
    do.call(stat_summary, c(list(fun = median, geom = "point"), median_pt)) +
    stat_compare_means(method = "wilcox.test",
                       label = "p.signif",
                       comparisons = comparisons) +
    #stat_pvalue_manual(lab_df, label = "label", tip.length = 0) +
    scale_y_continuous(labels = scales::comma) +
    labs(x = NULL, y = y_col, title = title) +
    theme_bw() +
    theme(legend.position = "none",
          plot.title = element_text(hjust = 0, face = "bold"))
  
  # 보기 범위만 자르는 안전한 방식
  if (!is.null(ylim)) {
    p <- p + coord_cartesian(ylim = ylim)
  }
  return(p)
}

#---------------------------------------------#
# 3) QC 전/후 한 번에 그려주는 래퍼
#---------------------------------------------#
plot_qc_before_after <- function(
    df_before, df_after,
    group_col = "group", y_col = "pct_counts_mt",
    comparisons = list(c("Quercetin","Vector")),
    metric = c("median","mean"),
    ylim_before = NULL, ylim_after = NULL,
    title_before = "Before QC", title_after = "After QC",
    order_levels = NULL
) {
  metric <- match.arg(metric)
  p1 <- plot_violin_with_stats(df_before, group_col, y_col,
                               comparisons = comparisons, metric = metric,
                               ylim = ylim_before, title = title_before, order_levels = order_levels)
  p2 <- plot_violin_with_stats(df_after,  group_col, y_col,
                               comparisons = comparisons, metric = metric,
                               ylim = ylim_after,  title = title_after, order_levels = order_levels)
  p1 / p2  # 위/아래로 배치 (patchwork)
}

#====================#
setwd("/home/ko/scRNA-analysis")
path <- getwd()
folder_path = paste0(path,"/data/GSE226225_RAW/QC_outputs")

# 또는 (권장: 운영체제 호환 안전)
df <- read.csv(file.path(folder_path, "dataDriven_ALL_files_QC_matrix_with_doublet.csv"))

order4 <- c("CTRL","ETO_day_0","ETO_day_1","ETO_day_2","ETO_day_4", "ETO_day_7", "ETO_day_10")

df2 <- df %>%
  mutate(
    # 파일명 안의 키워드로 그룹 매핑 (+ 오탈자/별칭 대응)
    group = case_when(
      str_detect(.data[["file"]], regex("CTRL", ignore_case = TRUE)) ~ "CTRL",
      str_detect(.data[["file"]], regex("ETO_day_10",      ignore_case = TRUE)) ~ "ETO_day_10",   # 'old'를 OIS로 매핑
      str_detect(.data[["file"]], regex("ETO_day_2",    ignore_case = TRUE)) ~ "ETO_day_2", # 'ddls' 오타 매핑
      str_detect(.data[["file"]], regex("ETO_day_4", ignore_case = TRUE)) ~ "ETO_day_4",
      str_detect(.data[["file"]], regex("ETO_day_7", ignore_case = TRUE)) ~ "ETO_day_7",
      str_detect(.data[["file"]], regex("ETO_day_1", ignore_case = TRUE)) ~ "ETO_day_1",
      str_detect(.data[["file"]], regex("ETO_day_0", ignore_case = TRUE)) ~ "ETO_day_0",
      TRUE ~ NA_character_
    )
  ) %>%
  filter(!is.na(group)) %>%                      # ‘Other’ 제거
  mutate(group = factor(group, levels = order4)) # 순서 고정

df <- df2
rm(df2)

is.factor(df$group)                 # TRUE ?
levels(df$group)                    # 원하는 순서?
unique(df$group)                    # 실제 등장하는 값?
anyNA(df$group)                     # NA 있으면 순서·매핑 깨짐

#데이터 분리
# 리스트로 나누기
df_split <- split(df, df$group)

# 리스트 이름 확인
names(df_split)

comparisons <-  lapply(order4[-1], function(g) c(g, "CTRL"))

out_prefix <- paste0("GSE226225", "data_driven")

#====================#


df_before <- df
df_after  <- df_before %>% filter(qc_status == "pass_all")


out_folder <- paste0(folder_path, "/figures/")

# 단일 플롯(예: QC 전만)
p_alone <- plot_violin_with_stats(
  df_before,
  group_col = "group",
  y_col = "pct_counts_mt",
  comparisons = comparisons,
  ylim = c(0, 0.2),
  metric = "median",
  title = "Before QC",
  order_levels = order4
)

p_alone

# 전/후 결합 플롯
p_pair <- plot_qc_before_after(
  df_before, df_after,
  group_col = "group", y_col = "pct_counts_mt",
  comparisons = comparisons,
  metric = "median",
  ylim_before = c(0, 0.2),
  title_before = "Before QC",
  title_after  = "After QC",
  order_levels = order4 #나열 순서
)
p_pair

ggsave(
  file.path(paste0(out_folder,out_prefix, "_beforeQC.jpg")),
  p_alone,
  width  = 6,
  height = 6,
  dpi    = 300
)

ggsave(
  file.path(paste0(out_folder ,out_prefix, "_changeQC.jpg")),
  p_pair,
  width  = 6,
  height = 14,
  dpi    = 600
)

