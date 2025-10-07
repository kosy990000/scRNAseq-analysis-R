suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(scales)
  library(stringr)
  library(scran)
  library(patchwork)
  library(ggpubr)
})

# 분위수 기반 선형 축 한계 계산 (mt_counts 포함)
.compute_limits_linear <- function(df, use_quantile = TRUE, q = c(0.01, 0.99)) {
  rng <- function(x) if (use_quantile) as.numeric(quantile(x, q, na.rm = TRUE)) else range(x, na.rm = TRUE)
  list(
    total_counts      = rng(df$total_counts),
    n_genes_by_counts = rng(df$n_genes_by_counts),
    pct_counts_mt     = rng(df$pct_counts_mt),
    mt_counts         = rng(df$mt_counts)
  )
}



qc_plots_AB_linear <- function(
    df,
    group_col,
    group_A,
    group_B,
    condition_col   = "condition",
    pct_is_fraction = TRUE,
    facet_nrow = 1,
    facet_ncol = NULL,
    use_quantile = TRUE,          # ★ 분위수 사용 여부
    q = c(0.01, 0.99)             # ★ 하위/상위 분위수 (기본 1%~99%)
) {
  # 전처리: 그룹 필터 + 파생변수
  df2 <- df %>%
    filter(.data[[group_col]] %in% c(group_A, group_B)) %>%
    mutate(
      mt_counts   = if (pct_is_fraction) pct_counts_mt * total_counts else (pct_counts_mt/100) * total_counts,
      facet_group = factor(.data[[group_col]], levels = c(group_A, group_B)),
      # FALSE -> Pass, TRUE -> Outlier
      cond        = factor(ifelse(.data[[condition_col]] == TRUE, "Outlier", "Pass"),
                           levels = c("Pass","Outlier")),
      x_all       = "ALL"
    )
  
  
  # 축 한계 계산 (분위수 또는 전체 범위)
  lims <- .compute_limits_linear(df2, use_quantile = use_quantile, q = q)
  
  # 색: Pass=회색, Outlier=빨강
  col_vals  <- c("Pass" = "#9E9E9E", "Outlier" = "red")
  fill_vals <- c("Pass" = "#BDBDBD", "Outlier" = "red")
  
#--------------------------------- total_counts ---------------------------------#
  # 1) total_counts (violin + jitter) — y 한계 적용
  p1 <- df2 %>%
    dplyr::filter(cond == "Pass") %>%
    ggplot(aes(x = facet_group, y = total_counts, fill = cond)) +
    geom_violin(trim = TRUE, width = 1, alpha = 0.6) +
    scale_y_continuous(labels = scales::comma) +
    scale_fill_manual(values = fill_vals) +
    labs(x = NULL, y = "total_counts", fill = "QC Status") +
    theme_bw() +
    theme(legend.position = "top")

#--------------------------------------------------------------------------------#  
#-------------------------- p2 total count frequency ----------------------------#
  
  # 3) total_counts 히스토그램 — strip 제거, x축 그룹
  p2 <- ggplot(df2, aes(x = total_counts, fill = cond)) +
    geom_histogram(bins = 100, color = "black", alpha = 0.45, position = "identity") +
    scale_x_continuous(labels = scales::comma) +
    scale_fill_manual(values = fill_vals) +
    coord_cartesian(xlim = lims$total_counts) +
    labs(x = "total_counts", y = "Frequency", fill = "QC Status") +
    theme_bw() + theme(legend.position = "top") +
    facet_grid(~ facet_group, switch = "x") +    # strip을 x축 밑으로
    theme(strip.placement = "outside",           # strip 위치 조정
          strip.background = element_blank())

  
#-------------------------- p3, p4 pct_MT proportion ----------------------------#  



  
  # (B) p3 본체: violin + 평균점/박스플롯 + 상관계수 텍스트
  p3 <- df2 %>%
    ggplot(aes(x = facet_group, y = pct_counts_mt)) +
    geom_violin(trim = TRUE, width = 1, fill = "red", alpha = 0.55) +
    # 평균점(화이트 다이아몬드) 추가
    # 박스플롯 오버레이(윤곽선만)
    geom_boxplot(width = 0.18, outlier.shape = NA,
                 fill = "skyblue", color = "black", alpha = 0.8) +
    scale_y_continuous(labels = scales::comma) +
    labs(x = NULL, y = "pct_counts_mt") +
    theme_bw() +
    theme(legend.position = "none") 
    # 그룹별 상단에 Spearman ρ 라벨
    #geom_text(
    #  data = pval_df_p3,
    # aes(x = facet_group, y = Inf, label = label),
    #  vjust = 1.2, size = 4.2, inherit.aes = FALSE
    #)
  
# ------------------------------------------------- p4 -------------------------------#

  # (B) p4 본체
  p4 <- df2 %>%
    filter(cond == "Pass") %>%
    ggplot(aes(x = facet_group, y = pct_counts_mt, fill = cond)) +
    geom_violin(trim = TRUE, width = 1, alpha = 0.55, fill="red") +
    # 박스플롯 오버레이
    geom_boxplot(width = 0.18, outlier.shape = NA,
                 fill = "skyblue", color = "black", alpha = 0.8) +
    scale_y_continuous(labels = scales::comma) +
    scale_fill_manual(values = fill_vals) +
    labs(x = NULL, y = "pct_counts_mt", fill = "QC Status") +
    theme_bw() +
    theme(legend.position = "top") 

  

#--------------------------------------------------------------------------------#  
#------------------------- p5 n_genes_by_counts vs total counts -----------------#
  
  # 2) Pass만 점도 + 상단 좌측에 상관계수 라벨
  corr_df <- df2 %>%
    filter(cond == "Pass") %>%
    group_by(facet_group) %>%
    summarise(r = cor(total_counts, n_genes_by_counts,
                      method = "spearman", use = "complete.obs"),
              .groups = "drop") 
  
    # 2) 그래프는 Pass + Outlier 모두 표시, 라벨은 Pass로 계산한 값만 각 패널 상단에 표시
  p5 <- ggplot(df2, aes(x = total_counts, y = n_genes_by_counts, color = cond)) +
    geom_point(alpha = 0.5, size = 0.5) +
    scale_x_continuous(labels = scales::comma) +
    scale_y_continuous(labels = scales::comma) +
    scale_color_manual(values = col_vals) +
    coord_cartesian(xlim = lims$total_counts, ylim = lims$n_genes_by_counts) +
    labs(x = "total_counts", y = "n_genes_by_counts", color = "QC Status") +
    theme_bw() + theme(legend.position = "top") +
    facet_grid(~ facet_group, switch = "x") +
    theme(strip.placement = "outside",
          strip.background = element_blank()) +
    # ★ 패널별 상단에 Pass 상관계수 텍스트
    geom_text(
      data = corr_df,
      aes(x = -Inf, y = Inf, label = label),
      hjust = -0.1, vjust = 1.2,
      inherit.aes = FALSE, size = 4.2
    )
  
    
# ------------------------------------------------------------------------------#  
# -------------------- p6 total_counts vs mt_counts ----------------------------# 
    
  # 6) total_counts vs mt_counts — strip 제거, x축 그룹
  p6 <- ggplot(df2, aes(x = total_counts, y = mt_counts, color = cond)) +
    geom_point(alpha = 0.5, size = 1) +
    scale_x_continuous(labels = scales::comma) +
    scale_y_continuous(labels = scales::comma) +
    scale_color_manual(values = col_vals) +
    coord_cartesian(xlim = lims$total_counts, ylim = lims$mt_counts) +
    labs(x = "total_counts", y = "mt_counts", color = "QC Status") +
    theme_bw() + theme(legend.position = "top") +
    facet_grid(~ facet_group, switch = "x") +
    facet_wrap(~ facet_group, nrow = 2) +
    theme(strip.placement = "outside",
          strip.background = element_blank()) +
    stat_function(fun = function(x) 0.1 * x,
                  aes(color = "10%"), linetype = "solid") +
    stat_function(fun = function(x) 0.2 * x,
                  aes(color = "20%"), linetype = "solid") +
    stat_function(fun = function(x) 0.3 * x,
                  aes(color = "30%"), linetype = "solid") +
    scale_color_manual(
      values = c("Pass"="#9E9E9E","Outlier"="red",
                 "10%"="#140e16","20%"="#0072B2","30%"="green"),
      breaks = c("Pass","Outlier","10%","20%","30%"),
    )
#------------------------------------------------------------------------------#
  
  list(p1 = p1, p2 = p2, p3 = p3, p4 = p4, p5 = p5, p6 = p6, data_used = df2, limits = lims)
}


## 입력 인자 받기
# 폴더 경로 지정 (예: 현재 작업 디렉토리라면 ".")





path <- getwd()
folder_path = paste0(path,"/data/GSE")

# 또는 (권장: 운영체제 호환 안전)
df <- read.csv(file.path(folder_path, "ALL_files_QC_matrix_with_doublet.csv"))



library(stringr)

df2 <- df %>%
  mutate(
    group = case_when(
      str_detect(file, "Q") ~ "Quercetin",
      str_detect(file, "V") ~ "Vector",
      TRUE ~ "Other"
    )
)
df <- df2
rm(df2)

#데이터 분리
# 리스트로 나누기
df_split <- split(df, df$group)

# 리스트 이름 확인
names(df_split)

# set group number
group_A = "Quercetin"
group_B = "Vector"

# 1. 리스트로 받기
plots_final  <- qc_plots_AB_linear(
  df          = df,
  group_col   = "group",       # 비교에 쓸 그룹 컬럼명
  group_A     = group_A,      # A
  group_B     = group_B,       # B
  condition_col   = "outlier_with_doublet",   # QC 통과/실패가 들어있는 열
  pct_is_fraction = TRUE
)




figure_final <- (plots_final$p3 | plots_final$p4)
# 3. 화면에 출력
figure_final

out_prefix <- paste0("figures/", group_A, "vs", group_B)

#--------------------------------------plot 출력 및 저장 ------------------------------#



# 1) total_count total cout frequency
figure_12 <- (plots_final$p1 | plots_final$p2)
ggsave(
  file.path(folder_path, paste0(out_prefix, "total_count_and_filtered_total_count.jpg")),
  figure_12,
  width  = 12,
  height = 6,
  dpi    = 300
)


# 2) pecentage of mt counts
figure_34 <- (plots_final$p3 | plots_final$p4)
ggsave(
  file.path(folder_path, paste0(out_prefix, "before_after_MT_proportion.jpg")),
  figure_34,
  width  = 12,
  height = 6,
  dpi    = 300
)

# 
figure_56 <- (plots_final$p5) | (plots_final$p6)
ggsave(
  file.path(folder_path, paste0(out_prefix, "total_countsvsn_genes_and_mtcounts_linear.jpg")),
  figure_56,
  width  = 12,
  height = 12,
  dpi    = 300
)

#------------------------------------------------------------------------------------------#
# 1) 방향 계산(중앙값 기준; 평균 쓰려면 mean으로 바꿔도 됨)
dir_tbl <- df %>%
  group_by(group) %>%
  summarise(med = median(pct_counts_mt, na.rm = TRUE), .groups = "drop")

g_hi   <- dir_tbl$group[which.max(dir_tbl$med)]
g_lo   <- dir_tbl$group[which.min(dir_tbl$med)]
ylab   <- max(df$pct_counts_mt, na.rm = TRUE)
lab_df <- data.frame(group1 = g_lo, group2 = g_hi,
                     y.position = ylab * 1.05,
                     label = paste0(g_hi, " > ", g_lo, " (by median)"))

# 2) 플롯 + 방향 라벨 수동 주석
conv1_dir <-
  ggplot(df, aes(x = group, y = pct_counts_mt)) +
  geom_violin(trim = TRUE, width = 1, fill = "red", alpha = 0.55) +
  geom_boxplot(width = 0.18, outlier.shape = NA,
               fill = "skyblue", color = "black", alpha = 0.8) +
  stat_summary(fun = mean,   geom = "point", shape = 16, size = 2.5) +  # 평균점
  stat_summary(fun = median, geom = "point", shape = 15, size = 3) +    # 중앙값점
  stat_compare_means(method = "wilcox.test",
                     label = "p.format",
                     comparisons = list(c("Quercetin","Vector"))) +
  stat_pvalue_manual(lab_df, label = "label", tip.length = 0) +
  scale_y_continuous(labels = scales::comma) +
  labs(x = NULL, y = "pct_counts_mt") +
  theme_bw() + theme(legend.position = "none")

conv1_dir
#---------------------------------------------------------------------------------------#
# 1) 방향 계산(중앙값 기준; 평균 쓰려면 mean으로 바꿔도 됨)
dir_tbl <- df %>%
  group_by(group) %>%
  summarise(med = median(pct_counts_mt, na.rm = TRUE), .groups = "drop")

g_hi   <- dir_tbl$group[which.max(dir_tbl$med)]
g_lo   <- dir_tbl$group[which.min(dir_tbl$med)]
ylab   <- max(df$pct_counts_mt, na.rm = TRUE)
lab_df <- data.frame(group1 = g_lo, group2 = g_hi,
                     y.position = ylab * 1.05,
                     label = paste0(g_hi, " > ", g_lo, " (by median)"))

# 2) 플롯 + 방향 라벨 수동 주석
conv2_dir <-
  ggplot(df, aes(x = group, y = pct_counts_mt)) +
  geom_violin(trim = TRUE, width = 1, fill = "red", alpha = 0.55) +
  geom_boxplot(width = 0.18, outlier.shape = NA,
               fill = "skyblue", color = "black", alpha = 0.8) +
  stat_summary(fun = mean,   geom = "point", shape = 16, size = 2.5) +  # 평균점
  stat_summary(fun = median, geom = "point", shape = 15, size = 3) +    # 중앙값점
  stat_compare_means(method = "wilcox.test",
                     label = "p.format",
                     comparisons = list(c("Quercetin","Vector"))) +
  stat_pvalue_manual(lab_df, label = "label", tip.length = 0) +
  scale_y_continuous(labels = scales::comma) +
  labs(x = NULL, y = "pct_counts_mt") +
  theme_bw() + theme(legend.position = "none")

conv2_dir
#---------------------------------------------------------------------------------------#
dir_tbl_2 <- df %>%
  filter(outlier_with_doublet == "FALSE") %>%
  group_by(group) %>%
  summarise(med = median(pct_counts_mt, na.rm = TRUE), .groups = "drop")

g_hi   <- dir_tbl$group[which.max(dir_tbl$med)]
g_lo   <- dir_tbl$group[which.min(dir_tbl$med)]
ylab   <- max(df$pct_counts_mt, na.rm = TRUE)
lab_df <- data.frame(group1 = g_lo, group2 = g_hi,
                     y.position = ylab * 1.05,
                     label = paste0(g_hi, " > ", g_lo, " (by median)"))

# 2) 플롯 + 방향 라벨 수동 주석
conv2 <- df %>%
  filter(outlier_with_doublet == "FALSE") %>%
  ggplot(aes(x = group, y = pct_counts_mt)) +
  geom_violin(trim = TRUE, width = 1, fill = "red", alpha = 0.55) +
  geom_boxplot(width = 0.18, outlier.shape = NA,
               fill = "skyblue", color = "black", alpha = 0.8) +
  stat_summary(fun = mean,   geom = "point", shape = 16, size = 2.5) +  # 평균점
  stat_summary(fun = median, geom = "point", shape = 15, size = 3) +    # 중앙값점
  stat_compare_means(method = "wilcox.test",
                     label = "p.format",
                     comparisons = list(c("Quercetin","Vector"))) +
  stat_pvalue_manual(lab_df, label = "label", tip.length = 0) +
  scale_y_continuous(labels = scales::comma) +
  labs(x = NULL, y = "pct_counts_mt") +
  theme_bw() + theme(legend.position = "none")

conv2


#------------------------------------------------------------------------------------------#
df$mt_counts <- df$pct_counts_mt * df$total_counts

# 2) 무한/NA 제외한 한계값 계산 (변수명: limits_ 로 충돌 회피)
limits_ <- list(
  total_counts = range(df$total_counts[is.finite(df$total_counts)], na.rm = TRUE),
  mt_counts    = range(df$mt_counts[is.finite(df$mt_counts)],    na.rm = TRUE)
)

# 3) 직접 df로 “직빵” 플롯 
conv <- ggplot(df, aes(x = total_counts, y = mt_counts, color = group)) +
  geom_point(alpha = 0.5, size = 0.5) +
  scale_x_continuous(labels = scales::comma) +
  scale_y_continuous(labels = scales::comma) +
  coord_cartesian(xlim = limits_$total_counts, ylim = limits_$mt_counts) +  # ← 여기!
  labs(x = "total_counts", y = "mt_counts", color = "group") +
  theme_bw() + theme(legend.position = "top") +
  theme(strip.placement = "outside", strip.background = element_blank()) +
  stat_function(fun = function(x) 0.1 * x, color = "black", linetype = "dashed") +
  stat_function(fun = function(x) 0.2 * x, color = "blue",  linetype = "dashed") +
  stat_function(fun = function(x) 0.3 * x, color = "red",   linetype = "dashed") +
  scale_color_manual(
    values = c(
      "Quercetin" = "red",       # Treatment A
      "Vector" = "gray" # Control
    )
  )

conv

# -----------------------------------------------------------------------------#
pval_df_1 <- df %>%
  summarise(
    test = list(wilcox.test(pct_counts_mt ~ group)),
    .groups = "drop"
  ) %>%
  rowwise() %>%
  mutate(
    p = test$p.value,
    signif = case_when(
      p > 0.05 ~ "ns",
      p <= 0.05 & p > 0.01 ~ "*",
      p <= 0.01 & p > 0.001 ~ "**",
      p <= 0.001 & p > 0.0001 ~ "***",
      p <= 0.0001 ~ "****"
    ),
    label = paste0("p = ", signif)
  ) %>%
  select(p, signif, label)


conv2 <- ggplot(df, aes(x = group, y = pct_counts_mt)) +
  geom_violin(trim = TRUE, width = 1, fill = "red", alpha = 0.55) +
  # 박스플롯 오버레이(윤곽선만)
  geom_boxplot(width = 0.18, outlier.shape = NA,
               fill = "skyblue", color = "black", alpha = 0.8) +
  scale_y_continuous(labels = scales::comma) +
  labs(x = NULL, y = "pct_counts_mt") +
  theme_bw() +
  theme(legend.position = "none") +
  stat_compare_means(method = "wilcox.test", 
                     label = "p.format",   # ns, *, **, ***
                     comparisons = list(c("Quercetin","Vector")))  # 그룹 이름 매칭

conv2
#------------------------------------------------------------------------------#


conv3 <- df %>%
  filter(outlier_with_doublet == "FALSE") %>%
  ggplot(aes(x = group, y = pct_counts_mt)) +
  geom_violin(trim = TRUE, width = 1, fill = "red", alpha = 0.55) +
  # 박스플롯 오버레이(윤곽선만)
  geom_boxplot(width = 0.18, outlier.shape = NA,
               fill = "skyblue", color = "black", alpha = 0.8) +
  scale_y_continuous(labels = scales::comma) +
  labs(x = NULL, y = "pct_counts_mt") +
  theme_bw() +
  theme(legend.position = "none") +
  stat_compare_means(method = "wilcox.test", 
                     label = "p.format",   # ns, *, **, ***
                     comparisons = list(c("Quercetin","Vector")))  # 그룹 이름 매칭


conv3


#------------------------------------------------------------------------------#

# "_outlier"로 끝나거나 "outlier"로 시작하는 컬럼 추출
outlier_cols <- grep("^outlier|_outlier$|_flag", colnames(df), value = TRUE)

# Age별 합산
df %>%
  group_by(Age) %>%
  summarise(across(all_of(outlier_cols), ~ sum(.x, na.rm = TRUE)))
.table(df$Age)
