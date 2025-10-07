suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(scales)
  library(stringr)
  library(scran)
  library(patchwork)
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
  col_vals  <- c("Pass" = "#9E9E9E", "Outlier" = "#fd6e1e")
  fill_vals <- c("Pass" = "#BDBDBD", "Outlier" = "#fd6e1e")
  
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
  
  p2 <- df2 %>%
    dplyr::filter(cond == "Pass") %>%
    ggplot(aes(x = facet_group, y = pct_counts_mt, fill = cond)) +
    geom_violin(trim = TRUE, width = 1, alpha = 0.55) +
    scale_y_continuous(labels = scales::comma) +
    scale_fill_manual(values = fill_vals) +
    labs(x = NULL, y = "pct_counts_mt", fill = "QC Status") +
    theme_bw() +
    theme(legend.position = "top")
  
  # 3) total_counts 히스토그램 — strip 제거, x축 그룹
  p3 <- ggplot(df2, aes(x = total_counts, fill = cond)) +
    geom_histogram(bins = 100, color = "black", alpha = 0.45, position = "identity") +
    scale_x_continuous(labels = scales::comma) +
    scale_fill_manual(values = fill_vals) +
    coord_cartesian(xlim = lims$total_counts) +
    labs(x = "total_counts", y = "Frequency", fill = "QC Status") +
    theme_bw() + theme(legend.position = "top") +
    facet_grid(~ facet_group, switch = "x") +    # strip을 x축 밑으로
    theme(strip.placement = "outside",           # strip 위치 조정
          strip.background = element_blank())
  
  corr_df <- df2 %>%
    filter(cond == "Pass") %>%
    group_by(facet_group) %>%
    summarise(r = cor(total_counts, n_genes_by_counts,
                      method = "spearman", use = "complete.obs"),
              .groups = "drop") %>%
    mutate(label = sprintf("ρ = %.2f", r))
  
  # 2) Pass만 점도 + 상단 좌측에 상관계수 라벨
  corr_df <- df2 %>%
    filter(cond == "Pass") %>%
    group_by(facet_group) %>%
    summarise(r = cor(total_counts, n_genes_by_counts,
                      method = "spearman", use = "complete.obs"),
              .groups = "drop") %>%
    mutate(label = sprintf("ρ(Pass) = %.2f", r))
  
  # 2) 그래프는 Pass + Outlier 모두 표시, 라벨은 Pass로 계산한 값만 각 패널 상단에 표시
  p4 <- ggplot(df2, aes(x = total_counts, y = n_genes_by_counts, color = cond)) +
    geom_point(alpha = 0.5, size = 0.6) +
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
  
  # 5) total_counts vs mt_counts — strip 제거, x축 그룹
  p5 <- ggplot(df2, aes(x = total_counts, y = mt_counts, color = cond)) +
    geom_point(alpha = 0.5, size = 0.6) +
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
      values = c("Pass"="#9E9E9E","Outlier"="#fd6e1e",
                 "10%"="#140e16","20%"="#0072B2","30%"="red"),
      breaks = c("Pass","Outlier","10%","20%","30%"),
    )
  
  list(p1 = p1, p2 = p2, p3 = p3, p4 = p4, p5 = p5, data_used = df2, limits = lims)
}


## 입력 인자 받기
# 폴더 경로 지정 (예: 현재 작업 디렉토리라면 ".")
folder_path <- "scRNA-analysis/"


# 또는 (권장: 운영체제 호환 안전)
df <- read.csv(file.path(folder_path, "ALL_files_QC_matrix_with_doublet.csv"))

group_mapping <- c(
  "M1" = "Old",
  "M2" = "Metformin",
  "M3" = "Rapamycin",
  "M4" = "Young"
)

df <- df %>%
  mutate(
    sample_group = str_extract(sample_tag, "M[0-9]+"),
    group = group_mapping[sample_group]
  ) 


#데이터 분리
# 리스트로 나누기
df_split <- split(df, df$group)

# 리스트 이름 확인
names(df_split)

# set group number
group_A = "Old"
group_B = "Young"

# 1. 리스트로 받기
plots_final  <- qc_plots_AB_linear(
  df          = df,
  group_col   = "group",       # 비교에 쓸 그룹 컬럼명
  group_A     = group_A,      # A
  group_B     = group_B,       # B
  condition_col   = "outlier_with_doublet",   # QC 통과/실패가 들어있는 열
  pct_is_fraction = TRUE
)




figure_final <- (plots_final$p1 | plots_final$p2) /
                (plots_final$p3 | plots_final$p4) /
                plots_final$p5

# 3. 화면에 출력
figure_final

out_prefix <- paste0("figures/", group_A, "vs", group_B)


# 1+2
figure_1234 <- (plots_final$p1 | plots_final$p2) /
                  plots_final$p3 | plots_final$p4
ggsave(
  file.path(folder_path, paste0(out_prefix, "_qc_panel_1234.jpg")),
  figure_1234,
  width  = 12,
  height = 6,
  dpi    = 300
)

figure_5 <- (plots_final$p5)
ggsave(
  file.path(folder_path, paste0(out_prefix, "_qc_panel_5.jpg")),
  figure_5,
  width  = 12,
  height = 12,
  dpi    = 300
)



