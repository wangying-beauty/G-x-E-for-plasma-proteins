# ================================================================
# Conditional Independence Analysis - Complete Pipeline
# ================================================================

# -------- 核心函数：逐步条件分析 --------

test_gxe <- function(exposure, bm, SNP, df, adj_exposures = c()) {
  bm_adj <- paste0(bm, "_adj")
  g <- df[[SNP]]
  e <- df[[exposure]]
  
  # 构建线性模型
  lm_str <- paste0(bm_adj, " ~ g * e")
  if (length(adj_exposures) > 0) {
    lm_str <- paste(lm_str, "+", paste0("g * ", adj_exposures, collapse = " + "))
  }
  
  # 提取交互项P值
  pval <- tryCatch({
    lm_res <- broom::tidy(lm(as.formula(lm_str), data = df))
    filter(lm_res, term == "g:e")$p.value
  }, error = function(e) as.numeric(NA))
  
  if (length(pval) == 1) pval else as.numeric(NA)
}

run_conditional_tests <- function(all_exposures, bm, SNP, p_threshold) {
  remaining_exposures <- all_exposures
  adj_exposures <- c()
  
  while (length(remaining_exposures) > 0) {
    # 测试剩余环境因子的条件独立性
    pvals <- map_dbl(remaining_exposures, test_gxe, bm, SNP, ca_df, adj_exposures)
    
    # 如果有显著信号，加入调整集
    if (any(pvals < p_threshold, na.rm = TRUE)) {
      adj_exposures <- c(adj_exposures, remaining_exposures[which.min(pvals)])
    }
    
    # 移除非显著的环境因子
    nonsig_exposures <- remaining_exposures[which(
      (pvals > p_threshold) | is.na(pvals)
    )]
    
    remaining_exposures <- setdiff(remaining_exposures, 
                                   c(adj_exposures, nonsig_exposures))
  }
  
  if (length(adj_exposures) > 0) adj_exposures else "none"
}

# -------- 运行条件分析 --------

sig_threshold <- 0.05 / 536 / 1357  # Bonferroni校正

indep_res <- sig_combos %>%
  rowwise() %>%
  mutate(
    independent = list(run_conditional_tests(
      str_split(sig_exposures, ",")[[1]], bm, SNP, sig_threshold
    )),
    independent_nom = list(run_conditional_tests(
      str_split(sig_exposures, ",")[[1]], bm, SNP, 0.05
    ))
  ) %>%
  ungroup() %>%
  mutate(
    n_sig = lengths(independent),
    n_sig_nom = lengths(independent_nom)
  )

saveRDS(indep_res, "independent_hits.rds")

# -------- 导出结果 --------

indep_res %>%
  mutate(
    indep_exposures = map_chr(independent, paste, collapse = "|"),
    indep_exposures_nom = map_chr(independent_nom, paste, collapse = "|")
  ) %>%
  select(bm, SNP, sig_exposures, 
         indep_exposures, n_sig, 
         indep_exposures_nom, n_sig_nom) %>%
  write_csv("independent_hits_results.csv")