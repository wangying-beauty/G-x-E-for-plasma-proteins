# -------- PROCESS中介分析 --------
# Model 4: 简单中介模型
# 路径：X (env1) → M (env2) → Y (蛋白质)

p1 <- PROCESS(
  data = merged_data,
  y = PHE,              # 因变量（蛋白质）
  x = envv,             # 自变量（env1）
  meds = med,           # 中介变量（env2）
  covs = covariates,    # 协变量（PC1-40, age, sex等）
  cov.path = 'a',       # 协变量作用于 X→M 路径（a路径）
  ci = "boot",          # Bootstrap置信区间
  nsim = 1000,          # Bootstrap重复1000次
  seed = 1234           # 随机种子（保证可重复）
)

# 显著性判断：
# - 校正阈值: 0.05/25 = 0.002