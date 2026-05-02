library(tidyverse)
library(data.table)

# -------- 1. 6MAD离群值处理 --------

# 公式: Outlier if |x - median| > 6 × MAD
# MAD = 1.4826 × median(|x - median|)

remove_outliers_6mad <- function(x) {
  med <- median(x, na.rm = TRUE)
  mad_val <- mad(x, na.rm = TRUE, constant = 1.4826)
  x[abs(x - med) > 6 * mad_val] <- NA
  return(x)
}



# -------- 2. PHESANT处理 --------

# 准备变量信息文件
variable_info <- tibble(
  FieldID = names(cleaned_data)[-1],  # 排除ID列
  Path = "Continuous",  # 根据实际调整
  TRAIT_OF_INTEREST = "YES"
)
fwrite(variable_info, "ewis_variable_info.tsv", sep = "\t")

# 运行PHESANT 
system("
Rscript phenomeScan.r \\
  --phenofile='input_phenos_536.csv' \\
  --variablelistfile='ewis_variable_info.tsv' \\
  --datacodingfile='data-coding-ordinal-info-nov2019-update.txt' \\
  --resDir='phesant_output/' \\
  --userId='id' \\
  --partIdx=1 \\
  --numParts=2 \\
  --out='ewis_phenotypes536'
")

# -------- 3. 零方差检验 --------

phesant_data <- fread("phesant_combined_536.csv")

# 计算方差并移除零方差变量
variance_check <- phesant_data %>%
  select(-userId) %>%
  summarise(across(everything(), ~ var(.x, na.rm = TRUE))) %>%
  pivot_longer(everything(), names_to = "variable", values_to = "variance")

zero_var_vars <- variance_check %>%
  filter(variance == 0 | is.na(variance)) %>%
  pull(variable)

final_data <- phesant_data %>%
  select(-all_of(zero_var_vars))

fwrite(final_data, "final_exposures_536.csv")