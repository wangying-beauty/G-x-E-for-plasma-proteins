#!/bin/bash

# ================================================================
# Core Logic for SNP × Environment Interaction Analysis
# ================================================================

# 基本参数设置
SNP_PROTEIN_PAIRS="list_snp_protein.txt"
ENV_DATA="environment_data.csv"
GENOTYPE_FILE="ukb_qc_merged"
P_THRESHOLD=6.87e-08

# 提取所有环境因子变量名 (跳过第一列和第二列的 FID, IID)
ENV_VARS=($(head -1 "$ENV_DATA" | tr ',' '\n' | tail -n +3))

# 外层循环：逐行读取需要分析的 SNP 和 对应的蛋白
while IFS=$'\t' read -r snp protein; do
    [ "$snp" = "SNP" ] && continue # 跳过表头
    
    pheno_file="res_${protein}.txt"
    pheno_name="Res${protein}"
    safe_snp=$(echo "$snp" | tr ':/' '__') # 将 SNP 名称中的特殊字符替换以用作文件名
    
    # 内层循环：遍历所有环境因子
    for env_var in "${ENV_VARS[@]}"; do
        out_prefix="inter_${protein}_${env_var}_${safe_snp}"
        
        # 步骤 1: 为当前环境因子单独提取协变量文件 (PLINK 需要特定的格式)
        awk -v env="$env_var" 'BEGIN {FS=","; OFS="\t"}
            # 找到目标环境因子所在的列
            NR==1 { for(i=1;i<=NF;i++) if($i==env) col=i; print "FID", "IID", env }
            # 提取数据并忽略缺失值
            NR>1 && $col!="" { print $1, $2, $col }
        ' "$ENV_DATA" > "${out_prefix}_cov.txt"
        
        # 步骤 2: 运行 PLINK2 核心交互项模型 (GxE)
        plink2 \
            --bfile "$GENOTYPE_FILE" \
            --pheno "$pheno_file" \
            --pheno-name "$pheno_name" \
            --snp "$snp" \
            --covar "${out_prefix}_cov.txt" \
            --covar-name "$env_var" \
            --parameters 1-3 \
            --variance-standardize \
            --glm interaction \
            --out "$out_prefix"
            
        # 步骤 3: 提取显著的交互项结果
        result_file="${out_prefix}.${pheno_name}.glm.linear"
        if [ -f "$result_file" ]; then
            # 提取交互项结果
            sed -n '1~3p' "$result_file" | \
            awk -v p="$P_THRESHOLD" 'BEGIN{FS=OFS="\t"} NR==1 || ($15 != "NA" && $15 < p)' \
            > "${out_prefix}_significant.txt"
            
            # 删除原始巨大的输出文件以节省空间
            rm -f "$result_file"
        fi
        
        # 清理临时的协变量文件
        rm -f "${out_prefix}_cov.txt"
        
    done
done < "$SNP_PROTEIN_PAIRS"

echo "Core interaction analysis completed."