# ================================================================
# Meta-Analysis Script (METAL)
# ================================================================

# 创建METAL配置文件
cat > meta_config.txt << 'EOF'
# Meta-analysis weighted by standard error
SCHEME STDERR
SEPARATOR TAB

# Track allele frequencies
AVERAGEFREQ ON
MINMAXFREQ ON

# Column mappings
MARKER SNP_PROTEIN_combo
WEIGHT NMISS
ALLELELABELS A1 A2
FREQLABEL freq
EFFECT beta
STDERR se
PVAL P
CUSTOMVARIABLE NMISS
LABEL NMISS as NMISS

# Input files (5 ancestries)
PROCESS White.txt
PROCESS Other.txt
PROCESS Asian.txt
PROCESS Black.txt
PROCESS Mixed.txt

# Output
OUTFILE meta_result_ .tbl
ANALYZE
QUIT
EOF

# 运行METAL
system("metal meta_config.txt")
