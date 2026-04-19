import os
import glob
import pandas as pd

# ================= 配置区域 =================
# RSEM 结果主目录 (根据你的路径设置)
INPUT_DIR = "D:/Graduation Project/RNASeq/genome_alignment/03_rsem_quantification"
# 输出目录
OUTPUT_DIR = "D:/Graduation Project/RNASeq/expression_matrices"

if not os.path.exists(OUTPUT_DIR):
    os.makedirs(OUTPUT_DIR)

print(f"正在扫描目录: {INPUT_DIR}")

# 1. 查找所有 .genes.results 文件
search_pattern = os.path.join(INPUT_DIR, "*", "*.genes.results")
files = glob.glob(search_pattern)

if not files:
    print("错误: 未找到 .genes.results 文件！请检查路径。")
    exit(1)

print(f"找到 {len(files)} 个样本文件，开始合并...")

# 2. 准备列表存储数据
counts_list = []
tpm_list = []

for f in files:
    # --- 提取样本名称 ---
    # 文件名示例: E250098172.genes.results -> E250098172
    filename = os.path.basename(f)
    sample_name = filename.replace(".genes.results", "")

    print(f"读取样本: {sample_name}")

    # --- 读取数据 ---
    # RSEM 的 gene 文件第一列就是 gene_id，直接用作索引
    df = pd.read_csv(f, sep='\t', index_col='gene_id')

    # --- 提取 Expected Count (用于差异分析) ---
    # rename 将 Series 的名字改成样本名，这样合并时列名就是样本名
    counts_list.append(df['expected_count'].rename(sample_name))

    # --- 提取 TPM (用于热图/PCA) ---
    tpm_list.append(df['TPM'].rename(sample_name))

# ================= 合并与保存 =================

print("正在合并矩阵...")

# 1. 生成 Gene Count Matrix
# axis=1 表示横向拼接（按列对齐）
gene_counts_matrix = pd.concat(counts_list, axis=1)
gene_counts_matrix = gene_counts_matrix.fillna(0)

# 保存 Count 矩阵
counts_output_file = os.path.join(OUTPUT_DIR, "gene_counts_matrix.csv")
gene_counts_matrix.to_csv(counts_output_file)

# 2. 生成 Gene TPM Matrix
gene_tpm_matrix = pd.concat(tpm_list, axis=1)
gene_tpm_matrix = gene_tpm_matrix.fillna(0)

# 保存 TPM 矩阵
tpm_output_file = os.path.join(OUTPUT_DIR, "gene_tpm_matrix.csv")
gene_tpm_matrix.to_csv(tpm_output_file)

print("================ 完成 ================")
print(f"Gene Count 矩阵已保存: {counts_output_file}")
print(f"Gene TPM 矩阵已保存:   {tpm_output_file}")
print("-" * 30)
print("矩阵预览 (Count 矩阵前5行):")
print(gene_counts_matrix.iloc[:5, :3])  # 打印前5行，前3个样本