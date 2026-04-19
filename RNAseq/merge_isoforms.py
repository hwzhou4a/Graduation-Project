import os
import glob
import pandas as pd

# ================= 配置区域 =================
# RSEM 结果主目录
INPUT_DIR = "D:/Graduation Project/RNASeq/genome_alignment/03_rsem_quantification"
# 输出目录
OUTPUT_DIR = "D:/Graduation Project/RNASeq/expression_matrices"

if not os.path.exists(OUTPUT_DIR):
    os.makedirs(OUTPUT_DIR)

print(f"正在处理 Isoform 结果文件...")

# 1. 查找所有 .isoforms.results 文件
search_pattern = os.path.join(INPUT_DIR, "*", "*.isoforms.results")
files = glob.glob(search_pattern)

if not files:
    print("错误: 未找到 .isoforms.results 文件！")
    exit(1)

# 2. 准备数据容器
counts_dict = {}
tpm_dict = {}
gene_map = None  # 用于存储 transcript_id 到 gene_id 的映射关系

# 3. 遍历文件
for f in files:
    # 获取样本名
    filename = os.path.basename(f)
    sample_name = filename.replace(".isoforms.results", "")
    print(f"读取样本: {sample_name}")

    # 读取数据 (Tab分隔, 索引设为 transcript_id)
    df = pd.read_csv(f, sep='\t', index_col='transcript_id')

    # --- 关键步骤：提取 Gene ID 映射 ---
    # 只需要从第一个文件提取一次映射关系，因为所有样本的 reference 是一样的
    if gene_map is None:
        gene_map = df[['gene_id']].copy()

    # --- 提取数据 ---
    counts_dict[sample_name] = df['expected_count']
    tpm_dict[sample_name] = df['TPM']

print("正在合并矩阵...")

# 4. 合并数据列 (axis=1)
counts_matrix = pd.DataFrame(counts_dict)
tpm_matrix = pd.DataFrame(tpm_dict)

# 5. --- 关键步骤：将 gene_id 加回到矩阵的第一列 ---
# 使用 merge 将 gene_map 合并进来
# left_index=True 表示利用行索引 (transcript_id) 进行对齐
final_counts = pd.merge(gene_map, counts_matrix, left_index=True, right_index=True)
final_tpm = pd.merge(gene_map, tpm_matrix, left_index=True, right_index=True)

# 填充可能出现的 NaN (通常不会有)
final_counts = final_counts.fillna(0)
final_tpm = final_tpm.fillna(0)

# 6. 保存结果
counts_out = os.path.join(OUTPUT_DIR, "isoform_counts_matrix_with_gene_info.csv")
tpm_out = os.path.join(OUTPUT_DIR, "isoform_tpm_matrix_with_gene_info.csv")

final_counts.to_csv(counts_out)
final_tpm.to_csv(tpm_out)

print("================ 完成 ================")
print(f"Isoform Count 矩阵已生成: {counts_out}")
print(f"Isoform TPM 矩阵已生成: {tpm_out}")
print(f"矩阵预览 (前几列):")
print(final_counts.iloc[:3, :3])  # 打印前3行前3列看看