#!/usr/bin/env bash
set -euo pipefail

# 输入文件前缀（不带 .bed/.bim/.fam 后缀）

input_prefix="merged_filtered"

# 指定输出目录
output_dir="output/"
mkdir -p "$output_dir"
cd "$output_dir"
# 参数设置
K_MIN=8         # 最小 K
K_MAX=20        # 最大 K
THREADS=1       # 线程数
CV_FOLDS=10     # 交叉验证折数
BOOTSTRAPS=100  # bootstrap 次数
SEED=12345      # 随机种子（可选，保证可重复性）

# 逐 K 循环运行
for K in $(seq $K_MIN $K_MAX); do
  echo "=== Running ADMIXTURE for K=${K} ==="
  
  # 为每个 K 设置输出前缀
  out_prefix="${output_dir}/merged_filtered_K${K}"
  
  admixture \
    --cv=$CV_FOLDS \
    -j$THREADS \
    -B$BOOTSTRAPS \
    -s$SEED \
    "${input_prefix}.bed" $K \
    | tee "${out_prefix}.log"
done

echo "All done. 结果保存在："
echo "  - ${output_dir}/merged_filtered_K{${K_MIN}..${K_MAX}}.Q  （群体构成矩阵）"
echo "  - ${output_dir}/merged_filtered_K{${K_MIN}..${K_MAX}}.P  （等位基因频率矩阵）"
echo "  - ${output_dir}/merged_filtered_K{${K_MIN}..${K_MAX}}.log（运行日志，包括 CV error）"