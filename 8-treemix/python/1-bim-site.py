import os
import shutil
import pandas as pd

# ——配置区——
# 原始 bim 文件路径
orig_path = "/mnt/d/幽门螺旋杆菌/Script/分析结果/4-treemix/data/3219_maf99/filter_filtered.bim"
# 备份目录（可以自定义）
bak_dir = os.path.join(os.path.dirname(orig_path), "bak")
# ——————

# 1. 创建备份目录（如果不存在）
os.makedirs(bak_dir, exist_ok=True)

# 2. 备份原始文件到 bak 目录
bak_path = os.path.join(bak_dir, os.path.basename(orig_path))
shutil.copy2(orig_path, bak_path)
print(f"已备份原始文件到：{bak_path}")

# 3. 读取原始文件（空白分隔、无表头）
df = pd.read_csv(orig_path, sep=r"\s+", header=None)

# 4. 将第 2 列（索引1）替换为第 4 列（索引3）
df[1] = df[3]

# 5. 原地写回——覆盖原始文件，不输出行号和表头
df.to_csv(orig_path, sep="\t", header=False, index=False)
print(f"已修改并覆盖原始文件：{orig_path}")
