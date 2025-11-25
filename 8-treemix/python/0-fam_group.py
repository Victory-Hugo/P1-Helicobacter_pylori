import pandas as pd
import shutil
import sys
import os

def main(pop_name, id_name, meta_file, fam_file, sheet_name=None):
    # 读取meta信息
    if meta_file.endswith('.xlsx') or meta_file.endswith('.xls'):
        df_meta = pd.read_excel(meta_file, sheet_name=sheet_name or 0).loc[:, [pop_name, id_name]]
    else:
        df_meta = pd.read_csv(meta_file).loc[:, [pop_name, id_name]]
    # ID转为字符串
    df_meta[id_name] = df_meta[id_name].astype(str)

    # 读取fam文件
    df_fam = pd.read_csv(fam_file, sep=' ', header=None, names=['ID1', 'ID2', 'P1', 'P2', 'P3', 'P4'])
    df_fam_merge = df_fam.merge(df_meta, left_on='ID1', right_on=id_name, how='left')
    df_fam_merge['ID1'] = df_fam_merge[pop_name]
    out_file = f"{fam_file}_updated"
    df_fam_merge.iloc[:, 0:6].to_csv(out_file, sep=' ', header=False, index=False)

    # 备份原始fam文件
    shutil.copy(fam_file, f"{fam_file}.bak")
    # 替换原始fam文件
    shutil.move(out_file, fam_file)

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description="Update fam file with new ID info from meta file.")
    parser.add_argument('--pop_name', required=True, help='列名：新ID')
    parser.add_argument('--id_name', required=True, help='列名：原ID')
    parser.add_argument('--meta_file', required=True, help='meta文件路径（支持csv/xlsx）')
    parser.add_argument('--fam_file', required=True, help='fam文件路径')
    parser.add_argument('--sheet_name', default=None, help='meta为excel时的sheet名，可选')
    args = parser.parse_args()
    main(args.pop_name, args.id_name, args.meta_file, args.fam_file, args.sheet_name)
