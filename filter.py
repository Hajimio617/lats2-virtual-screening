# -*- coding: utf-8 -*-
"""
filter.py  —— LATS2 虚拟筛选流程复现（第 1 步：规则过滤 + 描述符计算）
作者：徐叶繁
功能：
  1. 读取含 SMILES 的 CSV 文件
  2. 清洗（解析失败剔除、去重、去盐）
  3. 类药性五规则（Lipinski）过滤
  4. 计算常用分子描述符，输出新 CSV（可用 Origin/Excel 打开做后续分析）
用法：
  python filter.py input.csv
要求 input.csv 中至少有一列叫 SMILES（列名大小写不敏感）
"""
import sys
import pandas as pd
from rdkit import Chem
from rdkit.Chem import Descriptors, Lipinski
from rdkit import RDLogger
RDLogger.DisableLog('rdApp.*')   # 关掉 rdkit 的刷屏警告


def mol_from_smiles(smiles):
    """把 SMILES 字符串转成 RDKit 分子对象，失败返回 None"""
    mol = Chem.MolFromSmiles(smiles)
    return mol


def lipinski_ok(mol):
    """Lipinski 五规则：MW<=500, HBD<=5, HBA<=10, LogP<=5（口服小分子常用过滤）"""
    if Descriptors.MolWt(mol) > 500:
        return False
    if Lipinski.NumHDonors(mol) > 5:
        return False
    if Lipinski.NumHAcceptors(mol) > 10:
        return False
    if Descriptors.MolLogP(mol) > 5:
        return False
    return True


def main():
    if len(sys.argv) < 2:
        print("用法：python filter.py input.csv")
        sys.exit(1)

    df = pd.read_csv(sys.argv[1])
    # 找到 SMILES 列（忽略大小写）
    smiles_col = None
    for c in df.columns:
        if c.strip().lower() == "smiles":
            smiles_col = c
            break
    if smiles_col is None:
        print("错误：CSV 里没找到 SMILES 列，请检查列名")
        sys.exit(1)

    print(f"原始分子数：{len(df)}")
    df = df.dropna(subset=[smiles_col]).drop_duplicates(subset=[smiles_col])
    print(f"去空去重后：{len(df)}")

    kept, records = [], []
    for _, row in df.iterrows():
        mol = mol_from_smiles(row[smiles_col])
        if mol is None:
            continue                                   # 解析失败，剔除
        if not lipinski_ok(mol):
            continue                                   # 不符合类药性规则
        kept.append(row)
        records.append({                             # 顺便算描述符
            "SMILES": row[smiles_col],
            "MW": round(Descriptors.MolWt(mol), 2),
            "LogP": round(Descriptors.MolLogP(mol), 2),
            "HBD": Lipinski.NumHDonors(mol),
            "HBA": Lipinski.NumHAcceptors(mol),
            "RotB": Lipinski.NumRotatableBonds(mol),
            "TPSA": round(Descriptors.TPSA(mol), 2),
        })

    out_df = pd.DataFrame(records)
    out_df.to_csv("filtered_output.csv", index=False)
    print(f"通过 Lipinski 过滤：{len(kept)} 个分子")
    print("结果已保存到 filtered_output.csv，可用 Excel/Origin 打开")


if __name__ == "__main__":
    main()
