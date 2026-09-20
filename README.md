# LATS2 Virtual Screening Pipeline (复现版)

作者：徐叶繁（苏州大学 药学专硕）

## 项目说明
本科毕业设计"基于计算机辅助药物设计的 LATS2 小分子抑制剂发现"中，
虚拟筛选流程的 Python 复现版本。原始流程在 Discovery Studio 中完成，
本项目使用 RDKit 实现相同逻辑，便于复用与扩展。

## 流程
原始库（约170万化合物，ChemDiv + Specs）
→ 标准化预处理（解析、去重）
→ 类药性五规则过滤（Lipinski：MW<=500, HBD<=5, HBA<=10, LogP<=5）
→ 分子描述符计算（MW, LogP, HBD, HBA, RotB, TPSA）
→ 输出候选列表（供分子对接）

## 文件
- `filter.py`：主脚本（约60行）
- `test.csv`：测试数据（5个分子：3个应通过Lipinski，2个应被过滤）
  - 通过：aspirin, ibuprofen, naphthalene
  - 过滤：C30-alkane（LogP>5，疏水性过高）、polyalanine-10（MW>500）
- `filtered_output.csv`：运行输出示例（通过3个分子）

## 环境与依赖
- Python 3.10
- RDKit（化学信息学核心库：SMILES解析、描述符计算、规则过滤）
- pandas（表格读写）
- 运行平台：百度AI Studio Notebook（云端，免安装）

## 运行方法
    pip install rdkit pandas
    python filter.py test.csv

## 已知局限与后续计划
- [ ] 盐处理：当前版本对含盐分子未做片段剥离，计划用 GetMolFrags 取最大片段
- [ ] 接入 PAINS / 反应性基团过滤（RDKit FilterCatalog）
- [ ] RDKit 指纹相似性聚类，挑选结构多样性候选
- [ ] AutoDock Vina 批量对接
