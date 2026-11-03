# NCBI 注释标准化与 WGDI 预处理

本目录使用 `scripts/preprocess_wgdi_v2.py`，把 NCBI 的 GFF3、CDS FASTA 和蛋白 FASTA
转换为带全局唯一 ID、可直接进入 BLAST/WGDI 的一基因一代表转录本数据。

原始下载目录不会被修改：

```text
data/downloads/
resources/download_inventory.tsv
```

## 运行

```powershell
cd C:\Users\27563\Desktop\HXC-Project\plant-genome-resources

# 测试
python -m pytest -q

# 单物种试跑
python scripts/preprocess_wgdi_v2.py --species 996

# 全部物种
python scripts/preprocess_wgdi_v2.py

# 只处理部分物种
python scripts/preprocess_wgdi_v2.py --species 5,225,996

# 严格验收当前 manifest
python scripts/validate_wgdi_dataset.py `
  --project-root C:/Users/27563/Desktop/HXC-Project `
  --expected-count 20
```

可选参数：

```text
--min-genes-per-sequence 50   # 碎片组装中一个序列至少保留的基因数
--coverage 0.90               # 碎片组装兜底保留的基因覆盖率
--max-sequences 200           # WGDI 分析集最大序列数
--no-overwrite                # 已存在输出时拒绝覆盖
```

## 核心规则

1. GFF3 中的 gene、mRNA、CDS 和 protein_id 建立映射。
2. 一个基因只选一个代表转录本：
   - 先选蛋白最长；
   - 蛋白等长时选 CDS 最长；
   - 仍相同时按转录本 ID 排序。
3. 只接受 CDS 与 PEP 能实际翻译对齐的候选。
4. 拒绝以下数据：
   - CDS 长度不是 3 的倍数；
   - 内部终止密码子；
   - 氨基酸非法字符；
   - 翻译长度或内容不一致；
   - 没有 CDS/PEP 的非编码转录本；
   - 细胞器基因；
   - 跨核基因组 taxon 污染；
   - 同一蛋白 ID 被不同基因复用；
   - 重复且序列冲突的 FASTA ID。
5. 原始 ID、是否代表转录本、翻译状态和拒绝原因都保留在映射/报告文件中。

## 输出

```text
clean_v2/{code}/
├── {code}.gene_map.tsv
├── {code}.chr_map.tsv
└── {code}.rejected.tsv

analysis_v2/{code}/
├── {code}.cds                 # WGDI 筛选集 CDS
├── {code}.pep                 # WGDI 筛选集蛋白
├── {code}.new.gff             # WGDI 七列位置文件
├── {code}.lens                # WGDI 染色体长度/基因数
├── {code}.all.cds.fa          # 全部合格代表基因
├── {code}.all.pep.fa
├── {code}.all.gff
├── {code}.all.lens
├── {code}.qc.json
├── {code}.gene_map.tsv
└── {code}.chr_map.tsv
```

`{code}.new.gff` 为七列：

```text
简化序列号  唯一基因ID  起点  终点  链方向  基因序号  原始基因ID
```

WGDI 使用前六列，第七列用于回溯 NCBI 原始基因名。

## 汇总文件

```text
resources/stage_preprocess_v2_summary.tsv
resources/wgdi_v2_index.tsv
resources/species_registry_v2.tsv
resources/wgdi_v2_manifest_sha256.tsv
```

## 项目目录

```text
plant-genome-resources/
├── data/downloads/      # 原始下载，只读
├── resources/           # manifest、QC、哈希
├── clean_v2/            # 映射和拒绝原因
├── analysis_v2/         # WGDI 四件套
├── results/
│   ├── blast/
│   ├── wgdi/
│   ├── synteny/
│   ├── tree/
│   └── wgd/
├── logs/
├── runs/
├── scripts/
├── src/
└── tests/
```

目录决策见 `docs/adr/0001-data-layout-and-validation-contract.md`。

## WGDI 环境

本机 `bioinfo_win` 的 WGDI 0.75 在 matplotlib 绘图层会原生崩溃。已验证
`pced` 环境可以完成 WGDI dotplot：

```powershell
D:\Anaconda3\envs\pced\python.exe -c "import sys; sys.argv=['wgdi','-d','total.conf']; from wgdi.run import main; main()"
```

或使用：

```powershell
D:\Anaconda3\envs\pced\Scripts\wgdi.exe -d total.conf
```

WGDI 配置中的 CDS、PEP、GFF 和 lens 应使用同一套筛选集：

```text
{code}.cds
{code}.pep
{code}.new.gff
{code}.lens
```

