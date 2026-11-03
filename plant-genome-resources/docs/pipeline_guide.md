# 基因注释标准化流程说明（旧版 v1，已被 v2 取代）

> **状态：superseded。** 本文记录早期的 `wgdi_prep.py` 和 `clean/`、
> `analysis/` 结果。当前正式流程是 `scripts/preprocess_wgdi_v2.py`、
> `clean_v2/`、`analysis_v2/`，目录决策见
> `docs/adr/0001-data-layout-and-validation-contract.md`，使用说明见
> 项目根目录 `README.md`。

这份文档解释 `scripts/wgdi_prep.py` 在做什么、为什么这样定规则、出错怎么查。
所有结论都基于 20 个物种的真实数据跑出来的，不是推测。

---

## 0. 一句话总览

NCBI 给了三份"同一批基因的不同切面"：GFF 说基因在哪，CDS 给核酸，PEP 给蛋白。
把它们对齐成 **一基因一条代表转录本**，输出 BLASTP 和 WGDI 能直接吃的五件套。

```
data/downloads/  原始 gz（只读，永不修改）
      ↓  parse
clean/{code}/    全转录本映射表（忠实整理，不丢信息）
      ↓  build
analysis/{code}/ 一比一分析用文件（pep / cds / wgdi.gff / lens / chr_map）
      ↓  validate
校验报告         ID 集合一致性 + 坐标 + 顺序 + 翻译一致性
```

---

## 1. 为什么四方文件能对上：protein_id 是唯一钥匙

同一个基因在三份文件里的标识完全不同：

| 位置 | 实际写法 |
|---|---|
| GFF `gene` 行 | `ID=gene-LOC17897069`, `gene=LOC17897069` |
| GFF `mRNA` 行 | `ID=rna-XM_006303302.2`, `Parent=gene-LOC17897069` |
| GFF `CDS` 行 | `ID=cds-XP_006303364.2`, `protein_id=XP_006303364.2` |
| CDS FASTA 标题 | `>lcl\|NW_006238923.1_cds_XP_006303364.2_1 [gene=LOC17897069] [protein_id=XP_006303364.2]` |
| PEP FASTA 标题 | `>XP_006303364.2 protein mago nashi homolog [Capsella rubella]` |

能同时出现在三处的只有 **`protein_id`（XP_ 开头的蛋白号）**。
所以流水线的连接顺序是：

```
gene --(mRNA 的 Parent)--> mRNA --(CDS 行的 Parent)--> protein_id --> PEP 序列 / CDS 序列
```

一旦这条链走通，剩下的都是格式整理。走不通的地方必须显式记录，不能悄悄丢掉。

**怎么验证配对真的对了？** 用翻译一致性：
正确配对的 CDS 与蛋白应满足 `len(CDS)/3 - 1 == len(蛋白)`（减掉终止密码子）。
实测 20 个物种该比例都在 87%–100%，说明配对无误。残留的 1%–13% 是 NCBI
注释本身的 alternative processing / partial model，不是流程错。

---

## 2. 细胞器基因组必须剔除（这是 20 物种里才暴露的问题）

NCBI 有些物种会把**叶绿体、线粒体**一起写进同一个 GFF。它们从 `region` 行的
`genome=` 字段可以识别：

```
genome=chromosome / genomic      核基因组
genome=chloroplast               叶绿体
genome=mitochondrion             线粒体
```

实测 6 个物种带细胞器：

| 物种 | 细胞器基因数 |
|---|---:|
| Beta vulgaris | 148 |
| Andrographis paniculata | 87 |
| Ipomoea triloba | 85 |
| Solanum stenotomum | 85 |
| Cryptomeria japonica | 82 |
| Phalaenopsis equestris | 73 |
| Cucurbita pepo | 38 |
| Monoraphidium neglectum | 21 |

**必须剔除的三个理由**：

1. 叶绿体基因名会重名。Beta vulgaris 的 `rps12` 在叶绿体上有两个拷贝，
   直接拼 `物种.rps12` 会撞 ID，导致 WGDI 输入里出现重复基因（第一版就踩到了）。
2. 细胞器蛋白是原核起源，BLASTP 会跟核基因产生大量假阳性 hit，污染共线性。
3. 细胞器是环状基因组，没有"染色体顺序"概念，放进 WGDI 无意义。

处理方式：gate 在 `compartment` 字段，analysis 层直接排除，clean 层完整保留
（要研究细胞器随时能取）。

---

## 3. 代表转录本怎么选

一个基因常有多个转录本（isoform），WGDI/BLASTP 要求一基因一条。规则按优先级：

1. 蛋白最长（去终止符）—— 最可能完整的编码产物
2. 蛋白等长时，CDS 最长
3. 仍相同，按转录本 ID 字典序最小

第 3 条是为了**可复现**：不依赖文件顺序、不依赖 Python 字典顺序，重跑结果一样。

不进入代表集的：非编码 RNA（lncRNA/tRNA/rRNA/snRNA/snoRNA）、假基因、
以及有坐标但拿不到蛋白序列的转录本。它们保留在 clean 层，不进 analysis。

实测代表集覆盖率（蛋白编码基因中有代表转录本的比例）20 个物种全部 ≥ 99.9%，
基本等于 100%，说明 NCBI 注释质量不错。

---

## 4. 染色体筛选：50 基因阈值不够用

第一版按"至少 50 个基因才进 WGDI 集"过滤，结果暴露问题：

- Monoraphidium neglectum：6287 条序列（16734 个基因，平均 2.7 个基因/scaffold），一条都不满 50 个基因 → WGDI 集是空的
- Prosopis alba：4973 条序列，按阈值只剩 20.7% 的基因
- Phalaenopsis equestris：4363 条序列，只剩 41.9%

所以改成**两级自适应**：

```
第一级 threshold：保留基因数 >= 50 的序列
   → 若这样能覆盖 >= 50% 的基因，就用它（说明组装基本完整）

第二级 coverage：否则说明组装碎，按基因数从多到少累加，
   直到覆盖 90% 的基因 / 最多 200 条序列
```

同时给每个物种一个 `wgdi_ready` 判定：

| 值 | 含义 | 建议 |
|---|---|---|
| yes | 保留 ≥70% 基因且序列数 ≤200 | 直接进 WGDI |
| marginal | 保留 50%–70% | 能做，但共线性图会偏碎，需人工看 |
| no | 保留 <50% | 不适合染色体级共线性，考虑换组装版本或只做局部 |

实测 20 个物种：**17 个 yes、1 个 marginal（Phalaenopsis）、2 个 no（Monoraphidium、Prosopis）**。

---

## 5. 输出文件格式

### `{code}.wgdi.gff` — WGDI 六列

```
s996.NW_006238923.1   s996.LOC17897069   1864   3158   -   1
```
染色体 / 代表基因 ID / 起点 / 终点 / 链方向 / 该染色体上第几个基因（1 起，按起点排序连续编号）

注意：**不是 GFF3 九列**，没有 feature 类型、score、phase、attribute。
直接用 NCBI GFF3 改名喂给 WGDI 一定失败。

### `{code}.lens` — 三列

```
s996.NW_006238923.1   19566579   4325
```
染色体 / 长度 / 该染色体上的基因数。长度取自 GFF `##sequence-region`（最权威）。

### ID 命名方案

```
物种代码 s{show_id}         如 s996
基因     s996.LOC17897069   物种前缀 + 原基因名（稳定、可回溯）
转录本   s996.XM_006303302.2
染色体   s996.NW_006238923.1
```

带物种前缀是为了 20 个物种混在一起做 BLAST 时不会撞 ID。
原始 ID 全部保留在 `gene_map.tsv` 里，随时能反查。

---

## 6. 校验查了什么

`validate` 是硬门禁，不通过不进入 BLAST：

1. gff 列数 = 6
2. gff 里每条染色体都出现在 lens 中
3. 基因终点不超过染色体长度
4. start ≤ end
5. 链方向只能是 `+ - .`
6. 每条染色体的 `order` 是 1..N 连续无跳号
7. pep / cds / gff 三者 ID 集合完全相等
8. 三个文件内部各自无重复 ID
9. CDS 长度与蛋白长度一致率 ≥ 85%

实测 20 个物种全部 0 problem。

---

## 7. 怎么用

```bash
cd plant-genome-resources

# 扫描物种，生成清单
python scripts/wgdi_prep.py discover

# 单物种全流程（建议新物种先这样试）
python scripts/wgdi_prep.py all --species 996

# 全量 20 物种
python scripts/wgdi_prep.py batch

# 调整严格度
python scripts/wgdi_prep.py batch --min-genes-per-chr 100 --coverage 0.95 --max-chr 100
```

产出汇总：

| 文件 | 内容 |
|---|---|
| `resources/species_registry.tsv` | 20 物种清单 + accession + 文件名 |
| `resources/stage3_preprocess_summary.tsv` | 全量 QC 总表（每物种一行，30+ 指标） |
| `resources/wgdi_file_index.tsv` | 每物种 6 个文件的相对路径 + 是否可用 |
| `resources/stage3_analysis_manifest_sha256.tsv` | 所有输出文件的 SHA-256 |
| `clean/{code}/{code}.gene_map.tsv` | 全转录本映射表（34126 行级别） |
| `clean/{code}/{code}.chr_map.tsv` | 染色体表（含 genome_type / is_organelle） |
| `analysis/{code}/*` | BLAST/WGDI 输入五件套 |

---

## 8. 两个刻意的设计选择

**clean 层不复制序列。** 原始 CDS/PEP 是全转录本的，clean 层如果也存一份
就只是 4 GB 的纯重复。所以 clean 层只放**映射表**（谁对应谁），序列仍然
从 `data/downloads` 的 gz 直接读。原始文件永远是唯一真源。

**原始文件只读。** 整个流程没有一次写入 `data/downloads`。标准化结果全部
落到 `clean/` 和 `analysis/`，两者都可以随时删掉重跑。

---

## 9. 下一步（还没做）

1. **BLASTP**：环境里有 `blastp`、`makeblastdb`、`diamond`（`D:/ProgramFiles/bioinfo/bin`）。
   20 物种两两比对共 190 对，建议先用 DIAMOND 建库跑 `--very-sensitive`，
   再把结果转成 WGDI 需要的 12 列格式。
2. **total.conf**：没有自动生成。WGDI 各版本的 conf section 名称不完全一致，
   猜错会让整条流程静默失败，不如按你本机 WGDI 版本的手册填一次模板。
   `wgdi_file_index.tsv` 已经把 6 个文件的路径排好了，直接往里贴即可。
3. **三个不适合的物种**：Monoraphidium（藻类，2.7 基因/scaffold）、Prosopis、
   Phalaenopsis 若要做共线性，建议换更高连续度的组装版本，或改用基因家族层面分析。
