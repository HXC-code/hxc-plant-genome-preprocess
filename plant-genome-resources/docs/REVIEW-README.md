# HXC Preprocess Review Package

## 审阅重点

请重点审阅以下四件事：

1. `0check_data.sh` 是否应只作为 Bash 调度/门禁层。
2. `wgdi_v2_index.tsv` 作为唯一 manifest 是否合理。
3. 一基因一代表转录本、翻译一致性和细胞器隔离规则是否合理。
4. 目录是否应保持 `data/downloads`、`clean_v2`、`analysis_v2`、
   `results` 分层。

## 包内内容

```text
0check_data.sh
docs/data-validation.md
plant-genome-resources/
├── README.md
├── pyproject.toml
├── scripts/
│   ├── preprocess_wgdi_v2.py
│   └── validate_wgdi_dataset.py
├── src/pgci_preprocess/
├── tests/
├── docs/adr/
└── docs/REVIEW-README.md
```

包内不包含：

- `data/downloads` 原始基因组数据。
- `analysis_v2` 和 `clean_v2` 派生数据。
- BLAST/WGDI 结果。
- Python 缓存和测试缓存。

## 运行说明

`0check_data.sh` 不是独立完成全部标准化的单文件脚本。它是 Bash 包装器：

```text
0check_data.sh
→ Python preprocess_wgdi_v2.py
→ Python validate_wgdi_dataset.py
→ data-validation-report.tsv
```

因此不能只单独发送 `.sh` 文件。必须保留包内的 `scripts/` 和 `src/`。

测试：

```powershell
cd plant-genome-resources
python -m pytest -q
ruff check src tests scripts
mypy src --ignore-missing-imports
```

预处理与验收：

```powershell
D:\Git\bin\bash.exe ..\0check_data.sh `
  --preprocess `
  --refresh-links `
  --expected-count 20
```

只验收已有结果：

```powershell
D:\Git\bin\bash.exe ..\0check_data.sh --expected-count 20
```

## 当前验证状态

- Python tests：29 passed
- Ruff：passed
- Mypy：passed
- Bash syntax：passed
- 清空派生数据后：缺少 manifest 时返回退出码 2，不假通过
- 之前 20 物种完整预处理和验收：20 PASS / 0 FAIL
- WGDI 0.75 dotplot 冒烟测试：通过

## 尚未完成

- `1run_blast.sh` 仍依赖根目录 `s*` 兼容链接。
- 3,000 物种并行分片后的 manifest 合并和任务调度尚未实现。
- 当前验收检查数据结构和哈希，不替代 BUSCO、污染检测等生物学质量评估。

