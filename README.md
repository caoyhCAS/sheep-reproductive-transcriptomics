# Sheep reproductive transcriptomics

**方法重建版，非原始代码 / Method reconstruction, not original author code.**

本项目依据 Yang et al. 的正文、补充流程图和公开测序记录重建分析步骤：

> Comparative mRNA and miRNA expression in European mouflon (Ovis musimon)
> and sheep (Ovis aries) provides novel insights into the genetic mechanisms
> for female reproductive success. *Heredity* 122, 172–186 (2019),
> online 2018. [DOI:10.1038/s41437-018-0090-1](https://doi.org/10.1038/s41437-018-0090-1).

核对的补充材料没有作者自编 pipeline/script 压缩包。本仓库包含新编排的
工作流、数据清单、统计包装脚本和后处理工具，不代表作者原始实现，也不
保证复现原文数值或图形。没有在此项目中执行全量 FASTQ 分析。

## 覆盖范围

| 分支 | 已实现 | 仍需提供或人工完成 |
|---|---|---|
| 数据 | SRP142554 / PRJNA451237 的 31 个 run 清单，FASTQ URL、大小和 MD5；样本标签核对 | 大文件下载、MD5 验证、配置实际路径 |
| mRNA | FastQC → Cutadapt → TopHat → Cufflinks/Cuffmerge/Cuffdiff 命令计划和显式执行；原生 Cuffdiff 统计筛选 | 匹配的基因组/GTF、索引、adapter、文库方向和旧软件 |
| de novo | Trinity、RSEM 命令计划，计数矩阵拼接、edgeR R 包装脚本，BLAST 导出过滤 | 历史依赖、基因/转录本聚合选择、nt 数据库、注释映射和 PANTHER |
| miRNA | FastQC/Cutadapt/PRINSEQ 预处理；前体结果筛选；DESeq **1.22.0** 包装脚本 | miRDeep2/Bowtie 冲突处理、miRBase21 分类及计数导出 |
| 整合 | DE 阈值、TargetScan 导出结果过滤、相反表达方向网络和核心网络 TSV | TargetScan7 预测所需 UTR/模型输入、历史 DAVID、Cytoscape 布局、miRTarBase 查证 |

这里的“部署”是发布可审查、可配置的代码，不是在 GitHub Actions 中进行
大规模生物信息计算。CI 只检查脚本、输入约束和合成数据；历史统计包的
数值结果以及论文真实数据不属于 CI 验证范围。

## 快速检查（无需生物信息软件）

需要 Python 3.10+；生产辅助脚本只使用标准库。

```bash
python -m venv .venv
source .venv/bin/activate
python -m pip install -r requirements-dev.txt
python -m pytest -q
python scripts/demo.py --output results/demo
python scripts/workflow.py --config config/paper.json --samples config/samples.tsv --stage all --output results/plan.json
```

`demo.py` 使用明确标识的合成数据，验证 2 条相反方向的网络边和 1 条核心
网络边。它不使用论文表达量，不生成虚构的研究结论。运行计划默认成功
输出尚未解决的参数与阻断原因，**不会开始测序数据下载或调用旧版软件**。

## 使用真实数据

1. 阅读 [方法对应表](docs/METHODS.md)、[参数冲突与缺失项](docs/PARAMETERS.md)
   和 [数据清单说明](docs/DATA.md)，核对 `config/runs.tsv` 的 run、文库和个体。
2. 准备旧版环境、参考 FASTA/GTF/索引，下载并验证 FASTQ，将路径写入
   `config/samples.tsv`。默认路径是占位路径，未附带原始数据。
3. 在 `config/paper.json` 填入真实 adapter、方向、参考文件版本和问题解决记录，
   查看 [工作流文档](docs/WORKFLOW.md)。仅在阻断项清除后使用 `--execute`。
4. 按 [R 环境说明](docs/LEGACY_ENVIRONMENTS.md) 运行 edgeR/DESeq 包装脚本，
   按 [后处理接口](docs/POSTPROCESSING.md) 导出并整合结果。

所有比较使用 `log2(numerator / denominator)`，整合时 mRNA 和 miRNA 的
方向必须完全一致。工作流方向写在计划中；示例为演示而选择 Finnsheep /
Mouflon，不能不经核对直接套用到其他比较。

## 影响解释的关键限制

- Ensembl83 的羊注释对应 Oar v3.1，而正文 mRNA 方法写 Oar v4.0；必须
  提供坐标一致的参考文件并记录解决方案。
- 正文 miRNA 写 Oar v3.1、结果段又概括为 Oar v4.0；“Bowtie v2.2.1”
  与 Bowtie1 名称/引用不一致。不能悄悄替换为另一种比对器。
- 6 份盘羊卵巢文库来自 **3 只动物**，两份盘羊子宫内膜来自 **同一只 M2**。
  M2-EA 没有 miRNA 文库。文库数不是独立生物学重复数；历史比较不等价于
  有充分生物学重复的总体推断。
- miRNA 前体过滤的 AND 含义、BLAST 的“100% matching ratio”、TargetScan
  输入及核心 miRNA 排序等未完整报告。相关脚本要求显式选择或明确标注重建规则。
- 对向表达和预测靶标属于关联证据；网络输出本身不能证明因果调控关系。

## 文件索引

- `config/`：样本、公开 run 清单、未填参数的工作流配置。
- `scripts/workflow.py`：计划生成、执行前检查、顺序调用与日志。
- `scripts/edger_denovo.R`、`scripts/deseq_mirna.R`：历史方法的显式统计包装。
- `scripts/expression_analysis.py`：Cuffdiff 标准化、DE 筛选、网络整合。
- `scripts/assemble_counts.py`：按 feature ID 拼接 RSEM expected counts。
- `scripts/filter_evidence.py`：前体/BLAST/已校正富集结果过滤。
- `scripts/fetch_run_manifest.py`：获取 ENA 清单，不自动下载 FASTQ。
- `docs/`：方法、来源、缺失参数、历史环境及数据接口。
- `tests/`、`.github/workflows/ci.yml`：离线合成数据和接口检查。

新编写代码使用 MIT 许可。论文为 CC BY4.0；外部软件、数据和数据库各自
保留原许可，未打包进本仓库。引用原文并记录所用仓库 commit，见
[CITATION.cff](CITATION.cff) 和 [来源记录](docs/SOURCES.md)。
