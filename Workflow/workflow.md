17:20 2025/7/12
下游常规分析流程
细胞聚类

方法：t‑SNE / UMAP

目的：将高维表达数据投射到 2D/3D 空间，直观展示细胞群体结构

聚类结果注释

筛选各 cluster 的 Marker Gene（差异表达基因）

常见标记：

CD3E → T 细胞

MS4A1 → B 细胞

…

根据已知 marker gene 列表，注释细胞亚群

功能富集分析

差异表达分析

选取感兴趣的对比组或 cluster

得到 DEG（Differentially Expressed Genes）

GO/KEGG 富集

输入 DEG 列表

输出显著富集的生物过程 / 通路

GSEA（Gene Set Enrichment Analysis）

对全基因表达排序（按 fold change 或统计量）

分析通路在基因集矩阵中的整体激活趋势

GSVA（Gene Set Variation Analysis）

在细胞 / 样本水平对通路活性打分

代谢通路推断

使用专门针对代谢途径的基因集或工具（如 scMetabolism）

伪时序分析

常用工具：Monocle、Slingshot、PAGA（Scanpy）

目的：重现细胞状态转变路径

例如：干细胞 → 中间状态 → 分化细胞

细胞–细胞通讯分析

预测不同 cell cluster 之间的信号通路互作

常用工具：CellPhoneDB、CellChat、NATMI 等

VDJ（TCR/BCR）分析

克隆型分析

统计每个样本中的所有 clonotype

按细胞数或 UMI 数从高到低排序

识别哪些样本或哪些细胞亚群富集特定 clonotype

VDJ 基因片段使用

分析 V、D、J 各片段在样本／亚群中的过度使用情况

CDR3 长度多样性

比较不同 CDR3 长度分布

评估对亲和力和特异性的潜在影响

