
# 单细胞转录组入门简明教程

一直想写一个单细胞转录组入门的简明教程希望可以帮助更多的人来了解单细胞转录组的数据。
可以先看看基本分析流程，对处理过程有一个基本的印象：
本教程旨在帮助初学者了解单细胞转录组分析的基本流程和操作，使用的是 Seurat 包和 PBMC3k 数据集。

---

## 🔁 基本分析流程

1. 原始数据处理  
2. 质量控制与过滤  
3. 数据标准化  
4. 高变基因识别  
5. 降维  
6. 聚类  
7. Marker 基因识别  
8. 细胞类型注释  
9. 下游分析  


在入手单细胞转录组分析之前，强烈建议学习一下这个Seurat官方的经典教程https://satijalab.org/seurat/articles/pbmc3k_tutorial.html

本次分析数据来源是外周血(PBMC)，也是免疫学经典的分析对象。使用的是10X Genomics的Illumina NextSeq 500进行测序，总共有2700个单细胞sample。
#注意：不同测序技术分析的过程会有所不同，比如Smart-seq 和 10x Genomics（简称 10x） 就是两种完全不同的单细胞 RNA 测序技术。

接下来直接进入code环节：

## 🔧 A. 环境准备

1.下载安装R语言和Rstuido这两个软件。

2.创建好工作目录以及R包的安装目录： 这两者是R语言里很重要的两个文件夹，在R里面输入getwd()和.libPaths()可以看到。
# 本次使用的工作目录是"D:/RGuide/Workdir"。

3.安装Seurat包： 是R语言中用于分析单细胞测序很重要的一个R包(即函数)，他的核心作用是将表达矩阵，还存储降维、聚类、注释、元数据等信息到集成到一个Seurat 对象上，可以理解为是一个“细胞数据的容器”对象。之后的处理都紧紧围绕着Seurat这个对象。就好比是一个人拥有很多特征眼睛，鼻子，耳朵...

下载这个Seurat包可以使用BiocManager::install("Seurat"):通常来说 R 有两个主流的软件包来源：CRAN和生信专用的Bioconductor，下载方式分别是install.packages()和BiocManager::install()。

下载后记得library()加载!!!

4.下载本次分析数据：
链接来源(https://cf.10xgenomics.com/samples/cell/pbmc3k/pbmc3k_filtered_gene_bc_matrices.tar.gz)下载后解压到本次分析使用的工作目录。


## 📥 B. 创建 Seurat 对象
Seurat函数提供了很多方法来读取细胞稀疏数据矩阵，读取的文件有四种格式, 这里就是第一种10x官方提供的matrix.mtx, genes.tsv, barcodes.tsv，当然还有四种h5，rds，txt后缀的格式

对Seurat对象中的meta.data数据框做个解释，他是包含着重要的讯息，用于存储每个单细胞的元数据，包括每个细胞的总表达量，每个细胞的基因数，线粒体/核糖体基因百分比，聚类标签、样本分组、时间点等
后续分析中的各种计算结果（如细胞类型注释）等等。


本次数据中meta.data每一行是一个细胞（由条形码 ID 标识），每一列是这个细胞的属性（元数据）：每一行是如AAACATACAACCAC-1 是叫cell barcode类似于细胞的身份证id，orig.ident是初始项目名称，这里只有一个批次，所以不考虑其他批次。nCount_RNA，nFeature_RNA代表每个细胞的总表达量(因为基因往往会有一对多的转录本)和每个细胞中有表达的基因个数（特征数）。

## 🧹 C. 质量控制与过滤
主要是用三个指标来过滤独特基因的数量(nFeature_RNA), 每个细胞的总表达量(nCount_RNA), 线粒体基因比例(percent.mt)
使用三个核心指标进行细胞过滤：

- `nFeature_RNA`: 检测到的基因数（<200 或 >2500 可能为死细胞或双细胞）
- `nCount_RNA`: 总 UMI 数（用于辅助判断）
- `percent.mt`: 线粒体基因占比（>5% 说明细胞可能死亡

过滤后细胞数目减少到2638个细胞


## 📊 D. 标准化与降维

主要是对数据进行标准化，寻找前2000个高度可变的基因（HVG，Highly Variable Genes），中心化以及PCA降维处理
1. **NormalizeData**：消除测序深度差异  
2. **FindVariableFeatures**：识别前 2000 个高变基因（HVG）  
3. **ScaleData**：对基因进行标准化（均值为 0，标准差为 1）  
4. **RunPCA**：主成分分析（PCA）
单细胞数据中，不同细胞的测序深度（总转录本数）不一样，导致表达量不可直接比较，所以必须先进行标准化处理。接着按每个基因做中心化+方差缩放，使其均值为 0，标准差为 1，因为PCA 对数据敏感，需要均值变 0 方差变 1” ScaleData() 这一步。

之后可以用ElbowPlot图评估选择多少个主成分（PC）进行用于后续聚类。弯折处 → 通常代表信息量下降显著的位置


## 🔗 E. 聚类与可视化

根据前面ElbowPlot弯折处选取10个PCA成分做聚类。
可视化的方法有很多种，常见的可以使用FeaturePlot和DotPlot


## 🔍 F. Marker 基因与细胞注释

1. 查找每个簇的 Marker 基因：
markers <- FindAllMarkers(pbmc)
2. 注释细胞类型（根据 Marker 基因手动或自动注释）
3. 保存 Seurat 对象：
saveRDS(pbmc, file = "pbmc3k_annotated.rds")

---

📌 **至此，你已经完成了一个单细胞转录组的基本分析流程！**


