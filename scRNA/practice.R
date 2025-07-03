# 10:04 2025/7/3 tested by wky
# 该脚本用于对scRNA进行入门级的简单处理

library(dplyr)
library(Seurat)
library(patchwork)
#——————————————————————————————————————————————————————————————————————————————————————————————————————————————————————————————
# Step1：读取数据，对数据进行过滤
#——————————————————————————————————————————————————————————————————————————————————————————————————————————————————————————————
# 1.读取scRNA-Seq数据

count.data <- Read10X(data.dir = "/opt/exam/")

# 2.建立Seurat数据结构，此时的数据结构中共有多少个单细胞？

scdata <- CreateSeuratObject(counts = count.data$`Gene Expression`, project = "exam", min.cells = 3, min.features = 200)
scdata

# 结果2：共有7867个单细胞
# An object of class Seurat 
# 17872 features across 7867 samples within 1 assay 
# Active assay: RNA (17872 features, 0 variable features)

# 3.计算线粒体基因表达量占每个细胞内基因总表达的比例

scdata[["percent.mt"]] <- PercentageFeatureSet(scdata, pattern = "^MT-")
head(scdata@meta.data)

# 结果3：得到新的一列percent.mt，该列包含线粒体基因表达量占每个细胞内基因总表达的比例
# orig.ident nCount_RNA nFeature_RNA percent.mt
# AAACCTGAGATGCCAG-1       exam       8824         2322   2.878513
# AAACCTGAGGAATTAC-1       exam      20440         4244   3.444227
# AAACCTGAGGACACCA-1       exam       2479         1145   2.258975
# ...

# 4.绘制单细胞检测基因数、基因表达总量和线粒体基因表达占比的小提琴图，截图放置于此

VlnPlot(scdata, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# 5.根据质量过滤细胞

scdata <- subset(scdata, subset = nFeature_RNA > 200 & nFeature_RNA < 3000 & percent.mt < 10)
scdata

# 结果5：过滤后还有7413个细胞
# An object of class Seurat 
# 17872 features across 7413 samples within 1 assay 
# Active assay: RNA (17872 features, 0 variable features)

# 6. 上一步数据预处理的目的：
    # 去除检测基因过少的细胞（nFeature_RNA ≤ 200），这类细胞可能是空液滴（Empty droplets）或死亡细胞
    # 去除检测基因过多的细胞（nFeature_RNA ≥ 3000），这类“细胞”可能实际上是多个细胞混合在一起的“双细胞/多细胞”（Doublets or Multiplets）
    # 去除线粒体基因表达比例高的细胞（percent.mt ≥ 10）,线粒体基因表达占比高的细胞，通常是处于应激状态、细胞膜破裂或即将死亡的细胞。


#——————————————————————————————————————————————————————————————————————————————————————————————————————————————————————————————
# Step2：对数据聚类找cluster
#——————————————————————————————————————————————————————————————————————————————————————————————————————————————————————————————

# 7. 找出前2000个高变基因（VST算法）

scdata <- FindVariableFeatures(
  object = scdata,
  selection.method = "vst",
  nfeatures = 2000
)
# 直接绘图显示高变基因散点图
VariableFeaturePlot(scdata) + 
  ggtitle("Top 2000 Highly Variable Genes (VST)")

# 打印排名前十的高变基因
top10_hvg <- head(VariableFeatures(scdata), 10)
top10_hvg
 
 # 结果7：展示变异分数（variance stabilized score）排名前10的hvg
 # [1] "HBA2"     "HBG2"     "HBB"      "IGLV2-14" "HBA1"     "IGKV3-20" "IGHV1-2"  "PPBP"    
 # [9] "HBM"      "ALAS2"  

 # 8. 对数据进行缩放和中心化，默认使用高变基因

scdata <- ScaleData(scdata, features = VariableFeatures(scdata))
 # 运行PCA，使用高变基因进行降维
scdata <- RunPCA(scdata, features = VariableFeatures(scdata))

# 绘制ElbowPlot，帮助判断选取多少主成分
ElbowPlot(scdata)

### 9. 根据前10个主成分进行聚类分析

scdata <- FindNeighbors(scdata, dims = 1:10)  # 构建邻接图
scdata <- FindClusters(scdata, resolution = 0.5)  # Louvain 聚类

# 查看共获得多少个细胞群
table(scdata$seurat_clusters)  # 输出每个 cluster 的细胞数量
length(unique(scdata$seurat_clusters))  # 输出 cluster 的总数

# 结果9：一共获得13个细胞群
# [1] 13
#——————————————————————————————————————————————————————————————————————————————————————————————————————————————————————————————
# Step3： 执行 UMAP 降维
#——————————————————————————————————————————————————————————————————————————————————————————————————————————————————————————————

scdata <- RunUMAP(scdata, dims = 1:10)
# 绘制 UMAP 图
DimPlot(scdata, reduction = "umap", label = TRUE) +
  ggtitle("UMAP plot with clusters (resolution = 0.5)")

#——————————————————————————————————————————————————————————————————————————————————————————————————————————————————————————————
# Step4：找差异基因中的marker基因
#——————————————————————————————————————————————————————————————————————————————————————————————————————————————————————————————

# 11. 计算 Marker 基因：使用默认差异表达方法（Wilcoxon）
markers <- FindAllMarkers(
  scdata,
  only.pos = TRUE,           # 只保留正向 marker（上调基因）
  min.pct = 0.25,            # 至少在25%的细胞中表达
  logfc.threshold = 0.25     # log2(FC) > 0.25（大约是fold change > 1.19）
)
# 每个 cluster 取前5个 marker 基因（按 avg_log2FC 排序）
top5_markers <- markers %>% 
  group_by(cluster) %>% 
  top_n(n = 5, wt = avg_log2FC)
top5_markers
top5_markers %>% count(cluster)
top5_markers %>%
  group_by(cluster) %>%
  summarise(marker_genes = paste(gene, collapse = ", "))

# 打印结果11，结果包含所有每个细胞群（cluster）的前5个 marker 基因，但我发现总数目不对，根据筛选条件有些cluster只取到3个满足不了阈值，所以我
# A tibble: 13 × 2
#   cluster marker_genes                          
#   <fct>   <chr>                                 
# 1 0       LEF1, IL7R, SARAF, CCR7, TCF7         
# 2 1       IL7R, RPS18, LEF1, RPS6, CD8B         
# 3 2       S100A9, LYZ, THBS1, IL1B, CST3        
# 4 3       PIK3IP1, LINC00861, TXNIP             
# 5 4       NKG7, GNLY, CST7, PRF1, CCL5          
# 6 5       S100A8, CXCL8, S100A12, S100A6, FTH1  
# 7 6       CD79A, CD74, HLA-DRA, TCL1A, HLA-DPB1 
# 8 7       KLRB1, ZFP36L2, IL32, GZMK, CTSW      
# 9 8       PPBP, PF4, TUBB1, MYL9, SDPR          
#10 9       RPS6, RPS18, IL7R, CD8B, RP11-291B21.2
#11 10      THBS1, LYZ, S100A9, G0S2, FTH1        
#12 11      HBB, HBA2, HBA1, ALAS2, HBG2          
#13 12      JCHAIN, PTGDS, PLD4, GZMB, ITM2C 

# 筛选标准：
    # min.pct = 0.25	只考虑在该细胞群和其他细胞群中至少 25% 细胞表达的基因，避免选到极少数细胞表达的噪声基因。
    # logfc.threshold = 0.25	只保留log2表达倍数差异大于0.25的基因，即表达量至少有约1.19倍差异，确保差异有生物学意义。
    # only.pos = TRUE	只选在该细胞群中上调的基因（正向marker），这些基因更能代表细胞群的特征。

library(dplyr)
#——————————————————————————————————————————————————————————————————————————————————————————————————————————————————————————————
# Step5：找聚类后cluster里面的特征信息，比如最大群，最小群，哪一群是NK细胞等等
#——————————————————————————————————————————————————————————————————————————————————————————————————————————————————————————————
1. 统计每个细胞群大小
cluster_sizes <- table(scdata$seurat_clusters)
print(cluster_sizes)

# 找最大和最小群
max_cluster <- names(cluster_sizes)[which.max(cluster_sizes)]
min_cluster <- names(cluster_sizes)[which.min(cluster_sizes)]

cat("最大细胞群：Cluster", max_cluster, "，细胞数 =", max(cluster_sizes), "\n")
cat("最小细胞群：Cluster", min_cluster, "，细胞数 =", min(cluster_sizes), "\n")

# 最大细胞群：Cluster 0 ，细胞数 = 1423 
# 最小细胞群：Cluster 12 ，细胞数 = 26 

# 2. 找最大群的前5个marker基因
max_cluster_markers <- markers %>%
  filter(cluster == max_cluster) %>%
  arrange(desc(avg_log2FC)) %>%
  head(5)

print(max_cluster_markers$gene)

# 最大群的前5个marker基因（按照avg_log2FC排列）:[1] "IL7R"  "SARAF" "LEF1"  "TCF7"  "CCR7"

# 3. 判断最小群是否有marker基因
min_cluster_markers <- markers %>%
  filter(cluster == min_cluster)

min_cluster_markers$gene

# 最小群有marker基因：
    # [1] "LILRA4"        "RASD1"         "CLEC4C"        "LINC00996"     "SCT"           "MYBL2"      
    # [7] "LRRC26"        "PTPRS"         "LAMP5"         "KRT5"          "PACSIN1"       "LCNL1"
    # ... 
# 找不到潜在原因：细胞数量过少，该群细胞异质性低或无明显特征，技术噪声或数据质量问题，聚类结果过度细分或不合理，筛选参数过于严格。

# 3. NK细胞经典marker基因列表
nk_markers <- c("NCAM1", "NKG7", "GNLY", "KLRD1", "KLRF1", "PRF1", "GZMB")

# 画DotPlot展示这些基因在各cluster的表达情况
DotPlot(scdata, features = nk_markers) + RotatedAxis() + 
  ggtitle("Expression of NK cell markers across clusters")

# 结果：细胞群4非常有可能是NK细胞群   
#————————————————————————————————————————————————————————————————————————————————
### 最后. 提取每个cluster前5个marker基因，看一下这些基因来自什么细胞，来自什么组织
#——————————————————————————————————————————————————————————————————————————————————
top5_markers <- markers %>%
  group_by(cluster) %>%
  top_n(n = 5, wt = avg_log2FC)
top5_markers

# 根据对各个细胞群 marker 基因的分析，样本中包含丰富的 T 细胞、B 细胞、NK 细胞、单核细胞和血小板等群体
# 组成符合外周血单个核细胞（PBMC）的特征，因此推测该单细胞测序样本最可能来源于人类外周血。
# 可选择的看一下PBMC各个细胞的marker，可视化 marker gene 表达（例如免疫相关）
# common_markers <- c("CD3D", "CD8A", "CD4", "MS4A1", "NKG7", "GNLY", "LYZ", "CD14", "FCGR3A", "PPBP")
# DotPlot(scdata, features = common_markers) + RotatedAxis() + 
    # ggtitle("Key immune cell marker expression across clusters")
