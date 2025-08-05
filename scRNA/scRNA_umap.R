getwd()
.libPaths()
# setwd("./Desktop/scRNAworkdir")
rm(list =ls())
# 0. 导入库

BiocManager::install("Seurat") 
library(BiocManager)
library(Seurat)
library(dplyr)
library(patchwork)
packageVersion("Seurat")
#——----------------------------------------------------------------
# 1. 导入数据####
#——----------------------------------------------------------------
# 有四种格式, 这里就是第一种10x官方提供的matrix.mtx,genes.tsv,barcodes.tsv
# 当然还有h5，rds，txt格式

pbmc.data <- Read10X(data.dir = "C:/Users/tommyw/Desktop/scRNAworkdir/pbmc3k_filtered_gene_bc_matrices/filtered_gene_bc_matrices/hg19")
pbmc.data
# 初始化seurat数据，得到一个seurat对象
pbmc<-CreateSeuratObject(counts = pbmc.data,project = "pbmc3k",min.cells = 3,min.features = 200)
pbmc

# 特殊处理，若有多个样本，考虑行名称以及融合matrix 
# 假设已有 Seurat 对象 pbmc，且 orig.ident 为 "pbmc3k"
pbmc_label <- unique(pbmc$orig.ident)  # 获取样本标签，一般是 "pbmc3k"

# 给细胞名添加前缀（样本名）
new_cell_names <- paste0(pbmc_label, "_", colnames(pbmc))
colnames(pbmc) <- new_cell_names
rownames(pbmc@meta.data) <- new_cell_names  # 同步更新 meta.data 的行名

# 检查更新后的效果
head(pbmc@meta.data)


# 获得表达矩阵和注释信息，注意要注意指定active_assays的表达矩阵，即更改默认矩阵
# assays是表达矩阵，标准化，归一化，SCT处理后会增加多一层又一层layer
# pbmc@assays 存储表达矩阵的
# pbmc@active.assay 默认表达矩阵
# pbmc@meta.data$orig.ident存储细胞注释信息的细胞标识
#——----------------------------------------------------------------
# 2. 数据质检####
#——----------------------------------------------------------------
# QC的指标：在每个细胞中检测到的独特基因的数量(nFeature_RNA), 总的转录本数(测序深度), mt基因
# 其实还有红细胞以及核糖体基因但教程没有提到
# 使用PercentageFeatureSet函数计算线粒体QC指标
pbmc[['percent.mt']]<-PercentageFeatureSet(pbmc,pattern = "^MT-")
head(pbmc@meta.data,5)
# 使用violin plot可视化  QC指标，并使用这些指标过滤单元格
VlnPlot(pbmc,features = c("nFeature_RNA","nCount_RNA","percent.mt"),ncol = 3)

# FeatureScatter 通常用于可视化两个特征之间的关系，大概看一下数据，找过滤条件。
plot1<-FeatureScatter(pbmc,feature1 = "nCount_RNA",feature2 = "percent.mt")
plot2<-FeatureScatter(pbmc,feature1 = "nCount_RNA",feature2 = "nFeature_RNA")
plot1+plot2

# 将QC指标可视化，并使用这些指标过滤单元格
    # 某些细胞 nFeature_RNA > 5000，而大部分集中在 <2500，可能是双细胞
    # 某些细胞 percent.mt > 10，表示线粒体污染严重，考虑剔除
# 过滤具有2500或少于200的独特特征计数的单元格，过滤线粒体计数>5%的细胞
pbmc<-subset(pbmc,subset=nFeature_RNA>200 & nFeature_RNA<2500 & percent.mt<5)


#——----------------------------------------------------------------
# 3. 数据标准化####
#——----------------------------------------------------------------
# 标准化数据：pbmc对象增加一层layer：data
# 等价pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
# 原始表达 = 200, 细胞总 UMI = 20000, scale.factor = 10000,
# 归一化值 = (200 / 20000) * 10000 = 100, LogNormalize = log1p(100) ≈ 4.615
pbmc <- NormalizeData(pbmc)

#——----------------------------------------------------------------
# 4. PCA####
#——----------------------------------------------------------------
# 计算在数据集中表现出高细胞间变异的特征子集，识别高度可变的特征（特征选择）
pbmc<-FindVariableFeatures(pbmc,selection.method = "vst",nfeatures = 2000)
# 确定高表达的前十个基因
top10<-head(VariableFeatures(pbmc),10)
# 画高变基因散点图
plot1 <- VariableFeaturePlot(pbmc)
plot1
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE, xnudge = 0, ynudge = 0)
plot2

# 缩放数据,因为PCA对数据敏感，把数据均值变成0，方差变成1
all.genes<-rownames(pbmc)
pbmc<-ScaleData(pbmc,features = all.genes)

# 线性降维，看一下每个主成分的top基因
pbmc<-RunPCA(pbmc,features = VariableFeatures(object=pbmc))
print(pbmc[["pca"]],dims = 1:5,nfeatures = 5)

# Seurat提供可视化细胞和定义PCA
# 基因载荷loadings代表：每个基因对主成分的贡献大小
VizDimLoadings(pbmc,dims = 1:2,reduction = "pca")

DimPlot(pbmc,reduction = "pca")
#head看一下
pbmc[['pca']]@cell.embeddings[1:5,1:2]
# 热图用于确定哪些主成分对区分细胞亚群重要，以及看看标签
DimHeatmap(pbmc,dims = 1,cells = 500,balanced = TRUE)
# 前15个主成分（PC1 ~ PC15） 
DimHeatmap(pbmc,dims = 1:15,cells=500,balanced = TRUE)
# 自动判断保留前多少个主成分（PC）用于聚类、UMAP 是合理的，确定数据集的维度。
pbmc<-JackStraw(pbmc,num.replicate = 100)
pbmc<-ScoreJackStraw(pbmc,dims = 1:20)
# 可视化处理
JackStrawPlot(pbmc,dims = 1:15)

ElbowPlot(pbmc)
# 以上就是为了选择保留多少主成分分析，千万不要用太少（比如只用 5 个 PC）
# 会显著削弱聚类和分群的效果；宁愿选多一点 也别选太少（这是最安全的做法）。
#——----------------------------------------------------------------
# 5. 聚类细胞####
#——----------------------------------------------------------------
# 建立KNN图，并基于其局部领域中的共享重叠细化任意两个单元之间的边权重
pbmc<-FindNeighbors(pbmc,dims = 1:10)
# 对细胞进行聚类，使得模块化最大化，Louvain算法（默认）或SLM
# resolution的值可以尝试几个典型值（如 0.2、0.5、1.0）
# 看看后续UMAP聚类以及marker分布是否合理
pbmc<-FindClusters(pbmc,resolution = 0.5)
head(Idents(pbmc),5)

#——----------------------------------------------------------------
# 6. 运行非线性降维（UMAP/tSNE）####
#——----------------------------------------------------------------
pbmc<-RunUMAP(pbmc,dims = 1:10)
DimPlot(pbmc,reduction = "umap")

#——----------------------------------------------------------------
# 7. 寻找差异表达的特征（marker）####
#——----------------------------------------------------------------
# findmarkers为所有集群自动执行此过程，也可以测试集群组之间的对比，或针对所有单元格进行测试
# 默认情况下，ident.1与所有其他细胞相比，他识别单个簇的阳性和阴性标记。

# 寻找cluster2的所有markers(重要)
cluster2.markers<-FindMarkers(pbmc,ident.1 = 2,min.pct = 0.25)
head(cluster2.markers,n=5)

# 寻找cluster5中与cluster0和cluster3n不同的所有markers
cluster5.markers<-FindMarkers(pbmc,ident.1 = 5,ident.2=c(0,3),min.pct = 0.25)
head(cluster5.markers,n=5)

# 找出每个细胞簇的标记物，与所有剩余的细胞进行比较，只报告阳性细胞(重要) 
pbmc.markers<-FindAllMarkers(pbmc,only.pos = TRUE,min.pct = 0.25,logfc.threshold =0.25 )

# 从所有marker中，选出每个簇中表达差异最显著的前2个基因(重要)
pbmc.markers %>%
  group_by(cluster) %>%
  slice_max(n=2,order_by = avg_log2FC)

#——----------------------------------------------------------------
# 8.可视化差异基因的办法主要是五种
#——----------------------------------------------------------------
    # Vlnplot(显示跨集群的表达概率分布)
    # FeaturePlot()(在tSNE或PCA图上可视化特征表达)是最常用的可视化
VlnPlot(pbmc,features = c("MS4A1","CD79A"))
VlnPlot(pbmc,features = c("NKG7","PF4"),slot = "counts",log = TRUE)

FeaturePlot(pbmc,features = c("MS4A1", "GNLY", "CD3E", "CD14", "FCER1A", "FCGR3A", "LYZ", "PPBP","CD8A"))
# 使用不多的有3种
# 密度曲线，点阵图最常见，
RidgePlot(pbmc, features = c("CD3E", "MS4A1"))

DotPlot(pbmc, features = c("CD3E", "MS4A1", "GNLY")) + RotatedAxis()

DoHeatmap(pbmc,features= c("CD3E", "MS4A1"),size=1)

#——----------------------------------------------------------------
# 9. 将细胞类型标识分配给集群####
#——----------------------------------------------------------------
new.cluster.ids<-c("Naive CD4 T", "CD14+ Mono", "Memory CD4 T", "B", "CD8 T", "FCGR3A+ Mono",
                   "NK", "DC", "Platelet")
names(new.cluster.ids)<-levels(pbmc)
pbmc<-RenameIdents(pbmc,new.cluster.ids)
DimPlot(pbmc,reduction = "umap",label = TRUE,pt.size = 0.5)+NoLegend()

# 另一种方法：再meta.data添加结果构建细胞亚群注释表
celltype <- data.frame(
  ClusterID = 0:8,
  celltype = "unknown",
  stringsAsFactors = FALSE
)

# 对各个cluster赋予细胞类型标签
celltype[celltype$ClusterID == 0, "celltype"] <- "Naive CD4 T"
celltype[celltype$ClusterID == 1, "celltype"] <- "CD14+ Mono"
celltype[celltype$ClusterID == 2, "celltype"] <- "Memory CD4 T"
celltype[celltype$ClusterID == 3, "celltype"] <- "B"
celltype[celltype$ClusterID == 4, "celltype"] <- "CD8 T"
celltype[celltype$ClusterID == 5, "celltype"] <- "FCGR3A+ Mono"
celltype[celltype$ClusterID == 6, "celltype"] <- "NK"
celltype[celltype$ClusterID == 7, "celltype"] <- "DC"
celltype[celltype$ClusterID == 8, "celltype"] <- "Platelet"

# 查看注释统计
table(celltype$celltype)

# 将 celltype 信息映射到 pbmc 对象中
sce.all <- pbmc  # 复制一份对象（可选）
sce.all@meta.data$celltype <- "NA"

for (i in seq_len(nrow(celltype))) {
  cluster_id <- celltype$ClusterID[i]
  cluster_label <- celltype$celltype[i]
  sce.all@meta.data[sce.all@meta.data$seurat_clusters == cluster_id, "celltype"] <- cluster_label
}

# 查看分布情况
table(sce.all@meta.data$celltype)
table(sce.all@meta.data$celltype, sce.all@meta.data$orig.ident)
table(sce.all@meta.data$celltype, sce.all@meta.data$seurat_clusters)

# 画图：UMAP 按 cluster 和 celltype 展示
DimPlot(sce.all, reduction = "umap", label = TRUE, group.by = "seurat_clusters", pt.size = 0.5) + NoLegend()
DimPlot(sce.all, reduction = "umap", label = TRUE, group.by = "celltype", pt.size = 0.5) + NoLegend()

# 将 celltype 设置为默认的分组标识
Idents(sce.all) <- sce.all$celltype

# 保存了一个 R 对象的快照”
saveRDS(pbmc,file = "./pbmc3k_final.rds")

