# 批次效应的处理过程：

# 1. 创建环境，下载必要的数据集和包
# ————————————————————————————————————————————————————————
# conda create -n seurat_env r-base=4.2 r-essentials -c conda-forge
# 激活conda环境
# 设置 R 包安装路径
.libPaths("/data5/GPT/Wuky/Rlib")
# 设置工作目录
setwd("/data5/GPT/Wuky/Workdir/")
getwd()
.libPaths()
rm(list =ls())


# 数据集：提前下载ifnb.SeuratData_3.1.0.tar.gz
install.packages("/data5/GPT/Wuky/ifnb.SeuratData_3.1.0.tar.gz", repos = NULL, type = "source", lib = "/data5/GPT/Wuky/Rlib")

# 包的安装：下载R包SeuratData
if (!requireNamespace("remotes", quietly = TRUE)) {
  install.packages("remotes")
}

remotes::install_github("satijalab/seurat-data")

library(ggplot2)
library(Seurat)
library(patchwork)
library(SeuratData)

# 2. 加载数据，加载入数据集对象####
# ————————————————————————————————————————————————————————
# 识别两个数据集中存在的细胞亚群
# 获得在对照细胞和刺激细胞中均保守的细胞类型标记物
# 比较数据集以发现细胞类型对刺激的特异性反应

# Get Started！

# load dataset
ifnb <- LoadData("ifnb")

# 将多个样本的数据储存进多个Seurat Object中，并用一个list变量组织这些Seurat Object
ifnb[["RNA"]] <- split(ifnb[["RNA"]], f = ifnb$stim)
ifnb


# 核心概念理解
# counts：原始的基因表达计数矩阵（稀疏矩阵），行是基因，列是细胞。
# data：对 counts 进行归一化／对数转换后的矩阵（NormalizeData() 生成）。
# scale.data：在 data 基础上进一步中心化、标准化（ScaleData() 生成），用于 PCA 等下游分析。
# layers：用户自定义的其他矩阵，比如你刚才将 CTRL/STIM 拆分后的 counts.CTRL、counts.STIM 等。如counts.CTRL, counts.STIM, data.CTRL, data.STIM 

# 3. 在不集成的情况下执行分析####
# 预处理（必须为整合做准备）
ifnb <- NormalizeData(ifnb)
ifnb <- FindVariableFeatures(ifnb)
ifnb <- ScaleData(ifnb)
ifnb <- RunPCA(ifnb)

# 绘图
ifnb <- FindNeighbors(ifnb, dims = 1:30, reduction = "pca")
ifnb <- FindClusters(ifnb, resolution = 2, cluster.name = "unintegrated_clusters")
ifnb <- RunUMAP(ifnb, dims = 1:30, reduction = "pca", reduction.name = "umap.unintegrated")
DimPlot(ifnb, reduction = "umap.unintegrated", group.by = c("stim", "seurat_clusters"))

# 3. 整合多样本数据集####
# 消除批量效应、进行联合分析的基础，是找到数据集之间的“锚点”，并以锚点为基准来整合多个样本的数据集，利用锚点整合的一般方法为典型关联分析（CCA）方法。
# 锚点：不同数据集之间，由一组共同的分子特征定义的两个细胞（每个数据集一个），将对应关系称为锚定，将共同分子称为锚点。可以类比映射中的“双射”概念。
# 集成分析 则是在确认“条件效应”会干扰你对“细胞类型”比较时，专门把条件偏差校正掉。
# 整合CCA就是用数学方法“消除批次差异”，然后基于真实的细胞类型特征对细胞进行聚类和分析。


ifnb <- IntegrateLayers(object = ifnb, method = CCAIntegration, orig.reduction = "pca", new.reduction = "integrated.cca",
                        verbose = FALSE)

# re-join layers after integration
ifnb[["RNA"]] <- JoinLayers(ifnb[["RNA"]])

ifnb <- FindNeighbors(ifnb, reduction = "integrated.cca", dims = 1:30)
ifnb <- FindClusters(ifnb, resolution = 1)
ifnb <- RunUMAP(ifnb, dims = 1:30, reduction = "integrated.cca")
# Visualization
DimPlot(ifnb, reduction = "umap", group.by = c("stim", "seurat_annotations"))
DimPlot(ifnb, reduction = "umap", split.by = "stim")