# 核心思路：找出跨批次的HVG，分batch, PCA, harmony校正，
setwd("./scRNAworkdir")
rm(list =ls())
.libPaths()
install.packages("harmony")
library(data.table)
library(tidyverse)
library(ggthemes)
library(ggrepel)
library(harmony)
library(patchwork)
library(tidyr)
library(ggplot2)
library(Seurat)

## 一. 把内置的两个 PBMC 单细胞表达矩阵（刺激 vs 对照）拼接成一个 Seurat 对象
data("pbmc_stim")
pbmc <- CreateSeuratObject(counts = cbind(pbmc.stim, pbmc.ctrl), 
                           project = "PBMC", min.cells = 5)

## 二. 添加batch信息：Seurat 对象 pbmc 的元数据（meta.data）里新建一列，名叫 batch
pbmc@meta.data$batch <- c(rep("STIM", ncol(pbmc.stim)), rep("CTRL", ncol(pbmc.ctrl)))

## 三. 分别在 STIM 批次 和 CTRL 批次各找出 2000 个高变基因（HVG）
# 然后把两批结果合并、去重，最后把 “跨批次共有高变基因” 设为 Seurat 对象的 VariableFeatures
pbmc <- pbmc %>%
  NormalizeData(verbose = FALSE)
VariableFeatures(pbmc) <- split(row.names(pbmc@meta.data), pbmc@meta.data$batch) %>% lapply(function(cells_use) {
  pbmc[,cells_use] %>%
    FindVariableFeatures(selection.method = "vst", nfeatures = 2000) %>% 
    VariableFeatures()
}) %>% unlist %>% unique

## 四. PCA降维
pbmc <- pbmc %>% 
  ScaleData(verbose = FALSE) %>% 
  RunPCA(features = VariableFeatures(pbmc), npcs = 20, verbose = FALSE)

## 五. 去除批次效应前的可视化
## 未进行批次效应校正，反映存在批次效应，PC的散点图如果分的很开说明批次效应强
library(cowplot)
options(repr.plot.height = 5, repr.plot.width = 12)
p1 <- DimPlot(object = pbmc, reduction = "pca", pt.size = .1, group.by = "batch")
p2 <- VlnPlot(object = pbmc, features = "PC_1", group.by = "batch", pt.size = .1)
plot_grid(p1,p2)

## 六. 使用 Harmony 算法（Harmony 之前一定要先做 PCA）
pbmc <- pbmc %>% 
  RunHarmony("batch", plot_convergence = TRUE, 
             nclust = 50, max_iter = 10, early_stop = T)

# RunHarmony的运行结果在pbmc@reductions$harmony中。
# 把 Harmony 校正后的 20 维嵌入矩阵拿出来，瞄一眼前 4×4 的数值
harmony_embeddings <- Embeddings(pbmc, 'harmony')
harmony_embeddings[1:4, 1:4]

## 七. 去除批次效应后的可视化
options(repr.plot.height = 5, repr.plot.width = 12)
p1 <- DimPlot(object = pbmc, reduction = "harmony", pt.size = .1, 
              group.by = "batch")
p2 <- VlnPlot(object = pbmc, features = "harmony_1", 
              group.by = "batch", pt.size = .1)
plot_grid(p1,p2)
