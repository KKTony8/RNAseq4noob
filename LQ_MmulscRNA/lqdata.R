library(Seurat)
library(dplyr)

# 设置包含.h5文件的目录
setwd("/data5/GPT/Wuky/1LQ_scRNA/LQdata/h5d/")

#——————————————————————————————————————————————————————————
# 一.读取所有.h5文件####
#——————————————————————————————————————————————————————————
h5_files <- list.files(pattern = "\\.h5$")
print(paste("找到", length(h5_files), "个H5文件"))

# 1.初始化列表存储Seurat对象
seurat_list <- list()

# 2.循环读取每个文件
for (file in h5_files) {
  # 从文件名提取样本ID（去掉扩展名）
  sample_id <- tools::file_path_sans_ext(file)
  
  # 读取10x H5数据
  counts <- Read10X_h5(filename = file)
  
  # 创建Seurat对象
  seurat_obj <- CreateSeuratObject(
    counts = counts,
    project = sample_id,
    min.cells = 3,      # 基因至少在3个细胞中表达
    min.features = 200  # 细胞至少检测到200个基因
  )
  
  # 添加样本ID到元数据
  seurat_obj$sample <- sample_id
  
  # 存储到列表
  seurat_list[[sample_id]] <- seurat_obj
}

# 3.合并所有Seurat对象
combined_seurat <- merge(
  x = seurat_list[[1]],
  y = unlist(seurat_list[-1]),
  add.cell.ids = names(seurat_list)
)

# 查看合并后的对象
print(combined_seurat)

#——————————————————————————————————————————————————————————
# 二.假设你的全基因 Seurat 对象叫做 combined_seurat####
#——————————————————————————————————————————————————————————
seu_all <- combined_seurat

# 1. 查找所有以 Mmul10- 开头的猴子基因
monkey_genes <- grep("^Mmul10-", rownames(seu_all), value = TRUE)
message("找到 ", length(monkey_genes), " 个猴子基因")

# 2. 从 Seurat 对象中子集出这些基因
seu_monkey <- subset(seu_all, features = monkey_genes)

# 3. 查看结果
print(seu_monkey)

library(ggplot2)
library(patchwork)
library(Matrix)

#——————————————————————————————————————————————————————————
# 三.质控，猴子中线粒体基因命名和人不一样需要查找分析####
#——————————————————————————————————————————————————————————

# 假设你的合并对象叫 seu_monkey，active assay 为 "RNA"
# mt.genes 列表：
mt.genes <- c(
  "Mmul10-ND1","Mmul10-ND2","Mmul10-ND3","Mmul10-ND4","Mmul10-ND4L",
  "Mmul10-ND5","Mmul10-ND6","Mmul10-COX1","Mmul10-COX2","Mmul10-COX3",
  "Mmul10-CYTB","Mmul10-ATP6","Mmul10-ATP8"
)

# 1. 取出 layer 名
rna_assay   <- seu_monkey[["RNA"]]
layer_names <- names(rna_assay@layers)
message("检测到 ", length(layer_names), " 个 layer")

# 为所有细胞初始化一个向量
all_cells  <- colnames(seu_monkey)
percent_mt <- setNames(numeric(length(all_cells)), all_cells)

# 2. 分 layer 计算
for (L in layer_names) {
  message("计算 layer ", L, " …")
  mat <- GetAssayData(seu_monkey, assay = "RNA", layer = L)
  
  # 只保留本 layer 中存在的 mt.genes
  genes_in_layer <- intersect(mt.genes, rownames(mat))
  if (length(genes_in_layer) == 0) {
    # 如果这个 layer 完全没有 mt.genes，就跳过
    next
  }
  
  # 计算 mt counts 和 total counts
  mt_counts    <- colSums(mat[genes_in_layer, , drop = FALSE])
  total_counts <- colSums(mat)
  
  # 写入对应细胞
  cells_this_layer    <- colnames(mat)
  percent_mt[cells_this_layer] <- mt_counts[cells_this_layer] / total_counts[cells_this_layer] * 100
}

# 3. 写回 meta.data
seu_monkey$percent.mt <- percent_mt


#——————————————————————————————————————————————————————————
# 四.QC分别绘图观测####
#——————————————————————————————————————————————————————————
# 1. 小提琴图
vln_plot <- VlnPlot(
  seu_monkey,
  features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
  ncol     = 3,
  pt.size  = 0.1
) + theme(legend.position = "none")
ggsave("VlnPlot_qc_metrics.png", plot = vln_plot, width = 10, height = 4, dpi = 300)

# 2. nCount_RNA vs percent.mt 散点图
scatter1 <- FeatureScatter(seu_monkey, "nCount_RNA", "percent.mt") +
  ggtitle("nCount_RNA vs percent.mt")
ggsave("Scatter_nCount_vs_percentMT.png", plot = scatter1, width = 6, height = 5, dpi = 300)

# 3. nCount_RNA vs nFeature_RNA 散点图
scatter2 <- FeatureScatter(seu_monkey, "nCount_RNA", "nFeature_RNA") +
  ggtitle("nCount_RNA vs nFeature_RNA")
ggsave("Scatter_nCount_vs_nFeature.png", plot = scatter2, width = 6, height = 5, dpi = 300)

#汇总观测
library(patchwork)

p1 <- VlnPlot(seu_monkey, features = "percent.mt", pt.size = 0.1) + NoLegend()
p2 <- VlnPlot(seu_monkey, features = "nFeature_RNA", pt.size = 0.1) + NoLegend()
p3 <- VlnPlot(seu_monkey, features = "nCount_RNA", pt.size = 0.1) + NoLegend()

png("QC_plots/QC_ViolinPlots_combined.png", width = 1800, height = 600)
p1 + p2 + p3
dev.off()
#——————————————————————————————————————————————————————————
# 五.根据前面的质控图Vln_plot，规定一个过滤条件，开始过滤。####
#——————————————————————————————————————————————————————————
seu_qc <- subset(
  seu_monkey,
  subset = nFeature_RNA > 200 &
    nFeature_RNA < 6000 &
    percent.mt  < 5)

cat("过滤前细胞数：", ncol(seu_monkey), "\n")
cat("过滤后细胞数：", ncol(seu_qc), "\n")



# 1. 归一化 (Normalization)
seu_qc <- NormalizeData(
  object = seu_qc,
  normalization.method = "LogNormalize",
  scale.factor = 10000
)


# 2. 找到高变基因 (Feature Selection)
seu_qc <- FindVariableFeatures(
  object = seu_qc,
  selection.method = "vst",
  nfeatures = 2000
)
# 获取top10高变基因
top10 <- head(VariableFeatures(seu_qc), 10)
# 绘制高变基因图并添加标注
plot1 <- VariableFeaturePlot(seu_qc)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
# 保存图片到当前目录
ggsave("highly_variable_genes.png", plot2, width = 8, height = 6, dpi = 300)


# 3. 不必全基因缩放, 对hvg做scale就行(Scaling)
#all.genes <- rownames(seu_qc)
#seu_qc <- ScaleData(
  #object = seu_qc,
  #features = all.genes
#)
seu_qc <- ScaleData(seu_qc, features = VariableFeatures(seu_qc))


# 4. 主成分分析 (PCA)
seu_qc <- RunPCA(
  object = seu_qc,
  features = VariableFeatures(seu_qc),
  npcs = 30,
  verbose = FALSE
)
# 查看前几个主成分
print(seu_qc[["pca"]], dims = 1:5, nfeatures = 5)####
VizDimLoadings(seu_qc, dims = 1:2, reduction = "pca")
DimPlot(seu_qc, reduction = "pca")
# 画图
# 保存 VizDimLoadings（PC1、PC2 载荷条形图）
png("PCA_Loadings_PC1_2.png", width = 1000, height = 800, res = 150)
VizDimLoadings(seu_qc, dims = 1:2, reduction = "pca")
dev.off()
# 保存 PC1 vs PC2 的 DimPlot（关闭 raster，保留矢量点）
png("PCA_DimPlot_PC1_vs_PC2.png", width = 1000, height = 800, res = 150)
DimPlot(seu_qc, reduction = "pca", raster = FALSE) +
  ggtitle("PCA: PC1 vs PC2")
dev.off()
cat("✅ 已保存：\n",
    " - PCA_Loadings_PC1_2.png\n",
    " - PCA_DimPlot_PC1_vs_PC2.png\n")


# 5. 选择用于下游分析的主成分 (PC Elbow Plot / JackStraw)
ElbowPlot(seu_qc, ndims = 30)
# 绘图打开 PNG 设备
png("ElbowPlot_seu_qc.png", width = 1000, height = 800, res = 150)
# 绘制 ElbowPlot
ElbowPlot(seu_qc, ndims = 30)
# 关闭设备
dev.off()
cat("✅ 已保存 ElbowPlot 到文件：ElbowPlot_seu_qc.png\n")


# 6. 构建邻近图 (Find Neighbors)
# 这里假设选用前20个 PC
seu_qc <- FindNeighbors(
  object = seu_qc,
  dims   = 1:20
)


# 7. 聚类 (Clustering)
seu_qc <- FindClusters(
  object = seu_qc,
  resolution = 0.5
)
# 查看聚类标签
head(Idents(seu_qc), 10)


# 8. UMAP 可视化 (非线性降维)

seu_qc <- RunUMAP(
  object = seu_qc,
  dims   = 1:20
)


#——————————————————————————————————————————————————————————
# 六.去批次效应####
#——————————————————————————————————————————————————————————
# 不对劲，去批次,这是对所有的样本做去批次
library(harmony)

# 1. 运行 Harmony 校正，按 sample 分组，起码之前要跑完pca
seu_qc <- RunHarmony(
  object         = seu_qc,
  group.by.vars  = "sample",      # 用 sample 列区分批次
  reduction      = "pca",         # 在 PCA 空间上运行
  reduction.save = "harmony",     # 校正后的嵌入存到 seu_qc[["harmony"]]
  assay.use      = "RNA",
  project.dim    = FALSE
)

# 2. 基于 harmony 嵌入构建邻近图（前 20 个 PC）
seu_qc <- FindNeighbors(
  object    = seu_qc,
  reduction = "harmony",
  dims      = 1:20
)

# 3. 聚类
seu_qc <- FindClusters(
  object     = seu_qc,
  resolution = 0.5
)

# 4. 基于 harmony 嵌入运行 UMAP
seu_qc <- RunUMAP(
  object    = seu_qc,
  reduction = "harmony",
  dims      = 1:20
)

# 5. 可视化并保存结果
# 按 sample 着色
p1 <- DimPlot(
  seu_qc,
  reduction = "umap",
  group.by  = "sample"
) + ggtitle("Harmony-corrected UMAP by Sample")

ggsave("UMAP_Harmony_by_sample.png", plot = p1, width = 8, height = 6, dpi = 300)

# 按聚类着色并加标签
p2 <- DimPlot(
  seu_qc,
  reduction = "umap",
  label     = TRUE
) + ggtitle("Harmony-corrected UMAP Clusters")

ggsave("UMAP_Harmony_Clusters.png", plot = p2, width = 8, height = 6, dpi = 300)

cat("✅ 已完成 Harmony 校正并保存 UMAP 图：\n",
    " - UMAP_Harmony_by_sample.png\n",
    " - UMAP_Harmony_Clusters.png\n")
saveRDS(seu_qc, file = "sample_umap.rds")


#——————————————————————————————————————————————————————————
# 七. 差异表达marker####
#——————————————————————————————————————————————————————————
# 0. 确保数据层已合并（关键步骤）
setwd("/data5/GPT/Wuky/1LQ_scRNA/LQdata/h5d/")

seu_qc <- JoinLayers(seu_qc)

# 1. 寻找每个簇的标记基因
Idents(seu_qc) <- "seurat_clusters"  # 设置默认标识为聚类结果

cluster_markers <- FindAllMarkers(
  object = seu_qc,
  only.pos = TRUE,      # 只保留在簇中高表达的基因
  min.pct = 0.25,       # 基因至少在25%的细胞中表达
  logfc.threshold = 0.5, # 最小log2FC=0.5
  test.use = "wilcox",  # 非参数检验
  random.seed = 42,     # 确保可重复性
  verbose = TRUE        # 显示进度
)

write.csv(cluster_markers, "raw_cluster_markers.csv", row.names = FALSE)
cat("原始标记基因结果已保存到 raw_cluster_markers.csv\n")

# 筛选显著标记基因
library(dplyr)
top_markers <- cluster_markers %>%
  group_by(cluster) %>%
  filter(p_val_adj < 0.05) %>%
  slice_max(avg_log2FC, n = 10)  # 每个簇取前10个标记基因

# 保存结果
write.csv(top_markers, "cluster_top_markers.csv", row.names = FALSE)


# 4.2 保存注释后的Seurat对象
saveRDS(seu_qc, "annotated_seurat_object.rds")


#——————————————————————————————————————————————————————————
# 九. 利用top_marker对每一簇做注释，可以先做大注释####
#——————————————————————————————————————————————————————————
library(Seurat)
library(dplyr)
library(ggplot2)

setwd("/data5/GPT/Wuky/1LQ_scRNA/LQdata/h5d/")

# 1. 读取Seurat对象（含注释）
seu_qc <- readRDS("annotated_seurat_object.rds")

# 2. 读取标记基因结果
cluster_markers <- read.csv("raw_cluster_markers.csv")
top_markers <- read.csv("cluster_top_markers.csv")

library(Seurat)
library(ggplot2)

gene_dict <- c(
  "HBB"="Erythrocytes","HBA"="Erythrocytes",
  "MMP9"="Neutrophils","ICOS"="T_follicular_helper",
  "CD19"="B_cells","MS4A1"="B_cells",
  "CCR7"="Naive_CD4_T",
  "CD8B"="CD8_T","GZMH"="CD8_T",            # ← CD8 专属
  "CD4"="CD4_T","IL7R"="CD4_T",             # ← CD4 专属
  "GZMK"="Effector_memory_CD8_T",
  "CD4" ="Effector_memory_CD4_T",
  "RETN"="Monocytes","MS4A7"="Macrophages",
  "CCL5"="NK_cells","KLRC3"="Activated_NK_cells",
  "CD163"="M2_Macrophages","BANK1"="Germinal_center_B_cells",
  "S100A8"="Immature_Neutrophils","TCL1A"="Plasmablasts",
  "GP9"="Platelets","CXCL8"="Activated_Neutrophils",
  "UBE2C"="Cycling_cells",
  "AZU1"="Tissue_resident_Macrophages",
  "KLRG1"="Effector_memory_CD8_T",
  "LILRA4"="pDC","PF4V1"="Activated_Platelets",
  "CPA3"="Mast_cells","IL1RAP"="cDC","PAX5"="Follicular_B_cells"
)
## 2. 读入文件（base R）
top <- read.csv("cluster_top_markers.csv", stringsAsFactors = FALSE)

## 3. 去掉前缀并匹配字典
top$gene_clean <- sub("^Mmul10-", "", top$gene)
keep_idx <- top$gene_clean %in% names(gene_dict)
plot_df <- top[keep_idx, ]
plot_df$class <- gene_dict[plot_df$gene_clean]
plot_df$cluster <- factor(plot_df$cluster,
                          levels = sort(unique(plot_df$cluster)))

## 4. 画图
p <- ggplot(plot_df,
            aes(x = class, y = cluster,
                size = pct.1, color = avg_log2FC)) +
  geom_point(alpha = .9) +
  scale_color_viridis_c(option = "plasma") +
  scale_size_continuous(range = c(0.5, 3)) +
  labs(title = "Representative markers across PBMC clusters",
       x = "Cell class", y = "Cluster",
       size = "Pct1", color = "Avg log2FC") +
  theme_minimal(base_size = 14) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

## 5. 保存
png("bubble_PBMC_classes_representative.png",
    width = 3000, height = 1600, res = 300)
print(p)
dev.off()



## 0) 确保 cluster 列存在且为字符
cell_clusters <- as.character(seu_qc$seurat_clusters)

## 1) 构建完整映射表（0–28 全覆盖）
cluster_annotation <- setNames(
  c("Erythrocytes",
    "CD8_T_Eff",            # 1
    "Neutrophils",          # 2
    "CD4_T",                # 3
    "B_cells",              # 4
    "CD4_T",                # 5
    "CD8_T_Eff",            # 6
    "Macrophages",          # 7
    "Macrophages",          # 8
    "NK_cells",             # 9
    "CD8_T_Eff",            # 10
    "Macrophages",          # 11
    "B_cells",              # 12
    "Neutrophils",          # 13
    "B_cells",              # 14
    "CD4_T",                # 15
    "Platelets",            # 16
    "Classical monocytes",  # 17
    "B_cells",              # 18
    "Cycling_cells",        # 19
    "Non-classical monocytes", # 20
    "CD8_T_Eff",            # 21
    "Erythrocytes",         # 22
    "DCs",                  # 23
    "Platelets",            # 24
    "Macrophages",          # 25
    "DCs",                  # 26
    "DCs",                  # 27
    "B_cells"),             # 28
  nm = as.character(0:28)
)

## 2) 生成注释向量
cell_type_vec <- factor(cluster_annotation[cell_clusters], levels = unique(cluster_annotation))

## 3) 构建与 Seurat 行名完全一致的数据框
anno_df <- data.frame(
  cell_type = cell_type_vec,
  row.names = colnames(seu_qc)       # 关键：行名对齐 barcode
)

## 4) 用 AddMetaData 合并（永远不会报 overlap 错）
seu_qc <- AddMetaData(seu_qc, metadata = anno_df)

## 5) 校验
stopifnot(sum(is.na(seu_qc$cell_type)) == 0)
table(seu_qc$cell_type)   # 正常显示各类细胞数

saveRDS(seu_qc, file = "seu_qc_annotated.rds")

#——————————————————————————————————————————————————————————
# 十. 下次开始分析直接读取注释的seurat对象####
#——————————————————————————————————————————————————————————
library(Seurat)
library(dplyr)
library(ggplot2)

setwd("/data5/GPT/Wuky/1LQ_scRNA/LQdata/h5d/")

# 1. 读取Seurat对象（含注释）
seu_qc <- readRDS("seu_qc_annotated.rds")

# 2) 输出路径
dir.create("results", showWarnings = FALSE)
png_out <- "results/umap_cell_type_small.png"

# 3) 画图 & 保存
png(filename = png_out,
    width    = 1000,
    height   = 800,
    res      = 150)

DimPlot(seu_qc,
        reduction = "umap",
        group.by  = "cell_type",
        label     = TRUE,
        label.size = 3,
        pt.size   = 0.2,   # 点更小，防止重叠
        raster    = FALSE) +  # 不栅格化，更清晰
  NoAxes() +
  ggtitle("UMAP — cell type (simplified)") +
  theme(plot.title = element_text(hjust = 0.5, size = 14))

dev.off()

message("✅ UMAP 图已保存：", png_out)

#—————————————————————————————————————————————————————————
# 十一. 分组做umap：读取注释的seurat对象####
# 10:03 2025/8/2 tested by wky
#—————————————————————————————————————————————————————————
library(Seurat)
library(dplyr)
library(ggplot2)

setwd("/data5/GPT/Wuky/1LQ_scRNA/LQdata/h5d/")

# 1. 读取Seurat对象（含注释）
seu_qc <- readRDS("seu_qc_annotated.rds")
# 检查
colnames(seu_qc@meta.data)
unique(seu_qc$orig.ident)

# 1. 创建分组列
seu_qc$condition_group <- factor(
  case_when(
    grepl("T0_C", seu_qc$orig.ident) ~ "T0C",
    grepl("T20_C", seu_qc$orig.ident) ~ "T20C",
    grepl("T0_GD", seu_qc$orig.ident) ~ "T0GDT",
    grepl("T20_GD", seu_qc$orig.ident) ~ "T20GDT"
  ),
  levels = c("T0C", "T20C", "T0GDT", "T20GDT")
)
# 设置condition_group的顺序
seu_qc$condition_group <- factor(seu_qc$condition_group,
                                 levels = c("T0C", "T20C", "T0GDT", "T20GDT"))

# 定义颜色方案（使用Seurat默认颜色）
celltype_colors <- scales::hue_pal()(length(unique(seu_qc$cell_type)))
names(celltype_colors) <- unique(seu_qc$cell_type)

# 获取整个数据集的UMAP坐标范围
umap_coords <- Embeddings(seu_qc, "umap")
x_range <- range(umap_coords[, 1])
y_range <- range(umap_coords[, 2])

# 设置绘图主题
my_theme <- theme(
  plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
  legend.position = "none",  # 隐藏单个图中的图例
  axis.title = element_blank(),
  axis.text = element_blank(),
  axis.ticks = element_blank(),
  panel.background = element_blank(),
  panel.border = element_rect(colour = "black", fill = NA, size = 0.5)
)

# 创建存储图的列表
plot_list <- list()

# 为每个分组创建UMAP图
for (group in levels(seu_qc$condition_group)) {
  # 创建子集
  group_cells <- colnames(seu_qc)[seu_qc$condition_group == group]
  group_seu <- subset(seu_qc, cells = group_cells)
  
  # 绘制UMAP
  p <- DimPlot(group_seu, reduction = 'umap', 
               group.by = "cell_type",
               cols = celltype_colors,
               pt.size = 0.3) +
    xlim(x_range) + ylim(y_range) +  # 统一坐标范围
    ggtitle(group) +
    my_theme
  
  plot_list[[group]] <- p
}

# 提取共享图例
legend_plot <- DimPlot(seu_qc, reduction = 'umap', 
                       group.by = "cell_type",
                       cols = celltype_colors,
                       pt.size = 0.3) +
  theme(legend.position = "bottom",
        legend.text = element_text(size = 10),
        legend.title = element_blank())

legend <- cowplot::get_legend(legend_plot)

# 组合四个分组图
combined_plot <- (plot_list[["T0C"]] + plot_list[["T20C"]]) / 
  (plot_list[["T0GDT"]] + plot_list[["T20GDT"]]) / 
  legend +
  plot_layout(heights = c(1, 1, 0.1))  # 调整图例高度

# 添加主标题
combined_plot <- combined_plot + 
  plot_annotation(title = "UMAP Visualization by Experimental Group",
                  theme = theme(plot.title = element_text(size = 16, hjust = 0.5, face = "bold")))

# 保存结果
dir.create("results", showWarnings = FALSE)

# 保存组合图
ggsave("results/four_group_umap_plots.pdf", combined_plot, 
       width = 12, height = 10, dpi = 300)
ggsave("results/four_group_umap_plots.png", combined_plot, 
       width = 12, height = 10, dpi = 300)

# 单独保存每个分组图
for (group in levels(seu_qc$condition_group)) {
  # 创建包含图例的单独图
  p_single <- plot_list[[group]] + 
    theme(legend.position = "right") +  # 在右侧添加图例
    guides(color = guide_legend(override.aes = list(size = 3)))  # 增大图例点大小
  
  ggsave(paste0("results/umap_", group, ".png"), 
         p_single,
         width = 7, height = 5, dpi = 300)
}

message("✅ 四分组UMAP图已保存到results目录")

#—————————————————————————————————————————————————————————
# 十二。绘制所有样本的umap图以及条形图####
#—————————————————————————————————————————————————————————
library(Seurat)
library(ggplot2)
library(patchwork)
library(dplyr)

# 1. 准备样本信息
# 获取样本列表和细胞数量
sample_counts <- table(seu_qc$sample)
sample_levels <- names(sample_counts)

# 创建样本信息数据框
sample_info <- data.frame(
  sample = names(sample_counts),
  cell_count = as.numeric(sample_counts),
  time = ifelse(grepl("T0", names(sample_counts)), "T0", "T20"),
  treatment = ifelse(grepl("_C", names(sample_counts)), "Control", "GDT")
)

# 添加样本分组信息
sample_info$group <- paste0(sample_info$time, sample_info$treatment)

# 2. 设置绘图参数
# 定义颜色方案
group_colors <- c(
  "T0Control" = "#1f77b4",
  "T20Control" = "#aec7e8",
  "T0GDT" = "#ff7f0e",
  "T20GDT" = "#ffbb78"
)

# 获取UMAP坐标范围
umap_coords <- Embeddings(seu_qc, "umap")
x_range <- range(umap_coords[, 1])
y_range <- range(umap_coords[, 2])

# 设置绘图主题
my_theme <- theme(
  plot.title = element_text(hjust = 0.5, size = 10, face = "bold"),
  legend.position = "none",
  axis.title = element_blank(),
  axis.text = element_blank(),
  axis.ticks = element_blank(),
  panel.background = element_blank(),
  panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.5)
)

# 3. 绘制每个样本的UMAP图
plot_list <- list()

for (sample in sample_levels) {
  # 获取样本信息
  sample_info_row <- sample_info[sample_info$sample == sample, ]
  group <- paste0(sample_info_row$time, sample_info_row$treatment)
  
  # 创建子集
  sample_cells <- colnames(seu_qc)[seu_qc$sample == sample]
  sample_seu <- subset(seu_qc, cells = sample_cells)
  
  # 绘制UMAP
  p <- DimPlot(sample_seu, reduction = 'umap', 
               group.by = "cell_type",
               pt.size = 0.1) +
    xlim(x_range) + ylim(y_range) +
    ggtitle(paste0(sample, "\n(n = ", sample_info_row$cell_count, ")")) +
    my_theme +
    theme(plot.title = element_text(color = group_colors[group]))
  
  plot_list[[sample]] <- p
}

# 4. 组合所有样本图
# 创建4x4网格布局
combined_plot <- wrap_plots(plot_list, ncol = 4) +
  plot_annotation(title = "UMAP Visualization by Sample",
                  subtitle = "Color indicates group: T0Control (blue), T20Control (light blue), T0GDT (orange), T20GDT (light orange)",
                  theme = theme(plot.title = element_text(size = 16, hjust = 0.5, face = "bold"),
                                plot.subtitle = element_text(size = 12, hjust = 0.5)))

# 5. 按时间点和处理组分组比较

# 创建分组绘图函数
plot_group_comparison <- function(group_var, title) {
  p <- DimPlot(seu_qc, reduction = 'umap', 
               group.by = group_var,
               split.by = "sample",
               ncol = 4,
               pt.size = 0.05) +
    ggtitle(title) +
    theme(plot.title = element_text(size = 14, hjust = 0.5, face = "bold"))
  return(p)
}

# 按细胞类型分组比较
p_celltype <- plot_group_comparison("cell_type", "Cell Type Distribution by Sample")

# 按时间点分组比较
seu_qc$time_point <- ifelse(grepl("T0", seu_qc$orig.ident), "T0", "T20")
p_time <- plot_group_comparison("time_point", "Time Point Distribution by Sample")

# 按处理组分组比较
seu_qc$treatment <- ifelse(grepl("_C", seu_qc$orig.ident), "Control", "GDT")
p_treatment <- plot_group_comparison("treatment", "Treatment Group Distribution by Sample")

# 6. 保存结果
dir.create("results", showWarnings = FALSE)

# 保存样本级UMAP图
ggsave("results/sample_level_umap_plots.pdf", combined_plot, 
       width = 18, height = 16, dpi = 300)
ggsave("results/sample_level_umap_plots.png", combined_plot, 
       width = 18, height = 16, dpi = 300)

# 8. 生成样本信息可视化
# 样本细胞数量条形图
p_count <- ggplot(sample_info, aes(x = sample, y = cell_count, fill = group)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = cell_count), vjust = -0.3, size = 3) +
  scale_fill_manual(values = group_colors) +
  labs(title = "Cell Count by Sample", x = "Sample", y = "Number of Cells") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title = element_text(hjust = 0.5, size = 14, face = "bold"))

# 按组汇总细胞数量
group_summary <- sample_info %>%
  group_by(group) %>%
  summarise(total_cells = sum(cell_count),
            mean_cells = mean(cell_count),
            sd_cells = sd(cell_count))

# 组间比较条形图
p_group <- ggplot(group_summary, aes(x = group, y = total_cells, fill = group)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = total_cells), vjust = -0.3, size = 5) +
  scale_fill_manual(values = group_colors) +
  labs(title = "Total Cell Count by Experimental Group", 
       x = "Group", y = "Total Number of Cells") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"))

# 保存样本信息图
ggsave("results/sample_cell_counts.png", p_count, width = 10, height = 6, dpi = 300)
ggsave("results/group_cell_counts.png", p_group, width = 8, height = 6, dpi = 300)

message("✅ 所有分析结果已保存到results目录")

#—————————————————————————————————————————————————————————
# 十三.细胞亚群比例####
#—————————————————————————————————————————————————————————
library(Seurat)
library(ggplot2)
library(dplyr)
library(tidyr)
library(RColorBrewer)

# 1. 创建分组变量
seu_qc$time_point <- ifelse(grepl("T0", seu_qc$orig.ident), "T0", "T20")
seu_qc$treatment <- ifelse(grepl("_C", seu_qc$orig.ident), "Control", "GDT")
seu_qc$exp_group <- paste0(seu_qc$time_point, "_", seu_qc$treatment)

# 2. 计算细胞类型比例
# 样本级别比例计算
sample_prop <- seu_qc@meta.data %>%
  group_by(sample, cell_type) %>%
  summarise(count = n(), .groups = "drop") %>%  # 明确指定分组输出
  group_by(sample) %>%
  mutate(proportion = count / sum(count)) %>%
  ungroup()

# 添加样本信息
sample_info <- seu_qc@meta.data %>%
  distinct(sample, time_point, treatment, exp_group)

sample_prop <- left_join(sample_prop, sample_info, by = "sample")

# 分组级别比例计算
group_prop <- seu_qc@meta.data %>%
  group_by(exp_group, cell_type) %>%
  summarise(count = n(), .groups = "drop") %>%  # 明确指定分组输出
  group_by(exp_group) %>%
  mutate(proportion = count / sum(count)) %>%
  ungroup()

# 3. 创建自定义颜色方案
# 获取所有细胞类型
cell_types <- unique(seu_qc$cell_type)
n_cell_types <- length(cell_types)

# 创建扩展颜色方案
get_palette <- colorRampPalette(brewer.pal(9, "Set1"))
cell_colors <- get_palette(n_cell_types)
names(cell_colors) <- cell_types

# 4. 绘制16个样本的细胞亚群比例图
p_sample <- ggplot(sample_prop, 
                   aes(x = sample, y = proportion, fill = cell_type)) +
  geom_bar(stat = "identity", position = "stack", width = 0.8) +
  scale_fill_manual(values = cell_colors) +
  labs(title = "Cell Type Proportions by Sample",
       x = "Sample", 
       y = "Proportion of Cells",
       fill = "Cell Type") +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
    axis.title = element_text(size = 14),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
    axis.text.y = element_text(size = 10),
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 10),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank()
  ) +
  geom_text(
    data = sample_prop %>% distinct(sample, time_point, treatment),
    aes(x = sample, y = 1.05, label = paste("T:", time_point, "\nTx:", treatment)),
    size = 3, vjust = 0.5, hjust = 0.5, inherit.aes = FALSE
  ) +
  geom_text(
    data = sample_prop %>% group_by(sample) %>% summarise(total = sum(count), .groups = "drop"),
    aes(x = sample, y = -0.05, label = paste("n =", total)),
    size = 3, vjust = 1, hjust = 0.5, inherit.aes = FALSE
  ) +
  scale_y_continuous(expand = expansion(mult = c(0.1, 0.1)))

# 5. 绘制4个分组的细胞亚群比例图
p_group <- ggplot(group_prop, 
                  aes(x = exp_group, y = proportion, fill = cell_type)) +
  geom_bar(stat = "identity", position = "stack", width = 0.6) +
  scale_fill_manual(values = cell_colors) +
  labs(title = "Cell Type Proportions by Experimental Group",
       x = "Experimental Group", 
       y = "Proportion of Cells",
       fill = "Cell Type") +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
    axis.title = element_text(size = 14),
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 10),
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 10),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank()
  ) +
  geom_text(
    data = group_prop %>% group_by(exp_group) %>% summarise(total = sum(count), .groups = "drop"),
    aes(x = exp_group, y = -0.05, label = paste("n =", total)),
    size = 4, vjust = 1, hjust = 0.5, inherit.aes = FALSE
  ) +
  geom_text(
    aes(label = ifelse(proportion > 0.05, paste0(round(proportion*100, 1), "%"), "")),
    position = position_stack(vjust = 0.5),
    size = 3, color = "white", fontface = "bold"
  ) +
  scale_y_continuous(expand = expansion(mult = c(0.1, 0.1)))

# 6. 保存结果
dir.create("results", showWarnings = FALSE)

# 保存样本比例图
ggsave("results/cell_proportions_by_sample.png", p_sample, 
       width = 16, height = 8, dpi = 300)
ggsave("results/cell_proportions_by_sample.pdf", p_sample, 
       width = 16, height = 8, dpi = 300)

# 保存分组比例图
ggsave("results/cell_proportions_by_group.png", p_group, 
       width = 10, height = 8, dpi = 300)
ggsave("results/cell_proportions_by_group.pdf", p_group, 
       width = 10, height = 8, dpi = 300)

# 7. 保存数据表格
write.csv(sample_prop, "results/cell_proportions_by_sample.csv", row.names = FALSE)
write.csv(group_prop, "results/cell_proportions_by_group.csv", row.names = FALSE)

message("✅ 分析完成！所有结果已保存到results目录")

#—————————————————————————————————————————————————————————
# 十四，选取感兴趣的细胞如中性粒细胞的差异分析####
#—————————————————————————————————————————————————————————
# 加载Seurat包
library(Seurat)

# 步骤1: 提取中性粒细胞子集
neutrophils <- subset(seu_qc, subset = cell_type == "Neutrophils")

# 步骤2: 创建分组变量
# 从sample列提取处理组和时间信息
sample_info <- neutrophils@meta.data$sample
time_point <- sapply(strsplit(sample_info, "_"), `[`, 1)  # 提取T0/T20
treatment <- substr(sapply(strsplit(sample_info, "_"), `[`, 2), 1, 1)  # 提取C/G

# 创建四组标签
neutrophils@meta.data$group <- paste0(time_point, "_", ifelse(treatment == "G", "GDT", "C"))
Idents(neutrophils) <- "group"  # 设置分组标识

# 检查各组细胞数量
table(neutrophils@meta.data$group)

# 加载必要的包
library(Seurat)
library(ggplot2)
library(dplyr)
library(ggrepel)
library(viridis)

# 步骤3: 差异分析（核心比较）
# 主比较：T20_GDT vs T20_C
markers_T20vsC <- FindMarkers(
  neutrophils,
  ident.1 = "T20_GDT",
  ident.2 = "T20_C",
  min.pct = 0.1,
  logfc.threshold = 0.25,
  test.use = "wilcox"
)

# 辅助比较
markers_GDT_time <- FindMarkers(
  neutrophils,
  ident.1 = "T20_GDT",
  ident.2 = "T0_GDT",
  min.pct = 0.1
)

markers_C_time <- FindMarkers(
  neutrophils,
  ident.1 = "T20_C",
  ident.2 = "T0_C",
  min.pct = 0.1
)

# 步骤4: 结果筛选与注释
sig_genes <- markers_T20vsC %>%
  filter(p_val_adj < 0.05 & abs(avg_log2FC) > 0.5) %>%
  arrange(desc(avg_log2FC))

# 定义中性粒细胞标记基因
neutrophil_markers <- c(
  "Mmul10-S100A8", "Mmul10-S100A9", 
  "Mmul10-FCGR3B", "Mmul10-CSF3R", 
  "Mmul10-CEACAM8"
)

# 标记关键基因
sig_genes$is_neut_marker <- rownames(sig_genes) %in% neutrophil_markers
sig_genes$gene_symbol <- gsub("Mmul10-", "", rownames(sig_genes))

# 步骤5: 可视化并保存

## 1. 专业火山图
volcano_plot <- ggplot(sig_genes, aes(x = avg_log2FC, y = -log10(p_val_adj))) +
  geom_point(aes(color = is_neut_marker, size = is_neut_marker), alpha = 0.7) +
  scale_color_manual(
    name = "Gene Type",
    values = c("FALSE" = "gray60", "TRUE" = "red"),
    labels = c("Other DEGs", "Neutrophil Markers")
  ) +
  scale_size_manual(
    values = c(2, 3),
    guide = "none"
  ) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "blue") +
  geom_vline(xintercept = c(-0.5, 0.5), linetype = "dashed", color = "blue") +
  ggrepel::geom_text_repel(
    data = subset(sig_genes, is_neut_marker | abs(avg_log2FC) > 2),
    aes(label = gene_symbol),
    size = 3,
    max.overlaps = 20,
    box.padding = 0.5,
    segment.color = "grey50"
  ) +
  labs(
    title = "Differential Expression: T20_GDT vs T20_C",
    subtitle = "Neutrophil Population",
    x = "Log2 Fold Change",
    y = "-Log10(Adjusted p-value)",
    caption = paste("Total DEGs:", nrow(sig_genes), 
                    "| Neutrophil Markers:", sum(sig_genes$is_neut_marker))
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5),
    legend.position = "top"
  )

# 保存火山图
ggsave("Neutrophils_Volcano.pdf", plot = volcano_plot, 
       width = 8, height = 6, device = cairo_pdf)
ggsave("Neutrophils_Volcano.png", plot = volcano_plot, 
       width = 8, height = 6, dpi = 300, bg = "white")

## 2. 热图（关键差异基因表达）
# 准备数据
top_genes <- head(rownames(sig_genes), 20)
neutrophils <- ScaleData(neutrophils, features = top_genes)

# 创建分组颜色
group_colors <- c(
  "T0_C" = "#1f77b4", 
  "T0_GDT" = "#ff7f0e", 
  "T20_C" = "#2ca02c", 
  "T20_GDT" = "#d62728"
)

# 创建热图
heatmap_plot <- DoHeatmap(
  neutrophils,
  features = top_genes,
  group.by = "group",
  slot = "scale.data",
  disp.min = -2,
  disp.max = 2,
  angle = 45,
  size = 3.5,
  raster = FALSE
) +
  scale_fill_gradientn(
    colors = c("navy", "white", "firebrick"),
    breaks = c(-2, 0, 2),
    labels = c("Low", "Medium", "High"),
    name = "Expression\n(Z-score)"
  ) +
  labs(title = "Top 20 Differentially Expressed Genes",
       subtitle = "T20_GDT vs T20_C comparison") +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5),
    axis.text.y = element_text(face = "italic")
  )

# 保存热图
ggsave("Neutrophils_Heatmap.pdf", plot = heatmap_plot, 
       width = 10, height = 8, device = cairo_pdf)
ggsave("Neutrophils_Heatmap.png", plot = heatmap_plot, 
       width = 10, height = 8, dpi = 300, bg = "white")

## 3. 关键基因表达验证图
# 创建组合图
key_genes <- neutrophil_markers[neutrophil_markers %in% rownames(neutrophils)]
if(length(key_genes) > 0) {
  vln_plots <- VlnPlot(
    neutrophils,
    features = key_genes,
    group.by = "group",
    pt.size = 0,
    ncol = min(3, length(key_genes)),
    cols = group_colors
  ) +
    theme(plot.title = element_text(face = "bold", hjust = 0.5))
  
  # 保存关键基因图
  ggsave("Neutrophils_KeyMarkers.pdf", plot = vln_plots, 
         width = 8, height = 6, device = cairo_pdf)
  ggsave("Neutrophils_KeyMarkers.png", plot = vln_plots, 
         width = 8, height = 6, dpi = 300, bg = "white")
} else {
  message("No key neutrophil markers found in the dataset.")
}

## 4. 保存结果到CSV
write.csv(sig_genes, "Neutrophils_DEG_Results.csv", row.names = TRUE)
