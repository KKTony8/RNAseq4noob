library(monocle3)
# 1.加载数据，创建cds对象####
expression_matrix <- readRDS(url("https://depts.washington.edu:/trapnell-lab/software/monocle3/celegans/data/packer_embryo_expression.rds"))
cell_metadata <- readRDS(url("https://depts.washington.edu:/trapnell-lab/software/monocle3/celegans/data/packer_embryo_colData.rds"))
gene_annotation <- readRDS(url("https://depts.washington.edu:/trapnell-lab/software/monocle3/celegans/data/packer_embryo_rowData.rds"))

# 创建CDS对象
cds <- new_cell_data_set(expression_matrix,
                         cell_metadata = cell_metadata,
                         gene_metadata = gene_annotation)


# 2.但是更多时候我们是从已经处理好的 Seurat 对象构建CDS对象####
# srt是之前注释好的Seurat 对象
cds <- monocle3::new_cell_data_set(
  expression_data = SeuratObject::GetAssayData(srt, assay = "RNA", layer = "counts"),
  cell_metadata = srt@meta.data,
  gene_metadata = data.frame(gene_short_name = rownames(srt), row.names = rownames(srt))
)

# 3.数据预处理对数据进行归一化（normalization）####
# PCA 降维（这里保留前 50 个主成分）
cds <- preprocess_cds(cds, num_dim = 50)
# 批次校正
cds <- align_cds(cds, alignment_group = "batch")

# 4.降维，强烈建议使用UMAP方法（默认方法）####
cds <- reduce_dimension(cds)

#如果seurat对象已经做完umap降维可以使用该代码处理
#——————————————————————————————————————————————————————————————
cds.embed <- cds@int_colData$reducedDims$UMAP
int.embed <- Embeddings(srt, reduction = "umap")
int.embed <- int.embed[rownames(cds.embed), ]
cds@int_colData$reducedDims$UMAP <- int.embed
#————————————————————————————————————————————————————————————————
# 5.按分区可视化####
cds <- reduce_dimension(cds, reduction_method = "UMAP")  # 降维
cds <- cluster_cells(cds, reduction_method = "UMAP")     # ✅ 聚类
plot_cells(cds, color_cells_by = "partition")

# 学习轨迹图
cds <- learn_graph(cds) 
plot_cells(cds, label_groups_by_cluster = FALSE, color_cells_by = "cell.type")


# 6. 按伪时间排序细胞####
# cds <- order_cells(cds)
# 按伪时间可视化
# plot_cells(cds, color_cells_by = "pseudotime", ...)
# 编程方式指定根节点
# ------------------------------------------
# 🧬 定义函数：获取轨迹起点节点（根据早期时间点）
# ------------------------------------------
get_earliest_principal_node <- function(cds, time_bin = "130-170") {
  # 找出早期时间点的细胞（比如 embryo.time.bin == "130-170"）
  early_cells <- which(colData(cds)$embryo.time.bin == time_bin)
  
  # 获取每个细胞最接近的轨迹图顶点（principal graph node）
  closest_vertex <- cds@principal_graph_aux[["UMAP"]]$pr_graph_cell_proj_closest_vertex
  closest_vertex <- as.matrix(closest_vertex[colnames(cds), ])
  
  # 在早期细胞中，找出最常对应的轨迹节点作为起点
  root_pr_node <- igraph::V(principal_graph(cds)[["UMAP"]])$name[
    as.numeric(
      names(which.max(table(closest_vertex[early_cells, ])))
    )
  ]
  
  return(root_pr_node)
}

# ------------------------------------------
# 🚀 设置轨迹起点，并排序伪时间
# ------------------------------------------
cds <- order_cells(cds, root_pr_nodes = get_earliest_principal_node(cds, time_bin = "130-170"))

# ------------------------------------------
# 🎨 可视化拟时间轨迹
# ------------------------------------------
plot_cells(cds, 
           color_cells_by = "pseudotime", 
           label_groups_by_cluster = FALSE,
           label_leaves = FALSE,
           label_branch_points = FALSE)

# 7. 按伪时间排序细胞####
genes <- monocle3::graph_test(cds, neighbor_graph = "principal_graph", reduction_method = "UMAP", cores = 32)
top10 <- genes %>%
  top_n(n = 10, morans_I) %>%
  pull(gene_short_name) %>%
  as.character()

# 排序伪时间（只做一次）
root_node <- get_earliest_principal_node(cds, time_bin = "130-170")
cds <- order_cells(cds, root_pr_nodes = root_node)

# 提取 AFD 基因与细胞
AFD_genes <- c("gcy-8", "dac-1", "oig-8")
gene_idx <- which(!is.na(rowData(cds)$gene_short_name) & rowData(cds)$gene_short_name %in% AFD_genes)
cell_idx <- which(!is.na(colData(cds)$cell.type) & colData(cds)$cell.type == "AFD")
AFD_lineage_cds <- cds[gene_idx, cell_idx]

# 绘图
plot_genes_in_pseudotime(AFD_lineage_cds,
                         color_cells_by = "embryo.time.bin",
                         min_expr = 0.5)
library(dplyr)

# 8.找伪时间相关基因####
genes <- graph_test(cds, neighbor_graph = "principal_graph", reduction_method = "UMAP")

# 提取 top 50
top50 <- genes %>%
  top_n(n = 50, morans_I) %>%
  pull(gene_short_name) %>%
  as.character()

# 提取表达矩阵
mat <- exprs(cds[top50, ])  # 或 log-normalized matrix

# 聚类
ck <- clusterData(mat, cluster.method = "kmeans", cluster.num = 5)

# 可视化
pdf('monocle3.pdf', height = 10, width = 8)
visCluster(object = ck,
           plot.type = "both",
           add.sampleanno = FALSE,
           markGenes = sample(rownames(mat), 30, replace = FALSE))
dev.off()