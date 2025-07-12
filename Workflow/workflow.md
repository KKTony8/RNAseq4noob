17:20 2025/7/12 wky

# 下游常规分析流程

## 1. 细胞聚类  
```r
# 使用 Seurat 举例
sce <- RunPCA(sce, features = VariableFeatures(sce))
sce <- RunUMAP(sce, dims = 1:20)           # 或 RunTSNE()
DimPlot(sce, reduction = "umap", label = TRUE)
```

## 2. 聚类结果注释  
```r
# 找 marker genes
markers <- FindAllMarkers(sce, only.pos = TRUE, min.pct = 0.25)
# 查看常见标记
library(dplyr)
markers %>%
  filter(gene %in% c("CD3E", "MS4A1")) %>%
  arrange(cluster, desc(avg_log2FC))
# 根据 marker 列表给每个 cluster 命名
new.labels <- c("0"="T cell", "1"="B cell", "2"="Myeloid", "3"="NK cell")
sce$celltype <- plyr::mapvalues(sce$seurat_clusters,
                                from = names(new.labels),
                                to   = new.labels)
DimPlot(sce, group.by = "celltype", label = TRUE)
```

## 3. 功能富集分析  

### 3.1 差异表达分析  
```r
# 以 cluster 0 vs 1 为例
deg <- FindMarkers(sce, ident.1 = 0, ident.2 = 1)
```

### 3.2 GO/KEGG 富集  
```r
library(clusterProfiler)
ego <- enrichGO(gene          = rownames(deg)[deg$p_val_adj < 0.05],
                OrgDb         = org.Hs.eg.db,
                ont           = "BP",
                pAdjustMethod = "BH")
dotplot(ego)
```

### 3.3 GSEA  
```r
library(clusterProfiler)
geneList   <- deg$avg_log2FC
names(geneList) <- rownames(deg)
gsea.res   <- GSEA(geneList, TERM2GENE = kegg_sets)
gseaplot2(gsea.res, geneSetID = "hsa04010")
```

### 3.4 GSVA  
```r
library(GSVA)
expr       <- as.matrix(GetAssayData(sce, slot = "data"))
gsva.score <- gsva(expr, kegg_sets, method = "gsva")
pheatmap(gsva.score)
```

### 3.5 代谢通路推断  
```r
# 使用 scMetabolism 举例
library(scMetabolism)
met.res <- sc_metabolism(expr, method = "VISION")
```

## 4. 伪时序分析  
```r
library(monocle3)
cds <- new_cell_data_set(expr,
                         cell_metadata = sce@meta.data,
                         gene_metadata = rowData(sce))
cds <- preprocess_cds(cds)
cds <- learn_graph(cds)
plot_cells(cds, color_cells_by = "pseudotime")
```

## 5. 细胞–细胞通讯  
```r
library(CellChat)
cellchat <- createCellChat(object = sce, group.by = "celltype")
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- computeCommunProb(cellchat)
netVisual_circle(cellchat@net$weight)
```

## 6. VDJ 分析  

### 6.1 克隆型分析  
```r
library(scRepertoire)
contig     <- read.csv("vdj_contig.csv")
clonotype  <- combineTCR(contig, samples = "sample1")
clonalHomeostasis(clonotype)
```

### 6.2 V/D/J 片段使用  
```r
usage <- segmentUsage(contig, chain = "TRB")
with(usage, barplot(frequency ~ segment, las = 2))
```

### 6.3 CDR3 长度多样性  
```r
diversity <- alphaDiversity(clonotype,
                            group  = "celltype",
                            method = "shannon")
boxplot(diversity ~ celltype, data = diversity)
```
