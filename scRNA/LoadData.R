rm(list =ls())
# 学会读取GEO中多个样本，并整理成seurat对象
## 已经下载好https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE231755放在D:/test
#————————————————————————————————————————————————————————————————
# 1.整理文件夹
#————————————————————————————————————————————————————————————————
library(stringr)
# 读取文件
files <- list.files("D:/test", full.names = TRUE)
files

# 遍历文件并移动到对应的样本文件夹
for (file in files) {
  # 提取样本名（去掉 _barcodes/feature/matrix 后缀）
  sample_name <- gsub("_barcodes.tsv.gz|_features.tsv.gz|_matrix.mtx.gz", "", basename(file))
  
  # 创建样本文件夹（如果不存在）
  sample_folder <- file.path("D:/test", sample_name)
  if (!dir.exists(sample_folder)) {
    dir.create(sample_folder, showWarnings = FALSE)
  }
  
  # 构造目标文件路径：去掉 sample_name_ 前缀
  target_file <- file.path(sample_folder, gsub(paste0(sample_name, "_"), "", basename(file)))
  
  # 移动文件到对应文件夹
  file.rename(file, target_file)
}

# 可视化检查目录结构，fs是个对文件操作的R包
library(fs)
dir_tree("D:/test")

# install.packages("qs")
library(Seurat)
library(Matrix)
library(qs)

#----------------------------------------------------------------------
# 2.读取成seurat对象####
#————————————————————————————————————————————————————————————————————--
# Step 1: 列出样本子文件夹
data_dir <- "D:/test"
samples <- list.dirs(data_dir, recursive = FALSE, full.names = FALSE)
samples

# Step 2: 逐个样本读取数据，把seurat对象存储在一个list当中
scRNAlist <- lapply(samples, function(pro) {
  cat("正在读取：", pro, "\n")
  folder <- file.path(data_dir, pro)
  counts <- Read10X(data.dir = folder, gene.column = 2)
  sce <- CreateSeuratObject(counts, project = pro, min.cells = 3,min.features = 300)
  return(sce)
})
names(scRNAlist) <- samples

# Step 3: 合并所有样本
sce.all <- merge(scRNAlist[[1]], y = scRNAlist[-1], add.cell.ids = samples)
sce.all <- JoinLayers(sce.all)
sce.all
# Step 4: 检查数据和保存
table(sce.all$orig.ident)
head(sce.all@meta.data, 10)

# —— 查看 counts 矩阵的前几行 —— 
# Seurat v5 推荐使用 GetAssayData 访问 counts、data、scale.data 等
mat10 <- GetAssayData(sce.all, assay = "RNA", slot = "counts")
as.data.frame(mat10[1:10, 1:2])

# Step 5: 保存合并对象
dir.create("D:/test/output", showWarnings = FALSE)

# .qs 格式的读写速度更快
qsave(sce.all, file = "D:/test/output/sce.all.qs")
