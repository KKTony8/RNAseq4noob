getwd()
.libPaths()
setwd("C:/Users/tommyw/Desktop/scRNAworkdir")
rm(list =ls())

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("scRepertoire")
library(scRepertoire)
# 格式说明
S1 <- read.csv(".../Sample1/outs/filtered_contig_annotations.csv")
S2 <- read.csv(".../Sample2/outs/filtered_contig_annotations.csv")
S3 <- read.csv(".../Sample3/outs/filtered_contig_annotations.csv")
S4 <- read.csv(".../Sample4/outs/filtered_contig_annotations.csv")

contig_list <- list(S1, S2, S3, S4)
# —————————————————————————————————————————————————————————————————————————————————————————-
# 一.载入数据：####
# 例子1：TCR
# —————————————————————————————————————————————————————————————————————————————————————————-
# 四名急性呼吸窘迫患者的 T 细胞的数据集样品
# 成对的外周血 （B） 和支气管肺泡灌洗液 （L） 组成
# 1.总共8个contig annotation，查看第1个

data("contig_list") #the data built into scRepertoire
head(contig_list[[1]])

# 2.按 barcode 和 CDR3 分组，每个细胞条形码里各个克隆（clonotype）的汇总信息
# 四种模式contig划分策略：归到同一个克隆（clonotype）
    # “gene”：同一套基因片段
    # “nt”：CDR3 的核苷酸序列
    # “aa”：CDR3 的氨基酸序列；
    # CTstrict 基因片段（V/D/J/C）和 CDR3 核苷酸，最严格。
combined.TCR <- combineTCR(contig_list, 
                           samples = c("P17B", "P17L", "P18B", "P18L", 
                                       "P19B","P19L", "P20B", "P20L"),
                           removeNA = FALSE, 
                           removeMulti = FALSE, 
                           filterMulti = FALSE)

head(combined.TCR[[1]])
# —————————————————————————————————————————————————————————————————————————————————————————-
# 例子2：BCR
# —————————————————————————————————————————————————————————————————————————————————————————-
# 1. 从网络或本地读取 BCR 拼装结果表格（filtered_contig_annotations.csv 或同类格式）
BCR.contigs <- read.csv("https://www.borch.dev/uploads/contigs/b_contigs.csv")
# 2. 调用 combineBCR() 合并同条形码下的 contig，并根据序列相似性定义克隆
combined.BCR <- combineBCR(
  BCR.contigs,    # 上面读入的 contig 表
  samples   = "P1",    # 指定样本名字（将用于给条形码加前缀防重复）
  threshold = 0.85     # 定义在 0–1 之间，0.85 表示只有序列相似性 ≥85% 的 CDR3 才会归为同一克隆
)
# 3. 查看合并后第 1 个样本（这里就是 "P1"）的数据框前 6 行
head(combined.BCR[[1]])

# —————————————————————————————————————————————————————————————————————————————————————————-
# 二.可视化Basic Clonal Visualizations####
# —————————————————————————————————————————————————————————————————————————————————————————-
# 1.例子1的TCRlist接下来往下分析####
# 克隆库（clonal repertoire）中，有多少个独特的克隆
# 目的看看是否存在：
    # 克隆扩张（Clonal Expansion）：在免疫应答（感染、疫苗、肿瘤等）过程中，特异性 T 或 B 细胞克隆会成批增殖。
    # 就能定量比较哪种处理或哪段时间免疫应答更强、更聚焦。

# 统计每个样本中，按严格模式(strict)定义的 TCR α+β 克隆的相对比例
clonalQuant(
  combined.TCR, 
  cloneCall = "strict", 
  chain     = "both", 
  scale     = TRUE
)

# 2.特殊：分组分析独特克隆型####
combined.TCR_grp <- lapply(combined.TCR, function(df){
  # 从 sample 名称的最后一个字符提取 B 或 L
  df$Group <- ifelse(grepl("B$", df$sample), "B", "L")
  df
})
clonalQuant(
  combined.TCR_grp,
  cloneCall = "gene",
  chain     = "both",
  scale     = TRUE,
  group.by  = "Group"
)

# 3.“克隆丰度分布”clonalAbundance####
# 横坐标：也就是统计每个 clonotype 在一个样本里覆盖了多少个细胞。
# “大多数克隆只占 1–2 个细胞”（高多样性），还是有少数克隆占据 10、20、甚至更多细胞（寡头扩张）。
clonalAbundance(combined.TCR, 
                cloneCall = "gene", 
                scale = FALSE)

clonalAbundance(combined.TCR, 
                cloneCall = "gene", 
                scale = TRUE)

# 4.克隆长度clonalLength####

clonalLength(combined.TCR, 
             cloneCall="aa", 
             chain = "both") 
#可以见到四个峰，四条 TCR 链（TRA, TRB, TRD, TRG）的 CDR3 氨基酸长度：
# 第一个峰 (~12–15 aa)：TCR α 链 (TRA)
# 第二个峰 (~25–30 aa)：TCR β 链 (TRB)
# 第三个峰 (~40–45 aa)：TCR δ 链 (TRD)
# 第四个峰 (~55–60 aa)：TCR γ 链 (TRG)

clonalLength(combined.TCR, 
             cloneCall="aa", 
             chain = "TRA", 
             scale = TRUE) 

# 5. 克隆比较clonalCompare####
clonalCompare(combined.TCR, 
              top.clones = 10, 
              samples = c("P17B", "P17L"), 
              cloneCall="aa", 
              graph = "alluvial")

# 突出展示clonetype展示相对比例变化
clonalCompare(combined.TCR, 
              top.clones = 10,
              highlight.clones = c("CVVSDNTGGFKTIF_CASSVRRERANTGELFF", "NA_CASSVRRERANTGELFF"),
              relabel.clones = TRUE,
              samples = c("P17B", "P17L"), 
              cloneCall="aa", 
              graph = "alluvial")

# 选取一部分clonetype展示相对比例变化
clonalCompare(combined.TCR, 
              clones = c("CVVSDNTGGFKTIF_CASSVRRERANTGELFF", "NA_CASSVRRERANTGELFF"),
              relabel.clones = TRUE,
              samples = c("P17B", "P17L"), 
              cloneCall="aa", 
              graph = "alluvial")


# 6. 克隆散射clonalScatter####
# 若一个点位于对角线（y=x）上，表示该克隆在两个组织中分布类似
# 若点偏离对角线，例如 x 很高、y 很低，则说明该克隆在 P18B 中丰度高，在 P18L 中几乎没有（组织特异性）；
clonalScatter(combined.TCR, 
              cloneCall ="gene", 
              x.axis = "P18B", 
              y.axis = "P18L",
              dot.size = "total",
              graph = "proportion")

# 7.克隆稳态clonalHomeostasis####
# 克隆等级划分（cloneSize）
clonalHomeostasis(combined.TCR, 
                  cloneCall = "gene")
# 重新定义划分等级
clonalHomeostasis(combined.TCR, 
                  cloneCall = "gene",
                  cloneSize = c(Rare = 0.001, Small = 0.01, Medium = 0.1, Large = 0.3, Hyperexpanded =
                                  1))
combined.TCR <- addVariable(combined.TCR, 
                            variable.name = "Type", 
                            variables = rep(c("B", "L"), 4))

clonalHomeostasis(combined.TCR, 
                  group.by = "Type",
                  cloneCall = "gene")

# 按排名top多少分箱来绘制比例图
clonalProportion(combined.TCR, 
                 cloneCall = "gene") 

clonalProportion(combined.TCR, 
                 cloneCall = "nt",
                 clonalSplit = c(1, 5, 10, 100, 1000, 10000)) 


