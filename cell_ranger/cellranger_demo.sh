# 23:13 2025/7/10
# tested by wky
#————————————————————————————————————————————————————————————————————————————
# 1.下载安装cellranger
#————————————————————————————————————————————————————————————————————————————
wget -O cellranger-9.0.1.tar.gz "https://cf.10xgenomics.com/releases/cell-exp/cellranger-9.0.1.tar.gz?Expires=1752161309&Key-Pair-Id=APKAI7S6A5RYOXBWRPDA&Signature=ddbLH-QFRnEk1M~dBq9QzwGbe~IggXYbyyS9ywCY1xOwaqtV1v5NtQe0m1-4yVmzCeqMvi3rtbVifN8ChFwJt1KM~2P4LNa10GoRTihc~rvp4XxKJfFUMHDwvbMaPaKiV5gOhAaDLOtXyqf6ux2bc3DxltjYcrillYsvpornxglbp52TORvCSHKWByu7l4YT7geTKDCswBLkAxBVifrEDTrTIGNPR76URpFrMwYoJA73Ue4MJZXq7SWmtPbsNVjTQgZuwwGQqpsPGAURVHwKgKTOoiZU3xlpwN06Xi2MHXGInCE8YNi30fLYMvsZxMzdRhY9JtRpUKlFTcCjOojD3A__"
tar -xzvf cellranger-9.0.1.tar.gz
# 添加入环境
echo 'export PATH=/data5/GPT/Wuky/Tool/cellranger-9.0.1:$PATH' >> ~/.bashrc
source ~/.bashrc
# 检查安装是否成功
which cellranger
#————————————————————————————————————————————————————————————————————————————
# 2.下载fastq数据和参考基因组
#————————————————————————————————————————————————————————————————————————————
mkdir data5/GPT/Wuky/yard/run_cellranger_coun
cd data5/GPT/Wuky/yard/run_cellranger_coun

wget https://cf.10xgenomics.com/samples/cell-exp/3.0.0/pbmc_1k_v3/pbmc_1k_v3_fastqs.tar
tar -xvf pbmc_1k_v3_fastqs.tar

wget https://cf.10xgenomics.com/supp/cell-exp/refdata-gex-GRCh38-2020-A.tar.gz
tar -zxvf refdata-gex-GRCh38-2020-A.tar.gz


#————————————————————————————————————————————————————————————————————————————
# 3.使用cellranger count参数(这里只做演示，不要bam文件)
#————————————————————————————————————————————————————————————————————————————
cellranger count \
  --id=pbmc_1k_count \
  --fastqs=/data5/GPT/Wuky/yard/run_cellranger_coun/pbmc_1k_v3_fastqs \
  --sample=pbmc_1k_v3 \
  --transcriptome=/data5/GPT/Wuky/yard/run_cellranger_coun/refdata-gex-GRCh38-2020-A \
  --create-bam false 

#————————————————————————————————————————————————————————————————————————————
# 4.分析后得到3个重要文件： 
    # 质控：web_summary.html
    # 预分析 cloupe.cloupe
    # 核心下游分析：filtered_feature_bc_matrix.h5
#————————————————————————————————————————————————————————————————————————————

# .cloupe文件：可以查看参考B站上https://www.bilibili.com/opus/1010626218740940809

# 接下来我想查看filtered_feature_bc_matrix.h5
conda create -n scanpy_env python=3.9 -y
conda activate scanpy_env
mamba install -c conda-forge scanpy h5py -y
# python
    # import scanpy as sc
    # import pandas as pd

    # 读取 Cell Ranger 的 h5 矩阵文件
    # adata = sc.read_10x_h5("/data5/GPT/Wuky/yard/run_cellranger_coun/pbmc_1k_count/outs/filtered_feature_bc_matrix.h5")

    # 将稀疏矩阵转为稠密的 DataFrame（慎用：内存要求高！）
    # df = pd.DataFrame(adata.X.toarray(), index=adata.obs_names, columns=adata.var_names)
    # df.head(10).to_csv("pbmc_1k_matrix_preview.csv")

22:02 2025/7/10
# 测试一下Cell Ranger vdj
# 5′ ── ATG [V-segment] GGCATTAC [D-segment] TGTGAC [J-segment] AGCCTGA [C-region] AAA…AAA poly‑A ── 3′
cd /data5/GPT/Wuky/yard/dataset-vdj-practice/
tar -xf sc5p_v2_hs_B_1k_multi_5gex_b_Multiplex_fastqs.tar

curl -O https://cf.10xgenomics.com/supp/cell-vdj/refdata-cellranger-vdj-GRCh38-alts-ensembl-5.0.0.tar.gz
tar -xf refdata-cellranger-vdj-GRCh38-alts-ensembl-5.0.0.tar.gz

cd /data5/GPT/Wuky/yard/runs/
cellranger vdj \
  --id=HumanB_Cell \
  --reference=/data5/GPT/Wuky/yard/refdata-cellranger-vdj-GRCh38-alts-ensembl-5.0.0 \
  --fastqs=/data5/GPT/Wuky/yard/dataset-vdj-practice/sc5p_v2_hs_B_1k_multi_5gex_b_fastqs \
  --sample=sc5p_v2_hs_B_1k_b \
  --localcores=8 \
  --localmem=64

#————————————————————————————————————————————————————————————————————————————
# 5.分析后得到3个重要文件： 
    # 质控：web_summary.html
    # 预分析 vloupe.vloupe
    # 克隆型：clonotypes.csv，filtered_contig_annotations.csv，
#————————————————————————————————————————————————————————————————————————————

# 下载Loupe-VDJ-Browser-5.2.0查看vloupe.vloupe








