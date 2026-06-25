
library(SpatialExperiment)
library(SCORPIUS)
library(dplyr)
library(tidyr)
library(ggplot2)
library(pheatmap)
library(viridis)
library(SingleCellExperiment)
library(Rtsne)
library(ggplot2)
library(scater)
setwd('/home/yzx46/IMC_YZX/steinbock_ND_PD_T2D/review_answer/output/pesudotiom_trajectory')
islet_all_info=data.frame(read.table('Step1_tSNE_only medium and large.csv',sep=',',header=1))

#islet_features matrix
ome_exp_col=colnames(islet_all_info)[grep('ome_tiff_dist_1',colnames(islet_all_info))]
select_col=c( ome_exp_col,'HoleRatio')
islet_features=islet_all_info[,select_col]
islet_features <- apply(islet_features, 2, as.numeric)
islet_features[is.na(islet_features)] <- 0
islet_features=as.matrix(t(islet_features))
colnames(islet_features)=islet_all_info$donor_ROI_islet_id

metadata=islet_meta
sce <- SingleCellExperiment(
  assays  = list(counts = islet_features, exprs = islet_features),
  colData =data.frame(metadata)
)
tsne_input <- t(assay(sce, "exprs"))
sce_meta=sce@metadata
if (ncol(sce) > 30000) {
  pca <- prcomp(tsne_input, rank. = 20, center = FALSE, scale. = FALSE)
  tsne_input_pca <- pca$x
} else {
  tsne_input_pca <- tsne_input
}
set.seed(42)
tsne_res <- Rtsne(
  tsne_input_pca,
dims             = 2,         # NumDimensions = 2
perplexity       = 30,        # MATLAB 默认 perplexity = 30
max_iter         = 1000,      # MATLAB 默认 1000
theta            = 0.5,       # Barnes-Hut 近似
check_duplicates = FALSE,
pca              = FALSE # 若 data_norm 已降维,跳过内部 PCA
)
# 写回 SCE
reducedDim(sce, "TSNE") <- tsne_res$Y
colnames(reducedDim(sce, "TSNE")) <- c("tSNE_1", "tSNE_2")


library(monocle3)
# 构建 cell_data_set
# monocle3 要求 rows = features (markers), cols = cells
gene_metadata <- data.frame(
  gene_short_name = rownames(sce),
  row.names       = rownames(sce)
)

cds <- new_cell_data_set(
  expression_data = assay(sce, "exprs"),    
  cell_metadata   = as.data.frame(colData(sce)),
  gene_metadata   = gene_metadata
)

cds <- preprocess_cds(
  cds,
  num_dim     = 3,                           # IMC marker 数少,10-15 足够
  method      = "PCA",
  norm_method = "none",                       # 关键:跳过对数归一
  scaling     = F                          # 仍保留 z-score
)

reducedDims(cds)[["PCA"]] <- t(assay(sce, "exprs"))
cds <- reduce_dimension(
  cds,
  reduction_method = "UMAP",
  preprocess_method = "PCA",
  umap.min_dist     = 0.3,                   
  umap.n_neighbors  = 30
)
# 聚类(给 trajectory 提供 cluster prior)
cds <- cluster_cells(cds, resolution =( 1e-2+1e-3)/2 )
# 学习 principal graph(trajectory backbone)
cds <- learn_graph(
  cds,
  use_partition       = FALSE,                 # 多 partition 时分开学
  learn_graph_control = list(minimal_branch_len = 10)
)

cluster_type <- data.frame(
  donor_ROI_islet_id  = colnames(cds),
  cluster = as.character(clusters(cds)),       # 强制 character,剥离 factor + attr
  type    = as.character(colData(cds)$Type)
)

cluster_type_islet=merge(cluster_type,islet_all_info,by='donor_ROI_islet_id')
#root_cells 定义 principal graph 上的一个节点集合
tsne_root=cluster_type_islet[(cluster_type_islet$labels ==4) & (cluster_type_islet$cluster ==3) & (cluster_type_islet$type.x =='ND'),]
cds <- order_cells(cds, root_cells =tsne_root$donor_ROI_islet_id) 
pseudotime_vec <- pseudotime(cds)
sce$pseudotime <- pseudotime_vec[colnames(sce)]
cluster_vec <- as.character(clusters(cds))
names(cluster_vec) <- colnames(cds)
sce$cluster_type <- cluster_vec[colnames(sce)]
reducedDim(sce, "UMAP") <- reducedDims(cds)[["UMAP"]][colnames(sce), ]

metadata=colData(cds)[,c('donor_ROI_islet_id' ,     'type'     ,   'Type', 'Size_Factor')]
metadata_label_tsne=merge(metadata,islet_all_info[,c('donor_ROI_islet_id','label_paper')],by='donor_ROI_islet_id')
metadata_label_tsne=metadata_label_tsne[match(metadata$donor_ROI_islet_id,metadata_label_tsne$donor_ROI_islet_id),]
rownames(metadata_label_tsne)=metadata_label_tsne$donor_ROI_islet_id
colData(cds)=metadata_label_tsne

#monocle3绘图
monocle_cluster=plot_cells(cds,
           color_cells_by              = "label_paper",
           group_cells_by              = "label_paper",
           label_cell_groups           = TRUE,
           label_groups_by_cluster     = TRUE,
           group_label_size            = 3,
           label_leaves                = F,   # 叶节点数字
           label_branch_points         = F,   # 分支节点数字
           label_roots                 = T,   # 根节点标识
           label_principal_points      = F,   # 所有 principal graph 节点编号
           show_trajectory_graph       = T,   # 完全隐藏轨迹线(可选)
           cell_size                   = 0.3,
           trajectory_graph_color      = "black",
           trajectory_graph_segment_size = 0.5)
ggsave(
  filename = glue("{output_path}/monocle_cluster_tsnelabel.pdf"),
  plot     = monocle_cluster,
  width    = 4,
  height   = 3.5,
  units    = "in",
  device   = cairo_pdf,
  dpi      = 300
)
