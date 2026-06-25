%% ===================== 1. 数据准备（替换成你的胰岛数据） =====================
clear; clc; close all;

% 你的数据：行=细胞/ROI，列=蛋白 (Insulin/Glucagon/Somatostatin...)
load('Step1_tSNE.mat');   % 必须包含 data_mat (N细胞 × M指标)

data_mat = data{:,I_feature};
labels   = data{:,'labels'};
type     = data{:,'type'};

data_norm = data_mat;

%% ===================== 2. tSNE 降维 =====================
rng(59353);  % 固定结果
tsne_coords = tsne(data_norm,'NumDimensions',2);

% tsne_coords = tsne(data_norm, ...
%     'NumDimensions', 2, ...
%     'Perplexity', 20, ...
%     'Exaggeration', 12);

% 绘制 tSNE
figure;
scatter(tsne_coords(:,1), tsne_coords(:,2), 12,type, 'filled');
title('Islet tSNE');
xlabel('tSNE1'); ylabel('tSNE2');


figure;
scatter(tsne_coords(:,1), tsne_coords(:,2), 12,labels, 'filled');
title('Islet tSNE');
xlabel('tSNE1'); ylabel('tSNE2');
