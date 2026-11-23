%% t-SNE embedding and k-medoids clustering of islet-level features
% 输入：同目录下 Step1_tSNE.mat 中的 Islet.Data、k_tsne 和 k_cluster。
% 核心处理：检查胰岛级特征矩阵，计算二维 t-SNE，并在嵌入空间执行五类 k-medoids。
% 输出：tSNE_cluster_results.mat 和 tSNE_cluster_kmedoids.png。
% 验证：检查维度、有限值、常量列、特征量级、空簇及输出文件。
clearvars
clc
close all
scriptDir = fileparts(mfilename('fullpath'));
if isempty(scriptDir)
    scriptDir = pwd;
end
inputFile = fullfile(scriptDir, 'Step1_tSNE.mat');
outputMat = fullfile(scriptDir, 'tSNE_cluster_results.mat');
outputFigure = fullfile(scriptDir, 'tSNE_cluster_kmedoids.png');
nClusters = 5;
clusterDistance = 'cityblock';
fprintf('Loading input: %s\n', inputFile);
if ~isfile(inputFile)
    error('Input file not found: %s', inputFile);
end
loaded = load(inputFile);
requiredVariables = {'Islet', 'k_tsne', 'k_cluster'};
missingVariables = requiredVariables(~isfield(loaded, requiredVariables));
if ~isempty(missingVariables)
    error('Missing required MAT variables: %s', strjoin(missingVariables, ', '));
end
if ~isstruct(loaded.Islet) || ~isfield(loaded.Islet, 'Data')
    error('Islet.Data is required and must contain the islet-level feature matrix.');
end
dataMatrix = double(loaded.Islet.Data);
if ~ismatrix(dataMatrix) || isempty(dataMatrix) || ~isnumeric(dataMatrix)
    error('Islet.Data must be a non-empty numeric 2-D matrix.');
end
[nIslets, nFeatures] = size(dataMatrix);
if nIslets <= nClusters
    error('The number of islets (%d) must exceed the number of clusters (%d).', nIslets, nClusters);
end
nonFiniteMask = ~isfinite(dataMatrix);
if any(nonFiniteMask, 'all')
    badIndices = find(nonFiniteMask, 5, 'first');
    error('Islet.Data contains %d non-finite values. First linear indices: %s', sum(nonFiniteMask, 'all'), mat2str(badIndices'));
end
constantFeatures = find(range(dataMatrix, 1) == 0);
if ~isempty(constantFeatures)
    error('Constant feature columns detected: %s', mat2str(constantFeatures));
end
zeroCounts = sum(dataMatrix == 0, 1);
magnitudeRatios = nan(1, nFeatures);
for featureIndex = 1:nFeatures
    nonZeroValues = abs(dataMatrix(dataMatrix(:, featureIndex) ~= 0, featureIndex));
    if ~isempty(nonZeroValues)
        magnitudeRatios(featureIndex) = max(nonZeroValues) / min(nonZeroValues);
    end
end
largeRangeFeatures = find(magnitudeRatios >= 1e3);
fprintf('Input matrix: %d islets x %d features\n', nIslets, nFeatures);
fprintf('Zero values by feature: %s\n', mat2str(zeroCounts));
if isempty(largeRangeFeatures)
    fprintf('Magnitude check: no feature has a non-zero max/min ratio >= 1e3.\n');
else
    fprintf('Magnitude warning: feature columns with non-zero max/min ratio >= 1e3: %s\n', mat2str(largeRangeFeatures));
    fprintf('Magnitude ratios for flagged columns: %s\n', mat2str(magnitudeRatios(largeRangeFeatures), 5));
end
tsneSeed = double(loaded.k_tsne(1));
clusterSeed = double(loaded.k_cluster(1));
if ~isfinite(tsneSeed) || ~isfinite(clusterSeed)
    error('Random seeds k_tsne and k_cluster must be finite numeric values.');
end
% Preserve the prepared feature scale; no additional standardization is applied.
rng(tsneSeed, 'twister');
tsneCoords = tsne(dataMatrix, 'NumDimensions', 2);
if ~isequal(size(tsneCoords), [nIslets, 2]) || any(~isfinite(tsneCoords), 'all')
    error('t-SNE output must be a finite %d-by-2 matrix.', nIslets);
end
% Cluster the two-dimensional embedding with deterministic initialization.
rng(clusterSeed, 'twister');
[clusterLabels, clusterMedoids, sumd, distances, medoidIndices, clusterInfo] = kmedoids(tsneCoords, nClusters, 'Distance', clusterDistance, 'Replicates', 1);
clusterCounts = accumarray(clusterLabels, 1, [nClusters, 1]);
if numel(clusterLabels) ~= nIslets || any(~isfinite(clusterLabels)) || any(clusterCounts == 0)
    error('k-medoids returned invalid labels or one or more empty clusters.');
end
parameters = struct;
parameters.inputFile = inputFile;
parameters.nClusters = nClusters;
parameters.clusterDistance = clusterDistance;
parameters.tsneSeed = tsneSeed;
parameters.clusterSeed = clusterSeed;
parameters.tsneNumDimensions = 2;
parameters.additionalScaling = 'none';
parameters.matlabVersion = version;
parameters.matlabRelease = version('-release');
qc = struct;
qc.nIslets = nIslets;
qc.nFeatures = nFeatures;
qc.nonFiniteCount = sum(nonFiniteMask, 'all');
qc.zeroCounts = zeroCounts;
qc.constantFeatures = constantFeatures;
qc.magnitudeRatios = magnitudeRatios;
qc.largeRangeFeatures = largeRangeFeatures;
qc.clusterCounts = clusterCounts;
save(outputMat, 'tsneCoords', 'clusterLabels', 'clusterMedoids', 'sumd', 'distances', 'medoidIndices', 'clusterInfo', 'parameters', 'qc', '-v7');
figureHandle = figure('Color', 'w', 'Position', [100, 100, 760, 620], 'Visible', 'off');
gscatter(tsneCoords(:, 1), tsneCoords(:, 2), clusterLabels, lines(nClusters), '.', 16);
xlabel('t-SNE 1');
ylabel('t-SNE 2');
title('Islet-level t-SNE and k-medoids populations');
legend(compose('Population %d', 1:nClusters), 'Location', 'bestoutside');
axis tight
box off
exportgraphics(figureHandle, outputFigure, 'Resolution', 300);
close(figureHandle);
if ~isfile(outputMat) || ~isfile(outputFigure)
    error('Expected output files were not created.');
end
fprintf('Cluster counts: %s\n', mat2str(clusterCounts'));
fprintf('Saved MAT result: %s\n', outputMat);
fprintf('Saved cluster figure: %s\n', outputFigure);
fprintf('Validation completed with %d magnitude warning(s); execution was not blocked.\n', numel(largeRangeFeatures));
