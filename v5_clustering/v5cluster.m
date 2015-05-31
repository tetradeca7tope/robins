% A clustering application

close all;
clear all;

addpath ../kde/
addpath ../estimators/

numPCADims = 4;

load clusData;

% idxs = [1:5 51:55 101:105];
idxs = 1:150;

numIdxs = numel(idxs);
Xtr = cell(numIdxs, 1);
for j = 1:numIdxs
  Xtr{j} = train_data{idxs(j)};
end

size(Xtr),

[pred_labels, A] = distributionClustering(Xtr, 3, 0);

fprintf('Avg_diag of confusion matrix: %f\n', ...
  confusion_matrix(pred_labels, train_labels(idxs)));

I = mat2gray(A, [0.5 1]); imshow(I),
