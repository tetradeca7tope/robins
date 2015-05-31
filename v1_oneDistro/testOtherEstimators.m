% A script to test other estimators

addpath ../estimators/
addpath ../kde/
rmpath ../../if-estimators/

% addpath ../../ITE/code/
% ITE_add_to_path('/usr0/home/kkandasa/projects/Robins/ITE/code'); % install ITE
% ITE_install('/usr0/home/kkandasa/projects/Robins/ITE/code'); % install ITE

dim = 1;
dim = 2;

% Generate Data
N = 1000;
gamma = 10;
Z = rand(N, 1+gamma); B = double(rand(N, 1) < 0.5);
X = B.* Z(:,1) + (1-B).*max(Z(:,2:end), [], 2);
if dim == 2
  X = [X rand(N, 1)];
end
trueDensity = @(t) 0.5 + 0.5*gamma* t.^(gamma-1);
params.estLowerBound = 0.4;

% Compute the true Entropy
entropyFunc = @(t) trueDensity(t) .* log( trueDensity(t) );
t = linspace(0,1,1000);
trueEntropy = -mean(entropyFunc(t));

% First our estimator
params = struct;
params.doAsympAnalysis = false;
params.estLowerBound = 0.01;
[estDS, asympAnylsis] = estimateOneDistroFunctionals...
  (X, 'shannonEntropy', struct(), params);
errDS = abs(trueEntropy - estDS);
fprintf('Ours : %.4f,  Err : %.4f \n', estDS, errDS);

% Entropy kNN estimator 
co.k = N/100;
co.kNNmethod = 'knnFP1';
est = HShannon_kNN_k_estimation(X', co);
fprintf('kNN (%d): %.4f, Err: %.4f\n', co.k, est, abs(est-trueEntropy));

% kD partitioning
est = HShannon_KDP_estimation(X');
fprintf('kDP:  %.4f, Err: %.4f\n', est, abs(est-trueEntropy));

% VKDE
est = HShannon_Voronoi_estimation(X');
fprintf('Voronoi:  %.4f, Err: %.4f\n', est, abs(est-trueEntropy));

% Edgeworth expansion
est = -HShannon_Edgeworth_estimation(X');
fprintf('Edgeworth :  %.4f, Err: %.4f\n', est, abs(est-trueEntropy));
