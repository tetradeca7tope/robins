% Compare plugin and linear estimators

addpath ../kde/
addpath ../estimators/
addpath ../plugin/
% addpath ~/libs/kky-matlab/utils/
% addpath ~/libs/kky-matlab/ancillary/
clear all;
close all;

% First specify the functional
functional = 'shannonEntropy';
% And the distribution
distributionIndex = 2;

% Prelims
numExperiments = 30; nCands = [50 100 200:200:1000 1300 1600 2000]';
numExperiments = 30; nCands = [30:10:100 150:50:300 400:100:800 1000]';
% numExperiments = 5; nCands = [40 80]'; % Debug
maxN = max(nCands);
numNCands = numel(nCands);
functionalParams = struct;

% Set some of the parameters
params.doAsympAnalysis = false;

% Now for datasplit vs non-datasplit estimators
paramsds = params;
paramsds.numPartitions = 2;
paramsnds = params;
paramsnds.numPartitions = 1;

% Some parameters for generating the distribution
switch distributionIndex

  case 1
    trueVal = 0;

  case 2
    gamma = 10;
    trueDensity = @(t) 0.5 + 0.5*gamma* t.^(gamma-1);
    % Compute the true Entropy
    entropyFunc = @(t) trueDensity(t) .* log( trueDensity(t) );
    t = linspace(0,1,10000);
    trueVal = -mean(entropyFunc(t));
end

% True Value
fprintf('Functional: %s, Truth: %f\n', functional, trueVal);

% To store the results
plErrors = zeros(numExperiments, numNCands);
plEstimates = zeros(numExperiments, numNCands);
dsErrors = zeros(numExperiments, numNCands);
dsEstimates = zeros(numExperiments, numNCands);
ndsErrors = zeros(numExperiments, numNCands);
ndsEstimates = zeros(numExperiments, numNCands);


for expIter = 1:numExperiments

  fprintf('Experiment Iter : %d\n=================================\n', expIter);
  % Generate the Data
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  switch distributionIndex
    case 1
    % Uniform
      XX = rand(maxN, 1);
    case 2
    % Conv-1D
      Z = rand(maxN, 1+gamma); B = double(rand(maxN, 1) < 0.5);
      XX = B.* Z(:,1) + (1-B).*max(Z(:,2:end), [], 2);
  end

  for nIter = 1:numNCands

    n = nCands(nIter);
    X = XX(1:n, :);
    fprintf('n = %d\n', n);

    % First estimate the plug in
    plEstim = estimateOneDistroFunctionalsPlugin(X, functional, ...
      functionalParams, params);
    plEstimates(expIter, nIter) = plEstim;
    currErr = abs(trueVal - plEstim);
    plErrors(expIter, nIter) = currErr; 
    fprintf('  Plugin Error: %0.4f\n', currErr);

    % Now the linear estimator
    dsEstim = estimateOneDistroFunctionals(X, functional, functionalParams, paramsds);
    dsEstimates(expIter, nIter) = dsEstim;
    currErr = abs(trueVal - dsEstim);
    dsErrors(expIter, nIter) = currErr;
    fprintf('  DS Error: %0.4f\n', currErr);

    % Now the non-data split version
    ndsEstim = estimateOneDistroFunctionals(X, functional, functionalParams, paramsnds);
    ndsEstimates(expIter, nIter) = ndsEstim;
    currErr = abs(trueVal - ndsEstim);
    ndsErrors(expIter, nIter) = currErr;
    fprintf('  NDS Error: %0.4f\n', currErr);

  end

end



