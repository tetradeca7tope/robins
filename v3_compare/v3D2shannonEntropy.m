% Compare plugin and linear estimators

addpath ../kde/
addpath ../estimators/
addpath ../plugin/
% addpath ~/libs/kky-matlab/utils/
% addpath ~/libs/kky-matlab/ancillary/

% Add path to ITE & Copy and paste this
% addpath ../../ITE/code/
% ITE_add_to_path('/usr0/home/kkandasa/projects/Robins/ITE/code'); % install ITE

clear all;
close all;

% First specify the functional
functional = 'shannonEntropy';
functionalName = 'Shannon Entropy 2D';
% And the distribution
distributionIndex = 2;
numDims = 2;

% Methods
estimMethods = {'Plug-in', 'DS', 'LOO', 'kNN', 'KDP', 'Voronoi'};

% Prelims
% numExperiments = 30; nCands = [50 100 200:200:1000 1300 1600 2000]';
% numExperiments = 30; nCands = [30:10:100 150:50:300 400:100:800 1000]';
% numExperiments = 50; nCands = [30:20:100 150:50:400 500:100:1200 1400:200:2000]';
numExperiments = 50; maxN = 2000; numNCands = 15;
nCands = floor(logspace(log10(30), log10(maxN), numNCands));
% numExperiments = 5; nCands = [40 80]; % Debug

maxN = max(nCands);
numNCands = numel(nCands);
functionalParams = struct;

% Now for datasplit vs loo estimators
params.doAsympAnalysis = false;
params.doBoundaryCorrection = true;
% params.smoothness = 'gauss';
params.smoothness = 2;
% Data split params
paramsds = params;
paramsds.numPartitions = 2;
paramsds.averageAll = false;
% LOO params
paramsloo = params;
paramsloo.numPartitions = 'loo';
paramsloo.averageAll = true;
% Plug in params
paramspl = params;

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
numMethods = numel(estimMethods);
methErrors = cell(numMethods, 1);
methEstimates = cell(numMethods, 1);
for i = 1:numMethods
  methErrors{i} = zeros(numExperiments, numNCands);
  methEstimates{i} = zeros(numExperiments, numNCands);
end

for expIter = 1:numExperiments

  fprintf('Experiment Iter : %d/%d\n=========================================\n', ...
    expIter, numExperiments);

  % Generate the Data
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  switch distributionIndex
    case 1
    % Uniform
      XX = rand(maxN, 2);

    case 2
    % Conv-1D
      Z = rand(maxN, 1+gamma); B = double(rand(maxN, 1) < 0.5);
      ZZ = B.* Z(:,1) + (1-B).*max(Z(:,2:end), [], 2);
      XX = [ZZ rand(maxN, 1)];
  end

  for nIter = 1:numNCands

    n = nCands(nIter);
    X = XX(1:n, :);
    fprintf('n = %d\n', n);

    % Now estimate using each method
    for methIter = 1:numMethods
      
      % Now choose the experiment based on the label
      switch estimMethods{methIter}

        case 'Plug-in'
          estim = estimateOneDistroFunctionalsPlugin(X, functional, ...
            functionalParams, paramspl);

        case 'DS'
          estim = estimateOneDistroFunctionals(X, functional, functionalParams, ...
            paramsds);

        case 'LOO'
          estim = estimateOneDistroFunctionals(X, functional, functionalParams, ...
            paramsloo);

        case 'kNN'
          co.k = 1; %ceil(n/100);
          co.kNNmethod = 'knnFP1';
          estim = HShannon_kNN_k_estimation(X', co);

        case 'KDP'
          estim = HShannon_KDP_estimation(X');

        case 'Voronoi'
          estim = HShannon_Voronoi_estimation(X');

        case 'Vasicek-KDE'
          estim = HShannon_spacing_VKDE_estimation(X');

      end % switch estimMethods{methIter}

      % Now store the results
      currErr = abs(trueVal - estim);
      methEstimates{methIter}(expIter, nIter) = estim;      
      methErrors{methIter}(expIter, nIter) = currErr;
      fprintf('  %s: Estim: %0.4f, Error: %0.4f\n', estimMethods{methIter}, ...
        estim, currErr);

    end % for methIter 

  end % for nIter

end % for expIter

% Save results
saveFileName = sprintf('results/%s-D%s-%s.mat', ...
  functional, numDims, datestr(now, 'mmdd-HHMMSS'));
save(saveFileName, 'estimMethods', 'functional', 'methErrors', 'nCands', ...
  'numExperiments', 'numMethods', 'numNCands', 'functionalName');

% plot v3 results
plotV3results;
