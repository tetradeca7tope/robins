% Test for the Renyi Divergence

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
functional = 'renyiDiv';
functionalName = 'Renyi-0.75 Divergence';
% And the distribution
distributionIndex = 2;
numDims = 2;

% Methods
estimMethods = {'Plug-in', 'DS', 'LOO'};
% estimMethods = {'Plug-in', 'DS', 'LOO', 'kNN', 'Szego'};
estimMethods = {'Plug-in', 'DS', 'LOO', 'kNN'};

% Prelims
% numExperiments = 30; nCands = [50 100 200:200:1000 1300 1600 2000]';
% numExperiments = 30; nCands = [30:10:100 150:50:300 400:100:800 1000]';
% numExperiments = 50; nCands = [30:20:100 150:50:400 500:100:1200 1400:200:2000]';
numExperiments = 50; maxN = 2000; numNCands = 15;
nCands = floor(logspace(log10(30), log10(maxN), numNCands));
% numExperiments = 5; nCands = [40 80]'; % Debug

maxN = max(nCands);
numNCands = numel(nCands);
functionalParams = struct;
functionalParams.alpha = 0.75;

% Now for datasplit vs loo estimators
params.doAsympAnalysis = false;
params.doBoundaryCorrection = true;
params.smoothness = 2;
% params.smoothness = 'gauss';
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
    renyiDFunc = @(t) trueDensity(t).^functionalParams.alpha;
    t = linspace(0,1,10000);
    trueVal = 1/(functionalParams.alpha-1) * log( mean(renyiDFunc(t)) );
    % Since Y is uniform, the KL will be the negative entropy of q.

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
      XX = rand(maxN, numDims);
      YY = rand(maxN, numDims);
    case 2
      Z = rand(maxN, 1+gamma); B = double(rand(maxN, 1) < 0.5);
      ZZ = B.* Z(:,1) + (1-B).*max(Z(:,2:end), [], 2);
      XX = [ZZ rand(maxN, numDims-1)];
      YY = rand(maxN, numDims);
  end

  for nIter = 1:numNCands

    n = nCands(nIter);
    X = XX(1:n, :);
    Y = YY(1:n, :);
    fprintf('n = %d\n', n);

    % Now estimate using each method
    for methIter = 1:numMethods
      
      % Now choose the experiment based on the label
      switch estimMethods{methIter}

        case 'Plug-in'
          estim = estimateTwoDistroFunctionalsPlugin(X, Y, functional, ...
            functionalParams, paramspl);

        case 'DS'
          estim = estimateTwoDistroFunctionals(X,Y, functional, functionalParams, ...
            paramsds);

        case 'LOO'
          estim = estimateTwoDistroFunctionals(X,Y, functional, functionalParams, ...
            paramsloo);

        case 'kNN'
%           co.k = 1; %ceil(n/100);
%           co.kNNmethod = 'knnFP1';
          co = DRenyi_kNN_k_initialization(1, {});
          co.alpha = 0.75;
          estim = DRenyi_kNN_k_estimation(X', Y', co);

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
saveFileName = sprintf('results/%s-%s.mat', functional, datestr(now, 'mmdd-HHMMSS'));
save(saveFileName, 'estimMethods', 'functional', 'methErrors', 'nCands', ...
  'numExperiments', 'numMethods', 'numNCands', 'functionalName');

% plot v3 results
plotV3results;

