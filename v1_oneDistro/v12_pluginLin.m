% Compare plugin and linear estimators

addpath ../kde/
addpath ../estimators/
addpath ../plugin/
addpath ~/libs/kky-matlab/utils/
addpath ~/libs/kky-matlab/ancillary/
clear all;
close all;

plotFunc = @semilogy;
% First specify the functional
functional = 'entropy';
% And the distribution
distributionIndex = 2;

% Prelims
numExperiments = 30;
nCands = [50 100 200:200:2000]';
maxN = max(nCands);
numNCands = numel(nCands);
functionalParams = struct;

% Set some of the parameters
params.alpha = 0.05; % confidence level
params.doAsympAnalysis = false;
params.estLowerBound = 0.4; % lower bounds for density estimates
params.smoothness = 2;
% params.bandwidth = @(arg) 1.06 * std(arg) * size(arg,1)^(-1/5);
params.bandwidth = []; %@(arg) 1.06 * std(arg) * size(arg,1)^(-1/5);
% For the noise terms
params.numPartitions = 2;

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
errors = zeros(numExperiments, numNCands);
estimates = zeros(numExperiments, numNCands);


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
    % First estimate the plug in
    plEstim = estimateOneDistroFunctionalsPlugin(X, functional, ...
      functionalParams, params);
    plEstimates(expIter, nIter) = plEstim;
    plErrors(expIter, nIter) = abs(trueVal - plEstim);
    % Now the linear estimator
    estim = estimateOneDistroFunctionals(X, functional, functionalParams, params);
    estimates(expIter, nIter) = estim;
    errors(expIter, nIter) = abs(trueVal - estim);

  end

end


% Now plot the results out
figure;
meanPl = mean(plErrors);
stdPl = std(plErrors)/sqrt(numExperiments);
meanLin = mean(errors);
stdLin = std(errors)/sqrt(numExperiments);

plotFunc(nCands, meanPl, 'r'); hold on,
plotFunc(nCands, meanLin, 'b');
legend('Plug In', 'Linear');
errorbar(nCands, meanPl, stdPl, 'r');
errorbar(nCands, meanLin, stdLin, 'b');

saveFileName = sprintf('results/results-%s.mat', datestr(now, 'ddmm-HHMM'));
save(saveFileName, 'estimates', 'errors', 'plEstimates', 'plErrors');


