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
functional = 'hellingerDiv';
% And the distribution
distributionIndex = 3;

% Prelims
numDims = 2;
numExperiments = 30;
nCands = [50 100 200:200:2000]';
% nCands = [50 100 200]';
maxN = max(nCands);
numNCands = numel(nCands);
functionalParams = struct;

% Set some of the parameters
params.alpha = 0.05; % confidence level
params.doAsympAnalysis = false;
params.estLowerBound = 0.4; % lower bounds for density estimates
params.smoothness = 2;
% params.bandwidth = @(arg) 1.06 * std(arg) * size(arg,1)^(-1/5);
params.bandwidthX = []; 
params.bandwidthY = []; 
% For the noise terms
params.numPartitions = 2;

% Some parameters for generating the distribution
switch distributionIndex

  case 1
    trueVal = 0;
    trueYDensity = @(t) ones(size(t,1));

  case 2
    gamma = 10;
    trueYDensity = @(t) 0.5 + 0.5*gamma* t.^(gamma-1);
    trueVal = 1-0.921261;

  case 3
    Balpha = 20; Bbeta = 20;
    trueYDensity = @(t) 0.5 + ...
      0.5 * t.^(Balpha-1) .* (1-t).^(Bbeta-1) / beta(Balpha, Bbeta);
end

% The hellinger funtion handle
hellingerFuncHandle = @(t, dens) 1 - sqrt(dens(t));
th = linspace(0,1,10000); trueVal = mean(hellingerFuncHandle(th, trueYDensity));

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
      XX = rand(maxN, numDims);
      YY = rand(maxN, numDims);

    case 2
    % Conv-1D
      XX = rand(maxN, numDims);
      YY = rand(maxN, numDims);
      Z = rand(maxN, 1+gamma); B = double(rand(maxN, 1) < 0.5);
      YY(:,1) = B.* Z(:,1) + (1-B).*max(Z(:,2:end), [], 2);

    case 3
      XX = rand(maxN, numDims);
      YY = rand(maxN, numDims);
      Z = dirichlet_sample([Balpha, Bbeta], maxN);
      B = double(rand(maxN, 1) < 0.5);
      YY(:,1) = B.*YY(:,1) + (1-B).*Z(:,1);

  end

  for nIter = 1:numNCands

    n = nCands(nIter);
    X = XX(1:n, :);
    Y = YY(1:n, :);
    % First estimate the plug in
    plEstim = estimateTwoDistroFunctionalsPlugin(X, Y, functional, ...
      functionalParams, params);
    plEstimates(expIter, nIter) = plEstim;
    plErrors(expIter, nIter) = abs(trueVal - plEstim);
    % Now the linear estimator
    estim = estimateTwoDistroFunctionals(X, Y, functional, functionalParams, params);
    estimates(expIter, nIter) = estim;
    errors(expIter, nIter) = abs(trueVal - estim);
    fprintf('Plug in Err: %0.4f, Linear Err: %.4f\n', ...
      plErrors(expIter, nIter), errors(expIter, nIter) );

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


