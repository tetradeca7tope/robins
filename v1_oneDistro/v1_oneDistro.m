% Tests the asymptotic confidence intervals for the entropy.

% addpath ../../if-estimators/
% rmpath ../estimators/

addpath ../estimators/
addpath ../kde/
rmpath ../../if-estimators/

addpath ../plugin/
addpath ~/libs/kky-matlab/utils/
addpath ~/libs/kky-matlab/ancillary/
clear all;
close all;


% First specify the functional
functional = 'shannonEntropy';

% Define a function handle
switch functional

  case 'shannonEntropy'
    estFunctional = @shannonEntropy;

end
% And the distribution
distributionIndex = 2;

% Prelims
numTests = 200;
numData = 2000;
functionalParams = struct;

% Set some of the parameters
params.alpha = 0.05; % confidence level
params.doAsympAnalysis = true;
params.estLowerBound = 0.4; % lower bounds for density estimates
params.smoothness = 2;
params.bandwidth = @(arg) 1.06 * std(arg) * size(arg,1)^(-1/5);
% params.bandwidth = []; %@(arg) 1.06 * std(arg) * size(arg,1)^(-1/5);
params.numPartitions = 1;

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

% Store the results
estimates = zeros(numTests, 1);
confIntervals = zeros(numTests, 2);
trappedTruth = zeros(numTests, 1);
errors = zeros(numTests, 1);
% For the plugin estimator
plEstimates = zeros(numTests, 1);
plErrors = zeros(numTests, 1);

for testIdx = 1:numTests

  % Generate the Data
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  switch distributionIndex

    case 1
    % Uniform
      X = rand(numData, 1);

    case 2
    % Conv-1D
      N = numData;
      Z = rand(N, 1+gamma); B = double(rand(N, 1) < 0.5);
      X = B.* Z(:,1) + (1-B).*max(Z(:,2:end), [], 2);

  end

  % Estimate the quantity 
%   [estim, asympAnalysis] = estFunctional(X, ...
%     functionalParams, params);
  [estim, asympAnalysis] = estimateOneDistroFunctionals(X, functional, ...
    functionalParams, params);

  % Check for trapping and record values
  estimates(testIdx) = estim;
  confIntervals(testIdx, :) = asympAnalysis.confInterval;
  trappedTruth(testIdx) = (trueVal > asympAnalysis.confInterval(1)) & ...
                          (trueVal < asympAnalysis.confInterval(2));
  errors(testIdx) = abs(trueVal - estim);

  % print results
  fprintf('#%d: Estimate: %0.5f, Err: %0.5f, trapped: %d\n', testIdx, ...
    estim, errors(testIdx), trappedTruth(testIdx));

  % Also obtain the plugin estimator
  [plEstim] = estimateOneDistroFunctionalsPlugin(X, functional, ...
    functionalParams, params);
  plEstimates(testIdx) = plEstim;
  plErrors(testIdx) = abs(trueVal - plEstim);

end

% Now print the results
fprintf('\n');
fprintf('Avg Plugin Error: %f\n', mean(plErrors));
fprintf('Avg Error: %f\n', mean(errors));
fprintf('Trap Percentage: %f\n', mean(trappedTruth));

% Plot the results 
plotRange = [-6, 6];
N = numData;
asymDist = sqrt(N) * (estimates - trueVal) / sqrt(asympAnalysis.asympVar);
plot(asymDist, 0.02*rand(numTests, 1), 'kx'); hold on, % plot the pts
[~, asymDistEst] = kde(asymDist);
plot1DFunction(asymDistEst, plotRange, 'b', '--');
plot1DFunction( @(t) normpdf(t), plotRange, 'r', '-');
% axis([plotRange 0 0.45]);
legend('n^{-1/2}(That -T)/sigmahat', 'KDE on Pts', 'N(0,1)');
titlestr = sprintf('%s\nn=%d, m=%d, trap-perc:%.4f', ...
  functional, numData, numTests, mean(trappedTruth) );
title(titlestr);

% QQ plot
figure;
qqplot(asymDist);
