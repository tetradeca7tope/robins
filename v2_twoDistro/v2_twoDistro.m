% Tests the asymptotic confidence intervals for the entropy.

addpath ../kde/
addpath ../estimators/
addpath ~/libs/kky-matlab/utils/
addpath ~/libs/kky-matlab/ancillary/
clear all;
close all;

% First specify the functional
functional = 'hellingerDiv';
% And the distribution
distributionIndex = 3;
numDims = 2;

% Prelims
numTests = 200;
numDataX = 4000;
numDataY = 4000;
functionalParams = struct;

% Set some of the parameters
params.alpha = 0.05; % confidence level
params.doAsympAnalysis = true;
params.estLowerBound = 0.4; % lower bounds for density estimates
params.smoothness = 2;
% params.bandwidthX = @(arg) 1.06 * norm(std(arg))^(-1/(2*params.smoothness + numDims));
params.bandwidthX = [];
params.bandwidthY = params.bandwidthX;
params.doBoundaryCorrection = true;

params.bandwidthX = @(arg) 1.51;
params.bandwidthY = @(arg) 0.14;

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

% Store the results
estimates = zeros(numTests, 1);
confIntervals = zeros(numTests, 2);
trappedTruth = zeros(numTests, 1);
errors = zeros(numTests, 1);

for testIdx = 1:numTests

  % Generate the Data
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  switch distributionIndex

    case 1
    % Uniform and Uniform
      X = rand(numDataX, numDims);
      Y = rand(numDataY, numDims);

    case 2
    % Conv-1D
      X = rand(numDataX, numDims);
      Y = rand(numDataY, numDims);
      % Now change the first column of Y
      Z = rand(numDataY, 1+gamma); B = double(rand(numDataY, 1) < 0.5);
      Y(:,1) = B.* Z(:,1) + (1-B).*max(Z(:,2:end), [], 2);

    case 3
    % Conv-1D
      X = rand(numDataX, numDims);
      Y = rand(numDataY, numDims);
      % Now change the first column of Y
      Z = dirichlet_sample([Balpha Bbeta], numDataY);
      B = double(rand(numDataY, 1) < 0.5);
      Y(:,1) = B.* Y(:,1) + (1-B).*Z(:,1);

  end

  % Estimate the quantity 
  [estim, asympAnalysis, bwX, bwY] = estimateTwoDistroFunctionals...
    (X, Y, functional, functionalParams, params);

  % Check for trapping and record values
  estimates(testIdx) = estim;
  confIntervals(testIdx, :) = asympAnalysis.confInterval;
  trappedTruth(testIdx) = (trueVal > asympAnalysis.confInterval(1)) & ...
                          (trueVal < asympAnalysis.confInterval(2));
  errors(testIdx) = abs((trueVal - estim));

  % print results
  fprintf('#%d: Estimate: %0.5f, Err: %0.5f, trapped: %d bws=(%.4f, %.4f)\n',...
    testIdx, estim, errors(testIdx), trappedTruth(testIdx), bwX, bwY);

end

% Now print the results
fprintf('\n');
fprintf('Avg Error: %f\n', mean(errors));
fprintf('Trap Percentage: %f\n', mean(trappedTruth));

% Plot the results 
plotRange = [-6, 6];
N = numDataX + numDataY;
asymDist = sqrt(N) * (estimates - trueVal) / sqrt(asympAnalysis.asympVar);
plot(asymDist, 0.02*rand(numTests, 1), 'kx'); hold on, % plot the pts
[~, asymDistEst] = kde(asymDist);
plot1DFunction(asymDistEst, plotRange, 'b', '--');
plot1DFunction( @(t) normpdf(t), plotRange, 'r', '-');
% axis([plotRange 0 0.45]);
legend('n^{-1/2}(That -T)/sigmahat', 'KDE on Pts', 'N(0,1)');
titlestr = sprintf('%s\nn=%d, (n,m)=(%d,%d), trap-perc:%.4f ', ...
  functional, numDataX,numDataY, numTests, mean(trappedTruth) );
title(titlestr);

% QQ plot
figure;
qqplot(asymDist);


