function [optBW, kdeFuncH] = kdePickBW(X, smoothness, params, bwLogBounds)
% This picks a bandwidth for the KDE. We use k-fold cross validation in the
% range specified by bwLogBounds. The optimization is done via DiRect.
% If params.getKdeFuncH is True, then it also returns a function handle for the
% kde with the optimal bandiwidth.

  % prelims
  numData = size(X, 1);

  % Shuffle the data
  shuffleOrder = randperm(numData);
  X = X(shuffleOrder, :);

  % Obtain the Standard Deviation of X
  if numData == 1
    stdX = 1;
  else
    stdX = norm( std(X) );
  end
  
  % Set default parameter values
  if ~exist('params', 'var')
    params = struct;
  end
  if ~isfield(params, 'numPartsKFCV')
    params.numPartsKFCV = 5;
  end
  if ~isfield(params, 'numDiRectFunEvals')
    params.numDiRectFunEvals = 20;
  end
  if ~isfield(params, 'getKdeFuncH')
    params.getKdeFuncH = true;
  end
  if ~exist('bwLogBounds', 'var')
    bwLogBounds = log( [1e-4 10] * stdX );
  end

  % Now Set things up for DiRect
  diRectBounds = bwLogBounds;
  options.maxevals = params.numDiRectFunEvals;
  kFoldFunc = @(t) kdeKFoldCV(t, X, smoothness, numParts);
  [~, maxPt] = diRectWrap(kFoldFunc, diRectBounds, options);
  optBW = exp(maxPt);
  
  % Return a function handle
  if params.getKdeFuncH
    kdeFuncH = kdeGivenBW(X, optBW, smoothness);
  else
    kdeFuncH = [];
  end

end

function avgLogLikl = kdeKFoldCV(logBW, X, smoothness, numPartsKFCV) 

  h = exp(logBW);
  logLikls = zeros(numPartsKFCV, 1);
  numData = size(X, 1);

  for kFoldIter = 1:numPartsKFoldCV
    % Set the partition up
    testStartIdx = round( (kFoldIter-1)*numData/numPartsKFCV + 1 );
    testEndIdx = round( kFoldIter*numData/numPartsKFCV );
    trainIndices = [1:(testStartIdx-1), (testEndIdx+1):numData]';
    testIndices = [testStartIdx: testEndIdx]';
    numTestData = testEndIdx - testStartIdx;
    numTrainData = numData - numTestData;
    % Separate Training and Validation sets
    Xtr = X(trainIndices, :);
    Xte = X(testIndices, :);
    % Now Obtain the kde using Xtr
    kdeTr = kdeGivenBW(Xte, h, smoothness);
    % Compute Log Likelihood
    Pte = kdeTr(Xte);
    logLikls(kFoldIter) = mean(log(Pte));
  end

  avgLogLikl = mean(logLikls);

end
