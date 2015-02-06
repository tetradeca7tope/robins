function [estim, asympAnylsis] = estimateOneDistroFunctionals...
  (X, functional, functionalParams, params)
% X is data from densities f .
% functional is a string specifying which functional to estimate. See below.
% functionalParams specifies any parameters specific to the functional (viz.
%   alpha for the alpha-divergences)
% params is a struct with other parameters (viz. whether to data split or not)

  % prelims
  numX = size(X, 1);
  numDims = size(X, 2);

  if ~exist('params', 'var')
    params = struct;
  end
  if ~exist('functionalParams', 'var')
    functionalParams = struct;
  end

  % Whether to do Asymptotic Analysis or not
  if ~isfield(params, 'doAsympAnalysis')
    params.doAsympAnalysis = false;
  end
  % The smoothness of the function for the KDE
  if ~isfield(params, 'smoothness')
    params.smoothness = ceil(numDims/2);
  end
  % Number of partitions to split the data into
  if ~isfield(params, 'numPartitions')
    params.numPartitions = 2;
  end

  switch functional

    case 'entropy'
      [estim, asympAnylsis] = ...
        entropy(X, functionalParams, params);

    case 'pAlpha'
      [estim, asympAnylsis] = ...
        pAlpha(X, functionalParams, params);
      
  end

end

