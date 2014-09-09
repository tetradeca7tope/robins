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

  % If not specified, split the data
  if ~isfield(params, 'dataSplit')
    params.dataSplit = true;
  end
  % Whether to do Asymptotic Analysis or not
  if ~isfield(params, 'doAsympAnalysis')
    params.doAsympAnalysis = false;
  end
  % The smoothness of the function for the KDE
  if ~isfield(params, 'smoothness')
    params.smoothness = ceil(numDims/2);
  end
  % The noise level
  if ~isfield(params, 'normalNoiseLevel')
    params.normalNoiseLevel = 0;
  end

  % Split the data if needed.
  if params.dataSplit
    numX1 = ceil(numX/2);
    X1 = X(1:numX1, :);
    X2 = X( (numX1+1):end, :);
  else
    X1 = X;
    X2 = X;
  end


  switch functional

    case 'pAlpha'
      [estim, asympAnylsis] = ...
        pAlpha(X1, X2, functionalParams, params);

    case 'entropy'
      [estim, asympAnylsis] = ...
        entropy(X1, X2, functionalParams, params);
      
  end

end
