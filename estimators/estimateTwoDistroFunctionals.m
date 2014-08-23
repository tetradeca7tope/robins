function [estim, asympAnylsis] = estimateTwoDistroFunctionals...
  (X, Y, functional, functionalParams, params)
% X, Y are data from densities f and g respectively.
% functional is a string specifying which functional to estimate. See below.
% functionalParams specifies any parameters specific to the functional (viz.
%   alpha for the alpha-divergences)
% params is a struct with other parameters (viz. whether to data split or not)

  % prelims
  numX = size(X, 1);
  numY = size(Y, 1);
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
    smoothness = ceil(numDims/2);
  end

  % Split the data if needed.
  if params.dataSplit
    numX1 = ceil(numX/2);
    numY1 = ceil(numY/2);
    X1 = X(1:numX1, :);
    X2 = X( (numX1+1):end, :);
    Y1 = Y(1:numY1, :);
    Y2 = Y( (numY1+1):end, :);
  else
    X1 = X;
    X2 = X;
    Y1 = Y;
    Y2 = Y;
  end


  switch functional

    case 'condTsallisDiv'
      [estim, asympAnylsis] = ...
        conditionalTsallisDivergence(X1, X2, Y1, Y2, functionalParams, params);

    case 'condTsallisMI'

    case 'condRenyiDiv'

    case 'condRenyiMI'

  end

end
