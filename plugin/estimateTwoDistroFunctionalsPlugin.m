  function [estim] = estimateOneDistroFunctionalsPlugin...
    (X, Y, functional, functionalParams, params)
  % X is data from densities f .
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

  % The smoothness of the function for the KDE
  if ~isfield(params, 'smoothness')
    params.smoothness = ceil(numDims/2);
  end

  if isstr(functional)
    switch functional

      case 'hellingerDiv'
        func = @(t, densX, densY) hellingerFunctional(t, densX, densY);

    end
  else
    func = functional;
  end

  % First get a density estimate
    if ~isfield(params, 'bandwidth') | isempty(params.bandwidth)
      bwX = kdePickBW(X, params.smoothness, params);
    else
      bwX = params.bandwidth(X);
    end
    % for Y
    if ~isfield(params, 'bandwidth') | isempty(params.bandwidth)
      bwY = kdePickBW(Y, params.smoothness, params);
    else
      bwY = params.bandwidth(Y);
    end
    
  % Obtain the KDE
  densEstX = kdeGivenBW(X, bwX, params.smoothness, params);
  densEstY = kdeGivenBW(Y, bwY, params.smoothness, params);
  
  % Now obtain the estimate
  estim = numIntFuncCompute(func, densEstX, densEstY, numDims);

end


function val = numIntFuncCompute(func, densX, densY, numDims)
  if numDims == 1
    t = linspace(0, 1, 1002)'; t = t(2:end);
    val = mean( func(t, densX, densY) );

  elseif numDims == 2
    t = linspace(0, 1, 102); t = t(2:end);
    [r, s] = meshgrid(t,t);
    T = [r(:) s(:)];
    val = mean( func(T, densX, densY) );
  end
end


function hellingerFuncVals = hellingerFunctional(t, densX, densY)
  densXatT = densX(t);
  densYatT = densY(t);
  hellingerFuncVals = 1 - sqrt(densXatT .* densYatT);
end
