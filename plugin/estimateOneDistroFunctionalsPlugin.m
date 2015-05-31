  function [estim] = estimateOneDistroFunctionalsPlugin...
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

  % The smoothness of the function for the KDE
  if ~isfield(params, 'smoothness')
    params.smoothness = max(ceil(numDims/2), 2);
  end

  if isstr(functional)
    switch functional

      case 'shannonEntropy'
        func = @(t, dens) entropyFunctional(t, dens);

    end
  else
    func = functional;
  end

  % First get a density estimate
  if ~isfield(params, 'bandwidth') | isempty(params.bandwidth)
    bw = kdePickBW(X, params.smoothness, params);
  else
    bw = params.bandwidth(X);
  end
  % Obtain the KDE
  densEst = kdeGivenBW(X, bw, params.smoothness, params);
  
  % Now obtain the estimate
  estim = numIntFuncCompute(func, densEst, numDims);

end


function val = numIntFuncCompute(func, dens, numDims)
  if numDims == 1
    t = linspace(0, 1, 1002)'; t = t(2:end);
    val = nanmean( func(t, dens) );

  elseif numDims == 2
    t = linspace(0, 1, 52); t = t(2:end);
    [r, s] = meshgrid(t,t);
    T = [r(:) s(:)];
    val = nanmean( func(T, dens) );
  end
end


function entropyFuncVals = entropyFunctional(t, dens)
  densT = dens(t);
  entropyFuncVals = - densT .* log(densT);
end
