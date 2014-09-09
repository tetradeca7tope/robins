function [estim, asympAnalysis] = entropy(X1, X2, functionalParams, params)
% This estimates the shannon entropy -\int plog(p)

  % First Estimate the density
  if ~isfield(params, 'bandwidth') | isempty(params.bandwidth)
%     fprintf('First finding optimal Bandwidth for KDE\n');
    [optBW, kde] = kdePickBW(X1, params.smoothness, params);
    fprintf('Opt BW: %0.4f\n', optBW);
  else
%     fprintf('Using Given BW: %.5f\n');
    bw = params.bandwidth(X1);
    kde = kdeGivenBW(X1, bw, params.smoothness, params);
  end

  % Now estimate the quantity
  n2 = size(X2, 1);
  phatX2 = kde(X2);
  logEsts = log(phatX2);
%   estim = - mean( logEsts(isfinite(logEsts)) );
  estimNoNoise = - mean(logEsts);
  estim = estimNoNoise + params.normalNoiseTerm * mean( randn(n2, 1) );;

  if params.doAsympAnalysis
    w1 = norminv(1 - params.alpha/2);
    asympAnalysis = struct;
    asympAnalysis.asympVarEst = mean( log(phatX2).^2 ) - estimNoNoise^2 + ...
      params.normalNoiseTerm^2;
    asympStd = sqrt(asympAnalysis.asympVarEst); % + params.normalNoiseTerm;
    asympAnalysis.confInterval(1) = estim - w1 * asympStd/sqrt(n2);
    asympAnalysis.confInterval(2) = estim + w1 * asympStd/sqrt(n2);
  else
    asympAnalysis = [];
  end

end

