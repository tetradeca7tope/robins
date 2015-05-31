function [estim, asympAnalysis] = ...
  conditionalTsallisDivergence(V1, V2, W1, W2, functionalParams, params)
% Estimates the conditional Tsallis Divergence between p_{X|Z} and p_{Y|Z}
% conditioned on Z. 
% Inputs
%   The samples V1, V2, W1, W2. Where (X, Z)_i = Vi ~ p_{XZ} and
%     (Y, Z)_i = Wi ~ p_{YZ} for i = 1, 2. V1, W1 will be used for estimating
%     the densities and V2, W2 will be used for correcting the first order bias.
%   functionalParams: a sctructure which should specify the value of alpha
%   params: Some Run time parameters
% Outputs
%   estim: an Estimate of conditional Tsallis divergence
%   asymptAnalysis: returns some parameters/ function handles for the asymptotic
%     anlysis

  % Prelims
  numDims = size(V, 2);
  numV1 = size(V, 1);
  numV2 = size(V, 2);
  numW1 = size(W, 1);
  numW2 = size(W, 2);
  alpha = functionalParams.alpha;
  smoothness = params.smoothness;

  % First construct the KDE
  if ~isfield(params, 'bandwidth') | isempty(params.bandwidth)
  % Use cross validation to pick the bandwidth
    [optBW, kdeXZ] = kdePickBW(V1, smoothness, params);
    [optBW, kdeYZ] = kdePickBW(W1, smoothness, params);
    fprintf('Opt BW: %0.4f\n', optBW);
  else
    kdeXZ = kdeGivenBW(V1, bw, params.smoothness, params);
    kdeYZ = kdeGivenBW(W1, bw, params.smoothness, params);
  end

  % Now estimate the quantity
  XZhatAtXZ = kdeXZ(V2);
  XZhatAtYZ = kdeXZ(W2);
  YZhatAtXZ = kdeYZ(V2);
  YZhatAtYZ = kdeYZ(W2);
  estimNoNoise = 1/(1-alpha) + ...
    alpha/(alpha-1) * mean( (XZhatAtXZ ./ YZhatAtXZ).^(alpha -1) ) + ...
    - mean( (XZhatAtYZ ./ YZhatAtYZ).^alpha );
  estim = estimNoNoise + ...
    params.normalNoiseTerm * ( mean(randn(n2, 1)) + mean(randn(n2, 1)) );

  % asymptotic analysis
  if params.doAsympAnalysis

    % First precompute zeta, harmonic mean etc
    zeta = numW2/(numV2 + numW2);
    w1 = norminv(1 - params.alpha/2);
    asympAnalysis.asympVarEst = asympVarEst;
    asympStd = sqrt(asympAnalysis.asympVarEst);
    asympAnalysis.confInterval(1) = estim - w1 * asympStd/

  else
    asympAnalysis = struct();
  end
  
end

