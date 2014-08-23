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
  numV = size(V, 1);
  numW = size(W, 1);
  alpha = functionalParams.alpha;
  smoothness = params.smoothness;

  % First construct the KDE
  pXZhat = 
  
  asympAnalysis = struct();
  
end

