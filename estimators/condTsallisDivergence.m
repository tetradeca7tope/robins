function [estim, asympAnalysis, bwXZ, bwYZ] = condTsallisDivergence(X, ZX, Y, Z, ...
  functionalParams, params)
% Estimates the Conditional Tsallis Divergence between X and Y given Z.
  XZ = [X ZX];
  YZ = [Y ZY];
  [estim, asympAnalysis, bwXZ, bwYZ] = tsallisDivergence(XZ, YZ, ...
    functionalParams, params);
end

