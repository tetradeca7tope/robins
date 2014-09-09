function [X, props] = d1UniMaxConv(numData, params)
% This function generates data in the following manner
% 50% of the time the data just comes from a uniform density
% The other 50% of the time, we sample params.gamma times from a uniform
% distribution and take the maximum.

  N = numData;
  if ~isfield(params, 'gamma') | isempty(params.gamma)
    params.gamma = 10;
  end
  gamma = ceil(params.gamma);

  % Now generate the data
  Z = rand(N, 1+gamma);
  B = double(rand(N, 1) < 0.5);
  X = B.* Z(:,1) + (1-B).*max(Z(:,2:end), [], 2);

  % Now generate some properties of the distribution
  props.gamma = gamma
  props.density = @(t) 0.5 + 0.5* gamma * t.^(gamma-1);

end

