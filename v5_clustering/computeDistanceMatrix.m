function [D] = computeL2DistanceMatrix(data, bandwidth)
% data is a cell array containing num_dists matrices. Each matrix of size
% nixnum_dims gives the points sampled from each distribution.

  % Prelims
  num_dists = numel(data);
  num_dists,

  params.doBoundaryCorrection = false;
  params.smoothness = 2;
  params.doAsympAnalysis = false;
  params.bandwidthX = [];
  params.bandwidthY = [];
  
  % First compute the similarity matrix
  D = zeros(num_dists);

  for i = 1:num_dists
    if (num_dists > 40) && (mod(i,10) == 0)
      % then report progress
      fprintf('Computing L2 between %dth distribution and others\n', i);
    end
    fprintf('Computing L2 between %dth distribution and others\n', i);
    for j = (i+1):num_dists

      [l2ij, ~, bwX, bwY] = ...
        estimateTwoDistroFunctionals(data{i}, data{j}, 'hellingerDiv', ...
        struct(), params);

        %  
%       l2ij = norm(data{i} - data{j})^2/numel(data{i});

      % Plugin 
%       co = DHellinger_kNN_k_initialization(1, {});
%       l2ij = 0.5*DHellinger_kNN_k_estimation(data{i}', data{j}', co);

      params.bandwidthX = @(t) bwX;
      params.bandwidthY = @(t) bwY;
      if imag(l2ij) > 0, l2ij = 0;
      end
      D(i,j) = l2ij;
      D(j,i) = l2ij;
    end
  end
end
