function kde = kdeGivenBW(X, h, smoothness, params)
% Implements Kernel Density Estimator with kernels of order floor(smoothness)
% for the given bandwidth. You should cross validate h externally.
% Inputs
%   X: the nxd data matrix
%   h: bandwidth
%   smoothness: smoothness of the density. (Default value = 2)
% Outputs
%   kde: a function handle to estimate the density. kde takes in N points in a
%     Nxd matrix and outputs an Nx1 vector.

  % prelims
  numDims = size(X, 2);
  numPts = size(X, 1);

  if ~exist('params', 'var')
    params = struct;
  end
  if ~isfield(params, 'doBoundaryCorrection')
    params.doBoundaryCorrection = true;
  end
  if ~isfield(params, 'estLowerBound')
    params.estLowerBound = 0;
  end
  if ~isfield(params, 'estUpperBound')
    params.estUpperBound = Inf;
  end

  if ~params.doBoundaryCorrection
    augX = X;
  else
    % First augment the dataset by mirroring the points close to the boundaries.
    augX = zeros(0, numDims);
    % Our augmented space as 3^d regions. The centre region is the actual space
    % but all others are in the boundary. We iterate through them as follows
    for regionIdx = 0:(3^numDims -1)

      dimRegions = dec2base(regionIdx, 3);
      dimRegions = [repmat('0', 1, numDims - numel(dimRegions)) dimRegions];
      % Now dimRegions is a string of dimRegions characters with each character
      % corrsponding to each dimension. If the character is 0, we look on the
      % lower boundary of the dimension and if 2 we look at the higher boundary
      % of the dimension.

      % Now check for points within h of the bounary and add them to the dataset
      toReplicate = ones(numPts, 1);
      replicX = X;
      for d = 1:numDims
        if dimRegions(d) == '0'
          replicX(:,d) = -replicX(:,d);
          toReplicate = toReplicate .* double( X(:,d) < h );
        elseif dimRegions(d) == '2'
          replicX(:,d) = 2 - replicX(:,d);
          toReplicate = toReplicate .* double( 1 - X(:,d) < h ); 
        end
      end
      replicatedPts = replicX( logical(toReplicate), :);
      augX = [augX; replicatedPts];
      % Note that when dimRegions = '11...1', we will add the original X to augX

    end
    numAugPts = size(augX, 1);
%     fprintf('numPts = %d, numAugPts = %d\n', numPts, numAugPts);
  end % ~params.doBoundaryCorrection

  % Now return the function handle
  kde = @(arg) kdeIterative(arg, augX, h, smoothness, params, numPts);
end


% A function which estimates the KDE at pts. We use this to construct the
% function handle which will be returned.
function ests = kdeIterative(pts, augX, h, smoothness, params, numX)

  numPts = size(pts, 1);
  numData = size(augX, 1);
  maxNumPts = max(1e7, numData);
  ptsPerPartition = min( numPts, ceil(maxNumPts/numData) );

  ests = zeros(numPts, 1);
  % Now iterate through each 'partition' and obtain the relevant kernels
  cumNumPts = 0;
  while cumNumPts < numPts
    currNumPts = min(ptsPerPartition, numPts - cumNumPts);
    if isstr(smoothness) & strcmp(lower(smoothness(1:5)), 'gauss')
      K = kdeGaussKernel(pts(cumNumPts+1: cumNumPts+currNumPts, :), augX, h);
    else
      K = kdeLegendreKernel( pts(cumNumPts+1: cumNumPts+currNumPts, :), ...
            augX, h, smoothness);
    end
    ests(cumNumPts+1 : cumNumPts + currNumPts) = sum(K,2)/numX;
    cumNumPts = cumNumPts + currNumPts;
  end

  % Now truncate those values below and above the bounds
  ests = max(ests, params.estLowerBound);
  ests = min(ests, params.estUpperBound);
end


%   % Debug
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   % Plot the replicated points
%   if numDims == 1
%     plot(X, rand(numPts,1), 'b.'); hold on,
%     plot(augX, rand(numAugPts,1), 'r.');
%   elseif numDims == 2
%     plot(augX(:,1), augX(:,2), 'r.'); hold on,
%     plot(X(:,1), X(:,2), 'b.');
%   elseif numDims == 3
%     plot3(augX(:,1), augX(:,2), augX(:,3), 'r.'); hold on,
%     plot3(X(:,1), X(:,2), X(:,3), 'b.');
%   end
%   pause;

