% Script to plot v3 errors
% The following should be loaded: estimMethods, functional, methErrors,
% nCands, numExperiments, 

plotColours = {'r', 'k', 'b', 'g', 'm', 'c', [255 128 0]/255, ...
  [76, 0, 153]/253, [102 102 0]/255, 'y'};
plotShapes = {'+', '*','o',  'x', 's', 'd', '^', 'p', '>', 'v'};
plotShapesDot = {'+.', '*.','o.',  'x.', 's.', 'd.', '^.', 'p.', '>.', 'v.'};
plotShapesDash = {'+-', '*-','o-',  'x-', 's-', 'd-', '^-', 'p-', '>-', 'v-'};
plotFunc = @semilogx;
plotFunc = @loglog;
% plotFunc = @semilogy;

% Removing outliers
% outlierStart = 0.2; outlierEnd = 0.8;
outlierStart = 0; outlierEnd = 1;

% TODO: Set markersize and linewidth.
MS = 10;
LW = 2;

figure;

numMethods = numel(estimMethods);
numNCands = numel(nCands);

% Obtain the means and stds
methErrMeans = zeros(numMethods, numNCands);
methErrStds = zeros(numMethods, numNCands);
for i = 1:numMethods
  currErrs = sort(methErrors{i});
  currErrs = currErrs( round(outlierStart*numExperiments)+1: ...
    round(outlierEnd*numExperiments), : );
  sumIdxs = isfinite(currErrs(:,end));
  methErrMeans(i,:) = mean(currErrs(sumIdxs, :));
  methErrStds(i,:) = std(currErrs(sumIdxs, :))/sqrt(numExperiments);
end

% First plot the bullets
for i = 1:numMethods
  plotFunc(nCands, methErrMeans(i,:), plotShapesDash{i}, 'Color', plotColours{i});
  hold on,
end
legend(estimMethods);

% Now plot error bars
for i = 1:numMethods
  errorbar(nCands, methErrMeans(i,:), methErrStds(i,:), 'Color', plotColours{i});
end

% Limits
allVals = methErrMeans(:);
maxVal = max(allVals);
minVal = min(allVals);
yMaxVal = 1.10* maxVal; %maxVal + 0.05* (maxVal - minVal);
yMinVal = 0.90*minVal; %min(minVal - 0.2* (maxVal - minVal), 0);
% Set X and Y limits 
xlim([0 nCands(end)]*1.05);
ylim([yMinVal, yMaxVal]);
title(functionalName);
xlabel('n');
ylabel('Error');

set(0,'defaultAxesFontName', 'Dejavu Sans')
  set(findall(gca, '-property', 'FontSize'), 'FontSize', 18, ...
    'fontWeight', 'bold');

