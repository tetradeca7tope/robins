
plotFunc = @semilogy;
% plotFunc = @loglog;
% plotFunc = @plot;

% Now plot the results out
figure;
meanPl = mean(plErrors);
stdPl = std(plErrors)/sqrt(numExperiments);
meanDS = mean(dsErrors);
stdDS = std(dsErrors)/sqrt(numExperiments);
meanNDS = mean(ndsErrors);
stdNDS = std(ndsErrors)/sqrt(numExperiments);

plotFunc(nCands, meanPl, 'r'); hold on,
plotFunc(nCands, meanDS, 'b');
plotFunc(nCands, meanNDS, 'g');
legend('Plug In', 'DS', 'LOO');
errorbar(nCands, meanPl, stdPl, 'r');
errorbar(nCands, meanDS, stdDS, 'b');
errorbar(nCands, meanNDS, stdNDS, 'g');

saveFileName = sprintf('results/results-%s.mat', datestr(now, 'ddmm-HHMM'));
save(saveFileName, 'dsEstimates', 'dsErrors', 'plEstimates', 'plErrors', ...
  'ndsEstimates', 'ndsErrors');

allVals = [meanPl, meanDS, meanNDS];
maxVal = max(allVals);
minVal = min(allVals);
yMaxVal = maxVal + 0.1* (maxVal - minVal);
yMinVal = minVal - 0.1* (maxVal - minVal);

% Limits
xlim([0 nCands(end)]);
ylim([yMinVal, yMaxVal]);

