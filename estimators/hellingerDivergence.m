function [estim, asympAnalysis, bwX, bwY] = hellingerDivergence(X, Y, ...
  functionalParams, params)
% Estimate for the Hellinger Divergence
  [estim, asympAnalysis, bwX, bwY] = getTwoDistroInfFunAvgs(X, Y, ...
    @hellingerInfFunX, @hellingerInfFunY, @hellingerAsympVar, params);
end

function infFunVals = hellingerInfFunX(densXatX, densYatX)
% densXatX and densYatX are the densities of X and Y respectively (or their
% estimates) at the X points.
  infFunVals = 1 - 0.5 * sqrt(densYatX)./sqrt(densXatX);
end

function infFunVals = hellingerInfFunY(densXatY, densYatY)
% densXatY and densYatY are the densities of X and Y respectively (or their
% estimates) at the Y points.
  infFunVals = - 0.5 * sqrt(densXatY)./sqrt(densYatY);
end

function [asympVarX, asympVarY] = hellingerAsympVar(densXatX, densXatY, ...
  densYatX, densYatY)
% asympVarX and asympVarY are the asymptotic variances of the X and Y influence
% functions.
  asympVarX = 0.5 - ...
    0.25 * (mean(sqrt(densYatX)./sqrt(densXatX)) + ...
             mean(sqrt(densXatY)./sqrt(densYatY)) );
  asympVarY = asympVarX; 
end
