function dError = ComputeCircleEqErrFromCoeffs(dxCoords, dyCoords, dConicEqCoeffs)%#codegen
%% SIGNATURE
% dError = ComputeCircleEqErrFromCoeffs(dxCoords, dyCoords, dConicEqCoeffs)%#codegen
% -------------------------------------------------------------------------------------------------------------
%% DESCRIPTION
% Compute point-wise residuals of a conic-parameterized circle: distance of each point
% from the fitted center minus the fitted radius.
% -------------------------------------------------------------------------------------------------------------
%% INPUTS
% dxCoords        [1xN]     x coordinates of the data points.
% dyCoords        [1xN]     y coordinates of the data points.
% dConicEqCoeffs  [1x6]     Conic coefficients of the circle (A,B,C,D,E,F).
% -------------------------------------------------------------------------------------------------------------
%% OUTPUTS
% dError          [1xN]     Residuals for each data point (positive outside, negative inside).
% -------------------------------------------------------------------------------------------------------------
%% CHANGELOG
% 02-01-2026    Pietro Califano     Implement from legacy code
% -------------------------------------------------------------------------------------------------------------
%% DEPENDENCIES
% ComputeCircleParameterPolar
% -------------------------------------------------------------------------------------------------------------

%% Function code
dParametersCircle = ComputeCircleParameterPolar(dConicEqCoeffs);

dxCentre = dParametersCircle(1);
dyCentre = dParametersCircle(2);
dRadius  = dParametersCircle(3);

dError = sqrt((dxCoords - dxCentre).^2 + (dyCoords - dyCentre).^2) - dRadius;


end
