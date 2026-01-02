function bIsValidSet = CheckIfNonColinearSamples(dxCoords, dyCoords)%#codegen
arguments
    dxCoords (:,1) double
    dyCoords (:,1) double
end
%%SIGNATURE
% bIsValidSet = CheckIfNonColinearSamples(dxGuess, dyGuess)%#codegen
% -------------------------------------------------------------------------------------------------------------
%% DESCRIPTION
% Check that the sampled points span a full-rank line/plane system.
% -------------------------------------------------------------------------------------------------------------
%% INPUTS
% dxCoords (:,1) double  x coordinates of the sample points.
% dyCoords (:,1) double  y coordinates of the sample points.
% -------------------------------------------------------------------------------------------------------------
%% OUTPUTS
% bIsValidSet     [1,1] True if points are not colinear.
% -------------------------------------------------------------------------------------------------------------
%% CHANGELOG
% 02-01-2026    Pietro Califano     Implement from legacy code
% -------------------------------------------------------------------------------------------------------------

%% Function code
dxCoords = dxCoords(:);
dyCoords = dyCoords(:);
dPointsMatrix = [dxCoords, dyCoords, ones(size(dxCoords))];

% Evaluate rank of matrix
bIsValidSet = rank( dPointsMatrix ) == 3;

end