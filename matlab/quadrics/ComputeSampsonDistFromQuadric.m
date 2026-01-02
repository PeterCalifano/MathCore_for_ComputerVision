function dSampsDist = ComputeSampsonDistFromQuadric(dxCoords, dyCoords, dQuadricEqCoeffs)%#codegen
arguments
    dxCoords         (1,:) double
    dyCoords         (1,:) double
    dQuadricEqCoeffs (1,6) double 
end
%% SIGNATURE
% dSampsDist = ComputeSampsonDistFromQuadric(dxCoords, dyCoords, dQuadricEqCoeffs)
% -------------------------------------------------------------------------------------------------------------
%% DESCRIPTION
% Compute the Sampson distance of points (x, y) to the quadratic form defined by
% conic parameters A*x^2 + B*x*y + C*y^2 + D*x + E*y + F = 0.
% -------------------------------------------------------------------------------------------------------------
%% INPUTS
% dxCoords         (1,:) double    x coordinates of the data points.
% dyCoords         (1,:) double    y coordinates of the data points.
% dQuadricEqCoeffs (1,6) double    Conic parameters [A B C D E F].
% -------------------------------------------------------------------------------------------------------------
%% OUTPUTS
% dSampsDist       (1,:) Sampson distance for each data point.
% -------------------------------------------------------------------------------------------------------------
%% CHANGELOG
% 02-01-2026    Pietro Califano     Implement function from legacy code.
% -------------------------------------------------------------------------------------------------------------

%% Function code
% Parameters
A = dQuadricEqCoeffs(1);
B = dQuadricEqCoeffs(2);
C = dQuadricEqCoeffs(3);
D = dQuadricEqCoeffs(4);
E = dQuadricEqCoeffs(5);
F = dQuadricEqCoeffs(6);

% Conic equation at each (x, y)
dQuadricEqResidual = A*dxCoords.^2 + B.*dxCoords.*dyCoords + C*dyCoords.^2 + D.*dxCoords + E.*dyCoords + F;

% Gradients in x and y directions
dFx = 2*A*dxCoords + B*dyCoords + D;
dFy = 2*C*dyCoords + B*dxCoords + E;

dSquaredGradSum = dFx.^2 + dFy.^2;
% Set machine precision if zero detected
dSquaredGradSum( dSquaredGradSum < eps) = 1-25 * eps;

% Sampson error value
dSampsDist = (dQuadricEqResidual.^2) ./ dSquaredGradSum;

end
