function [strConicData, i8ReturnId] = ComputeConicDataFromConicMatrix(dConicMatrix, ...
                                                                       bALLOW_OPEN_CONICS)%#codegen
arguments
    dConicMatrix        (3,3) double {mustBeNumeric, mustBeFinite}
    bALLOW_OPEN_CONICS  (1,1) logical = true
end
%% SIGNATURE
% [strConicData, i8ReturnId] = ComputeConicDataFromConicMatrix(dConicMatrix, ...
%                                                              bALLOW_OPEN_CONICS)%#codegen
% -------------------------------------------------------------------------------------------------------------
%% DESCRIPTION
% What the function does
% 
% REFERENCE:
% [1] J. A. Christian, “A Tutorial on Horizon-Based Optical Navigation and Attitude Determination With Space
% Imaging Systems,” IEEE Access, vol. 9, pp. 19819–19853, 2021, doi: 10.1109/ACCESS.2021.3051914.
% [2] Wikipedia page: "Ellipse", section: General parameters representation
% [3] J. A. Christian, Optical Navigation Using Planet s Centroid and Apparent Diameter in Image, 2015
% -------------------------------------------------------------------------------------------------------------
%% INPUT
% in1 [dim] description
% -------------------------------------------------------------------------------------------------------------
%% OUTPUT
% out1 [dim] description
% -------------------------------------------------------------------------------------------------------------
%% CHANGELOG
% 01-01-2026    Pietro Califano     Implement general purpose version from legacy ComputeHorizonConic().
% -------------------------------------------------------------------------------------------------------------
%% DEPENDENCIES
% [-]
% -------------------------------------------------------------------------------------------------------------

%% Function code

% Set inline properties
coder.inline('default');

% Extract coefficients from matrix
dCoeffsABCDEF = zeros(1,6);
dCoeffsABCDEF(1) = dConicMatrix(1,1);
dCoeffsABCDEF(2) = 2 * dConicMatrix(1,2);
dCoeffsABCDEF(3) = dConicMatrix(2,2);
dCoeffsABCDEF(4) = 2 * dConicMatrix(1,3);
dCoeffsABCDEF(5) = 2 * dConicMatrix(2,3);
dCoeffsABCDEF(6) = dConicMatrix(3,3);
        
% Call implementation
[strConicData, i8ReturnId] = ComputeConicDataFromEqCoeffs(dCoeffsABCDEF, ...
                                                          bALLOW_OPEN_CONICS);

end