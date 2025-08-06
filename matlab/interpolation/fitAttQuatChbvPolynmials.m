function [dChbvCoeffs, dScaledInterpDomain, dswitchIntervals, strfitStats, ui32PtrToLastCoeff] = fitAttQuatChbvPolynmials(ui32PolyDeg, ...
                                                                                                                        dInterpDomain, ...
                                                                                                                        dDataMatrix, ...
                                                                                                                        dDomainLB, ...
                                                                                                                        dDomainUB, ...
                                                                                                                        bENABLE_FIT_CHECK, ...
                                                                                                                        ui32PolyMaxDeg, ...
                                                                                                                        ui32OutputSize) %#codegen
arguments
    ui32PolyDeg         (1, 1) uint32
    dInterpDomain       (:, 1) double
    dDataMatrix         (:, :) double
    dDomainLB           (1, 1) double
    dDomainUB           (1, 1) double
    bENABLE_FIT_CHECK   (1, 1) logical = true
    ui32PolyMaxDeg      (1, 1) uint32  = ui32PolyDeg;
    ui32OutputSize      (1, 1) uint32  = size(dDataMatrix, 1);
end
%% PROTOTYPE
% [dChbvCoeffs, dScaledInterpDomain, dswitchIntervals, strfitStats, ui32PtrToLastCoeff] = fitAttQuatChbvPolynmials(ui32PolyDeg, ...
%                                                                                                                   dInterpDomain, ...
%                                                                                                                   dDataMatrix, ...
%                                                                                                                   dDomainLB, ...
%                                                                                                                   dDomainUB, ...
%                                                                                                                   bENABLE_FIT_CHECK, ...
%                                                                                                                   ui32PolyMaxDeg, ...
%                                                                                                                   ui32OutputSize) %#codegen
% -------------------------------------------------------------------------------------------------------------
%% DESCRIPTION
% Function implementing the fitting of attitude quaternion interpolation polynomials based on Chebyshev
% polynomial basis. The function supports automatic fitting check over a uniform grid against data points.
% -------------------------------------------------------------------------------------------------------------
%% INPUT
% ui32PolyDeg         (1, 1) uint32
% dInterpDomain       (:, 1) double
% dDataMatrix         (:, :) double
% dDomainLB           (1, 1) double
% dDomainUB           (1, 1) double
% bENABLE_FIT_CHECK   (1, 1) logical = true
% ui32PolyMaxDeg      (1, 1) uint32  = ui32PolyDeg;
% ui32OutputSize      (1, 1) uint32  = size(dDataMatrix, 1);
% -------------------------------------------------------------------------------------------------------------
%% OUTPUT
% dChbvCoeffs
% dScaledInterpDomain
% dswitchIntervals
% strfitStats
% ui32PtrToLastCoeff
% -------------------------------------------------------------------------------------------------------------
%% CHANGELOG
% 08-05-2024        Pietro Califano     fitChbvPolynomials specification for Attitude quaternions, with error checks.
% 01-02-2025        Pietro Califano     Minor changes for better compatibility with fncs of toolbox
% 18-07-2025        Pietro Califano     Fix basis and fitting problem errors
% -------------------------------------------------------------------------------------------------------------
%% DEPENDENCIES
% [-]
% -------------------------------------------------------------------------------------------------------------
%% Function code

% Check input dimensions
% i_dDataMatrix: [L, N] where N is the number of points, L is the output vector size
assert(size(dDataMatrix, 2) == length(dInterpDomain));
assert(ui32PolyDeg >= 2);

assert(length(dInterpDomain) >= ui32PolyDeg +1);

% Allocate output matrix
ui32PtrToLastCoeff = ui32OutputSize * (ui32PolyMaxDeg+1);
dChbvCoeffs = zeros(ui32PtrToLastCoeff, 1);

if nargin < 3
    dDomainUB = max(dInterpDomain, [], 'all');
    dDomainLB = min(dInterpDomain, [], 'all');
end

% AUTOMATIC CHECK AND FIX OF DISCONTINUITIES
[dDataMatrix, bIsSignSwitched, ui8howManySwitches, ...
    bsignSwitchDetectionMask] = fixQuatSignDiscontinuity(transpose(dDataMatrix));

% Determine sign switch intervals
dswitchIntervals = zeros(ui8howManySwitches, 2);

if ui8howManySwitches > 0
    switchesIDs = find(bsignSwitchDetectionMask, ui8howManySwitches);
    dswitchIntervals(:, 1) = switchesIDs;

    howManySamples = length(bIsSignSwitched);

    for idC = 1:ui8howManySwitches
        idtmp = switchesIDs(idC);
        while idtmp < howManySamples && bIsSignSwitched(idtmp) == 1
            idtmp = idtmp + 1;
        end
        dswitchIntervals(idC, 2) = idtmp;
    end
end

% Compute scaled domain
dScaledInterpDomain = (2.*dInterpDomain - (dDomainUB+dDomainLB))./(dDomainUB-dDomainLB);

% Compute regressors matrix on scaled domain
dRegrMatrix = zeros(ui32PolyDeg+1, size(dDataMatrix, 2));

for idN = 1:size(dDataMatrix, 2)

    % Evaluate Chebyshev polynomial at scaled point
    dTmpChbvPoly = EvalRecursiveChbv(ui32PolyDeg, dScaledInterpDomain(idN));
    dRegrMatrix(:, idN) = dTmpChbvPoly;
    
end

% Compute fit coefficients matrix (transposed)
% Xmat = Cmat * Phi: [LxN] = [LxM]*[MxN] where M: poly degree, N: number of samples, L: output vector size
% ith Chebyshev polynomial i=1,...N, along each column
% jth element of ith sample has coefficients along each row of Cmat
dChbvCoeffs_matrixT = dRegrMatrix' \ dDataMatrix'; % Solve the transposed problem
% Flatten matrix to 1D vector
dChbvCoeffs(1:ui32PtrToLastCoeff) = dChbvCoeffs_matrixT(:);

%% Automatic error check
if bENABLE_FIT_CHECK == true
    [strfitStats] = checkFitChbvPoly(ui32PolyDeg, dInterpDomain, dChbvCoeffs, ...
        dDataMatrix, dDomainLB, dDomainUB, true, dswitchIntervals);
else
    strfitStats = struct();
end


end
