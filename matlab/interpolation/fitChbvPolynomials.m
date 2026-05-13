function [dChbvCoeffs, dScaledInterpDomain, strfitStats] = fitChbvPolynomials(ui32PolyDeg, ...
                                                                                dInterpDomain, ...
                                                                                dDataMatrix, ...
                                                                                dDomainLB, ...
                                                                                dDomainUB, ...
                                                                                bENABLE_FIT_CHECK, ...
                                                                                bEnableErrorThrow, ...
                                                                                dPercRelErrorTol) %#codegen
arguments
    ui32PolyDeg         (1,1) uint32
    dInterpDomain       (:,1) double
    dDataMatrix         (:,:) double
    dDomainLB           (1,1) double 
    dDomainUB           (1,1) double  
    bENABLE_FIT_CHECK   (1,1) logical = true
    bEnableErrorThrow   (1,1) logical = true
    dPercRelErrorTol    (1,1) double {mustBeNumeric} = 0.1
end
%% PROTOTYPE
% [dChbvCoeffs, dScaledInterpDomain, strfitStats] = fitChbvPolynomials(ui32PolyDeg, ...
%                                                                      dInterpDomain, ...
%                                                                      dDataMatrix, ...
%                                                                      dDomainLB, ...
%                                                                      dDomainUB, ...
%                                                                      bENABLE_FIT_CHECK, ...
%                                                                      bEnableErrorThrow, ...
%                                                                      dPercRelErrorTol) %#codegen
% -------------------------------------------------------------------------------------------------------------
%% DESCRIPTION
% Function for fitting interpolation coefficients of Chebyshev Polynomial up to the specified degree. The
% interpolant maps the input 1D domain to a N-dimensional domain as specified by the 1st dimension of the
% data matrix. However, note that each entry of the jth sample is interpolated by a different polynomial.
% The output is a 1D vector containing the coefficients for Chebyshev Polynomials from the 1st to PolyDeg-th
% degree. The scaling to [-1,1] domain is automatically handled.
% -------------------------------------------------------------------------------------------------------------
%% INPUT
% ui32PolyDeg         (1,1) uint32
% dInterpDomain       (:,1) double
% dDataMatrix         (:,:) double
% dDomainLB           (1,1) double
% dDomainUB           (1,1) double
% bENABLE_FIT_CHECK   (1,1) logical = true
% bEnableErrorThrow   (1,1) logical {islogical, isscalar} = true
% dPercRelErrorTol    (1,1) double {mustBeNumeric, isscalar} = 0.1
% -------------------------------------------------------------------------------------------------------------
%% OUTPUT
% dChbvCoeffs
% dScaledInterpDomain
% strfitStats
% -------------------------------------------------------------------------------------------------------------
%% CHANGELOG
% 07-04-2024        Pietro Califano     First version. Validated.
% 08-05-2024        Pietro Califano     Updated with error checks.
% 18-07-2025        Pietro Califano     Fix basis and fitting problem errors
% 06-08-2025        Pietro Califano     Update input options to set error tols
% -------------------------------------------------------------------------------------------------------------
%% DEPENDENCIES
% checkFitChbvPoly() if check is enabled
% EvalRecursiveChbv()
% -------------------------------------------------------------------------------------------------------------

%% Function code

% Check input dimensions
% dDataMatrix: [L, N] where N is the number of points, L is the output vector size
assert(size(dDataMatrix, 2) == length(dInterpDomain));
assert(ui32PolyDeg >= 2);

assert(length(dInterpDomain) >= ui32PolyDeg+1);

% Get size of the output vector
ui8OutputSize = size(dDataMatrix, 1);

% Allocate output matrix
dChbvCoeffs = zeros(ui8OutputSize*(ui32PolyDeg+1), 1);

if nargin < 3
    dDomainUB = max(dInterpDomain, [], 'all');
    dDomainLB = min(dInterpDomain, [], 'all');
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
dChbvCoeffs(1:end) = dChbvCoeffs_matrixT(:);

if bENABLE_FIT_CHECK == true
        [strfitStats] = checkFitChbvPoly(ui32PolyDeg, dInterpDomain, dChbvCoeffs, ...
            dDataMatrix, dDomainLB, dDomainUB, false, [], bEnableErrorThrow, dPercRelErrorTol);
else
    strfitStats = struct();
end

end
