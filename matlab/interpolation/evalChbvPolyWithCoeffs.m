function [dChbvInterpVector] = evalChbvPolyWithCoeffs(ui32PolyDeg, ...
                                                      ui32OutputSize, ...
                                                      dEvalPoint, ...
                                                      dChbvCoeffs, ...
                                                      dDomainLB, ...
                                                      dDomainUB, ...
                                                      ui32PtrToLastCoeff, ...
                                                      ui32PolyMaxDeg) %#codegen
arguments
    ui32PolyDeg         (1,1) uint32    % {isscalar, mustBeNumeric} % Commented for speed-up
    ui32OutputSize      (1,1) uint32    % {isscalar, mustBeNumeric}
    dEvalPoint          (1,1) double    % {isscalar, mustBeNumeric}
    dChbvCoeffs         (:,1) double    % {mustBeNumeric, ismatrix}
    dDomainLB           (1,1) double    % {isscalar, mustBeNumeric}
    dDomainUB           (1,1) double    % {isscalar, mustBeNumeric}
    ui32PtrToLastCoeff  (1,1) uint32    = ui32OutputSize * (ui32PolyDeg + 1)% {isscalar, mustBeNumeric} = length(dChbvCoeffs)
    ui32PolyMaxDeg      (1,1) uint32    = ui32PolyDeg% {isscalar, mustBeNumeric} = ui32PolyDeg
end
%% PROTOTYPE
% TODO: update doc
% [dChbvInterpVector] = evalChbvPolyWithCoeffs(ui8PolyDeg, ...
%                                              ui8OutputSize, ...
%                                              dEvalPoint, ...
%                                              dChbvCoeffs, ...
%                                              dDomainLB, ...
%                                              dDomainUB) %#codegen
% -------------------------------------------------------------------------------------------------------------
%% DESCRIPTION
% What the function does
% -------------------------------------------------------------------------------------------------------------
%% INPUT
% ui32PolyDeg         (1,1) uint32    % {isscalar, mustBeNumeric} % Commented for speed-
% ui32OutputSize      (1,1) uint32    % {isscalar, mustBeNumeric}
% dEvalPoint          (1,1) double    % {isscalar, mustBeNumeric}
% dChbvCoeffs         (:,1) double    % {mustBeNumeric, ismatrix}
% dDomainLB           (1,1) double    % {isscalar, mustBeNumeric}
% dDomainUB           (1,1) double    % {isscalar, mustBeNumeric}
% ui32PtrToLastCoeff  (1,1) uint32    % {isscalar, mustBeNumeric} = length(dChbvCoeffs)
% ui32PolyMaxDeg      (1,1) uint32    % {isscalar, mustBeNumeric} = ui32PolyDeg
% -------------------------------------------------------------------------------------------------------------
%% OUTPUT
% dChbvInterpVector
% -------------------------------------------------------------------------------------------------------------
%% CHANGELOG
% 07-04-2024    Pietro Califano     First version, verified in unit test.
% 01-02-2025    Pietro Califano     Upgrade of functions for codegen with static-sized arrays
% 18-07-2025    Pietro Califano     Fix basis and fitting problem errors
% -------------------------------------------------------------------------------------------------------------
%% DEPENDENCIES
% [-]
% -------------------------------------------------------------------------------------------------------------
%% Function code

if coder.target('MATLAB') || coder.target('MEX')
    assert(dEvalPoint >= dDomainLB && dEvalPoint <= dDomainUB, 'ERROR: invalid evaluation point. Out of interpolation bound.');
end
% Check length of used coefficients 
assert(ui32PtrToLastCoeff == (ui32PolyDeg+1)* ui32OutputSize, ...
    'Number of coefficients does not match the expected size!')

% Variables definition (static size)
dChbvPolynomial     = zeros(ui32PolyMaxDeg + 1, 1);
dChbvInterpVector   = zeros(ui32OutputSize, 1);

% Compute scaled evaluation point 
dScaledPoint = coder.nullcopy(0.0);
dScaledPoint(:) = (2 * dEvalPoint - (dDomainLB + dDomainUB)) / (dDomainUB - dDomainLB); % scalar

% Get evaluated Chebyshev polynomials at scaled point
dChbvPolynomial(1:ui32PolyDeg+1) = EvalRecursiveChbv(ui32PolyDeg, dScaledPoint, ui32PolyMaxDeg); % TODO 

% Compute interpolated output value by inner product with coefficients matrix
dChbvInterpVector(1:ui32OutputSize) = transpose( reshape(dChbvCoeffs, ...
                                                ui32PolyDeg + 1, ui32OutputSize) ) * ...
                                                dChbvPolynomial;


end
