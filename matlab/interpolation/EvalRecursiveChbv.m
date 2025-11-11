function dChbvPolynomial = EvalRecursiveChbv(ui32PolyDeg, dScaledPoint, ui32PolyMaxDeg) %#codegen
arguments
    ui32PolyDeg         (1,1) uint32 {isscalar, mustBeNumeric}
    dScaledPoint        (1,1) double {isscalar, mustBeNumeric}
    ui32PolyMaxDeg      (1,1) uint32 {isscalar, mustBeNumeric} = ui32PolyDeg
end
%% PROTOTYPE
% TODO: update doc
% dChbvPolynomial = EvalRecursiveChbv(ui32PolyDeg, dScaledPoint, ui32PolyMaxDeg)
% -------------------------------------------------------------------------------------------------------------
%% DESCRIPTION
% Function constructing the Chebyshev Polynomial evaluated at dScaledPoint point in [-1, 1]
% interval up to dPolyDeg degree. No coefficient is applied (assumed as ones).
% -------------------------------------------------------------------------------------------------------------
%% INPUT
% ui32PolyDeg         (1,1) uint32 {isscalar, mustBeNumeric}  Degree of the Chebyshev polynomial
% dScaledPoint        (1,1) double {isscalar, mustBeNumeric}  Chebyshev polynomial evaluation point (scalar only). Must be in [-1,1]
% ui32PolyMaxDeg      (1,1) uint32 {isscalar, mustBeNumeric} = ui32PolyDeg
% -------------------------------------------------------------------------------------------------------------
%% OUTPUT
% dChbvPolynomial: [ui32PolyMaxDeg+1, 1]  Vector of evaluation Chebyshev polynomials (no coefficients) 
%                                         defined over [-1,1] domain. Max size is determined by ui32PolyMaxDeg.
% -------------------------------------------------------------------------------------------------------------
%% CHANGELOG
% 05-04-2024    Pietro Califano         First version. Validated.
% 01-02-2025    Pietro Califano         Function modified for compatibility with static size requirements
% 18-07-2025    Pietro Califano     Fix basis and fitting problem errors
% -------------------------------------------------------------------------------------------------------------
%% DEPENDENCIES
% [-]
% -------------------------------------------------------------------------------------------------------------
%% Future upgrades
% [-]
% -------------------------------------------------------------------------------------------------------------
%% Function code

assert(ui32PolyDeg >= 2, 'Error: selected degree is too low!')
dChbvPolynomial = zeros(ui32PolyMaxDeg + 1, 1);

% Initialize recursion
dChbvPolynomial(1) = 1.0;
dChbvPolynomial(2) = dScaledPoint;

for idN = 3:ui32PolyDeg + 1
    dChbvPolynomial(idN) = 2.0 * dScaledPoint * dChbvPolynomial(idN-1) - dChbvPolynomial(idN-2);
end

end
