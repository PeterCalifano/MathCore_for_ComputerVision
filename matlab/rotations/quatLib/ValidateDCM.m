function [bValidity] = ValidateDCM(dDCM, dTol, bAssertOnInvalid)%#codegen
arguments
    dDCM              (3,3) {mustBeNumeric}
    dTol              (1,1) {mustBePositive} = eps('single')
    bAssertOnInvalid  (1,1) logical {coder.mustBeConst} = true
end
%% PROTOTYPE
% [bValidity] = ValidateDCM(dDCM, dTol, bAssertOnInvalid)%#codegen
% -------------------------------------------------------------------------------------------------------------
%% DESCRIPTION
% Check whether dDCM is a valid direction cosine matrix. A valid DCM is real-valued and finite, each column
% has unit norm within dTol (abs(dot(Ci, Ci) - 1) <= dTol), and the columns are orthogonal and right-handed
% within dTol (norm(cross(C1, C2) - C3) <= dTol).
% The function has two branches, selected by bAssertOnInvalid:
% - true (default), assert branch: in MATLAB and MEX, assert on the first failed condition, in the order
%   real-valued, finite, unit norm, orthonormal and right-handed. Generated code does not assert and only
%   returns bValidity.
% - false, output branch: return bValidity and never assert or error on an invalid matrix, in MATLAB, MEX
%   or generated code. Use it to test a matrix as a predicate.
% bAssertOnInvalid must be a compile-time constant, so generated code keeps only one branch.
% -------------------------------------------------------------------------------------------------------------
%% INPUT
% dDCM              (3,3) {mustBeNumeric}, matrix to check
% dTol              (1,1) {mustBePositive} = eps('single'), tolerance of the unit-norm and cross-product checks
% bAssertOnInvalid  (1,1) logical {coder.mustBeConst} = true, compile-time flag: true asserts on an invalid
%                   matrix in MATLAB and MEX, false only returns bValidity
% -------------------------------------------------------------------------------------------------------------
%% OUTPUT
% bValidity (1,1) logical, true if dDCM is real-valued, finite, has unit-norm columns and is orthonormal and
%           right-handed within dTol. In the assert branch, MATLAB and MEX error instead of returning false.
% -------------------------------------------------------------------------------------------------------------
%% CHANGELOG
% 08-06-2026    Pietro Califano     First implementation
% 08-11-2026    Pietro Califano, Claude Opus 5.5    Add assert and output branches (compile-time flag) and a
%                                                   real-valued check: complex input such as diag([1i, -1i, 1])
%                                                   is now invalid.
% -------------------------------------------------------------------------------------------------------------
%% DEPENDENCIES
% [-]
% -------------------------------------------------------------------------------------------------------------

%% Function code

C1 = dDCM(:,1);
C2 = dDCM(:,2);
C3 = dDCM(:,3);

% Evaluate conditions
bIsReal = isreal(dDCM);

bIsFinite = all(isfinite(dDCM), 'all');

bHasUnitNorm = abs(dot(C1,C1) - 1.0) <= dTol && ...
    abs(dot(C2,C2) - 1.0) <= dTol && ...
    abs(dot(C3,C3) - 1.0) <= dTol;

bIsOrthoRightHanded = norm(cross(C1,C2) - C3) <= dTol;

bValidity = bIsReal && bIsFinite && bHasUnitNorm && bIsOrthoRightHanded;

% Select the branch at compile time
if coder.const(bAssertOnInvalid)

    % Assert branch: assert if not codegen
    if coder.target('MATLAB') || coder.target('MEX')

        assert(bIsReal, ...
            'ERROR: invalid DCM. Must be real-valued.');

        assert(bIsFinite, ...
            'ERROR: invalid DCM. Contains NaN or Inf.');

        assert(bHasUnitNorm, ...
            'ERROR: invalid DCM. Columns must have unit norm.');

        assert(bIsOrthoRightHanded, ...
            'ERROR: invalid DCM. Columns must be orthonormal and right-handed.');
    end

else
    % Output branch: return bValidity without asserting, in every target
end

end

