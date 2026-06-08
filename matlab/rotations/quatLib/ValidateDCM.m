function [bValidity] = ValidateDCM(dDCM, dTol)%#codegen
arguments
    dDCM  (3,3) {mustBeNumeric}
    dTol  (1,1) {mustBePositive} = eps('single')
end
%% PROTOTYPE
% [bValidity] = ValidateDCM(dDCM, dTol)%#codegen
% -------------------------------------------------------------------------------------------------------------
%% DESCRIPTION
% 
% -------------------------------------------------------------------------------------------------------------
%% INPUT
% dDCM  (3,3) {mustBeNumeric}
% dTol  (1,1) {mustBePositive} = eps('single')
% -------------------------------------------------------------------------------------------------------------
%% OUTPUT
% bValidity (1,1) Bool indicating validity
% -------------------------------------------------------------------------------------------------------------
%% CHANGELOG
% 08-06-2026    Pietro Califano     First implementation
% -------------------------------------------------------------------------------------------------------------
%% DEPENDENCIES
% [-]
% -------------------------------------------------------------------------------------------------------------

%% Function code

C1 = dDCM(:,1);
C2 = dDCM(:,2);
C3 = dDCM(:,3);

% Evaluate conditions
bIsFinite = all(isfinite(dDCM), 'all');

bHasUnitNorm = abs(dot(C1,C1) - 1.0) <= dTol && ...
    abs(dot(C2,C2) - 1.0) <= dTol && ...
    abs(dot(C3,C3) - 1.0) <= dTol;

bIsOrthoRightHanded = norm(cross(C1,C2) - dC3) <= dTol;

bValidity = bIsFinite && bHasUnitNorm && bIsOrthoRightHanded;

% Assert if not codegen
if coder.target('MATLAB') || coder.target('MEX')

    assert(bIsFinite, ...
        'ERROR: invalid DCM. Contains NaN or Inf.');

    assert(bHasUnitNorm, ...
        'ERROR: invalid DCM. Columns must have unit norm.');

    assert(bIsOrthoRightHanded, ...
        'ERROR: invalid DCM. Columns must be orthonormal and right-handed.');
end

end

