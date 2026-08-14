function [dAxisDirectionError, dSignedTwist, bTwistSingular] = ...
        DecomposeSO3ErrorIntoAxisAngleError(dErrorDCM, dUnitTwistAxis) %#codegen
%% SIGNATURE
% [dAxisDirectionError, dSignedTwist, bTwistSingular] = ...
%     DecomposeSO3ErrorIntoAxisAngleError(dErrorDCM, dUnitTwistAxis)
% -------------------------------------------------------------------------------------------------------------
%% DESCRIPTION
% Decompose a truth-relative SO(3) error into spin-axis direction error and signed twist about the true axis. The
% quaternion projection is exact; a 180-degree swing about an orthogonal axis has undefined twist and is flagged.
% -------------------------------------------------------------------------------------------------------------
%% INPUT
% dErrorDCM             Truth-relative error direction cosine matrix.
% dUnitTwistAxis        True unit twist axis expressed in the error DCM coordinates.
% -------------------------------------------------------------------------------------------------------------
%% OUTPUT
% dAxisDirectionError   Angle between the true axis and its error-rotated direction [rad].
% dSignedTwist          Signed principal twist angle about the true axis [rad], or NaN when singular.
% bTwistSingular        True when no unique normalized twist quaternion exists.
% -------------------------------------------------------------------------------------------------------------
%% CHANGELOG
% 13-08-2026  Pietro Califano, Codex gpt-5.6     First implementation.
% -------------------------------------------------------------------------------------------------------------
%% DEPENDENCIES
% LogMap_SO3toR3.
% -------------------------------------------------------------------------------------------------------------

arguments (Input)
    dErrorDCM (3,3) double {mustBeFinite}
    dUnitTwistAxis (3,1) double {mustBeFinite}
end
arguments (Output)
    dAxisDirectionError (1,1) double
    dSignedTwist (1,1) double
    bTwistSingular (1,1) logical
end

dAxisNorm = norm(dUnitTwistAxis);
if dAxisNorm <= eps
    error('DecomposeSO3ErrorIntoAxisAngleError:InvalidAxis', ...
        'The twist axis must be a finite nonzero three-vector.');
end
dUnitTwistAxis = dUnitTwistAxis ./ dAxisNorm;

% Compute the angle between the true axis and its comparison direction
dRotatedAxis = dErrorDCM * dUnitTwistAxis;
dAxisDotProduct = min(1.0, max(-1.0, dot(dUnitTwistAxis, dRotatedAxis)));
dAxisDirectionError = acos(dAxisDotProduct);

% Convert the error DCM to a rotation vector and compute the quaternion representation
dErrorRotationVector = LogMap_SO3toR3(dErrorDCM);
dErrorTheta = norm(dErrorRotationVector);

if dErrorTheta < 1.0e-10
    dQuaternionVector = 0.5 .* dErrorRotationVector;
else
    dQuaternionVector = sin(0.5 * dErrorTheta) .* dErrorRotationVector ./ dErrorTheta;
end

dQuaternionScalar = cos(0.5 * dErrorTheta);
dProjectedQuaternionVector = dot(dUnitTwistAxis, dQuaternionVector);
dTwistQuaternionNorm = hypot(dQuaternionScalar, dProjectedQuaternionVector);

bTwistSingular = dTwistQuaternionNorm <= 1.0e-10;
if bTwistSingular
    dSignedTwist = NaN;
    return
end

% Compute the signed twist angle about the true axis using the quaternion representation
% NOTE: this is the rotation angle error about the true axis, removing the swing component.
dSignedTwist = 2.0 * atan2(dProjectedQuaternionVector / dTwistQuaternionNorm, ...
    dQuaternionScalar / dTwistQuaternionNorm);
dSignedTwist = atan2(sin(dSignedTwist), cos(dSignedTwist));

end
