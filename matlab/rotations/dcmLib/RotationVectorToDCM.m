function dDCM = RotationVectorToDCM(dRotationVector, dSmallAngleThreshold) %#codegen
arguments
    dRotationVector        (:,1) {mustBeNumeric}
    dSmallAngleThreshold   (1,1) {mustBeNumeric, mustBeNonnegative, coder.mustBeConst} = 1e-3 % [rad]
end

%% PROTOTYPE
% dDCM = RotationVectorToDCM(dRotationVector)
% -------------------------------------------------------------------------------------------------------------
%% DESCRIPTION
% Function computing the direction cosine matrix corresponding to a 3D
% rotation vector. The vector direction defines the rotation axis and its
% norm defines the rotation angle in radians. The implementation uses the
% exact exponential map on SO(3), not a first-order approximation.
% -------------------------------------------------------------------------------------------------------------
%% INPUT
% dRotationVector: [3, 1] Rotation vector. Norm is the rotation angle [rad].
% dSmallAngleThreshold: Scalar threshold for small angles (default: 1e-3 rad). For angles below this, a first-order approximation is used.
% -------------------------------------------------------------------------------------------------------------
%% OUTPUT
% dDCM: [3, 3] Direction cosine matrix associated with dRotationVector.
% -------------------------------------------------------------------------------------------------------------
%% CHANGELOG
% 14-06-2026    Pietro Califano     Add exact DCM construction from rotation vector.
% -------------------------------------------------------------------------------------------------------------
%% DEPENDENCIES
% skewSymm
% -------------------------------------------------------------------------------------------------------------

%% Function code

assert(numel(dRotationVector) == 3, ...
    'RotationVectorToDCM:InvalidInput', ...
    'Rotation vector must contain exactly 3 elements.');

dRotationVector = dRotationVector(:);
dRotationAngle = norm(dRotationVector);

% No meaningful rotation
if dRotationAngle <= eps
    dDCM = eye(3);
    return
end

% For very small angles, use first-order approximation
if dRotationAngle <= eps + coder.const(dSmallAngleThreshold)
   dDCM = eye(3) + skewSymm(dRotationVector);
   return
end

% Compute DCM using Rodrigues' rotation formula (exact exponential map on SO(3))
dRotationAxis = dRotationVector / dRotationAngle;
dSkewAxis = skewSymm(dRotationAxis);

dDCM = eye(3) + sin(dRotationAngle) * dSkewAxis + ...
    (1.0 - cos(dRotationAngle)) * (dSkewAxis * dSkewAxis);

end
