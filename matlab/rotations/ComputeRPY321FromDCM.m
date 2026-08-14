function [dRPY, bNearGimbalLock] = ComputeRPY321FromDCM(dDCM, dSingularityTolerance) %#codegen
%% SIGNATURE
% [dRPY, bNearGimbalLock] = ComputeRPY321FromDCM(dDCM, dSingularityTolerance)
% -------------------------------------------------------------------------------------------------------------
%% DESCRIPTION
% Decompose a direction cosine matrix as Rz(yaw) * Ry(pitch) * Rx(roll), returning [roll; pitch; yaw]. At gimbal lock
% yaw is set to zero and roll carries the observable coupled angle; the singularity flag must accompany diagnostics.
% -------------------------------------------------------------------------------------------------------------
%% INPUT
% dDCM                    Proper orthogonal direction cosine matrix.
% dSingularityTolerance   Threshold on abs(cos(pitch)) used to flag gimbal lock.
% -------------------------------------------------------------------------------------------------------------
%% OUTPUT
% dRPY                    Roll, pitch, and yaw angles [rad].
% bNearGimbalLock         True when the ZYX representation is singular or near singular.
% -------------------------------------------------------------------------------------------------------------
%% CHANGELOG
% 13-08-2026  Pietro Califano, Codex gpt-5.6     First implementation.
% -------------------------------------------------------------------------------------------------------------
%% DEPENDENCIES
% [-]
% -------------------------------------------------------------------------------------------------------------

arguments (Input)
    dDCM (3,3) double {mustBeFinite}
    dSingularityTolerance (1,1) double {mustBeNonnegative, mustBeFinite} = 1.0e-8
end
arguments (Output)
    dRPY (3,1) double
    bNearGimbalLock (1,1) logical
end

dHorizontalProjection = hypot(dDCM(1,1), dDCM(2,1));
dPitch = atan2(-dDCM(3,1), dHorizontalProjection);
bNearGimbalLock = dHorizontalProjection <= dSingularityTolerance;

if not(bNearGimbalLock)
    dRoll = atan2(dDCM(3,2), dDCM(3,3));
    dYaw = atan2(dDCM(2,1), dDCM(1,1));
else
    dYaw = 0.0;
    dRoll = atan2(sign(dPitch) * dDCM(1,2), dDCM(2,2));
end
dRPY = [dRoll; dPitch; dYaw];

end
