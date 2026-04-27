function dAngle_rad = ComputAngleBetweenVectors(dVec1, dVec2, bNearParallelRobust) %#codegen
arguments
    dVec1               (3,1) double  {mustBeNumeric}
    dVec2               (3,1) double  {mustBeNumeric}
    bNearParallelRobust (1,1) logical {coder.mustBeConst} = false
end
%% PROTOTYPE
% dAngle_rad = ComputAngleBetweenVectors(dVec1, dVec2, bNearParallelRobust)
% -------------------------------------------------------------------------------------------------------------
%% DESCRIPTION
% Computes the angle between two 3D vectors in radians. Two formulations are available:
%   - Default: acos(dot) — simple and sufficient for general use.
%   - Robust:  atan2(norm(cross), dot) — numerically stable near 0 and pi where acos loses precision.
% The formulation is selected via bNearParallelRobust (compile-time constant via coder.const).
% -------------------------------------------------------------------------------------------------------------
%% INPUT
% dVec1                 (3,1) First vector [any unit]
% dVec2                 (3,1) Second vector [any unit]
% bNearParallelRobust   (1,1) If true, use atan2 formulation for near-parallel robustness [default: false]
% -------------------------------------------------------------------------------------------------------------
%% OUTPUT
% dAngle_rad            (1,1) Angle between the two vectors [rad]
% -------------------------------------------------------------------------------------------------------------
%% CHANGELOG
% 12-04-2026    Pietro Califano     Modernized interface, imported to MathCore from RCS-1. Added dual formulation with bNearParallelRobust switch.
% -------------------------------------------------------------------------------------------------------------
%% DEPENDENCIES
% [-]
% -------------------------------------------------------------------------------------------------------------

%% Function code
if coder.const(bNearParallelRobust)
    % atan2 formulation: stable for angles near 0 and pi
    dVec1Unit = dVec1 / norm(dVec1);
    dVec2Unit = dVec2 / norm(dVec2);
    dAngle_rad = atan2(norm(cross(dVec1Unit, dVec2Unit)), dot(dVec1Unit, dVec2Unit));
else
    % acos formulation: simple, sufficient for general use
    dAngle_rad = acos(dot(dVec1, dVec2) / (norm(dVec1) * norm(dVec2)));
end

end
