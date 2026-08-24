function dNextQuaternion = IntegrateQuaternionRKMK4Step( ...
    dInitialQuaternion, dAngularVelocityStages, dStepSize) %#codegen
%% SIGNATURE
% dNextQuaternion = IntegrateQuaternionRKMK4Step(dInitialQuaternion, dAngularVelocityStages, dStepSize)
% -------------------------------------------------------------------------------------------------------------
%% DESCRIPTION
% Advance one scalar-first Hamilton quaternion through a fourth-order Runge-Kutta-Munthe-Kaas step. The angular
% velocity samples are expressed in the body frame at the classical RK4 nodes. The update is applied on the
% right, and the returned quaternion is normalized against floating-point drift.
% -------------------------------------------------------------------------------------------------------------
%% INPUT
% dInitialQuaternion       (4,1) double initial quaternion
% dAngularVelocityStages   (3,4) double body-frame angular velocities at RK4 stages
% dStepSize                (1,1) double positive integration step [s]
% -------------------------------------------------------------------------------------------------------------
%% OUTPUT
% dNextQuaternion          (4,1) double normalized propagated quaternion
% -------------------------------------------------------------------------------------------------------------
%% CHANGELOG
% 23-08-2026  Pietro Califano, Codex     Extract stateless code-generation-safe RKMK4 primitive.
% -------------------------------------------------------------------------------------------------------------
%% DEPENDENCIES
% qCross
% -------------------------------------------------------------------------------------------------------------

arguments (Input)
    dInitialQuaternion (4,1) double {mustBeFinite}
    dAngularVelocityStages (3,4) double {mustBeFinite}
    dStepSize (1,1) double {mustBePositive, mustBeFinite}
end
arguments (Output)
    dNextQuaternion (4,1) double
end

% Pull each noncommuting stage rate back to the Lie algebra associated with
% the existing right quaternion update.
dLieRateStages = 0.5 * dAngularVelocityStages;
dStageIncrement1 = dStepSize * dLieRateStages(:, 1);
dStageIncrement2 = dStepSize * ApplyInvRightDexp_( ...
    0.5 * dStageIncrement1, dLieRateStages(:, 2));
dStageIncrement3 = dStepSize * ApplyInvRightDexp_( ...
    0.5 * dStageIncrement2, dLieRateStages(:, 3));
dStageIncrement4 = dStepSize * ApplyInvRightDexp_( ...
    dStageIncrement3, dLieRateStages(:, 4));

% Compose the RK4 Lie-algebra increment once, then renormalize only for
% floating-point drift; the exact exponential update is norm preserving.
dRotationIncrement = (dStageIncrement1 + 2.0 * dStageIncrement2 + ...
    2.0 * dStageIncrement3 + dStageIncrement4) / 6.0;
dIncrementQuaternion = QuaternionExp_(dRotationIncrement);
dNextQuaternion = qCross(dInitialQuaternion, dIncrementQuaternion, false);
dQuaternionNorm = norm(dNextQuaternion);
assert(dQuaternionNorm > eps, ...
    'IntegrateQuaternionRKMK4Step:InvalidQuaternion', ...
    'Quaternion norm must remain strictly positive.');
dNextQuaternion = dNextQuaternion / dQuaternionNorm;

end

function dPulledBackRate = ApplyInvRightDexp_(dIncrement, dLieRate)
%% SIGNATURE
% dPulledBackRate = ApplyInvRightDexp_(dIncrement, dLieRate)
% -------------------------------------------------------------------------------------------------------------
%% DESCRIPTION
% Apply the fourth-order inverse right-trivialized exponential differential for pure Hamilton quaternions.
% -------------------------------------------------------------------------------------------------------------
%% INPUT
% dIncrement       (3,1) double current Lie-algebra stage increment
% dLieRate         (3,1) double stage Lie-algebra rate
% -------------------------------------------------------------------------------------------------------------
%% OUTPUT
% dPulledBackRate  (3,1) double pulled-back stage rate
% -------------------------------------------------------------------------------------------------------------
%% CHANGELOG
% 23-08-2026  Pietro Califano, Codex     First implementation.
% -------------------------------------------------------------------------------------------------------------
%% DEPENDENCIES
% None.
% -------------------------------------------------------------------------------------------------------------

arguments (Input)
    dIncrement (3,1) double {mustBeFinite}
    dLieRate (3,1) double {mustBeFinite}
end
arguments (Output)
    dPulledBackRate (3,1) double
end

% For the right update, dexp^-1_{-u}(v) has a positive first commutator.
dFirstBracket = 2.0 * cross(dIncrement, dLieRate);
dSecondBracket = 2.0 * cross(dIncrement, dFirstBracket);
dPulledBackRate = dLieRate + 0.5 * dFirstBracket + dSecondBracket / 12.0;

end

function dQuaternion = QuaternionExp_(dRotationIncrement)
%% SIGNATURE
% dQuaternion = QuaternionExp_(dRotationIncrement)
% -------------------------------------------------------------------------------------------------------------
%% DESCRIPTION
% Map one pure-quaternion Lie-algebra increment to a scalar-first Hamilton quaternion.
% -------------------------------------------------------------------------------------------------------------
%% INPUT
% dRotationIncrement  (3,1) double Lie-algebra increment
% -------------------------------------------------------------------------------------------------------------
%% OUTPUT
% dQuaternion         (4,1) double unit quaternion
% -------------------------------------------------------------------------------------------------------------
%% CHANGELOG
% 23-08-2026  Pietro Califano, Codex     First implementation.
% -------------------------------------------------------------------------------------------------------------
%% DEPENDENCIES
% None.
% -------------------------------------------------------------------------------------------------------------

arguments (Input)
    dRotationIncrement (3,1) double {mustBeFinite}
end
arguments (Output)
    dQuaternion (4,1) double
end

dIncrementNorm = norm(dRotationIncrement);
if dIncrementNorm < 1.0e-8
    % Preserve sub-threshold rotations with the leading terms of sinc(x).
    dIncrementNormSquared = dIncrementNorm * dIncrementNorm;
    dVectorScale = 1.0 - dIncrementNormSquared / 6.0;
else
    dVectorScale = sin(dIncrementNorm) / dIncrementNorm;
end

dQuaternion = [cos(dIncrementNorm); dVectorScale * dRotationIncrement];

end
