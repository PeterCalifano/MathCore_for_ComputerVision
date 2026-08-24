function [dNextQuaternion, dNextAngularVelocity] = ...
    IntegrateRigidBodyRK4RKMK4Step(dInertiaMatrix, dInitialQuaternion, ...
    dInitialAngularVelocity, dConstantTorque, dStepSize) %#codegen
%% SIGNATURE
% [dNextQuaternion, dNextAngularVelocity] = IntegrateRigidBodyRK4RKMK4Step(dInertiaMatrix, ...
%     dInitialQuaternion, dInitialAngularVelocity, dConstantTorque, dStepSize)
% -------------------------------------------------------------------------------------------------------------
%% DESCRIPTION
% Advance one rigid-body state with classical RK4 for body-frame angular velocity and RKMK4 for its associated
% scalar-first Hamilton/passive quaternion. Torque is constant in the body frame over this step. Matrix validity
% belongs to the owning history integrator so this performance-sensitive primitive does not repeat factorization.
% -------------------------------------------------------------------------------------------------------------
%% INPUT
% dInertiaMatrix          (3,3) double invertible inertia matrix
% dInitialQuaternion      (4,1) double initial quaternion
% dInitialAngularVelocity (3,1) double initial body-frame angular velocity [rad/s]
% dConstantTorque         (3,1) double constant body-frame torque over the step
% dStepSize               (1,1) double positive integration step [s]
% -------------------------------------------------------------------------------------------------------------
%% OUTPUT
% dNextQuaternion         (4,1) double propagated normalized quaternion
% dNextAngularVelocity    (3,1) double propagated body-frame angular velocity [rad/s]
% -------------------------------------------------------------------------------------------------------------
%% CHANGELOG
% 23-08-2026  Pietro Califano, Codex     Extract stateless constant-torque rigid-body step.
% -------------------------------------------------------------------------------------------------------------
%% DEPENDENCIES
% IntegrateQuaternionRKMK4Step
% -------------------------------------------------------------------------------------------------------------

arguments (Input)
    dInertiaMatrix (3,3) double {mustBeFinite}
    dInitialQuaternion (4,1) double {mustBeFinite}
    dInitialAngularVelocity (3,1) double {mustBeFinite}
    dConstantTorque (3,1) double {mustBeFinite}
    dStepSize (1,1) double {mustBePositive, mustBeFinite}
end
arguments (Output)
    dNextQuaternion (4,1) double
    dNextAngularVelocity (3,1) double
end

% Evaluate the four classical RK stages for Euler's body-frame dynamics.
dDerivative1 = EvaluateAngAccel_( ...
    dInertiaMatrix, dInitialAngularVelocity, dConstantTorque);
dAngularVelocity2 = dInitialAngularVelocity + 0.5 * dStepSize * dDerivative1;
dDerivative2 = EvaluateAngAccel_( ...
    dInertiaMatrix, dAngularVelocity2, dConstantTorque);
dAngularVelocity3 = dInitialAngularVelocity + 0.5 * dStepSize * dDerivative2;
dDerivative3 = EvaluateAngAccel_( ...
    dInertiaMatrix, dAngularVelocity3, dConstantTorque);
dAngularVelocity4 = dInitialAngularVelocity + dStepSize * dDerivative3;
dDerivative4 = EvaluateAngAccel_( ...
    dInertiaMatrix, dAngularVelocity4, dConstantTorque);

dNextAngularVelocity = dInitialAngularVelocity + dStepSize * ...
    (dDerivative1 + 2.0 * dDerivative2 + 2.0 * dDerivative3 + dDerivative4) / 6.0;
dAngularVelocityStages = [dInitialAngularVelocity, dAngularVelocity2, ...
                          dAngularVelocity3, dAngularVelocity4];
dNextQuaternion = IntegrateQuaternionRKMK4Step( ...
    dInitialQuaternion, dAngularVelocityStages, dStepSize);

end

function dAngularAcceleration = EvaluateAngAccel_( ...
    dInertiaMatrix, dAngularVelocity, dConstantTorque)
%% SIGNATURE
% dAngularAcceleration = EvaluateAngAccel_(dInertiaMatrix, dAngularVelocity, dConstantTorque)
% -------------------------------------------------------------------------------------------------------------
%% DESCRIPTION
% Evaluate Euler's rigid-body equation in the body frame.
% -------------------------------------------------------------------------------------------------------------
%% INPUT
% dInertiaMatrix      (3,3) double inertia matrix
% dAngularVelocity    (3,1) double body-frame angular velocity [rad/s]
% dConstantTorque     (3,1) double body-frame torque
% -------------------------------------------------------------------------------------------------------------
%% OUTPUT
% dAngularAcceleration (3,1) double body-frame angular acceleration [rad/s^2]
% -------------------------------------------------------------------------------------------------------------
%% CHANGELOG
% 23-08-2026  Pietro Califano, Codex     First implementation.
% -------------------------------------------------------------------------------------------------------------
%% DEPENDENCIES
% None.
% -------------------------------------------------------------------------------------------------------------

arguments (Input)
    dInertiaMatrix (3,3) double {mustBeFinite}
    dAngularVelocity (3,1) double {mustBeFinite}
    dConstantTorque (3,1) double {mustBeFinite}
end
arguments (Output)
    dAngularAcceleration (3,1) double
end

dAngularMomentum = dInertiaMatrix * dAngularVelocity;
dAngularAcceleration = dInertiaMatrix \ ...
    (dConstantTorque - cross(dAngularVelocity, dAngularMomentum));

end
