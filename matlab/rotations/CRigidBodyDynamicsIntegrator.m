classdef CRigidBodyDynamicsIntegrator < handle & matlab.mixin.Copyable
    %% SIGNATURE
    % objIntegrator = CRigidBodyDynamicsIntegrator(dInertiaMatrix)
    % -------------------------------------------------------------------------------------------------------------
    %% DESCRIPTION
    % Integrate passive scalar-first Hamilton attitude quaternions together with body-frame angular velocity.
    % The class preserves callback-based legacy schemes and delegates the common constant-body-torque RK4-RKMK4
    % step to a stateless code-generation-safe numerical primitive.
    % -------------------------------------------------------------------------------------------------------------
    %% INPUT
    % dInertiaMatrix       (3,3) double rigid-body inertia matrix expressed in the body frame
    % -------------------------------------------------------------------------------------------------------------
    %% OUTPUT
    % objIntegrator        (1,1) configured CRigidBodyDynamicsIntegrator
    % -------------------------------------------------------------------------------------------------------------
    %% CHANGELOG
    % 24-08-2026  Pietro Califano, Codex     Remove unused quaternion-integrator injection.
    % 23-08-2026  Pietro Califano, Codex     Document the frame contract and add the stateless RKMK4 fast path.
    % -------------------------------------------------------------------------------------------------------------
    %% DEPENDENCIES
    % CQuatKinematicsIntegrator, IntegrateQuaternionRKMK4Step, IntegrateRigidBodyRK4RKMK4Step
    % -------------------------------------------------------------------------------------------------------------

    properties
        dInertiaMatrix    (3,3) double {mustBeFinite} = zeros(3,3)
    end

    methods (Access = public)
        function self = CRigidBodyDynamicsIntegrator(dInertiaMatrix)
            %% SIGNATURE
            % self = CRigidBodyDynamicsIntegrator(dInertiaMatrix)
            % -----------------------------------------------------------------------------------------------------
            %% DESCRIPTION
            % Construct a rigid-body integrator with a body-frame inertia matrix.
            % -----------------------------------------------------------------------------------------------------
            %% INPUT
            % dInertiaMatrix  (3,3) double body-frame inertia matrix
            % -----------------------------------------------------------------------------------------------------
            %% OUTPUT
            % self            (1,1) configured CRigidBodyDynamicsIntegrator
            % -----------------------------------------------------------------------------------------------------
            %% CHANGELOG
            % 24-08-2026  Pietro Califano, Codex     Remove unused quaternion-integrator injection.
            % 23-08-2026  Pietro Califano, Codex     Consolidate the constructor contract.
            % -----------------------------------------------------------------------------------------------------
            %% DEPENDENCIES
            % None.
            % -----------------------------------------------------------------------------------------------------

            arguments
                dInertiaMatrix (3,3) double {mustBeFinite} = zeros(3,3)
            end

            self.dInertiaMatrix = dInertiaMatrix;
        end

        function [dQuatOutSeq, dOmegaOutSeq] = integrate(self, ...
                                                       dTimegrid, ...
                                                       dQuat0, ...
                                                       dOmega0, ...
                                                       varTorque, ...
                                                       dDeltaT, ...
                                                       enumMethod, ...
                                                       bDecoupledKinoDynamics, ...
                                                       dDefaultMaxDeltaT)
            %% SIGNATURE
            % [dQuatOutSeq, dOmegaOutSeq] = integrate(self, dTimegrid, dQuat0, dOmega0, varTorque, dDeltaT, ...
            %     enumMethod, bDecoupledKinoDynamics, dDefaultMaxDeltaT)
            % -----------------------------------------------------------------------------------------------------
            %% DESCRIPTION
            % Propagate passive attitude and body angular velocity at the requested output epochs. Internal steps
            % are bounded by dDeltaT or dDefaultMaxDeltaT according to the existing class integration contract.
            % -----------------------------------------------------------------------------------------------------
            %% INPUT
            % self                      (1,1) configured rigid-body integrator
            % dTimegrid                 (1,:) double requested output epochs [s]
            % dQuat0                    (4,1) double initial passive scalar-first Hamilton quaternion
            % dOmega0                   (3,1) double initial body-frame angular velocity [rad/s]
            % varTorque                 function handle or double torque model
            % dDeltaT                   (1,1) double explicit internal step [s], or zero to deduce it
            % enumMethod                char integration method
            % bDecoupledKinoDynamics    logical true when torque is independent of attitude
            % dDefaultMaxDeltaT         (1,1) double maximum deduced internal step [s]
            % -----------------------------------------------------------------------------------------------------
            %% OUTPUT
            % dQuatOutSeq               (4,:) double propagated passive quaternion history
            % dOmegaOutSeq              (3,:) double propagated body-frame angular-velocity history [rad/s]
            % -----------------------------------------------------------------------------------------------------
            %% CHANGELOG
            % 23-08-2026  Pietro Califano, Codex     Consolidate the public propagation contract.
            % -----------------------------------------------------------------------------------------------------
            %% DEPENDENCIES
            % CRigidBodyDynamicsIntegrator.integrateStep_
            % -----------------------------------------------------------------------------------------------------

            arguments
                self                    (1,1) CRigidBodyDynamicsIntegrator
                dTimegrid               (1,:) double {mustBeFinite, mustBeVector}
                dQuat0                  (4,1) double {mustBeFinite}
                dOmega0                 (3,1) double {mustBeFinite}
                varTorque               {mustBeA(varTorque, ["function_handle", "double"])}
                dDeltaT                 (1,1) double {mustBeScalarOrEmpty, mustBeGreaterThanOrEqual(dDeltaT, 0.0)} = 0.0;
                enumMethod              (1,:) char {coder.mustBeConst, mustBeMember(enumMethod, {'lie_euler', 'rk4_rkmk4', 'rk4'})} = 'rk4_rkmk4'
                bDecoupledKinoDynamics  (1,1) logical = true % True if torque does not depend on attitude
                dDefaultMaxDeltaT       (1,1) double {mustBeScalarOrEmpty, mustBeGreaterThanOrEqual(dDefaultMaxDeltaT, 0.0)} = 1.0;
            end 

            bDeduceDeltaT = dDeltaT == 0.0;
            dT0 = dTimegrid(1);
            dTf = dTimegrid(end);

            if size(dTimegrid,2) == 2
                assert(dDeltaT > 0.0, 'ERROR: timestep must be provided if timegrid is given as [t0,tf].')
                dTimegrid = dT0:dDeltaT:dTf;
                bDeduceDeltaT = false;
            end

            % Allocate output sequence array
            ui32NumTimes = size(dTimegrid, 2);

            dQuatOutSeq  = zeros(4,  ui32NumTimes);
            dOmegaOutSeq = zeros(3, ui32NumTimes);
            dTimegridOut = zeros(1, ui32NumTimes);

            % Select implementation of step
            switch coder.const(enumMethod)

                case {'lie_euler', 'rk4', 'rk4_rkmk4'}
                    fcnIntegrStep = @(dTimestamp, dQuat0, dOmega0, dStepSize, varTorque) CRigidBodyDynamicsIntegrator.integrateStep_(dTimestamp, ...
                                                                                                                    self.dInertiaMatrix, ...
                                                                                                                    dQuat0, ...
                                                                                                                    dOmega0, ...
                                                                                                                    dStepSize, ...
                                                                                                                    varTorque, ...
                                                                                                                    bDecoupledKinoDynamics, ...
                                                                                                                    enumMethod);
                otherwise
                    error('Unsupported rigid-body integration method: %s', enumMethod);
            end

            % Initialize variables at t0
            dTmpQuat = dQuat0;
            dTmpOmega = dOmega0;

            dCurrentTime = dT0;
            dQuatOutSeq(:,1)  = dQuat0;
            dOmegaOutSeq(:,1) = dOmega0;
            dTimegridOut(1)   = dCurrentTime;
            ui32CurrentIntegrStepIdx = uint32(2);

            % Integration cycle from t0 to tf
            while dCurrentTime < dTf

                % Determine internal grid step
                dInternalLoopTime_ = dCurrentTime;
                dNextTargetTime = dTimegrid(ui32CurrentIntegrStepIdx);

                % Internal loop (over single step between timegrid entries
                dAccumStepTime = 0.0;
                while (dNextTargetTime - dInternalLoopTime_) > 1.5 * eps
                
                    % Update integration time step
                    dInternalStepInterval_ = dTimegrid(ui32CurrentIntegrStepIdx) - dTimegrid(ui32CurrentIntegrStepIdx-1);

                    if ui32CurrentIntegrStepIdx < ui32NumTimes && bDeduceDeltaT
                        dDeltaT_ = min(dInternalStepInterval_, dDefaultMaxDeltaT);

                    elseif not(bDeduceDeltaT)
                        % Not constrained by default max delta T
                        dDeltaT_ = dDeltaT;
                    end

                    dTmpDeltaStep = min(dDeltaT_, dTimegrid(ui32CurrentIntegrStepIdx) - dInternalLoopTime_);
                    
                    % Integrate over dTmpDeltaStep time (actual integration step)
                    [dTmpQuat, dTmpOmega] = fcnIntegrStep(dInternalLoopTime_, ...
                                                          dTmpQuat, ...
                                                          dTmpOmega, ...
                                                          dTmpDeltaStep, ...
                                                          varTorque);

                    dInternalLoopTime_ = dInternalLoopTime_ + dTmpDeltaStep;
                    dAccumStepTime = dAccumStepTime + dTmpDeltaStep;
                end

                % Update current index and timegrid
                dCurrentTime = round(dInternalLoopTime_, 16); % dTimegrid(ui32CurrentIntegrStepIdx) + dAccumStepTime;

                % Store in sequence
                % dTmpQuatOut = QuatKinematicsIntegrator.NormalizeSeq(dTmpQuatOut);
                dQuatOutSeq(:, ui32CurrentIntegrStepIdx) = dTmpQuat;
                dOmegaOutSeq(:, ui32CurrentIntegrStepIdx) = dTmpOmega;
                dTimegridOut(ui32CurrentIntegrStepIdx) = dCurrentTime;

                ui32CurrentIntegrStepIdx = ui32CurrentIntegrStepIdx + uint32(1);
            end
        end

    end

    methods (Static)
        %%% RHS components
        % Angular acceleration RHS
        function dAngAccel = EvalRHS_AngAccel_(dInertiaMatrix, ...
                                               dOmegaAngVel, ...
                                               dExtTorque) %#codegen
            arguments
                dInertiaMatrix  (3,3) double {mustBeFinite}
                dOmegaAngVel    (3,1) double {mustBeFinite}
                dExtTorque      (3,1) double {mustBeFinite} = zeros(3,1);
            end
            arguments

            end

            % Compute angular acceleration alpha = I^-1 * (T - w x Iw)
            dAngAccel =  dInertiaMatrix \ (dExtTorque - cross(dOmegaAngVel, dInertiaMatrix * dOmegaAngVel));
        end

        function [dQuatOut, dOmegaOut] = integrateStep_(dTimestamp, ...
                                                        dInertiaMatrix, ...
                                                        dQuat0, ...
                                                        dOmega0, ...
                                                        dDeltaT, ...
                                                        varTorque, ...
                                                        bDecoupledKinoDynamics, ...
                                                        enumMethod) %#codegen
            %% SIGNATURE
            % [dQuatOut, dOmegaOut] = integrateStep_(dTimestamp, dInertiaMatrix, dQuat0, dOmega0, dDeltaT, ...
            %     varTorque, bDecoupledKinoDynamics, enumMethod)
            % -----------------------------------------------------------------------------------------------------
            %% DESCRIPTION
            % Advance one rigid-body step. Constant body torque with RK4-RKMK4 uses the stateless primitive;
            % callback and alternate-method requests retain the compatibility implementation.
            % -----------------------------------------------------------------------------------------------------
            %% INPUT
            % dTimestamp                 (1,1) double step start epoch [s]
            % dInertiaMatrix             (3,3) double body-frame inertia matrix
            % dQuat0                     (4,1) double initial passive quaternion
            % dOmega0                    (3,1) double initial body-frame angular velocity [rad/s]
            % dDeltaT                    (1,1) double positive integration step [s]
            % varTorque                  function handle or constant body-frame torque
            % bDecoupledKinoDynamics     logical true when torque is independent of attitude
            % enumMethod                 char integration method
            % -----------------------------------------------------------------------------------------------------
            %% OUTPUT
            % dQuatOut                   (4,1) double propagated passive quaternion
            % dOmegaOut                  (3,1) double propagated body-frame angular velocity [rad/s]
            % -----------------------------------------------------------------------------------------------------
            %% CHANGELOG
            % 23-08-2026  Pietro Califano, Codex     Add the stateless constant-torque RKMK4 path.
            % -----------------------------------------------------------------------------------------------------
            %% DEPENDENCIES
            % IntegrateRigidBodyRK4RKMK4Step
            % -----------------------------------------------------------------------------------------------------

            arguments
                dTimestamp              (1,1) double {mustBeGreaterThanOrEqual(dTimestamp, 0.0)}
                dInertiaMatrix          (3,3) double {mustBeNumeric}
                dQuat0                  (4,1) double {mustBeFinite}
                dOmega0                 (3,1) double {mustBeFinite}
                dDeltaT                 (1,1) double {mustBePositive}
                varTorque               {mustBeA(varTorque, ["function_handle", "double"])} = zeros(3,1)
                bDecoupledKinoDynamics  (1,1) logical {coder.mustBeConst} = true % True if torque does not depend on attitude
                enumMethod              (1,:) char {coder.mustBeConst, mustBeMember(enumMethod, {'lie_euler', 'rk4_rkmk4', 'rk4'})} = 'rk4_rkmk4'
            end
            % Combined integration step (decoupled or coupled)

            assert(all(diag(dInertiaMatrix) > 0.0), 'ERROR: invalid inertia matrix. Diagonal cannot be non-positive.')

            % Return before constructing callback closures for the common
            % decoupled constant-torque RKMK4 path.
            if coder.const(bDecoupledKinoDynamics && ...
                    strcmp(enumMethod, 'rk4_rkmk4') && isa(varTorque, "double"))
                [dQuatOut, dOmegaOut] = IntegrateRigidBodyRK4RKMK4Step( ...
                    dInertiaMatrix, dQuat0, dOmega0, varTorque, dDeltaT);
                return
            end

            if coder.const(isa(varTorque, "function_handle"))
                dvTmpTorque_ = varTorque;
            else
                dvTmpTorque_ = @(dTstamp, dOmega, dQuat0) varTorque; % Constant, defaults to zero
            end
            fcnEvalOmegaRHS = @(dTstamp, dOmega, dQuat0) CRigidBodyDynamicsIntegrator.EvalRHS_AngAccel_(dInertiaMatrix, ...
                                                                                                       dOmega, ...
                                                                                                       dvTmpTorque_(dTstamp, dOmega, dQuat0));

            if coder.const(bDecoupledKinoDynamics)
                switch coder.const(enumMethod)
                    case 'lie_euler'
                        % Lie group Euler step: integrate angular velocity with explicit Euler and update quaternion with Lie group Euler step using the initial angular velocity sample.
                        [dOmegaOut, dQuatOut] = CRigidBodyDynamicsIntegrator.IntegrStep_LieEuler(dOmega0, ...
                                                                                                  dQuat0, ...
                                                                                                  fcnEvalOmegaRHS, ...
                                                                                                  dTimestamp, ...
                                                                                                  dDeltaT);
                    case 'rk4'
                        % RK4 step: integrate angular velocity with RK4 and update quaternion with RK4 using the same angular velocity samples at the stages.
                        [dOmegaOut, dQuatOut] = CRigidBodyDynamicsIntegrator.IntegrStep_RK4(dOmega0, ...
                                                                                            dQuat0, ...
                                                                                            fcnEvalOmegaRHS, ...
                                                                                            dTimestamp, ...
                                                                                            dDeltaT);
                    case 'rk4_rkmk4'
                        [dOmegaOut, dQuatOut] = ...
                            CRigidBodyDynamicsIntegrator.IntegrStep_RK4_RKMK4( ...
                                dOmega0, dQuat0, fcnEvalOmegaRHS, ...
                                dTimestamp, dDeltaT);
                end
            else
                error('Not yet implemented')
            end
        end

        function [dNextOmega, dNextQuat] = IntegrStep_LieEuler(dOmega0, ...
                                                               dQuat0, ...
                                                               fcnEvalOmegaRHS, ...
                                                               dTstamp, ...
                                                               dDeltaTime)
            arguments
                dOmega0         (3,1) double {mustBeFinite}
                dQuat0          (4,1) double {mustBeFinite}
                fcnEvalOmegaRHS function_handle % Expected signature: fcnEvalOmegaRHS(dTstamp, dOmega, dQuat0)
                dTstamp         (1,1) double {mustBeFinite}
                dDeltaTime      (1,1) double {mustBePositive}
            end

            dTmpRHS = fcnEvalOmegaRHS(dTstamp, dOmega0, dQuat0);
            dNextOmega = dOmega0 + dDeltaTime * dTmpRHS;
            fcnOmega = @(~) dOmega0;
            dNextQuat = CQuatKinematicsIntegrator.IntegrStep_LieGroupEuler(dQuat0, fcnOmega, dTstamp, dDeltaTime);
        end

        function [dNextOmega, dNextQuat, dOmegaStagesVals, dStagesTimes] = IntegrStep_RK4(dOmega0, ...
                                                                                         dQuat0, ...
                                                                                         fcnEvalOmegaRHS, ...
                                                                                         dTstamp, ...
                                                                                         dDeltaTime)
            arguments
                dOmega0         (3,1) double {mustBeFinite}
                dQuat0          (4,1) double {mustBeFinite}
                fcnEvalOmegaRHS function_handle % Expected signature: fcnEvalOmegaRHS(dTstamp, dOmega, dQuat0)
                dTstamp         (1,1) double {mustBeFinite}
                dDeltaTime      (1,1) double {mustBePositive}
            end

            [dNextOmega, dOmegaStagesVals, dStagesTimes] = CRigidBodyDynamicsIntegrator.IntegrStep_RK4_OmegaStages_(dOmega0, ...
                                                                                                                     dQuat0, ...
                                                                                                                     fcnEvalOmegaRHS, ...
                                                                                                                     dTstamp, ...
                                                                                                                     dDeltaTime);
            dNextQuat = CQuatKinematicsIntegrator.IntegrStep_RK4_Samples(dQuat0, dOmegaStagesVals, dDeltaTime);
        end

        function [dNextOmega, dNextQuat, dOmegaStagesVals, dStagesTimes] = IntegrStep_RK4_RKMK4(dOmega0, ...
                                                                              dQuat0, ...
                                                                              fcnEvalOmegaRHS, ...
                                                                              dTstamp, ...
                                                                              dDeltaTime)
            %% SIGNATURE
            % [dNextOmega, dNextQuat, dOmegaStagesVals, dStagesTimes] = IntegrStep_RK4_RKMK4( ...
            %     dOmega0, dQuat0, fcnEvalOmegaRHS, dTstamp, dDeltaTime)
            % -----------------------------------------------------------------------------------------------------
            %% DESCRIPTION
            % Advance one callback-driven rigid-body step with classical RK4 angular-velocity stages and the
            % canonical RKMK4 quaternion primitive. This compatibility path supports acceleration callbacks;
            % constant body torque uses IntegrateRigidBodyRK4RKMK4Step through integrateStep_.
            % -----------------------------------------------------------------------------------------------------
            %% INPUT
            % dOmega0             (3,1) double initial body-frame angular velocity [rad/s]
            % dQuat0              (4,1) double initial passive scalar-first Hamilton quaternion
            % fcnEvalOmegaRHS     function handle returning body-frame angular acceleration [rad/s^2]
            % dTstamp             (1,1) double step start epoch [s]
            % dDeltaTime          (1,1) double positive integration step [s]
            % -----------------------------------------------------------------------------------------------------
            %% OUTPUT
            % dNextOmega          (3,1) double propagated body-frame angular velocity [rad/s]
            % dNextQuat           (4,1) double propagated normalized passive quaternion
            % dOmegaStagesVals    (3,4) double body-frame angular-velocity stage values [rad/s]
            % dStagesTimes        (1,4) double RK4 stage epochs [s]
            % -----------------------------------------------------------------------------------------------------
            %% CHANGELOG
            % 24-08-2026  Pietro Califano, Codex     Delegate quaternion propagation to the stateless primitive.
            % -----------------------------------------------------------------------------------------------------
            %% DEPENDENCIES
            % CRigidBodyDynamicsIntegrator.IntegrStep_RK4_OmegaStages_, IntegrateQuaternionRKMK4Step
            % -----------------------------------------------------------------------------------------------------

            arguments
                dOmega0         (3,1) double {mustBeFinite}
                dQuat0          (4,1) double {mustBeFinite}
                fcnEvalOmegaRHS  function_handle % Expected signature: fcnEvalOmegaRHS(dTstamp, dOmega, dQuat0)
                dTstamp         (1,1) double {mustBeFinite}
                dDeltaTime      (1,1) double {mustBePositive}
            end
            [dNextOmega, dOmegaStagesVals, dStagesTimes] = CRigidBodyDynamicsIntegrator.IntegrStep_RK4_OmegaStages_(dOmega0, ...
                                                                                                                     dQuat0, ...
                                                                                                                     fcnEvalOmegaRHS, ...
                                                                                                                     dTstamp, ...
                                                                                                                     dDeltaTime);
            dNextQuat = IntegrateQuaternionRKMK4Step(dQuat0, dOmegaStagesVals, dDeltaTime);
        end
    end

    methods (Static, Access = protected)
        function [dNextOmega, dOmegaStagesVals, dStagesTimes] = IntegrStep_RK4_OmegaStages_(dOmega0, ...
                                                                                           dQuat0, ...
                                                                                           fcnEvalOmegaRHS, ...
                                                                                           dTstamp, ...
                                                                                           dDeltaTime)
            arguments
                dOmega0         (3,1) double {mustBeFinite}
                dQuat0          (4,1) double {mustBeFinite}
                fcnEvalOmegaRHS function_handle
                dTstamp         (1,1) double {mustBeFinite}
                dDeltaTime      (1,1) double {mustBePositive}
            end

            dStagesTimes = dTstamp + [0.0, 0.5 * dDeltaTime, 0.5 * dDeltaTime, dDeltaTime];

            % Stage 1
            dTmpRHS_K1 = fcnEvalOmegaRHS(dStagesTimes(1), dOmega0, dQuat0);
            dOmega_K1 = dOmega0 + 0.5 * dDeltaTime * dTmpRHS_K1;
            dQuat_K1 = dQuat0;

            % Stage 2
            dTmpRHS_K2 = fcnEvalOmegaRHS(dStagesTimes(2), dOmega_K1, dQuat_K1); 
            dOmega_K2 = dOmega0 + 0.5 * dDeltaTime * dTmpRHS_K2;
            dQuat_K2 = dQuat0;

            % Stage 3
            dTmpRHS_K3 = fcnEvalOmegaRHS(dStagesTimes(3), dOmega_K2, dQuat_K2);
            dOmega_K3 = dOmega0 + dDeltaTime * dTmpRHS_K3;
            dQuat_K3 = dQuat0;

            % Stage 4
            dTmpRHS_K4 = fcnEvalOmegaRHS(dStagesTimes(4), dOmega_K3, dQuat_K3);
            
            % Define output
            dOmegaStagesVals = [dOmega0, dOmega_K1, dOmega_K2, dOmega_K3];

            % Update solution
            dOneOverSix = coder.const(1/6);
            dNextOmega =  dOmega0 + dDeltaTime * dOneOverSix * (dTmpRHS_K1 + 2*dTmpRHS_K2 + 2*dTmpRHS_K3 + dTmpRHS_K4);
        end
    end
end
