%% RigidBodyDynamicsIntegrator.m
classdef CRigidBodyDynamicsIntegrator < handle & matlab.mixin.Copyable
    properties
        dInertiaMatrix          (3,3) double {mustBeFinite}     = zeros(3,3)
        objQuatKinIntegrator    CQuatKinematicsIntegrator       = CQuatKinematicsIntegrator()
    end

    methods (Access = public)
         % CONSTRUCTOR
        function self = CRigidBodyDynamicsIntegrator(dInertiaMatrix, objQuatInt)
            arguments
                dInertiaMatrix (3,3) double {mustBeFinite} = zeros(3,3)
                objQuatInt (1,1) CQuatKinematicsIntegrator = CQuatKinematicsIntegrator()
            end
            % Assign attributes
            self.dInertiaMatrix = dInertiaMatrix;
            self.objQuatKinIntegrator = objQuatInt;
        end

        % PUBLIC METHODS
        function [dQuatOutSeq, dOmegaOutSeq] = integrate(self, ...
                                                       dTimegrid, ...
                                                       dQuat0, ...
                                                       dOmega0, ...
                                                       varTorque, ...
                                                       dDeltaT, ...
                                                       enumMethod, ...
                                                       bDecoupledKinoDynamics, ...
                                                       dDefaultMaxDeltaT)
            arguments
                self                    (1,1) CRigidBodyDynamicsIntegrator
                dTimegrid               (1,:) double {mustBeFinite, mustBeVector}
                dQuat0                  (4,1) double {mustBeFinite}
                dOmega0                 (3,1) double {mustBeFinite}
                varTorque               {mustBeA(varTorque, ["function_handle", "double"])}
                dDeltaT                 (1,1) double {mustBeScalarOrEmpty, mustBeGreaterThanOrEqual(dDeltaT, 0.0)} = 0.0;
                enumMethod              (1,:) char {mustBeMember(enumMethod, {'lie_euler', 'rk4_rkmk4', 'rk4'})} = 'rk4_rkmk4'
                bDecoupledKinoDynamics  (1,1) {islogical, isscalar} = true % True if torque does not depend on attitude
                dDefaultMaxDeltaT       (1,1) double {mustBeScalarOrEmpty, mustBeGreaterThanOrEqual(dDefaultMaxDeltaT, 0.0)} = 1.0;
            end 

            % TODO implement integration step from t0 to tf
            % Preliminary checks
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
            switch enumMethod

                case 'rkmk4'
                    error('Coupled case with rkmk4: not implemented yet')

                case 'rk4_rkmk4'
                    % [dQuatOut, dOmegaOut]
                    fcnIntegrStep = @(dTimestamp, dQuat0, dOmega0, varTorque) CRigidBodyDynamicsIntegrator.integrateStep_(dTimestamp, ...
                                                                                                                    self.dInertiaMatrix, ...
                                                                                                                    dQuat0, ...
                                                                                                                    dOmega0, ...
                                                                                                                    dDeltaT, ...
                                                                                                                    varTorque, ...
                                                                                                                    bDecoupledKinoDynamics);
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
                    [dTmpQuat, dTmpOmega] = fcnIntegrStep(dCurrentTime, ...
                                                            dTmpQuat, ...
                                                            dTmpOmega, ...
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
                                                        bDecoupledKinoDynamics) %#codegen
            arguments
                dTimestamp              (1,1) double {mustBeGreaterThanOrEqual(dTimestamp, 0.0)}
                dInertiaMatrix          (3,3) double {ismatrix, mustBeNumeric}
                dQuat0                  (4,1) double {mustBeFinite}
                dOmega0                 (3,1) double {mustBeFinite}
                dDeltaT                 (1,1) double {mustBePositive}
                varTorque               {mustBeA(varTorque, ["function_handle", "double"])} = zeros(3,1)
                bDecoupledKinoDynamics  (1,1) {islogical, isscalar} = true % True if torque does not depend on attitude
            end
            % Combined integration step (decoupled or coupled)
            % persistent fcnEvalOmegaRHS 

            if isa(varTorque, "function_handle")
                dvTmpTorque_ = varTorque;
            else
                dvTmpTorque_ = @(dTstamp, dOmega, dQuat0) varTorque; % Constant, defaults to zero
            end

            % if isempty(fcnEvalOmegaRHS)
            % TODO how to improve need of declaring function handle here? persistent seems a bad idea. There
            % seems no way of avoiding this method to be a class method...
            assert(all(diag(dInertiaMatrix) > 0.0), 'ERROR: invalid inertia matrix. Diagonal cannot be non-positive.')
            fcnEvalOmegaRHS = @(dTstamp, dOmega, dQuat0) CRigidBodyDynamicsIntegrator.EvalRHS_AngAccel_(dInertiaMatrix, ...
                                                                                                       dOmega, ...
                                                                                                       dvTmpTorque_(dTstamp, dOmega, dQuat0));
            % end


            if bDecoupledKinoDynamics
                
                % Integrate RK4 step for angular velocity
                [dOmegaOut, ~, dOmegaStages, dStagesTimes] = CRigidBodyDynamicsIntegrator.IntegrStep_RK4_RKMK4(dOmega0, ...
                                                                                                          dQuat0, ...
                                                                                                          fcnEvalOmegaRHS, ...
                                                                                                          dTimestamp, ...
                                                                                                          dDeltaT);

                % Integrate attitude quaternion using linearly interpolated angular velocity between integration times
                % dDiffs = diff(dOmegaStages, 1,2) ~= 0.0;
                % if any(dDiffs, 'all')
                % 
                %     for id = 1:3
                % 
                %     end
                % 
                %     dOmegaInterp_ = @(dTime) [interp1(dStagesTimes, dOmegaStages(1,:), dTime, 'linear'), ...
                %                                   interp1(dStagesTimes, dOmegaStages(2,:), dTime, 'linear'), ...
                %                                   interp1(dStagesTimes, dOmegaStages(3,:), dTime, 'linear')];
                % 
                % else
                %     dOmegaInterp_ = @(dTime) dOmegaStages(:,1);
                % end

                % Build omega interpolant for the stage
                fcnOmega1 = CRigidBodyDynamicsIntegrator.getUniqueOmegaInterp_(dStagesTimes, dOmegaStages(1,:));
                fcnOmega2 = CRigidBodyDynamicsIntegrator.getUniqueOmegaInterp_(dStagesTimes, dOmegaStages(2,:));
                fcnOmega3 = CRigidBodyDynamicsIntegrator.getUniqueOmegaInterp_(dStagesTimes, dOmegaStages(3,:));

                dOmegaInterp_ = @(dTime) [ ...
                                            fcnOmega1(dTime); ...
                                            fcnOmega2(dTime); ...
                                            fcnOmega3(dTime)  ...
                                          ];

                dQuatOut = CQuatKinematicsIntegrator.IntegrStep_RKMK4(dQuat0, ...
                                                                    dOmegaInterp_, ...
                                                                    dTimestamp, ...
                                                                    dDeltaT);
            else
                error('Not yet implemented')
            end
        end

        function [dNextOmega, dNextQuat, dOmegaStagesVals, dStagesTimes] = IntegrStep_RK4_RKMK4(dOmega0, ...
                                                                              dQuat0, ...
                                                                              fcnEvalOmegaRHS, ...
                                                                              dTstamp, ...
                                                                              dDeltaTime)
            arguments
                dOmega0         (3,1) double {mustBeFinite}
                dQuat0          (4,1) double {mustBeFinite}
                fcnEvalOmegaRHS  function_handle % Expected signatfcnEvalOmegaRHS(dTstamp, dOmega, dQuat0)
                dTstamp         (1,1) double {mustBeFinite}
                dDeltaTime      (1,1) double {mustBePositive}
            end
            %%% IntegrStep_RK4  Single RK4 step for attitude dynamics of a rigid body
            %
            % Butcher table:
            %    A = [0    0    0    0;
            %         1/2  0    0    0;
            %         0    1/2  0    0;
            %         0    0    1    0];
            %    b = [1/6; 1/3; 1/3; 1/6];
            %    c = [0; 1/2; 1/2; 1];

            % DEVNOTE: current implementation only supports torque-free dynamics integration
            dQuat0   = [1;0;0;0];
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
            dNextOmega = dOmega0 + dDeltaTime * (dTmpRHS_K1 + 2*dTmpRHS_K2 + 2*dTmpRHS_K3 + dTmpRHS_K4)/6;
            dNextQuat = dQuat0;
        end
    end

    methods (Static, Access = protected)
        function [fcnOmegaInterp] = getUniqueOmegaInterp_(dStagesTimes, dOmegaStages)
            arguments
                dStagesTimes (1,:) double
                dOmegaStages (1,:) double
            end

            [~, ui32UniqueIdx] = unique(dStagesTimes, 'stable');
            dUniqueTimes    = dStagesTimes(ui32UniqueIdx);
            dOmegaUnique    = dOmegaStages(:, ui32UniqueIdx);

            if numel(dUniqueTimes) == 1

                % Single point: constant for any query
                fcnOmegaInterp = @(dTime) dOmegaUnique(1,1);

            else
                % ≥2 points: plain linear interp
                fcnOmegaInterp = @(dTime) interp1(dUniqueTimes, dOmegaUnique, dTime, 'linear');
            end
        end

    end
end
