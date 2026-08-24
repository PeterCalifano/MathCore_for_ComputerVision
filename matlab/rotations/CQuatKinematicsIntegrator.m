classdef CQuatKinematicsIntegrator < handle & matlab.mixin.Copyable
    %% SIGNATURE
    % objIntegrator = CQuatKinematicsIntegrator()
    % -------------------------------------------------------------------------------------------------------------
    %% DESCRIPTION
    % Class containing methods to perform integration of quaternion kinematics (Hamilton quaternion
    % convention) using on-manifold or classical integration schemes. Additional utility methods are
    % provided (quaternion multiplication, normalization). Step methods are implemented as static methods.
    % -------------------------------------------------------------------------------------------------------------
    %% INPUT
    % None.
    % -------------------------------------------------------------------------------------------------------------
    %% OUTPUT
    % objIntegrator    (1,1) CQuatKinematicsIntegrator
    % -------------------------------------------------------------------------------------------------------------
    %% CHANGELOG
    % 14-07-2025    Pietro Califano, GPT o4-mini-high       First prototype implementation
    % 15-07-2025    Pietro Califano                         Complete implementation with RK4 and LieGroup methods; 
    %                                                       add capability to work with constant, discrete and fcn 
    %                                                       handle angular velocity.
    % 16-07-2025    Pietro Califano                         Update integration
    % 23-08-2026    Pietro Califano, Codex                  Extract corrected RKMK4 step into a stateless primitive.
    % 24-08-2026    Pietro Califano, Codex                  Remove the redundant sampled-step adapter.
    % -------------------------------------------------------------------------------------------------------------
    %% DEPENDENCIES
    % qCross, IntegrateQuaternionRKMK4Step
    % -------------------------------------------------------------------------------------------------------------

    properties

    end

    methods
        function self = CQuatKinematicsIntegrator()
            % Default constructor
        end

        function [dTmpQuatOut, dTimegridOut, dQuatOutSeq] = integrate(self, ...
                                                                    dQuat0, ...
                                                                    varOmegaAngVel, ...
                                                                    dTimegrid, ...
                                                                    enumMethod, ...
                                                                    dDeltaT, ...
                                                                    dDefaultMaxDeltaT, ...
                                                                    dAngVelTimegrid, ...
                                                                    enumInterpMethod)
            arguments
                self                (1,1) CQuatKinematicsIntegrator
                dQuat0              (4,1) double {mustBeFinite}
                varOmegaAngVel      {mustBeA(varOmegaAngVel, ["function_handle", "double"])}
                dTimegrid           (1,:) double {mustBeFinite, mustBeVector}
                enumMethod          (1,:) char {coder.mustBeConst, mustBeMember(enumMethod, {'lie_euler', 'rkmk4', 'rk4'})}
                dDeltaT             (1,1) double {mustBeScalarOrEmpty, mustBeGreaterThanOrEqual(dDeltaT, 0.0)} = 0.0;
                dDefaultMaxDeltaT   (1,1) double {mustBeScalarOrEmpty, mustBeGreaterThanOrEqual(dDefaultMaxDeltaT, 0.0)} = 1.0;
                dAngVelTimegrid     = [];
                enumInterpMethod    char {coder.mustBeConst, mustBeMember(enumInterpMethod, ["linear", "spline"])} = "linear";
            end

            bDeduceDeltaT = dDeltaT == 0.0;
            % Initialize quaternion and timestamps
            dTmpQuatOut = CQuatKinematicsIntegrator.NormalizeSeq(dQuat0);

            dT0 = dTimegrid(1);
            dTf = dTimegrid(end);

            if size(dTimegrid,2) == 2
                assert(dDeltaT > 0.0, 'ERROR: timestep must be provided if timegrid is given as [t0,tf].')
                dTimegrid = dT0:dDeltaT:dTf;
                bDeduceDeltaT = false;
            end

            % Allocate output sequence array
            ui32NumSteps = size(dTimegrid, 2);
            dQuatOutSeq = zeros(4, ui32NumSteps);
            dTimegridOut = zeros(1, ui32NumSteps);

            % Determine angular velocity mode
            if coder.const(isa(varOmegaAngVel, "double"))

                if size(varOmegaAngVel, 2) == 1
                    % Constant angular velocity case
                    varOmegaAngVel_ = @(dTstamp) varOmegaAngVel;
                else
                    % Tabulated angular velocity case. Interpolation function is defined based on input timegrid and method.
                    if isempty(dAngVelTimegrid)
                        dAngVelTimegrid = dTimegrid;
                        assert(size(varOmegaAngVel, 2) == ui32NumSteps, 'ERROR: angular velocity profile size must match size of timegrid.');
                    else
                        assert(size(varOmegaAngVel, 2) == length(dAngVelTimegrid), ...
                            'ERROR: angular velocity profile size must match size of input timegrid dAngVelTimegrid.');
                    end

                    dAngVelTimegrid = transpose(dAngVelTimegrid(:));
                    assert(all(diff(dAngVelTimegrid) > 0.0), ...
                        'ERROR: angular velocity timegrid must be strictly increasing.');

                    if numel(dAngVelTimegrid) == 1
                        % Single time point provided for tabulated angular velocity: assume constant
                        dOmegaConst = varOmegaAngVel(:,1);
                        varOmegaAngVel_ = @(~) dOmegaConst;
                    else
                        % Interpolation of angular velocity profile
                        objOmegaInterp = griddedInterpolant(dAngVelTimegrid(:), ...
                                                           transpose(varOmegaAngVel), ...
                                                           char(enumInterpMethod), ...
                                                           'none');
                        varOmegaAngVel_ = @(dTstamp) transpose(objOmegaInterp(dTstamp));
                    end

                end

            elseif coder.const(isa(varOmegaAngVel, "function_handle"))
                varOmegaAngVel_ = varOmegaAngVel;
            end

            % Select implementation of step
            switch coder.const(enumMethod)
                case 'lie_euler'
                    % Lie group Euler method: compute delta quaternion from angular velocity at the current time
                    fcnQuatIntegr = @(dQuatOut, dOmegaAngVel, dTstamp, dTmpStep) CQuatKinematicsIntegrator.IntegrStep_LieGroupEuler(dQuatOut, dOmegaAngVel, dTstamp, dTmpStep);

                case 'rkmk4'
                    % RKMK4 method: compute delta quaternion from a combination of angular velocity samples at the RK4 stages
                    fcnQuatIntegr = @(dQuatOut, dOmegaAngVel, dTstamp, dTmpStep) CQuatKinematicsIntegrator.IntegrStep_RKMK4(dQuatOut, ...
                                                                                                                dOmegaAngVel, ...
                                                                                                                dTstamp, ...
                                                                                                                dTmpStep);

                case 'rk4'
                    % Classical RK4 method: compute quaternion increment without on-manifold operations
                    fcnQuatIntegr = @(dQuatOut, dOmegaAngVel, dTstamp, dTmpStep) CQuatKinematicsIntegrator.IntegrStep_RK4(dQuatOut, dOmegaAngVel, dTstamp, dTmpStep);

            end

            % Initialize variables
            dCurrentTime = dT0;
            dQuatOutSeq(:,1) = dTmpQuatOut;
            dTimegridOut(1)  = dCurrentTime;
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

                    if ui32CurrentIntegrStepIdx < ui32NumSteps && bDeduceDeltaT
                        dDeltaT_ = min(dInternalStepInterval_, dDefaultMaxDeltaT);

                    elseif not(bDeduceDeltaT)
                        % Not constrained by default max delta T
                        dDeltaT_ = dDeltaT;
                    end

                    dTmpDeltaStep = min(dDeltaT_, dTimegrid(ui32CurrentIntegrStepIdx) - dInternalLoopTime_);
                    
                    % Integrate over dTmpDeltaStep time
                    dTmpQuatOut = fcnQuatIntegr(dTmpQuatOut, varOmegaAngVel_, dInternalLoopTime_, dTmpDeltaStep);

                    % Update time counters
                    dInternalLoopTime_ = dInternalLoopTime_ + dTmpDeltaStep;
                    dAccumStepTime = dAccumStepTime + dTmpDeltaStep;
                end
                

                % Update current index and timegrid
                dCurrentTime = round(dInternalLoopTime_, 16); % dTimegrid(ui32CurrentIntegrStepIdx) + dAccumStepTime;

                % Store in sequence
                % dTmpQuatOut = QuatKinematicsIntegrator.NormalizeSeq(dTmpQuatOut);
                dQuatOutSeq(:, ui32CurrentIntegrStepIdx) = dTmpQuatOut;
                dTimegridOut(ui32CurrentIntegrStepIdx) = dCurrentTime;

                ui32CurrentIntegrStepIdx = ui32CurrentIntegrStepIdx + uint32(1);
            end

        end
    end

    methods (Static)

        function dQuatOut = ExpMap(dDeltaAngle)
            arguments
                dDeltaAngle (3,:) double {mustBeFinite}
            end
            % Function applying delta rotation to a quaternion using ExpMap operation for multiple samples in parallel.
            
            dCurrentTimeTheta = vecnorm(dDeltaAngle, 2, 1);
            ui32NumSamples = size(dDeltaAngle, 2);
            dQuatOut = zeros(4, ui32NumSamples);

            % Handle small angle case with Taylor expansion to avoid numerical issues
            bSmallAngle = dCurrentTimeTheta < 1e-12;
            dQuatOut(1, bSmallAngle) = 1.0;

            % Handle non-small angle case with standard ExpMap formula
            for ui32Idx = find(~bSmallAngle)
                
                % Compute rotation axis and angle
                dvTmpAxis = dDeltaAngle(:, ui32Idx) / dCurrentTimeTheta(ui32Idx);
                dTmpHalfTheta = 0.5 * dCurrentTimeTheta(ui32Idx);

                % Apply ExpMap formula for quaternion update
                dQuatOut(:, ui32Idx) = [cos(dTmpHalfTheta); sin(dTmpHalfTheta) * dvTmpAxis];
            
            end
        
        end

        function dQuatOut = expMap(dDeltaAngle)
            arguments
                dDeltaAngle (3,1) double {mustBeFinite}
            end
            % Function applying delta rotation to a quaternion using ExpMap operation for a single sample.

            % Compute rotation angle 
            dCurrentDeltaTheta = norm(dDeltaAngle);

            if dCurrentDeltaTheta < 1e-12
                dQuatOut = [1.0; 0.0; 0.0; 0.0];
                return;
            end

            % Apply ExpMap formula for quaternion update
            dvTmpAxis = dDeltaAngle / dCurrentDeltaTheta;
            dQuatOut = [cos(dCurrentDeltaTheta); sin(dCurrentDeltaTheta) * dvTmpAxis];
        end


        function dQuatSeqOut = QuatSeqCross(dQuatSeq1, dQuatSeq2, bIS_VSRPplus)
            arguments
                dQuatSeq1 (4,:) double {mustBeFinite}
                dQuatSeq2 (4,:) double {mustBeFinite}
                bIS_VSRPplus (1,1) logical {isscalar, islogical} = false
            end
            % Quaternion convention definition
            % (SV) Scalar first, Vector last
            % (P) Passive
            % (R) Successive coordinate transformations have the unmodified quaternion chain on the Right side of
            %     the triple product.
            % (plus) Right-Handed Rule for the imaginary numbers i, j, k. (aka Hamilton)

            dQuatSeqOut = zeros(size(dQuatSeq1));

            for idQ = 1:size(dQuatSeq1, 2)
                dQuatSeqOut(:,idQ) = qCross(dQuatSeq1(:,idQ), ...
                                        dQuatSeq2(:,idQ), ...
                                        bIS_VSRPplus);
            end

        end

        function dQuatOut = NormalizeSeq(dQuat)
            arguments
                dQuat (4,:) double {mustBeFinite, mustBeNumeric}
            end
            dQuatOut = dQuat ./ vecnorm(dQuat);
        end

        function dQuatOut = IntegrStep_RK4(dQuat0, ...
                                        fcnEvalOmegaAngVel, ...
                                        dTstamp, ...
                                        dDeltaTime)
            arguments
                dQuat0          (4,1) double {mustBeFinite}
                fcnEvalOmegaAngVel   function_handle % Gives angular velocity at each time
                dTstamp         (1,1) double {mustBeFinite}
                dDeltaTime      (1,1) double {mustBePositive}
            end
            %%% IntegrStep_RK4  Single RK4 step for quaternion kinematics with constant omega
            %
            % Butcher table:
            %    A = [0    0    0    0;
            %         1/2  0    0    0;
            %         0    1/2  0    0;
            %         0    0    1    0];
            %    b = [1/6; 1/3; 1/3; 1/6];
            %    c = [0; 1/2; 1/2; 1];

            dOmegaStages = [ ...
                reshape(fcnEvalOmegaAngVel(dTstamp), 3, 1), ...
                reshape(fcnEvalOmegaAngVel(dTstamp + 0.5 * dDeltaTime), 3, 1), ...
                reshape(fcnEvalOmegaAngVel(dTstamp + 0.5 * dDeltaTime), 3, 1), ...
                reshape(fcnEvalOmegaAngVel(dTstamp + dDeltaTime), 3, 1)];

            dQuatOut = CQuatKinematicsIntegrator.IntegrStep_RK4_Samples(dQuat0, dOmegaStages, dDeltaTime);
        end

        function dQuatOut = IntegrStep_RK4_Samples(dQuat0, dOmegaStages, dDeltaTime)
            arguments
                dQuat0      (4,1) double {mustBeFinite}
                dOmegaStages (3,4) double {mustBeFinite}
                dDeltaTime  (1,1) double {mustBePositive}
            end

            fcnOmegaMat = @(dvOmegaAngVel) [0.0,              -dvOmegaAngVel(1), -dvOmegaAngVel(2), -dvOmegaAngVel(3); ...
                                            dvOmegaAngVel(1),  0.0,               dvOmegaAngVel(3), -dvOmegaAngVel(2); ...
                                            dvOmegaAngVel(2), -dvOmegaAngVel(3),  0.0,              dvOmegaAngVel(1); ...
                                            dvOmegaAngVel(3),  dvOmegaAngVel(2), -dvOmegaAngVel(1),  0.0];

            dvTmpK1 = 0.5 * fcnOmegaMat(dOmegaStages(:,1)) * dQuat0;
            dvTmpQuat = dQuat0 + 0.5 * dDeltaTime * dvTmpK1;

            dvTmpK2 = 0.5 * fcnOmegaMat(dOmegaStages(:,2)) * dvTmpQuat;
            dvTmpQuat = dQuat0 + 0.5 * dDeltaTime * dvTmpK2;

            dvTmpK3 = 0.5 * fcnOmegaMat(dOmegaStages(:,3)) * dvTmpQuat;
            dvTmpQuat = dQuat0 + dDeltaTime * dvTmpK3;

            dvTmpK4 = 0.5 * fcnOmegaMat(dOmegaStages(:,4)) * dvTmpQuat;

            dQuatOut = dQuat0 + dDeltaTime * (dvTmpK1 + 2*dvTmpK2 + 2*dvTmpK3 + dvTmpK4)/6;
            dQuatOut = CQuatKinematicsIntegrator.NormalizeSeq(dQuatOut);
        end
        
        function dQuatOut = IntegrStep_LieGroupEuler(dQuat0, ...
                                                    fcnEvalOmegaAngVel, ...
                                                    dTstamp, ...
                                                    dDeltaTime)
            arguments
                dQuat0          (4,1) double {mustBeFinite}
                fcnEvalOmegaAngVel   function_handle % Gives angular velocity at each time
                dTstamp         (1,1) double {mustBeFinite}
                dDeltaTime      (1,1) double {mustBePositive}
            end
            % Function performing an integration setp assuming constant angular velocity and computing the
            % delta quaternion accordingly over the dDeltaTime timestep. ExpMap operation is then used to
            % rotate the initial quaternion.
            
            % Compute delta quaternion over the time interval
            dTmpDeltaAngle = 0.5 * reshape(fcnEvalOmegaAngVel(dTstamp), 3, 1) * dDeltaTime;

            % Use qdot = 0.5 * Omega(omega) * q, equivalent to right-applying the incremental quaternion.
            dQuatOut = CQuatKinematicsIntegrator.QuatSeqCross(dQuat0, CQuatKinematicsIntegrator.expMap(dTmpDeltaAngle));

        end

        function dQuatOut = IntegrStep_RKMK4(dQuat0, ...
                                            fcnEvalOmegaAngVel, ...
                                            dTstamp, ...
                                            dDeltaTime)
            %% SIGNATURE
            % dQuatOut = IntegrStep_RKMK4(dQuat0, fcnEvalOmegaAngVel, dTstamp, dDeltaTime)
            % -----------------------------------------------------------------------------------------------------
            %% DESCRIPTION
            % Sample a prescribed body-frame angular-velocity function at the classical RK4 nodes and advance one
            % passive scalar-first Hamilton quaternion with the canonical stateless RKMK4 primitive.
            % -----------------------------------------------------------------------------------------------------
            %% INPUT
            % dQuat0                 (4,1) double initial passive quaternion
            % fcnEvalOmegaAngVel     function handle returning body-frame angular velocity [rad/s]
            % dTstamp                (1,1) double step start epoch [s]
            % dDeltaTime             (1,1) double positive integration step [s]
            % -----------------------------------------------------------------------------------------------------
            %% OUTPUT
            % dQuatOut               (4,1) double propagated normalized quaternion
            % -----------------------------------------------------------------------------------------------------
            %% CHANGELOG
            % 24-08-2026  Pietro Califano, Codex     Delegate directly to the stateless RKMK4 primitive.
            % -----------------------------------------------------------------------------------------------------
            %% DEPENDENCIES
            % IntegrateQuaternionRKMK4Step
            % -----------------------------------------------------------------------------------------------------

            arguments
                dQuat0          (4,1) double {mustBeFinite}
                fcnEvalOmegaAngVel   function_handle % Gives angular velocity at each time
                dTstamp         (1,1) double {mustBeFinite}
                dDeltaTime      (1,1) double {mustBePositive}
            end
            %%% RKMK4 butcher table
            % A-matrix (stage coefficients)
            % A = [ ...
            %     0,   0,   0,   0; ...
            %     1/2, 0,   0,   0; ...
            %     0,   1/2, 0,   0; ...
            %     0,   0,   1,   0  ...
            % ];

            % b-vector (weights for final combination)
            % B = [ ...
            %     1/6; ...
            %     1/3; ...
            %     1/3; ...
            %     1/6  ...
            % ];

            % c-vector (nodes / time fractions)
            % C = [ ...
            %     0; ...
            %     1/2; ...
            %     1/2; ...
            %     1    ...
            % ];

            % Compute angular velocity samples at the RK4 stages
            dOmegaStages = [ ...
                reshape(fcnEvalOmegaAngVel(dTstamp), 3, 1), ...
                reshape(fcnEvalOmegaAngVel(dTstamp + 0.5 * dDeltaTime), 3, 1), ...
                reshape(fcnEvalOmegaAngVel(dTstamp + 0.5 * dDeltaTime), 3, 1), ...
                reshape(fcnEvalOmegaAngVel(dTstamp + dDeltaTime), 3, 1)];

            assert(all(isfinite(dOmegaStages), 'all'), ...
                'ERROR: nan detected in integration step. Interpolant of angular velocity may have failed.');

            % Delegate sampled Lie-group integration to the canonical
            % code-generation-safe numerical primitive.
            dQuatOut = IntegrateQuaternionRKMK4Step(dQuat0, dOmegaStages, dDeltaTime);

        end
    end
end
