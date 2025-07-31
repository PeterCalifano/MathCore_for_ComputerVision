classdef CQuatKinematicsIntegrator < handle & matlab.mixin.Copyable
    %% DESCRIPTION
    % Class containing methods to perform integration of quaternion kinematics (Hamilton quaternion
    % convention) using on-manifold or classical integration schemes. Additional utility methods are
    % provided (quaternion multiplication, normalization). Step methods are implemented as static methods.
    % -------------------------------------------------------------------------------------------------------------
    %% CHANGELOG
    % 14-07-2025    Pietro Califano, GPT o4-mini-high       First prototype implementation
    % 15-07-2025    Pietro Califano                         Complete implementation with RK4 and LieGroup methods; 
    %                                                       add capability to work with constant, discrete and fcn 
    %                                                       handle angular velocity.
    % 16-07-2025    Pietro Califano                         Update integration
    % -------------------------------------------------------------------------------------------------------------
    %% DEPENDENCIES
    % qCross function
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
                enumMethod          (1,:) char {mustBeMember(enumMethod, {'lie_euler', 'rkmk4', 'rk4'})}
                dDeltaT             (1,1) double {mustBeScalarOrEmpty, mustBeGreaterThanOrEqual(dDeltaT, 0.0)} = 0.0;
                dDefaultMaxDeltaT   (1,1) double {mustBeScalarOrEmpty, mustBeGreaterThanOrEqual(dDefaultMaxDeltaT, 0.0)} = 1.0;
                dAngVelTimegrid     = [];
                enumInterpMethod    char {mustBeMember(enumInterpMethod, ["linear", "spline"])} = "linear";
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
            if isa(varOmegaAngVel, "double")
                if size(varOmegaAngVel, 2) == 1
                    varOmegaAngVel_ = @(dTstamp) varOmegaAngVel;
                else
                    % TODO add interpolation of omega angular velocity to perform step (use splining)
                    
                    if isempty(dAngVelTimegrid)
                        dAngVelTimegrid = dTimegrid;
                        assert(size(varOmegaAngVel, 2) == ui32NumSteps, 'ERROR: angular velocity profile size must match size of timegrid.');
                    else
                        assert(size(varOmegaAngVel, 2) == length(dAngVelTimegrid), ...
                            'ERROR: angular velocity profile size must match size of input timegrid dAngVelTimegrid.');
                    end

                    % Define angular velocity interpolant
                    % varOmegaAngVel_ = @(dTstamp) spline(dAngVelTimegrid, varOmegaAngVel, dTstamp);
                    varOmegaAngVel_ = @(dTstamp) transpose(interp1(dAngVelTimegrid - dAngVelTimegrid(1), varOmegaAngVel', dTstamp, enumInterpMethod));

                end

            elseif isa(varOmegaAngVel, "function_handle")
                varOmegaAngVel_ = varOmegaAngVel;
            end

            % Select implementation of step
            switch enumMethod
                case 'lie_euler'
                    fcnQuatIntegr = @(dQuatOut, dOmegaAngVel, dTstamp, dTmpStep) CQuatKinematicsIntegrator.IntegrStep_LieGroupEuler(dQuatOut, dOmegaAngVel, dTstamp, dTmpStep);

                case 'rkmk4'
                    fcnQuatIntegr = @(dQuatOut, dOmegaAngVel, dTstamp, dTmpStep) CQuatKinematicsIntegrator.IntegrStep_RKMK4(dQuatOut, ...
                                                                                                                dOmegaAngVel, ...
                                                                                                                dTstamp, ...
                                                                                                                dTmpStep);

                case 'rk4'
                    fcnQuatIntegr = @(dQuatOut, dOmegaAngVel, dTstamp, dTmpStep) CQuatKinematicsIntegrator.IntegrStep_RK4(dQuatOut, dOmegaAngVel, dTstamp, dTmpStep);

            end

            % Define function handle if needed
            % varOmegaAngVel_
            
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
                    dTmpQuatOut = fcnQuatIntegr(dTmpQuatOut, varOmegaAngVel_, dCurrentTime, dTmpDeltaStep);

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
                dDeltaAngle (3,1) double {mustBeFinite}
            end

            dCurrentTimeheta = norm(dDeltaAngle);

            if dCurrentTimeheta < 1e-12
                dQuatOut = [1.0; 0.0; 0.0; 0.0];
                return;
            end
            % Function applying delta rotation to a quaternion using ExpMap operation.
            % TODO: extend to support sequences of quaternions.

            dvTmpAxis = dDeltaAngle / dCurrentTimeheta;
            dTmpHalfTheta = 0.5 * dCurrentTimeheta;
            dQuatOut = [cos(dTmpHalfTheta); sin(dTmpHalfTheta) * dvTmpAxis];
        
        end

        function dQuatOut = expMap(dDeltaAngle)
            arguments
                dDeltaAngle (3,1) double {mustBeFinite}
            end

            dCurrentDeltaTheta = norm(dDeltaAngle);

            if dCurrentDeltaTheta < 1e-12
                dQuatOut = [1.0; 0.0; 0.0; 0.0];
                return;
            end
            % Function applying delta rotation to a quaternion using ExpMap operation.
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

            % Build Omega matrix
            fcnOmegaMat = @(dvOmegaAngVel) [0.0,               -dvOmegaAngVel(1), -dvOmegaAngVel(2), -dvOmegaAngVel(3); ...
                                            dvOmegaAngVel(1), 0.0,                dvOmegaAngVel(3), -dvOmegaAngVel(2); ...
                                            dvOmegaAngVel(2), -dvOmegaAngVel(3), 0.0,               dvOmegaAngVel(1); ...
                                            dvOmegaAngVel(3), dvOmegaAngVel(2), -dvOmegaAngVel(1), 0.0];

            % Stage 1
            dTmpOmegaVelK1 = fcnEvalOmegaAngVel(dTstamp);
            dvTmpK1 = 0.5 * fcnOmegaMat(dTmpOmegaVelK1) * dQuat0; % Eval RHS at t0
            dvTmpQuat = dQuat0 + 0.5 * dDeltaTime * dvTmpK1;
            
            % Stage 2
            dTmpOmegaVelK2 = fcnEvalOmegaAngVel(dTstamp + 0.5 * dDeltaTime); 
            dvTmpK2 = 0.5 * fcnOmegaMat(dTmpOmegaVelK2) * dvTmpQuat; % Eval RHS at t0 + 0.5*dt
            dvTmpQuat = dQuat0 + 0.5 * dDeltaTime * dvTmpK2;
            
            % Stage 3
            dTmpOmegaVelK3 = fcnEvalOmegaAngVel(dTstamp + 0.5 * dDeltaTime);
            dvTmpK3 = 0.5 * fcnOmegaMat(dTmpOmegaVelK3) * dvTmpQuat; % Eval RHS at t0 + 0.5*dt
            dvTmpQuat = dQuat0 + dDeltaTime * dvTmpK3;
            
            % Stage 4
            dTmpOmegaVelK4 = fcnEvalOmegaAngVel(dTstamp + dDeltaTime);
            dvTmpK4 = 0.5 * fcnOmegaMat(dTmpOmegaVelK4) * dvTmpQuat; % Eval RHS at t0 + dt
            
            % Update solution and enforce normalization
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
            dTmpDeltaAngle = 0.5 * fcnEvalOmegaAngVel(dTstamp) * dDeltaTime;
            % Update solution
            dQuatOut = CQuatKinematicsIntegrator.QuatSeqCross( CQuatKinematicsIntegrator.expMap(dTmpDeltaAngle), dQuat0 );

        end

        function dQuatOut = IntegrStep_RKMK4(dQuat0, ...
                                            fcnEvalOmegaAngVel, ...
                                            dTstamp, ...
                                            dDeltaTime)
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

            % Compute angular velocity at each timestamp of the stages
            dInvRightJacPsi = 0.5; % Inverse right jacobian for the quaternion --> 1/2 * I since manifold is linearized at identity rotation

            % Evaluate 4 stages
            dvTmpK1 = dInvRightJacPsi * fcnEvalOmegaAngVel(dTstamp);
            dvTmpK2 = dInvRightJacPsi * fcnEvalOmegaAngVel(dTstamp + 0.5 * dDeltaTime);
            dvTmpK3 = dInvRightJacPsi * fcnEvalOmegaAngVel(dTstamp + 0.5 * dDeltaTime);
            dvTmpK4 = dInvRightJacPsi * fcnEvalOmegaAngVel(dTstamp + dDeltaTime);

            assert(all(not(isnan(dvTmpK1))) , 'ERROR: nan detected in integration step. Interpolant of angular velocity may have failed.');

            % Compute delta quaternion and update current solution
            dTmpOmegaDelta  = (dvTmpK1 + 2*dvTmpK2 + 2*dvTmpK3 + dvTmpK4) * (dDeltaTime/6);
            dQuatOut        = CQuatKinematicsIntegrator.QuatSeqCross(CQuatKinematicsIntegrator.expMap(dTmpOmegaDelta), dQuat0);

        end
    end
end

