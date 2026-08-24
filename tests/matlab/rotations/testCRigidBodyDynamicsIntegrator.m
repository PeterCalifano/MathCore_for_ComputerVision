classdef testCRigidBodyDynamicsIntegrator < matlab.unittest.TestCase
    methods (Test)
        function testDefaultConstructor(testCase)
            %% SIGNATURE
            % testDefaultConstructor(testCase)
            % -------------------------------------------------------------------------------------------------
            %% DESCRIPTION
            % Verify the default rigid-body integrator owns only its default zero inertia matrix.
            % -------------------------------------------------------------------------------------------------
            %% INPUT
            % testCase    Active MATLAB unit-test instance.
            % -------------------------------------------------------------------------------------------------
            %% OUTPUT
            % None.
            % -------------------------------------------------------------------------------------------------
            %% CHANGELOG
            % 24-08-2026  Pietro Califano, Codex     Update for the simplified constructor contract.
            % -------------------------------------------------------------------------------------------------
            %% DEPENDENCIES
            % CRigidBodyDynamicsIntegrator
            % -------------------------------------------------------------------------------------------------

            objIntegrator = CRigidBodyDynamicsIntegrator();

            testCase.verifyEqual(objIntegrator.dInertiaMatrix, zeros(3,3));
        end

        function testCustomConstructor(testCase)
            %% SIGNATURE
            % testCustomConstructor(testCase)
            % -------------------------------------------------------------------------------------------------
            %% DESCRIPTION
            % Verify construction with a supplied body-frame inertia matrix.
            % -------------------------------------------------------------------------------------------------
            %% INPUT
            % testCase    Active MATLAB unit-test instance.
            % -------------------------------------------------------------------------------------------------
            %% OUTPUT
            % None.
            % -------------------------------------------------------------------------------------------------
            %% CHANGELOG
            % 24-08-2026  Pietro Califano, Codex     Update for the simplified constructor contract.
            % -------------------------------------------------------------------------------------------------
            %% DEPENDENCIES
            % CRigidBodyDynamicsIntegrator
            % -------------------------------------------------------------------------------------------------

            dInertia = diag([2,3,4]);
            objIntegrator = CRigidBodyDynamicsIntegrator(dInertia);

            testCase.verifyEqual(objIntegrator.dInertiaMatrix, dInertia);
        end

        function testEvalRHSAngAccel(testCase)

            % Test static RHS for angular acceleration (single RHS call)
            dInertia        = eye(3);
            dOmegaAngVel    = [0.1;0.2;0.3];
            dExtTorque      = [1;0;0];

            dAngAccel = CRigidBodyDynamicsIntegrator.EvalRHS_AngAccel_(dInertia, dOmegaAngVel, dExtTorque);
            dExpected = dInertia \ (dExtTorque - cross(dOmegaAngVel, dInertia*dOmegaAngVel));

            testCase.verifyEqual(dAngAccel, dExpected, 'AbsTol',1e-12);
        end

        function testIntegrStepRK4Signature(testCase)
            % Ensure RK4 step returns correct sizes
            dQuat0      = [0;1;0;0];
            dOmega0     = zeros(3,1);
            fcnEval     = @(t,w,q) zeros(3,1);

            dTimestamp  = 0.0;
            dDeltaTime  = 0.1;
            [dNextOmega, dNextQuat, dOmegaStagesVals, dStagesTimes] = CRigidBodyDynamicsIntegrator.IntegrStep_RK4_RKMK4(dOmega0, ...
                                                                                                    dQuat0, ...
                                                                                                    fcnEval, ...
                                                                                                    dTimestamp, ...
                                                                                                    dDeltaTime);
            
            testCase.verifyEqual(dOmegaStagesVals, repmat(dNextOmega, 1, 4)); % Mocked dynamics RHS = zero --> omega should remain constant
            testCase.verifySize(dNextOmega, [3 1]);
            testCase.verifySize(dNextQuat,  [4,1]);         % 4 stages × 3 components
            testCase.verifyEqual(dNextQuat, dQuat0, 'AbsTol', 1e-12);
            testCase.verifySize(dOmegaStagesVals, [3, 4]);  % 4 stages × 3 components
            testCase.verifySize(dStagesTimes, [1 4]);
        end

        function testIntegrateStepInterfaceAllConstant(testCase)
            % Test integrateStep_ interface and output sizes
            dInertia    = eye(3);
            dQuat0      = [1;0;0;0];
            dOmega0     = zeros(3,1);
            dTimestamp  = 0.0;
            dDeltaT     = 0.05;
            varTorque   = zeros(3,1);
            % Method signature
            % [dQuatOut, dOmegaOut] = integrateStep_(dTimestamp, ...
            %                                        dInertiaMatrix, ...
            %                                        dQuat0, ...
            %                                        dOmega0, ...
            %                                        dDeltaT, ...
            %                                        varTorque, ...
            %                                        bDecoupledKinoDynamics)

            [dQuatOut, dOmegaOut] = CRigidBodyDynamicsIntegrator.integrateStep_(dTimestamp, ...
                                                                                dInertia, ...
                                                                                dQuat0, ...
                                                                                dOmega0, ...
                                                                                dDeltaT);

            testCase.verifySize(dQuatOut,  [4 1]);
            testCase.verifySize(dOmegaOut, [3 1]);
            testCase.verifyEqual(dQuatOut, dQuat0); 
            testCase.verifyEqual(dOmegaOut, dOmega0);
        end

        function testIntegrateInterfaceAllConstant(testCase)
            % Test high-level integrate method interface
            objIntegrator = CRigidBodyDynamicsIntegrator();

            dTimegrid    = linspace(0,1,6);
            dQuat0       = [1;0;0;0];
            dOmega0      = zeros(3,1);
            varTorque    = zeros(3,1);
            dDeltaT      = 0.1;

            objIntegrator.dInertiaMatrix = eye(3);

            [dQuatSeq, dOmegaSeq] = objIntegrator.integrate(dTimegrid, dQuat0, dOmega0, varTorque, dDeltaT);
            
            testCase.verifySize(dQuatSeq,  [4 numel(dTimegrid)]);
            testCase.verifySize(dOmegaSeq, [3 numel(dTimegrid)]);
        end

        function testIntegrateSupportsLieEulerAndRK4ConstantRate(testCase)
            objIntegrator = CRigidBodyDynamicsIntegrator(eye(3));

            dTimegrid = 0:0.1:1.0;
            dQuat0 = [1;0;0;0];
            dOmega0 = [0;0;pi];
            varTorque = zeros(3,1);
            expectedQuat = [cos(pi/2); 0; 0; sin(pi/2)];
            cellMethods = {'lie_euler', 'rk4'};

            for ui32Idx = 1:numel(cellMethods)
                [dQuatSeq, dOmegaSeq] = objIntegrator.integrate(dTimegrid, ...
                                                                 dQuat0, ...
                                                                 dOmega0, ...
                                                                 varTorque, ...
                                                                 0.1, ...
                                                                 cellMethods{ui32Idx}, ...
                                                                 true, ...
                                                                 1.0);

                testCase.verifyEqual(dOmegaSeq, repmat(dOmega0, 1, numel(dTimegrid)), 'AbsTol', 1e-12);
                testCase.verifyEqual(dQuatSeq(:,end), expectedQuat, 'AbsTol', 1e-5);
            end
        end

        function testIntegrateUsesActualPartialInternalStep(testCase)
            objIntegrator = CRigidBodyDynamicsIntegrator(eye(3));

            dTimegrid = [0.0, 0.25, 0.5];
            dQuat0 = [1;0;0;0];
            dOmega0 = zeros(3,1);
            varTorque = @(dTstamp, dOmega, dQuat0) [1;0;0];

            [~, dOmegaSeq] = objIntegrator.integrate(dTimegrid, ...
                                                     dQuat0, ...
                                                     dOmega0, ...
                                                     varTorque, ...
                                                     0.2, ...
                                                     'rk4_rkmk4', ...
                                                     true, ...
                                                     1.0);

            expectedOmega = [0.0, 0.25, 0.5;
                             0.0, 0.0,  0.0;
                             0.0, 0.0,  0.0];
            testCase.verifyEqual(dOmegaSeq, expectedOmega, 'AbsTol', 1e-10);
        end

        function testIntegrateDeducesInternalStepFromDefaultLimit(testCase)
            objIntegrator = CRigidBodyDynamicsIntegrator(eye(3));

            dTimegrid = [0.0, 0.25, 0.5];
            dQuat0 = [1;0;0;0];
            dOmega0 = zeros(3,1);
            varTorque = @(dTstamp, dOmega, dQuat0) [1;0;0];

            [~, dOmegaSeq] = objIntegrator.integrate(dTimegrid, ...
                                                     dQuat0, ...
                                                     dOmega0, ...
                                                     varTorque, ...
                                                     0.0, ...
                                                     'rk4_rkmk4', ...
                                                     true, ...
                                                     0.1);

            expectedOmega = [0.0, 0.25, 0.5;
                             0.0, 0.0,  0.0;
                             0.0, 0.0,  0.0];
            testCase.verifyEqual(dOmegaSeq, expectedOmega, 'AbsTol', 1e-10);
        end


        function testZeroTorqueNonZeroInitialRate(testCase)

            % With zero torque, angular rate should vary keeping angular momentum constant
            objIntegrator = CRigidBodyDynamicsIntegrator();
            objIntegrator.dInertiaMatrix = eye(3);
            objIntegrator.dInertiaMatrix(1,1) = 5;

            dTimegrid = 0:0.1:1.0;
            dQuat0    = [1;0;0;0];
            dOmega0   = [0.5; -0.3; 0.2];
            varTorque = zeros(3,1);

            [dQuatSeq, dOmegaSeq] = objIntegrator.integrate(dTimegrid, ...
                                                        dQuat0, ...
                                                        dOmega0, ...
                                                        varTorque, ...
                                                        0.01, ...
                                                        'rk4_rkmk4', ...
                                                        true, ...
                                                        1.0);

            testCase.verifyEqual(dOmegaSeq(:,1), dOmega0, 'AbsTol',1e-12);
            for idk = 2:size(dOmegaSeq, 2)
                testCase.verifyNotEqual(dOmegaSeq(:,idk), dOmega0);
            end

            % Compute angular momentum at each time
            dAngMom = objIntegrator.dInertiaMatrix * dOmegaSeq;
            dNormAngMom = vecnorm(dAngMom, 2, 1);
            testCase.verifyLessThanOrEqual( abs(dNormAngMom - dNormAngMom(1))./dNormAngMom, 1e-6)

        end

        function testConstantTorqueNonZeroInitialRate(testCase)

            % With constant torque, angular momentum must not be constant
            I = diag([2,2,2]);
            objIntegrator = CRigidBodyDynamicsIntegrator(I);
            
            dTimegrid = 0:0.2:1.0;
            dQuat0    = [1;0;0;0];
            dOmega0   = [0.1; 0.0; -0.1];
            dTorque   = [0.2; 0.0; 0.0];
            varTorque = @(dTstamp, dOmega, dQuat0) dTorque;

            [~, dOmegaSeq] = objIntegrator.integrate(dTimegrid, ...
                dQuat0, dOmega0, varTorque, 0.2, 'rk4_rkmk4', true, 1.0);

            % Analytical angular acceleration: alpha = I^{-1} * torque = [0.1;0;0]
            alpha = I \ dTorque;

            for k=1:size(dOmegaSeq,2)
                t = dTimegrid(k);
                expectedOmega = dOmega0 + alpha * t;
                testCase.verifyEqual(dOmegaSeq(:,k), expectedOmega, 'AbsTol',1e-6);
            end

        end

        function testTorqueFreeMotionAnalytical(testCase)

            % Validate quaternion and omega against analytical torque-free solution
            objIntegrator = CRigidBodyDynamicsIntegrator();
            
            objIntegrator.dInertiaMatrix = eye(3);
            objIntegrator.dInertiaMatrix(3,3) = 10;

            dTimegrid = 0:0.1:1.0;
            omega0 = [0.0; 0.0; pi]; % rotation about z-axis at pi rad/s
            dQuat0 = [1; 0; 0; 0];
            varTorque = zeros(3,1);
            [dQuatSeq, dOmegaSeq] = objIntegrator.integrate(dTimegrid, dQuat0, omega0, varTorque, 0.1, 'rk4_rkmk4', true, 1.0);
            
            % Analytical: constant omega and quaternion q(t) = [cos(theta/2); 0; 0; sin(theta/2)]
            for k = 1:length(dTimegrid)
            
                t = dTimegrid(k);
                theta = norm(omega0) * t;

                qExp = [cos(theta/2); 0; 0; sin(theta/2)];
                qAct = dQuatSeq(:,k) / norm(dQuatSeq(:,k));
                
                % Account for sign ambiguity
                err1 = norm(qAct - qExp);
                err2 = norm(qAct + qExp);
                testCase.verifyLessThan(min(err1, err2), 1e-6);
                
                % Omega remains constant
                testCase.verifyEqual(dOmegaSeq(:,k), omega0, 'AbsTol',1e-12);
            end
        end

        function testRk4Rkmk4MatchesOdeReferenceForChangingBodyRate(testCase)
            dInertia = diag([1.5, 2.0, 3.5]);
            objIntegrator = CRigidBodyDynamicsIntegrator(dInertia);

            dTimegrid = 0:1:60;
            dQuat0 = [1;0;0;0];
            dOmega0 = [0.02; 0.01; 0.2];
            varTorque = zeros(3,1);

            [dQuatSeqRKMK4, dOmegaSeqRKMK4] = objIntegrator.integrate(dTimegrid, ...
                                                                      dQuat0, ...
                                                                      dOmega0, ...
                                                                      varTorque, ...
                                                                      0.1, ...
                                                                      'rk4_rkmk4', ...
                                                                      true, ...
                                                                      1.0);

            [dQuatRefSeq, dOmegaRefSeq] = testCRigidBodyDynamicsIntegrator.referenceRigidBodyOde_(dTimegrid, ...
                                                                                                  dQuat0, ...
                                                                                                  dOmega0, ...
                                                                                                  dInertia);
            dQuatErr = testCRigidBodyDynamicsIntegrator.quatSequenceError_(dQuatSeqRKMK4, dQuatRefSeq);
            dOmegaErr = vecnorm(dOmegaSeqRKMK4 - dOmegaRefSeq, 2, 1);

            testCase.verifyLessThan(max(dQuatErr), 2e-5);
            testCase.verifyLessThan(max(dOmegaErr), 5e-8);
        end

        function testRk4Rkmk4PreservesInertialAngularMomentumDirection(testCase)
            %% SIGNATURE
            % testRk4Rkmk4PreservesInertialAngularMomentumDirection(testCase)
            % -------------------------------------------------------------------------------------------------
            %% DESCRIPTION
            % Verify the passive-attitude/body-rate contract through conservation of the inertial momentum vector.
            % -------------------------------------------------------------------------------------------------
            %% INPUT
            % testCase    Active MATLAB unit-test instance.
            % -------------------------------------------------------------------------------------------------
            %% OUTPUT
            % None.
            % -------------------------------------------------------------------------------------------------
            %% CHANGELOG
            % 23-08-2026  Pietro Califano, Codex     First implementation.
            % -------------------------------------------------------------------------------------------------
            %% DEPENDENCIES
            % CRigidBodyDynamicsIntegrator, QuatSeq2DCM, RotationVectorToDCM
            % -------------------------------------------------------------------------------------------------

            dInertia = diag([1.8, 2.3, 3.1]);
            dInitialDCM_TBfromIN = RotationVectorToDCM([0.31; -0.18; 0.22]);
            dInitialQuat_TBfromIN = DCM2quatSeq(dInitialDCM_TBfromIN, false);
            dInitialAngVel_TB = [0.11; -0.07; 0.24];
            dTimegrid = 0.0:0.5:20.0;

            objIntegrator = CRigidBodyDynamicsIntegrator(dInertia);
            [dQuatHistory_TBfromIN, dAngVelHistory_TB] = objIntegrator.integrate( ...
                dTimegrid, dInitialQuat_TBfromIN, dInitialAngVel_TB, ...
                zeros(3,1), 0.02, 'rk4_rkmk4', true, 0.02);

            dDCMHistory_TBfromIN = QuatSeq2DCM(dQuatHistory_TBfromIN, false);
            dAngMomHistory_IN = zeros(3, numel(dTimegrid));
            for ui32TimeIdx = uint32(1):uint32(numel(dTimegrid))
                dAngMomHistory_IN(:, ui32TimeIdx) = ...
                    transpose(dDCMHistory_TBfromIN(:, :, ui32TimeIdx)) * ...
                    dInertia * dAngVelHistory_TB(:, ui32TimeIdx);
            end

            dMaxAngMomDirectionDrift = max(vecnorm( ...
                dAngMomHistory_IN - dAngMomHistory_IN(:, 1), 2, 1));
            testCase.verifyLessThan(dMaxAngMomDirectionDrift, 2.0e-10);
        end

        function testStatelessConstantTorqueStepMatchesCallbackPath(testCase)
            %% SIGNATURE
            % testStatelessConstantTorqueStepMatchesCallbackPath(testCase)
            % -------------------------------------------------------------------------------------------------
            %% DESCRIPTION
            % Compare the stateless constant-body-torque step with the independent callback-based compatibility path.
            % -------------------------------------------------------------------------------------------------
            %% INPUT
            % testCase    Active MATLAB unit-test instance.
            % -------------------------------------------------------------------------------------------------
            %% OUTPUT
            % None.
            % -------------------------------------------------------------------------------------------------
            %% CHANGELOG
            % 23-08-2026  Pietro Califano, Codex     First implementation.
            % -------------------------------------------------------------------------------------------------
            %% DEPENDENCIES
            % CRigidBodyDynamicsIntegrator.IntegrStep_RK4_RKMK4, IntegrateRigidBodyRK4RKMK4Step
            % -------------------------------------------------------------------------------------------------

            dInertia = [2.1, 0.08, -0.03; ...
                        0.08, 2.7, 0.05; ...
                       -0.03, 0.05, 3.4];
            dInitialQuat = CQuatKinematicsIntegrator.NormalizeSeq( ...
                [0.88; 0.19; -0.31; 0.23]);
            dInitialAngVel = [0.14; -0.09; 0.27];
            dTorque = [2.0e-3; -1.0e-3; 1.5e-3];
            dStepSize = 0.04;
            fcnAngAccel = @(~, dAngVel, ~) dInertia \ ...
                (dTorque - cross(dAngVel, dInertia * dAngVel));

            [dExpectedAngVel, dExpectedQuat] = ...
                CRigidBodyDynamicsIntegrator.IntegrStep_RK4_RKMK4( ...
                    dInitialAngVel, dInitialQuat, fcnAngAccel, 0.0, dStepSize);
            [dActualQuat, dActualAngVel] = IntegrateRigidBodyRK4RKMK4Step( ...
                dInertia, dInitialQuat, dInitialAngVel, dTorque, dStepSize);

            testCase.verifyEqual(dActualAngVel, dExpectedAngVel, 'AbsTol', 2.0e-15);
            testCase.verifyEqual(dActualQuat, dExpectedQuat, 'AbsTol', 2.0e-15);
        end

        function testTorqueFreeSymmetricTop(testCase)

            % Torque-free motion for symmetric top (I1=I2 != I3)
            dI1 = 2; 
            dI2 = 2; 
            dI3 = 1;
            
            dI = diag([dI1, dI2, dI3]);
            objIntegrator = CRigidBodyDynamicsIntegrator(dI);

            % Time grid
            dTimegrid = 0:0.1:2.0;

            % Initial angular velocity: small transverse and dominant spin about symmetry axis
            dOmega0 = [0.1; 0.0; 1.0];
            dQuat0 = [1; 0; 0; 0];
            
            varTorque = zeros(3,1);
            [dQuatSeq, dOmegaSeq] = objIntegrator.integrate(dTimegrid, ...
                    dQuat0, dOmega0, varTorque, 0.1, 'rk4_rkmk4', true, 1.0);
            
            % Analytical precession frequency: dOmega_p = (I3 - I1)/I1 * omega3
            dOmega_p = (dI3 - dI1)/dI1 * dOmega0(3);
            dA = dOmega0(1);
            
            for k = 1:length(dTimegrid)

                t = dTimegrid(k);
                dExp1 = dA * cos(dOmega_p * t);
                dExp2 = dA * sin(dOmega_p * t);
                
                % Compare transverse components
                testCase.verifyEqual(dOmegaSeq(1,k), dExp1, 'AbsTol',1e-4);
                testCase.verifyEqual(dOmegaSeq(2,k), dExp2, 'AbsTol',1e-4);
                
                % Spin component remains constant
                testCase.verifyEqual(dOmegaSeq(3,k), dOmega0(3), 'AbsTol',1e-12);
            end
        end

    end

    methods (Static, Access = private)
        function [dQuatRefSeq, dOmegaRefSeq] = referenceRigidBodyOde_(dTimegrid, dQuat0, dOmega0, dInertia)
            stOptions = odeset('RelTol', 1e-12, 'AbsTol', 1e-14);
            [~, dStateSeq] = ode113(@(dTstamp, dState) testCRigidBodyDynamicsIntegrator.rigidBodyRhs_(dTstamp, dState, dInertia), ...
                                    dTimegrid, ...
                                    [dQuat0; dOmega0], ...
                                    stOptions);
            dQuatRefSeq = CQuatKinematicsIntegrator.NormalizeSeq(transpose(dStateSeq(:,1:4)));
            dOmegaRefSeq = transpose(dStateSeq(:,5:7));
        end

        function dStateDot = rigidBodyRhs_(~, dState, dInertia)
            dQuat = dState(1:4);
            dQuat = dQuat / norm(dQuat);
            dOmega = dState(5:7);

            dOmegaMat = [0.0,       -dOmega(1), -dOmega(2), -dOmega(3); ...
                         dOmega(1),  0.0,        dOmega(3), -dOmega(2); ...
                         dOmega(2), -dOmega(3),  0.0,        dOmega(1); ...
                         dOmega(3),  dOmega(2), -dOmega(1),  0.0];
            dQuatDot = 0.5 * dOmegaMat * dQuat;
            dOmegaDot = CRigidBodyDynamicsIntegrator.EvalRHS_AngAccel_(dInertia, dOmega, zeros(3,1));
            dStateDot = [dQuatDot; dOmegaDot];
        end

        function dQuatErr = quatSequenceError_(dQuatSeq, dQuatRefSeq)
            dQuatErr = min(vecnorm(dQuatSeq - dQuatRefSeq, 2, 1), ...
                           vecnorm(dQuatSeq + dQuatRefSeq, 2, 1));
        end
    end
end
