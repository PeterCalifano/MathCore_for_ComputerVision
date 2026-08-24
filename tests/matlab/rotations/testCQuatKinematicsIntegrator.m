classdef testCQuatKinematicsIntegrator < matlab.unittest.TestCase
    % Unit tests for CQuatKinematicsIntegrator

    properties
        Integr  CQuatKinematicsIntegrator
        Tolerance = 1e-10
    end

    methods (TestMethodSetup)
        function createIntegrator(testCase)
            testCase.Integr = CQuatKinematicsIntegrator();
        end
    end

    methods (Test)

        function testNormalizeSeq(testCase)
            % Norm of each column should be 1
            Q = [2 0; 0 3; 0 4; 0 0];
            Qn = CQuatKinematicsIntegrator.NormalizeSeq(Q);
            norms = vecnorm(Qn);
            testCase.verifyLessThanOrEqual(abs(norms - 1), testCase.Tolerance);
        end

        function testExpMapIdentity(testCase)
            % Zero rotation should map to identity quaternion
            dq = [0;0;0];
            q = CQuatKinematicsIntegrator.ExpMap(dq);
            testCase.verifyEqual(q, [1;0;0;0], 'AbsTol', testCase.Tolerance);
        end

        function testExpMapPiRotation(testCase)
            % Rotation of pi about x-axis
            dq = [pi;0;0];
            q = CQuatKinematicsIntegrator.ExpMap(dq);
            expected = [cos(pi/2); sin(pi/2); 0; 0];
            testCase.verifyEqual(q, expected, 'AbsTol', testCase.Tolerance);
        end

        function testExpMapSequence(testCase)
            dq = [0,  pi, 0;
                  0,  0,  pi;
                  0,  0,  0];

            q = CQuatKinematicsIntegrator.ExpMap(dq);
            expected = [1,  cos(pi/2),  cos(pi/2);
                        0,  sin(pi/2),  0;
                        0,  0,          sin(pi/2);
                        0,  0,          0];

            testCase.verifyEqual(q, expected, 'AbsTol', testCase.Tolerance);
        end

        function testZeroOmegaRK4(testCase)
            % Integrating zero angular velocity yields constant quaternion
            q0   = [1;0;0;0];
            omega = [0;0;0];
            tgrid = linspace(0,1,11);

            % Integrate
            [qEnd, dTimegridOut, qSeq] = testCase.Integr.integrate(q0, omega, tgrid, 'rk4');
            % All outputs should equal the initial quaternion
            testCase.verifyEqual(qSeq, repmat(q0,1,numel(tgrid)), 'AbsTol', testCase.Tolerance);
            testCase.verifyEqual(qEnd, q0, 'AbsTol', testCase.Tolerance);
        end

        function testConstantOmegaRK4(testCase)
            % Constant rotation about Z at pi rad/s for 1 second
            q0   = [1;0;0;0];
            omega = [0;0;pi];
            dt    = 0.1;
            tgrid = 0:dt:1;

            % Integrate
            [qEnd, dTimegridOut, qSeq] = testCase.Integr.integrate(q0, omega, tgrid, 'rk4', dt);
            % Final rotation is by pi radians about Z:
            expected = [cos(pi/2); 0;0;sin(pi/2)];
            testCase.verifyEqual(qEnd, expected, 'AbsTol', 1e-3);
            testCase.verifyEqual(tgrid, dTimegridOut, 'AbsTol', 1e-6)
        end

        function testConstantOmegaRK4_CustomTimegrid(testCase)
            % Constant rotation about Z at pi rad/s for 1 second
            q0   = [1;0;0;0];
            omega = [0;0;pi];
            dt    = 0.06;
            tgrid = 0:0.1:1;

            % Integrate
            [qEnd, dTimegridOut, qSeq] = testCase.Integr.integrate(q0, omega, tgrid, 'rk4', dt);
            % Final rotation is by pi radians about Z:
            expected = [cos(pi/2); 0;0;sin(pi/2)];
            testCase.verifyEqual(qEnd, expected, 'AbsTol', 1e-3);
            testCase.verifyEqual(tgrid, dTimegridOut, 'AbsTol', 1e-6)
        end


        function testConstantOmegaLieGroupEuler(testCase)
            % Constant rotation about Z at pi rad/s for 1 second
            q0   = [1;0;0;0];
            omega = [0;0;pi];
            dt    = 0.1;
            tgrid = 0:dt:1;

            [qEnd, dTimegridOut, qSeq] = testCase.Integr.integrate(q0, omega, tgrid, 'lie_euler', dt);
            % Final rotation is by pi radians about Z:
            expected = [cos(pi/2); 0;0;sin(pi/2)];
            testCase.verifyEqual(qEnd, expected, 'AbsTol', 1e-3);
            testCase.verifyEqual(tgrid, dTimegridOut, 'AbsTol', 1e-6)
        end

        function testConstantOmegaRKMK4(testCase)
            % Constant rotation about Z at pi rad/s for 1 second
            q0   = [1;0;0;0];
            omega = [0;0;pi];
            dt    = 0.1;
            tgrid = 0:dt:1;

            % Integrate
            [qEnd, dTimegridOut, qSeq] = testCase.Integr.integrate(q0, omega, tgrid, 'rkmk4', dt);

            % Final rotation is by pi radians about Z:
            expected = [cos(pi/2); 0;0;sin(pi/2)];
            testCase.verifyEqual(qEnd, expected, 'AbsTol', 1e-3);
            testCase.verifyEqual(tgrid, dTimegridOut, 'AbsTol', 1e-6)
        end

        function testLieEulerVsRKMK4(testCase)
            % For small dt, Lie‐Euler and RKMK4 should agree to first order
            q0 = [1;0;0;0];
            omega = [0.1; 0.2; 0.3];
            dt = 1e-3;
            tgrid = [0 dt];

            % Integrate and compare
            [qLE, dTimegridOutEuler, qSeqEuler] = testCase.Integr.integrate(q0, omega, tgrid, 'lie_euler', dt);
            [qRK, dTimegridOutRK, qSeqRK] = testCase.Integr.integrate(q0, omega, tgrid, 'rkmk4',     dt);
            testCase.verifyLessThan(norm(qLE - qRK), 1e-3);
            testCase.verifyEqual(dTimegridOutEuler, dTimegridOutRK, 'AbsTol', 1e-6)
            testCase.verifyEqual(qSeqEuler, qSeqRK, 'AbsTol', 1e-4)
        end

        function testOmegaAngVelProfileRKMK4(testCase)
            q0 = [1;0;0;0];
            omegaProfile = [0, 0.5*pi, pi;
                            0, 0,      0;
                            0, 0,      0];
            tgrid = [10, 10.5, 11.0];

            [qEnd, dTimegridOut, qSeq] = testCase.Integr.integrate(q0, ...
                                                                   omegaProfile, ...
                                                                   tgrid, ...
                                                                   'rkmk4', ...
                                                                   0.25, ...
                                                                   1.0, ...
                                                                   tgrid, ...
                                                                   'linear');

            expected = [cos(pi/4); sin(pi/4); 0; 0];
            testCase.verifyEqual(qEnd, expected, 'AbsTol', 5e-4);
            testCase.verifyEqual(dTimegridOut, tgrid, 'AbsTol', 1e-10);
            testCase.verifyEqual(qSeq(:,1), q0, 'AbsTol', testCase.Tolerance);
        end

        function testOmegaFcnHandleRKMK4UsesInternalSubstepTime(testCase)
            q0 = [1;0;0;0];
            omega = @(dT) [0; 0; pi * (dT - 10)];
            tgrid = [10, 10.5, 11.0];

            [qEnd, dTimegridOut, qSeq] = testCase.Integr.integrate(q0, omega, tgrid, 'rkmk4', 0.25);

            expected = [cos(pi/4); 0; 0; sin(pi/4)];
            testCase.verifyEqual(qEnd, expected, 'AbsTol', 5e-4);
            testCase.verifyEqual(dTimegridOut, tgrid, 'AbsTol', 1e-10);
            testCase.verifyEqual(qSeq(:,1), q0, 'AbsTol', testCase.Tolerance);
        end

        function testRKMK4FunctionOmegaNonCommutingMatchesOdeReference(testCase)
            q0 = [1;0;0;0];
            omega = @(dT) [0.35 + 0.08 * sin(0.9 * dT); ...
                          -0.22 + 0.06 * cos(1.3 * dT); ...
                           0.17 + 0.05 * sin(1.7 * dT)];
            tgrid = 0:0.2:4.0;

            [~, dTimegridOut, qSeq] = testCase.Integr.integrate(q0, omega, tgrid, 'rkmk4', 0.02);
            qRefSeq = testCQuatKinematicsIntegrator.referenceQuatOde_(q0, omega, tgrid);

            dQuatErr = testCQuatKinematicsIntegrator.quatSequenceError_(qSeq, qRefSeq);
            testCase.verifyLessThan(max(dQuatErr), 5e-6);
            testCase.verifyEqual(dTimegridOut, tgrid, 'AbsTol', 1e-12);
        end

        function testRKMK4TabulatedOmegaNonCommutingMatchesOdeReference(testCase)
            q0 = [1;0;0;0];
            tgrid = 0:0.4:4.0;
            dProfileTimes = 0:0.05:4.0;
            dOmegaProfile = [0.35 + 0.08 * sin(0.9 * dProfileTimes); ...
                            -0.22 + 0.06 * cos(1.3 * dProfileTimes); ...
                             0.17 + 0.05 * sin(1.7 * dProfileTimes)];
            cellInterpMethods = {'linear', 'spline'};

            for ui32Idx = 1:numel(cellInterpMethods)
                enumInterpMethod = cellInterpMethods{ui32Idx};
                objOmegaInterp = griddedInterpolant(dProfileTimes(:), ...
                                                    transpose(dOmegaProfile), ...
                                                    enumInterpMethod, ...
                                                    'none');
                omegaReference = @(dT) transpose(objOmegaInterp(dT));

                [~, dTimegridOut, qSeq] = testCase.Integr.integrate(q0, ...
                                                                    dOmegaProfile, ...
                                                                    tgrid, ...
                                                                    'rkmk4', ...
                                                                    0.02, ...
                                                                    1.0, ...
                                                                    dProfileTimes, ...
                                                                    enumInterpMethod);
                qRefSeq = testCQuatKinematicsIntegrator.referenceQuatOde_(q0, omegaReference, tgrid);

                dQuatErr = testCQuatKinematicsIntegrator.quatSequenceError_(qSeq, qRefSeq);
                testCase.verifyLessThan(max(dQuatErr), 5e-6);
                testCase.verifyEqual(dTimegridOut, tgrid, 'AbsTol', 1e-12);
            end
        end

        function testOmegaFcnHandleRKMK4(testCase)
            
            % Constant rotation about Z at pi rad/s for 1 second
            q0   = [1;0;0;0];
            omega = @(dT) [0;0;pi];
            dt    = 0.1;
            tgrid = 0:dt:1;

            % Integrate
            [qEnd, dTimegridOut, qSeq] = testCase.Integr.integrate(q0, omega, tgrid, 'rkmk4', dt);

            % Final rotation is by pi radians about Z:
            expected = [cos(pi/2); 0;0;sin(pi/2)];
            testCase.verifyEqual(qEnd, expected, 'AbsTol', 1e-3);
            testCase.verifyEqual(tgrid, dTimegridOut, 'AbsTol', 1e-6)

        end

        function testRkmk4NoncommutingRateFourthOrder(testCase)
            %% SIGNATURE
            % testRkmk4NoncommutingRateFourthOrder(testCase)
            % -------------------------------------------------------------------------------------------------
            %% DESCRIPTION
            % Verify fourth-order convergence for a time-varying, noncommuting body angular-velocity profile.
            % -------------------------------------------------------------------------------------------------
            %% INPUT
            % testCase    Active MATLAB unit-test instance.
            % -------------------------------------------------------------------------------------------------
            %% OUTPUT
            % None.
            % -------------------------------------------------------------------------------------------------
            %% CHANGELOG
            % 24-08-2026  Pietro Califano, Codex     First implementation.
            % -------------------------------------------------------------------------------------------------
            %% DEPENDENCIES
            % IntegrateQuaternionRKMK4Step, testCQuatKinematicsIntegrator.referenceQuatOde_
            % -------------------------------------------------------------------------------------------------

            dInitialQuat = CQuatKinematicsIntegrator.NormalizeSeq( ...
                [0.91; -0.14; 0.27; 0.19]);
            fcnAngVel = @(dTstamp) [0.35 + 0.08 * sin(0.9 * dTstamp); ...
                                   -0.22 + 0.06 * cos(1.3 * dTstamp); ...
                                    0.17 + 0.05 * sin(1.7 * dTstamp)];
            dFinalTime = 4.0;
            dStepSizes = [0.4, 0.2, 0.1, 0.05, 0.025];
            dReferenceQuatSeq = testCQuatKinematicsIntegrator.referenceQuatOde_( ...
                dInitialQuat, fcnAngVel, [0.0, dFinalTime]);
            dReferenceQuat = dReferenceQuatSeq(:, end);
            dEndpointErrors = zeros(size(dStepSizes));

            % Propagate each resolution with the same RKMK4 stage sampling
            % so the measured slope isolates the canonical numerical step.
            for ui32ResolutionIdx = uint32(1):uint32(numel(dStepSizes))
                dStepSize = dStepSizes(ui32ResolutionIdx);
                ui32NumSteps = uint32(round(dFinalTime / dStepSize));
                dActualQuat = dInitialQuat;

                for ui32StepIdx = uint32(1):ui32NumSteps
                    dTstamp = (double(ui32StepIdx) - 1.0) * dStepSize;
                    dAngVelStages = [fcnAngVel(dTstamp), ...
                        fcnAngVel(dTstamp + 0.5 * dStepSize), ...
                        fcnAngVel(dTstamp + 0.5 * dStepSize), ...
                        fcnAngVel(dTstamp + dStepSize)];
                    dActualQuat = IntegrateQuaternionRKMK4Step( ...
                        dActualQuat, dAngVelStages, dStepSize);
                end

                dEndpointErrors(ui32ResolutionIdx) = min( ...
                    norm(dActualQuat - dReferenceQuat), ...
                    norm(dActualQuat + dReferenceQuat));
            end

            dObservedOrders = log2( ...
                dEndpointErrors(1:end-1) ./ dEndpointErrors(2:end));
            testCase.verifyGreaterThan(min(dObservedOrders), 3.9);
        end

        function testStatelessRkmk4StepPreservesTinyRotation(testCase)
            %% SIGNATURE
            % testStatelessRkmk4StepPreservesTinyRotation(testCase)
            % -------------------------------------------------------------------------------------------------
            %% DESCRIPTION
            % Verify that the exponential map preserves a nonzero rotation below the former small-angle cutoff.
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
            % IntegrateQuaternionRKMK4Step
            % -------------------------------------------------------------------------------------------------

            dInitialQuat = [1.0; 0.0; 0.0; 0.0];
            dAngVel = [1.0e-13; -0.6e-13; 0.4e-13];
            dAngVelStages = repmat(dAngVel, 1, 4);
            dStepSize = 1.0;

            dActualQuat = IntegrateQuaternionRKMK4Step( ...
                dInitialQuat, dAngVelStages, dStepSize);
            dLieIncrement = 0.5 * dStepSize * dAngVel;
            dIncrementNorm = norm(dLieIncrement);
            dExpectedQuat = [cos(dIncrementNorm); ...
                sin(dIncrementNorm) * dLieIncrement / dIncrementNorm];

            testCase.verifyEqual(dActualQuat, dExpectedQuat, 'AbsTol', 1.0e-28);
            testCase.verifyGreaterThan(norm(dActualQuat(2:4)), 0.0);
        end

    end

    methods (Static, Access = private)
        function qRefSeq = referenceQuatOde_(q0, omegaFcn, tgrid)
            stOptions = odeset('RelTol', 1e-12, 'AbsTol', 1e-14);
            [~, dStateSeq] = ode113(@(dTstamp, q) testCQuatKinematicsIntegrator.quatRhs_(dTstamp, q, omegaFcn), ...
                                    tgrid, ...
                                    q0, ...
                                    stOptions);
            qRefSeq = CQuatKinematicsIntegrator.NormalizeSeq(transpose(dStateSeq));
        end

        function qDot = quatRhs_(dTstamp, q, omegaFcn)
            q = q / norm(q);
            dOmega = reshape(omegaFcn(dTstamp), 3, 1);
            dOmegaMat = [0.0,       -dOmega(1), -dOmega(2), -dOmega(3); ...
                         dOmega(1),  0.0,        dOmega(3), -dOmega(2); ...
                         dOmega(2), -dOmega(3),  0.0,        dOmega(1); ...
                         dOmega(3),  dOmega(2), -dOmega(1),  0.0];
            qDot = 0.5 * dOmegaMat * q;
        end

        function dQuatErr = quatSequenceError_(qSeq, qRefSeq)
            dQuatErr = min(vecnorm(qSeq - qRefSeq, 2, 1), ...
                           vecnorm(qSeq + qRefSeq, 2, 1));
        end
    end
end
