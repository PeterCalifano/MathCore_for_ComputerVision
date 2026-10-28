classdef testPropagateKeplerianElems < matlab.unittest.TestCase
    % Test suite for PropagateKeplerianElems function. Besides the shape, period-return and error checks, the
    % suite validates the propagator against a two-body numerical integration (ode113 at 1e-13) on elliptic
    % and hyperbolic arcs started away from periapsis, where the anomaly conversions matter (the 28-10-2026
    % regression: Kepler's equation evaluated with the true anomaly gave 100 km errors within minutes).

    methods (Test)
        function TestZeroScalarHorizonReturnsInitialState(testCase)
            dxStateKeplInit = [7000; 0.1; 0; 0; 0; 1.0];
            dGravParam = 398600.4418;
            dScalarHorizonSec = 0.0;

            dxStatesOut = PropagateKeplerianElems(dxStateKeplInit, dGravParam, dScalarHorizonSec);
            testCase.verifySize(dxStatesOut, [6, 1]);
            testCase.verifyEqual(dxStatesOut(:, 1), dxStateKeplInit, 'AbsTol', 1e-12);
        end

        function TestPositiveScalarHorizonReturnsTwoSamples(testCase)
            % Circular orbit: true anomaly equals mean anomaly.
            dSma = 7000;
            dxStateKeplInit = [dSma; 0; 0; 0; 0; 0];
            dGravParam = 398600.4418;
            dMeanMotion = sqrt(dGravParam / (dSma^3));
            dScalarHorizonSec = (pi/2) / dMeanMotion; % Quarter period

            dxStatesOut = PropagateKeplerianElems(dxStateKeplInit, dGravParam, dScalarHorizonSec);
            dExpectedTrueAnomaly = mod(dMeanMotion * dScalarHorizonSec, 2*pi);

            testCase.verifySize(dxStatesOut, [6, 2]);
            testCase.verifyEqual(dxStatesOut(6, 2), dExpectedTrueAnomaly, 'AbsTol', 1e-9);
        end

        function TestCircularOnePeriodReturnsInitialState(testCase)
            dSma = 7000;
            dEcc = 0.0;
            dIncl = 0.25;
            dRaan = 0.30;
            dArgPer = 0.40;
            dTrueAnom = 1.10;
            dxStateKeplInit = [dSma; dEcc; dIncl; dRaan; dArgPer; dTrueAnom];

            dGravParam = 398600.4418;
            dPeriodSec = 2.0 * pi * sqrt(dSma^3 / dGravParam);
            dTimeGridSec = [0.0, dPeriodSec];

            dxStatesOut = PropagateKeplerianElems(dxStateKeplInit, dGravParam, dTimeGridSec);
            dxStateFinal = dxStatesOut(:, end);

            testCase.verifyEqual(dxStateFinal(1:5), dxStateKeplInit(1:5), 'AbsTol', 1e-8);

            dWrappedAnomalyErr = testPropagateKeplerianElems.ComputeWrappedAngleError_(dxStateFinal(6), dxStateKeplInit(6));
            testCase.verifyLessThan(dWrappedAnomalyErr, 1e-8);
        end

        function TestCircularOnePeriodReturnsInitialCartesianState(testCase)
            dSma = 7000;
            dEcc = 0.0;
            dIncl = 0.25;
            dRaan = 0.30;
            dArgPer = 0.40;
            dTrueAnom = 1.10;
            dxStateKeplInit = [dSma; dEcc; dIncl; dRaan; dArgPer; dTrueAnom];

            dGravParam = 398600.4418;
            dPeriodSec = 2.0 * pi * sqrt(dSma^3 / dGravParam);
            dTimeGridSec = [0.0, dPeriodSec];

            dxStateCartInit = kepl2rv(dxStateKeplInit, dGravParam, "rad");
            dxStatesCartOut = PropagateKeplerianElems(dxStateKeplInit, ...
                                                      dGravParam, ...
                                                      dTimeGridSec, ...
                                                      "bConvert2Cart", true);
            dxStateCartFinal = dxStatesCartOut(:, end);

            testCase.verifyEqual(dxStateCartFinal(1:3), dxStateCartInit(1:3), 'AbsTol', 1e-5);
            testCase.verifyEqual(dxStateCartFinal(4:6), dxStateCartInit(4:6), 'AbsTol', 1e-8);
        end

        function TestEllipticalOnePeriodReturnsInitialState(testCase)
            dSma = 12000;
            dEcc = 0.35;
            dIncl = 0.55;
            dRaan = 1.10;
            dArgPer = 0.70;
            dTrueAnom = 0.0;
            dxStateKeplInit = [dSma; dEcc; dIncl; dRaan; dArgPer; dTrueAnom];

            dGravParam = 398600.4418;
            dPeriodSec = 2.0 * pi * sqrt(dSma^3 / dGravParam);
            dTimeGridSec = [0.0, dPeriodSec];

            dxStatesOut = PropagateKeplerianElems(dxStateKeplInit, dGravParam, dTimeGridSec);
            dxStateFinal = dxStatesOut(:, end);

            testCase.verifyEqual(dxStateFinal(1:5), dxStateKeplInit(1:5), 'AbsTol', 1e-8);

            dWrappedAnomalyErr = testPropagateKeplerianElems.ComputeWrappedAngleError_(dxStateFinal(6), dxStateKeplInit(6));
            testCase.verifyLessThan(dWrappedAnomalyErr, 1e-8);
        end

        function TestEllipticalOnePeriodReturnsInitialCartesianState(testCase)
            dSma = 12000;
            dEcc = 0.35;
            dIncl = 0.55;
            dRaan = 1.10;
            dArgPer = 0.70;
            dTrueAnom = 0.0;
            dxStateKeplInit = [dSma; dEcc; dIncl; dRaan; dArgPer; dTrueAnom];

            dGravParam = 398600.4418;
            dPeriodSec = 2.0 * pi * sqrt(dSma^3 / dGravParam);
            dTimeGridSec = [0.0, dPeriodSec];

            dxStateCartInit = kepl2rv(dxStateKeplInit, dGravParam, "rad");
            dxStatesCartOut = PropagateKeplerianElems(dxStateKeplInit, ...
                                                      dGravParam, ...
                                                      dTimeGridSec, ...
                                                      "bConvert2Cart", true);
            dxStateCartFinal = dxStatesCartOut(:, end);

            testCase.verifyEqual(dxStateCartFinal(1:3), dxStateCartInit(1:3), 'AbsTol', 1e-4);
            testCase.verifyEqual(dxStateCartFinal(4:6), dxStateCartInit(4:6), 'AbsTol', 1e-7);
        end

        function TestAbsoluteTimeGridOutputSize(testCase)
            dxStateKeplInit = [7000; 0.1; 0; 0; 0; 0.2];
            dGravParam = 398600.4418;
            dTimeGridSec = [0, 100, 200, 300];

            dxStatesOut = PropagateKeplerianElems(dxStateKeplInit, dGravParam, dTimeGridSec);
            testCase.verifySize(dxStatesOut, [6, numel(dTimeGridSec)]);
        end

        function TestNonIncreasingTimeGridThrowsAssertion(testCase)
            dxStateKeplInit = [7000; 0.2; 0; 0; 0; 0.1];
            dGravParam = 398600.4418;
            dTimeGridSec = [0, 10, 10];

            try
                PropagateKeplerianElems(dxStateKeplInit, dGravParam, dTimeGridSec);
                testCase.verifyFail('Expected an assertion for non-increasing time grid.');
            catch ME
                testCase.verifyTrue(contains(ME.message, 'Time grid must be strictly increasing'));
            end
        end

        function TestHyperbolicBranchRunsWithoutUndefinedSymbols(testCase)
            dxStateKeplInit = [-10000; 1.5; 0; 0; 0; 0.2];
            dGravParam = 398600.4418;
            dScalarHorizonSec = 50.0;

            dxStatesOut = PropagateKeplerianElems(dxStateKeplInit, dGravParam, dScalarHorizonSec);
            testCase.verifySize(dxStatesOut, [6, 2]);
            testCase.verifyTrue(all(isfinite(dxStatesOut), 'all'));
        end

        function TestTrueAnomalyIsContinuousAtStart(testCase)
            % Regression of 28-10-2026: a propagation over a negligible interval must leave the true anomaly
            % unchanged. With the true anomaly fed to Kepler's equation the first step jumped by about e
            % (0.1 rad here) for an elliptic state at nu = 90 deg, and similarly on the hyperbola.
            dGravParam = 398600.4418;
            cellStates = {[7000; 0.1; 0.3; 0.2; 0.1; pi/2], [-10000; 1.5; 0.3; 0.2; 0.1; 0.8]};
            for idCase = 1:numel(cellStates)
                dxStateKeplInit = cellStates{idCase};
                dxStatesOut = PropagateKeplerianElems(dxStateKeplInit, dGravParam, [0, 1e-6]);
                dAnomalyErr = testPropagateKeplerianElems.ComputeWrappedAngleError_(dxStatesOut(6, end), dxStateKeplInit(6));
                testCase.verifyLessThan(dAnomalyErr, 1e-6, sprintf('case %d', idCase));
            end
        end

        function TestEllipticalArcMatchesTwoBodyIntegration(testCase)
            % Three periods of an e = 0.35 orbit started at nu = 0.9 rad, sampled every 60 s, against ode113.
            dSma = 12000; dEcc = 0.35;
            dxStateKeplInit = [dSma; dEcc; 0.55; 1.10; 0.70; 0.9];
            dGravParam = 398600.4418;
            dPeriodSec = 2.0 * pi * sqrt(dSma^3 / dGravParam);
            dTimeGridSec = 0:60:(3 * dPeriodSec);
            testPropagateKeplerianElems.VerifyAgainstTwoBodyIntegration_(testCase, dxStateKeplInit, dGravParam, dTimeGridSec);
        end

        function TestHyperbolicArcMatchesTwoBodyIntegration(testCase)
            % Inbound hyperbola (e = 1.5, nu0 = -0.8 rad) through periapsis, sampled every 10 s, against ode113.
            dxStateKeplInit = [-10000; 1.5; 0.3; 0.2; 0.1; -0.8];
            dGravParam = 398600.4418;
            dTimeGridSec = 0:10:4000;
            testPropagateKeplerianElems.VerifyAgainstTwoBodyIntegration_(testCase, dxStateKeplInit, dGravParam, dTimeGridSec);
        end

        function TestSmallBodyApproachArcMatchesTwoBodyIntegration(testCase)
            % Metre-scale approach of a 500 m body (GM 69.89 m^3/s^2, e about 44): the arc that exposed the
            % anomaly defect in the image-motion recording generator (radius 3015 m -> 1487 m after one second).
            dGravParam = 69.89;
            dxCartInit = [0; -300; 3000; 0.5; 0.2; -2.0];
            dxStateKeplInit = rv2kepl(dxCartInit, dGravParam, 'rad');
            dTimeGridSec = 0:1:120;
            testPropagateKeplerianElems.VerifyAgainstTwoBodyIntegration_(testCase, dxStateKeplInit, dGravParam, dTimeGridSec);
        end

        function TestSequentialStepsMatchSingleStep(testCase)
            % The state after one long step equals the state after many short steps (the true anomaly carried
            % between steps is converted consistently), elliptic and hyperbolic.
            dGravParam = 398600.4418;
            cellStates = {[9000; 0.4; 0.3; 0.2; 0.1; 2.0], [-10000; 1.5; 0.3; 0.2; 0.1; -0.8]};
            for idCase = 1:numel(cellStates)
                dxStateKeplInit = cellStates{idCase};
                dxSingle = PropagateKeplerianElems(dxStateKeplInit, dGravParam, [0, 3000], 'bConvert2Cart', true, 'dRadAngleTol', 1e-13);
                dxMany = PropagateKeplerianElems(dxStateKeplInit, dGravParam, 0:30:3000, 'bConvert2Cart', true, 'dRadAngleTol', 1e-13);
                % 100 sequential Newton solves at 1e-13 rad on a 1e4 km orbit: below 1e-6 km of position.
                testCase.verifyEqual(dxMany(1:3, end), dxSingle(1:3, end), 'AbsTol', 1e-6, sprintf('case %d', idCase));
                testCase.verifyEqual(dxMany(4:6, end), dxSingle(4:6, end), 'AbsTol', 1e-9, sprintf('case %d', idCase));
            end
        end

        function TestInvalidEccentricityError(testCase)
            dxStateKeplInit = [7000; -0.1; 0; 0; 0; 0];
            dGravParam = 398600.4418;
            testCase.verifyError(@() PropagateKeplerianElems(dxStateKeplInit, dGravParam, 0), ...
                                 'MATLAB:assertion:failed');
        end
    end

    methods (Static, Access = private)
        function dWrappedErr = ComputeWrappedAngleError_(dAngleA, dAngleB)
            dWrappedErr = abs(atan2(sin(dAngleA - dAngleB), cos(dAngleA - dAngleB)));
        end

        function VerifyAgainstTwoBodyIntegration_(testCase, dxStateKeplInit, dGravParam, dTimeGridSec)
            % Compare the Cartesian states on the grid with ode113 (relative and absolute tolerance 1e-13) of
            % r'' = -GM r / |r|^3 from the same initial state, and assert ten times the error budget:
            % propagator N_steps (|r|max + |a|) tol_anomaly (the state changes by at most about |r| + |a| per
            % radian of anomaly) plus integrator 1e4 RelTol |r|max (step count bound of the adaptive solver);
            % velocity likewise with |v|max + sqrt(GM/|a|). The budget is checked to be discriminative
            % (below 1e-4 of the smallest radius), so the anomaly defect of 28-10-2026 cannot pass.
            dAnomalyTol = 1e-12;
            dIntegratorRelTol = 1e-13;
            dxStates = PropagateKeplerianElems(dxStateKeplInit, dGravParam, dTimeGridSec, ...
                'bConvert2Cart', true, 'dRadAngleTol', dAnomalyTol);
            dxCartInit = kepl2rv(dxStateKeplInit, dGravParam, "rad");
            objOptions = odeset('RelTol', dIntegratorRelTol, 'AbsTol', dIntegratorRelTol);
            [~, dxReference] = ode113(@(dTime, dxState) [dxState(4:6); -dGravParam * dxState(1:3) / norm(dxState(1:3))^3], ...
                dTimeGridSec, dxCartInit, objOptions);
            dxReference = dxReference';
            dRadiusMax = max(vecnorm(dxReference(1:3, :))); dRadiusMin = min(vecnorm(dxReference(1:3, :)));
            dSpeedMax = max(vecnorm(dxReference(4:6, :)));
            dSmaAbs = abs(dxStateKeplInit(1));
            dNumSteps = numel(dTimeGridSec);
            dPositionBudget = dNumSteps * (dRadiusMax + dSmaAbs) * dAnomalyTol + 1e4 * dIntegratorRelTol * dRadiusMax;
            dVelocityBudget = dNumSteps * (dSpeedMax + sqrt(dGravParam / dSmaAbs)) * dAnomalyTol + 1e4 * dIntegratorRelTol * dSpeedMax;
            testCase.verifyLessThan(10 * dPositionBudget, 1e-4 * dRadiusMin, 'budget must be discriminative');
            dPositionErr = max(vecnorm(dxStates(1:3, :) - dxReference(1:3, :)));
            dVelocityErr = max(vecnorm(dxStates(4:6, :) - dxReference(4:6, :)));
            testCase.verifyLessThan(dPositionErr, 10 * dPositionBudget, sprintf('position error %.3e', dPositionErr));
            testCase.verifyLessThan(dVelocityErr, 10 * dVelocityBudget, sprintf('velocity error %.3e', dVelocityErr));
        end
    end
end
