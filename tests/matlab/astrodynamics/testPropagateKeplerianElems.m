classdef testPropagateKeplerianElems < matlab.unittest.TestCase
    % Test suite for PropagateKeplerianElems function

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
    end
end
