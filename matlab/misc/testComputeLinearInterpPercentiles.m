function tests = testComputeLinearInterpPercentiles
%% DESCRIPTION
% Unit tests of ComputeLinearInterpPercentiles. Reference values come from numpy.percentile
% (NumPy 1.26, default 'linear' method) or are derived by hand from the definition.
% -------------------------------------------------------------------------------------------------------------
%% CHANGELOG
% 15-11-2026  Pietro Califano, Claude Opus 5.5  Add tests of the linear-interpolation percentiles.
% -------------------------------------------------------------------------------------------------------------
tests = functiontests(localfunctions);
end

function testMatchesNumpyReference(testCase)
% numpy.percentile([3.2, -1.5, 7.7, 0.4, 2.2, 9.1, -4.0, 5.5], [5, 25, 50, 90, 95]).
dValues = [3.2, -1.5, 7.7, 0.4, 2.2, 9.1, -4.0, 5.5];
dExpected = [-3.125, -0.075, 2.7, 8.12, 8.61];
verifyEqual(testCase, ComputeLinearInterpPercentiles(dValues, [5, 25, 50, 90, 95]), dExpected, ...
    'AbsTol', 1e-12);
end

function testEndpointsAndSingleSample(testCase)
verifyEqual(testCase, ComputeLinearInterpPercentiles(1:5, [0, 50, 95, 100]), [1, 3, 4.8, 5], ...
    'AbsTol', 1e-12);
verifyEqual(testCase, ComputeLinearInterpPercentiles(7, [5, 95]), [7; 7]);
end

function testIntegerCountsDifferFromMidpointRule(testCase)
% 100 counts with x(95) = 4 and x(96) = 5: the linear rule gives 4.05, the midpoint rule 4.5.
dCounts = [zeros(1, 94), 4, 5 * ones(1, 5)];
verifyEqual(testCase, ComputeLinearInterpPercentiles(dCounts, 95), 4.05, 'AbsTol', 1e-12);
end

function testColumnPercentileVector(testCase)
verifyEqual(testCase, ComputeLinearInterpPercentiles(1:5, [5; 95]), [1.2, 4.8], 'AbsTol', 1e-12);
end

function testOmitsNanPerSlice(testCase)
dValues = [1, NaN; 2, NaN; NaN, NaN; 4, NaN];
dResult = ComputeLinearInterpPercentiles(dValues, [0, 50, 100], 1);
verifyEqual(testCase, dResult(:, 1), [1; 2; 4], 'AbsTol', 1e-12);
verifyTrue(testCase, all(isnan(dResult(:, 2))));
end

function testDimensionAndOutputSize(testCase)
dValues = reshape(1:24, 4, 3, 2);
verifySize(testCase, ComputeLinearInterpPercentiles(dValues, [5, 50, 95]), [3, 3, 2]);
verifySize(testCase, ComputeLinearInterpPercentiles(dValues, [5, 50, 95], 2), [4, 3, 2]);
verifySize(testCase, ComputeLinearInterpPercentiles(dValues, [50, 95], 3), [4, 3, 2]);
verifyEqual(testCase, ComputeLinearInterpPercentiles(dValues, 50, 2), median(dValues, 2), ...
    'AbsTol', 1e-12);
verifyEqual(testCase, ComputeLinearInterpPercentiles((1:5)', [0, 100]), [1; 5]);
end

function testEmptyVectorReturnsNan(testCase)
verifyEqual(testCase, ComputeLinearInterpPercentiles(zeros(0, 1), [5, 95]), [NaN; NaN]);
verifyEqual(testCase, ComputeLinearInterpPercentiles(zeros(1, 0), 50), NaN);
end

function testInfiniteSamples(testCase)
verifyEqual(testCase, ComputeLinearInterpPercentiles([1, Inf], [0, 50, 100]), [1, Inf, Inf]);
end

function testRejectsOutOfRangePercentiles(testCase)
verifyError(testCase, @() ComputeLinearInterpPercentiles(1:5, 101), ...
    'MATLAB:validators:mustBeInRange');
end
