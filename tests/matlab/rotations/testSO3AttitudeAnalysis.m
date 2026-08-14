function tests = testSO3AttitudeAnalysis
%% SIGNATURE
% tests = testSO3AttitudeAnalysis
% -------------------------------------------------------------------------------------------------------------
%% DESCRIPTION
% Verify truth-relative SO(3) error, ZYX RPY, right-local covariance transport, robust large-angle logarithms, and
% exact swing-twist diagnostics used by navigation attitude analysis.
% -------------------------------------------------------------------------------------------------------------
%% INPUT
% None.
% -------------------------------------------------------------------------------------------------------------
%% OUTPUT
% tests    MATLAB unit-test suite.
% -------------------------------------------------------------------------------------------------------------
%% CHANGELOG
% 13-08-2026  Pietro Califano, Codex gpt-5.6     First implementation.
% -------------------------------------------------------------------------------------------------------------
%% DEPENDENCIES
% ComputeSO3AttitudeError, ComputeRPY321FromDCM, MapRightLocalRotCovToErrorLog,
% DecomposeSO3ErrorIntoAxisAngleError, RotationVectorToDCM, LogMap_SO3toR3.
% -------------------------------------------------------------------------------------------------------------

tests = functiontests(localfunctions);
end

function setupOnce(testCase)
charTestFolder = fileparts(mfilename('fullpath'));
charRootFolder = fullfile(charTestFolder, '..', '..', '..');
charMatlabFolder = fullfile(charRootFolder, 'matlab');
testCase.TestData.charOriginalPath = path;
addpath(charTestFolder, charMatlabFolder, ...
    fullfile(charMatlabFolder, 'rotations'), ...
    fullfile(charMatlabFolder, 'rotations', 'dcmLib'), ...
    fullfile(charMatlabFolder, 'linearAlgebra'));
end

function teardownOnce(testCase)
path(testCase.TestData.charOriginalPath);
end

function testTruthRelativeErrorUsesEstimatedTimesTrueTranspose(testCase)
dTrueDCM = RotationVectorToDCM([0.22; -0.11; 0.08]);
dExpectedErrorVector = [0.12; -0.08; 0.04];
dEstimatedDCM = RotationVectorToDCM(dExpectedErrorVector) * dTrueDCM;

[dErrorVector, dTheta, dErrorDCM] = ComputeSO3AttitudeError(dEstimatedDCM, dTrueDCM);

verifyEqual(testCase, dErrorVector, dExpectedErrorVector, 'AbsTol', 2.0e-14);
verifyEqual(testCase, dTheta, norm(dExpectedErrorVector), 'AbsTol', 2.0e-14);
verifyEqual(testCase, dErrorDCM, dEstimatedDCM * transpose(dTrueDCM), 'AbsTol', 0.0);
end

function testZYXReturnsRollPitchYawAndFlagsGimbalLock(testCase)
dExpectedRPY = [0.31; -0.27; 0.42];
dDCM = RotationVectorToDCM([0; 0; dExpectedRPY(3)]) * ...
    RotationVectorToDCM([0; dExpectedRPY(2); 0]) * ...
    RotationVectorToDCM([dExpectedRPY(1); 0; 0]);

[dRPY, bNearGimbalLock] = ComputeRPY321FromDCM(dDCM);
[~, bAtGimbalLock] = ComputeRPY321FromDCM( ...
    RotationVectorToDCM([0; pi / 2; 0]));

verifyEqual(testCase, dRPY, dExpectedRPY, 'AbsTol', 2.0e-14);
verifyFalse(testCase, bNearGimbalLock);
verifyTrue(testCase, bAtGimbalLock);
end

function testRightLocalCovarianceMapMatchesFiniteDifference(testCase)
dTrueDCM = RotationVectorToDCM([0.41; -0.18; 0.23]);
dErrorDCM = RotationVectorToDCM([0.32; 0.17; -0.21]);
dEstimatedDCM = dErrorDCM * dTrueDCM;
dCovariance_RightLocal = [4.0e-4, 1.2e-4, -0.4e-4; ...
                          1.2e-4, 8.0e-4,  0.7e-4; ...
                         -0.4e-4, 0.7e-4,  2.0e-4];

[dMappedCovariance, dErrorJacobian] = MapRightLocalRotCovToErrorLog( ...
    dEstimatedDCM, dTrueDCM, dCovariance_RightLocal);
dFiniteDifferenceJacobian = ComputeErrorFiniteDifference_(dEstimatedDCM, dTrueDCM, 1.0e-7);

verifyEqual(testCase, dErrorJacobian, dFiniteDifferenceJacobian, 'AbsTol', 5.0e-8);
verifyEqual(testCase, dMappedCovariance, ...
    dFiniteDifferenceJacobian * dCovariance_RightLocal * transpose(dFiniteDifferenceJacobian), ...
    'AbsTol', 2.0e-10);
end

function testZeroErrorCovarianceMapRotatesNativeAxesThroughTruth(testCase)
dTrueDCM = RotationVectorToDCM([0.31; -0.16; 0.09]);
dCovariance_RightLocal = diag([1.0e-4, 4.0e-4, 9.0e-4]);

[dMappedCovariance, dErrorJacobian] = MapRightLocalRotCovToErrorLog( ...
    dTrueDCM, dTrueDCM, dCovariance_RightLocal);

verifyEqual(testCase, dErrorJacobian, dTrueDCM, 'AbsTol', 2.0e-14);
verifyEqual(testCase, dMappedCovariance, ...
    dTrueDCM * dCovariance_RightLocal * transpose(dTrueDCM), 'AbsTol', 2.0e-14);
end

function testSwingTwistSeparatesAxisDirectionAndSpin(testCase)
dTrueSpinAxis = [0; 0; 1];
dSwingAngle = 0.24;
dTwistAngle = -0.37;
dErrorDCM = RotationVectorToDCM([dSwingAngle; 0; 0]) * ...
    RotationVectorToDCM(dTrueSpinAxis .* dTwistAngle);

[dAxisDirectionError, dSignedTwist, bTwistSingular] = ...
    DecomposeSO3ErrorIntoAxisAngleError(dErrorDCM, dTrueSpinAxis);

verifyEqual(testCase, dAxisDirectionError, dSwingAngle, 'AbsTol', 2.0e-14);
verifyEqual(testCase, dSignedTwist, dTwistAngle, 'AbsTol', 2.0e-14);
verifyFalse(testCase, bTwistSingular);
end

function testNearPiLogarithmAndSingularTwistAreExplicit(testCase)
dRotationAxis = [1; -2; 3];
dRotationAxis = dRotationAxis ./ norm(dRotationAxis);
dExpectedAngle = pi - 1.0e-8;
dNearPiDCM = RotationVectorToDCM(dRotationAxis .* dExpectedAngle);

dRecoveredVector = LogMap_SO3toR3(dNearPiDCM);
[~, dSignedTwist, bTwistSingular] = DecomposeSO3ErrorIntoAxisAngleError( ...
    RotationVectorToDCM([pi; 0; 0]), [0; 0; 1]);

verifyEqual(testCase, RotationVectorToDCM(dRecoveredVector), dNearPiDCM, 'AbsTol', 2.0e-8);
verifyTrue(testCase, isnan(dSignedTwist));
verifyTrue(testCase, bTwistSingular);
end

function testMappedCovarianceMatchesNonlinearMonteCarlo(testCase)
dTrueDCM = RotationVectorToDCM([0.38; -0.24; 0.19]);
dEstimatedDCM = RotationVectorToDCM([0.28; 0.14; -0.17]) * dTrueDCM;
dCovariance_RightLocal = [3.0e-5, 0.8e-5, -0.3e-5; ...
                          0.8e-5, 6.0e-5,  0.4e-5; ...
                         -0.3e-5, 0.4e-5,  2.0e-5];
dMappedCovariance = MapRightLocalRotCovToErrorLog( ...
    dEstimatedDCM, dTrueDCM, dCovariance_RightLocal);

rng(6143, 'twister');
ui32NumSamples = uint32(20000);
dPerturbationSamples = chol(dCovariance_RightLocal, 'lower') * ...
    randn(3, double(ui32NumSamples));
dErrorSamples = zeros(3, double(ui32NumSamples));
for ui32SampleIndex = uint32(1):ui32NumSamples
    dErrorSamples(:, ui32SampleIndex) = LogMap_SO3toR3( ...
        dEstimatedDCM * RotationVectorToDCM(dPerturbationSamples(:, ui32SampleIndex)) * ...
        transpose(dTrueDCM));
end
dEmpiricalCovariance = cov(transpose(dErrorSamples), 1);

verifyLessThan(testCase, norm(dEmpiricalCovariance - dMappedCovariance, 'fro') / ...
    norm(dMappedCovariance, 'fro'), 0.035);
end

function dFiniteDifferenceJacobian = ComputeErrorFiniteDifference_(dEstimatedDCM, dTrueDCM, dStep)
dFiniteDifferenceJacobian = zeros(3, 3);
for ui32AxisIndex = uint32(1):uint32(3)
    dPerturbation = zeros(3, 1);
    dPerturbation(ui32AxisIndex) = dStep;
    dErrorPlus = LogMap_SO3toR3( ...
        dEstimatedDCM * RotationVectorToDCM(dPerturbation) * transpose(dTrueDCM));
    dErrorMinus = LogMap_SO3toR3( ...
        dEstimatedDCM * RotationVectorToDCM(-dPerturbation) * transpose(dTrueDCM));
    dFiniteDifferenceJacobian(:, ui32AxisIndex) = (dErrorPlus - dErrorMinus) ./ (2.0 * dStep);
end
end
