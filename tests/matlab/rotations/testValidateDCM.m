function tests = testValidateDCM
%% SIGNATURE
% tests = testValidateDCM
% -------------------------------------------------------------------------------------------------------------
%% DESCRIPTION
% Validate both branches of ValidateDCM: the assert branch errors on an invalid DCM, the output branch
% returns false without error. Cover real rotations, reflections, scaling, non-finite and complex entries,
% and the column-norm tolerance boundary.
% -------------------------------------------------------------------------------------------------------------
%% INPUT
% None.
% -------------------------------------------------------------------------------------------------------------
%% OUTPUT
% tests    Function-based MATLAB test suite.
% -------------------------------------------------------------------------------------------------------------
%% CHANGELOG
% 08-11-2026  Pietro Califano, Claude Opus 5.5    Test the assert and output branches of ValidateDCM.
% -------------------------------------------------------------------------------------------------------------
%% DEPENDENCIES
% ValidateDCM.
% -------------------------------------------------------------------------------------------------------------
tests = functiontests(localfunctions);
end

function setupOnce(testCase)
charTestFolder = fileparts(mfilename('fullpath'));
charMatlabFolder = fullfile(charTestFolder, '..', '..', '..', 'matlab');
charRotationsFolder = fullfile(charMatlabFolder, 'rotations');
testCase.TestData.charOriginalPath = path;
addpath(charMatlabFolder, charRotationsFolder, fullfile(charRotationsFolder, 'quatLib'));
end

function teardownOnce(testCase)
path(testCase.TestData.charOriginalPath);
end

function testValidRotationsInBothBranches(testCase)
dTol = eps('single');

% Build a random rotation from a seeded local stream and the matrix exponential
objStream = RandStream('mt19937ar', 'Seed', 20261108);
dRotVec = randn(objStream, 3, 1);
dSkew = [0.0,         -dRotVec(3), dRotVec(2);
         dRotVec(3),  0.0,         -dRotVec(1);
         -dRotVec(2), dRotVec(1),  0.0];
dRandomDCM = expm(dSkew);

cellValidDCMs = {eye(3), dRandomDCM};

for idDCM = 1:numel(cellValidDCMs)
    verifyTrue(testCase, ValidateDCM(cellValidDCMs{idDCM}, dTol, true));
    verifyTrue(testCase, ValidateDCM(cellValidDCMs{idDCM}, dTol, false));
end
end

function testOutputBranchReturnsFalseWithoutError(testCase)
dTol = eps('single');
cellInvalidDCMs = BuildInvalidDCMs_();

for idDCM = 1:numel(cellInvalidDCMs)
    % Fail on any error or warning, then check the returned flag
    bValidity = verifyWarningFree(testCase, @() ValidateDCM(cellInvalidDCMs{idDCM}, dTol, false));
    verifyClass(testCase, bValidity, 'logical');
    verifyFalse(testCase, bValidity);
end
end

function testAssertBranchErrorsOnInvalidMatrices(testCase)
dTol = eps('single');
cellInvalidDCMs = BuildInvalidDCMs_();

for idDCM = 1:numel(cellInvalidDCMs)
    verifyError(testCase, @() ValidateDCM(cellInvalidDCMs{idDCM}, dTol, true), ?MException);
end
end

function testDefaultBranchIsAssertBranch(testCase)
dTol = eps('single');
cellInvalidDCMs = BuildInvalidDCMs_();

for idDCM = 1:numel(cellInvalidDCMs)
    verifyError(testCase, @() ValidateDCM(cellInvalidDCMs{idDCM}, dTol), ?MException);
end
end

function testRealValuedAssertComesFirst(testCase)
dTol = eps('single');

% Complex and non-finite: the real-valued assert must fire before the finite one
dComplexNonFiniteDCM = diag([1i, -1i, NaN]);

charMessage = '';
try
    ValidateDCM(dComplexNonFiniteDCM, dTol, true);
catch objException
    charMessage = objException.message;
end

verifySubstring(testCase, charMessage, 'real-valued');
end

function testColumnNormToleranceBoundary(testCase)
dTol = 1.0e-6;

% Scale the first column so that abs(dot(C1, C1) - 1) is 0.9*dTol or 1.1*dTol. The cross-product error is
% about half of that, so only the column-norm check decides the result.
dBelowDCM = diag([sqrt(1.0 + 0.9*dTol), 1.0, 1.0]);
dAboveDCM = diag([sqrt(1.0 + 1.1*dTol), 1.0, 1.0]);

verifyLessThan(testCase, norm(cross(dAboveDCM(:,1), dAboveDCM(:,2)) - dAboveDCM(:,3)), dTol);

verifyTrue(testCase, ValidateDCM(dBelowDCM, dTol, false));
verifyTrue(testCase, ValidateDCM(dBelowDCM, dTol, true));

verifyFalse(testCase, ValidateDCM(dAboveDCM, dTol, false));
verifyError(testCase, @() ValidateDCM(dAboveDCM, dTol, true), ?MException);
end

function cellInvalidDCMs = BuildInvalidDCMs_()
% Build the invalid matrices shared by the branch tests
dNaNDCM = eye(3);
dNaNDCM(2,3) = NaN;

dInfDCM = eye(3);
dInfDCM(1,2) = Inf;

cellInvalidDCMs = {diag([1.0, 1.0, -1.0]), ...
                   2.0 * eye(3), ...
                   dNaNDCM, ...
                   dInfDCM, ...
                   diag([1i, -1i, 1.0])};
end
