function tests = testRotationVectorToDCM
tests = functiontests(localfunctions);
end

function setupOnce(testCase)
testFolder = fileparts(mfilename('fullpath'));
rootFolder = fullfile(testFolder, '..', '..', '..');
matlabFolder = fullfile(rootFolder, 'matlab');
rotationsFolder = fullfile(matlabFolder, 'rotations');
dcmFolder = fullfile(rotationsFolder, 'dcmLib');
quatFolder = fullfile(rotationsFolder, 'quatLib');
linearAlgebraFolder = fullfile(matlabFolder, 'linearAlgebra');
testCase.TestData.OriginalPath = path;
addpath(testFolder, matlabFolder, rotationsFolder, dcmFolder, quatFolder, linearAlgebraFolder);
end

function teardownOnce(testCase)
path(testCase.TestData.OriginalPath);
end

function testZeroRotationVectorReturnsIdentity(testCase)
dDCM = RotationVectorToDCM(zeros(3, 1));

verifyEqual(testCase, dDCM, eye(3), "AbsTol", 0.0);
verifyTrue(testCase, ValidateDCM(dDCM));
end

function testSmallAngleRotationVectorReturnsValidDCM(testCase)
dRotationVector = [1.0e-4; 2.0e-4; 3.0e-4];

dDCM = RotationVectorToDCM(dRotationVector);
dQuat = DCM2quat(dDCM, false);

verifyTrue(testCase, ValidateDCM(dDCM));
verifyLessThanOrEqual(testCase, max(abs(vecnorm(dDCM).^2 - 1.0)), eps('single'));
verifyLessThanOrEqual(testCase, norm(cross(dDCM(:, 1), dDCM(:, 2)) - dDCM(:, 3)), eps('single'));
verifyEqual(testCase, size(dQuat), [4, 1]);
verifyLessThanOrEqual(testCase, abs(norm(dQuat) - 1.0), double(eps('single')));
end
