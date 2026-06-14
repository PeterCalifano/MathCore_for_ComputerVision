function tests = testRodriguezRotation
tests = functiontests(localfunctions);
end

function setupOnce(testCase)
testFolder = fileparts(mfilename('fullpath'));
rootFolder = fullfile(testFolder, '..', '..', '..');
matlabFolder = fullfile(rootFolder, 'matlab');
rotationsFolder = fullfile(matlabFolder, 'rotations');
testCase.TestData.OriginalPath = path;
addpath(testFolder, matlabFolder, rotationsFolder);
end

function teardownOnce(testCase)
path(testCase.TestData.OriginalPath);
end

function testRot3dVecAboutDirFullVector(testCase)
dTestVector = [1, 1, 0;
               0, 0, 1;
               0, 0, 0];
dDirections = [0, 0, 0;
               0, 0, 0;
               1, 1, 1];
dAngles = deg2rad([90, 0, -90]);

dRotateVectors = Rot3dVecAboutDir(dDirections, dTestVector, dAngles);

verifyEqual(testCase, dRotateVectors(:, 1), [0; 1; 0], "AbsTol", 1.0e-12);
verifyEqual(testCase, dRotateVectors(:, 2), [1; 0; 0], "AbsTol", 1.0e-12);
verifyEqual(testCase, dRotateVectors(:, 3), [1; 0; 0], "AbsTol", 1.0e-12);
end

function testRot3dVecAboutDirConvenienceCases(testCase)
dTestVector = [1, 1, 0;
               0, 0, 1;
               0, 0, 0];
dDirection = [0; 0; 1];
dAngles = deg2rad(90);

dRotateVectors = Rot3dVecAboutDir(dDirection, dTestVector, dAngles);

verifyEqual(testCase, dRotateVectors(:, 1), [0; 1; 0], "AbsTol", 1.0e-12);
verifyEqual(testCase, dRotateVectors(:, 2), [0; 1; 0], "AbsTol", 1.0e-12);
verifyEqual(testCase, dRotateVectors(:, 3), [-1; 0; 0], "AbsTol", 1.0e-12);
end
