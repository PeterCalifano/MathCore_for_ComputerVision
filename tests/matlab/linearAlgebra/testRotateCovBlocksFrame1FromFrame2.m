function tests = testRotateCovBlocksFrame1FromFrame2
tests = functiontests(localfunctions);
end

function setupOnce(testCase)
testFolder = fileparts(mfilename('fullpath'));
utilsFolder = fullfile(testFolder, '..');
testCase.TestData.OriginalPath = path;
addpath(testFolder, utilsFolder);
end

function teardownOnce(testCase)
path(testCase.TestData.OriginalPath);
end

function testFullRotationDefaultMask(testCase)
[dCovMatrix_Frame2, dRot3_Frame1FromFrame2] = makeTwoBlockInputs_();

dCovMatrix_Frame1 = RotateCovBlocksFrame1FromFrame2(dCovMatrix_Frame2, dRot3_Frame1FromFrame2);

[A, B, C, D] = splitBlocks_(dCovMatrix_Frame2);
expected = [rotateBlock_(dRot3_Frame1FromFrame2, A), rotateBlock_(dRot3_Frame1FromFrame2, B); ...
            rotateBlock_(dRot3_Frame1FromFrame2, C), rotateBlock_(dRot3_Frame1FromFrame2, D)];

verifyEqual(testCase, dCovMatrix_Frame1, expected, "AbsTol", 1e-12);
end

function testMixedRotationMask(testCase)
[dCovMatrix_Frame2, dRot3_Frame1FromFrame2] = makeTwoBlockInputs_();
bRotateBlockMask = [true, false];

dCovMatrix_Frame1 = RotateCovBlocksFrame1FromFrame2(dCovMatrix_Frame2, dRot3_Frame1FromFrame2, bRotateBlockMask);

[A, B, C, D] = splitBlocks_(dCovMatrix_Frame2);
expected = [rotateBlock_(dRot3_Frame1FromFrame2, A), dRot3_Frame1FromFrame2 * B; ...
            C * transpose(dRot3_Frame1FromFrame2), D];

verifyEqual(testCase, dCovMatrix_Frame1, expected, "AbsTol", 1e-12);
end

function testNoRotationMask(testCase)
[dCovMatrix_Frame2, dRot3_Frame1FromFrame2] = makeTwoBlockInputs_();
bRotateBlockMask = [false, false];

dCovMatrix_Frame1 = RotateCovBlocksFrame1FromFrame2(dCovMatrix_Frame2, dRot3_Frame1FromFrame2, bRotateBlockMask);

verifyEqual(testCase, dCovMatrix_Frame1, dCovMatrix_Frame2, "AbsTol", 0);
end

function testSingleBlockRotation(testCase)
theta = pi / 6;
dRot3_Frame1FromFrame2 = [cos(theta), -sin(theta), 0; ...
                          sin(theta), cos(theta), 0; ...
                          0, 0, 1];
dCovMatrix_Frame2 = [3, 1, 2; 0, -1, 4; 2, 5, 6];
bRotateBlockMask = true;

dCovMatrix_Frame1 = RotateCovBlocksFrame1FromFrame2(dCovMatrix_Frame2, dRot3_Frame1FromFrame2, bRotateBlockMask);

expected = rotateBlock_(dRot3_Frame1FromFrame2, dCovMatrix_Frame2);
verifyEqual(testCase, dCovMatrix_Frame1, expected, "AbsTol", 1e-12);
end

function testMaskLengthMismatch(testCase)
[dCovMatrix_Frame2, dRot3_Frame1FromFrame2] = makeTwoBlockInputs_();
bRotateBlockMask = true;

verifyError(testCase, ...
    @() RotateCovBlocksFrame1FromFrame2(dCovMatrix_Frame2, dRot3_Frame1FromFrame2, bRotateBlockMask), ...
    "RotateCovBlocksFrame1FromFrame2:InvalidInput");
end

function testNonSquareMatrixError(testCase)
dCovMatrix_Frame2 = zeros(3, 6);
dRot3_Frame1FromFrame2 = eye(3);
bRotateBlockMask = true;

verifyError(testCase, ...
    @() RotateCovBlocksFrame1FromFrame2(dCovMatrix_Frame2, dRot3_Frame1FromFrame2, bRotateBlockMask), ...
    "RotateCovBlocksFrame1FromFrame2:InvalidInput");
end

function testSizeNotMultipleOfThreeError(testCase)
dCovMatrix_Frame2 = zeros(4, 4);
dRot3_Frame1FromFrame2 = eye(3);
bRotateBlockMask = true;

verifyError(testCase, ...
    @() RotateCovBlocksFrame1FromFrame2(dCovMatrix_Frame2, dRot3_Frame1FromFrame2, bRotateBlockMask), ...
    "RotateCovBlocksFrame1FromFrame2:InvalidInput");
end

function [dCovMatrix_Frame2, dRot3_Frame1FromFrame2] = makeTwoBlockInputs_()
theta = pi / 6;
dRot3_Frame1FromFrame2 = [cos(theta), -sin(theta), 0; ...
                          sin(theta), cos(theta), 0; ...
                          0, 0, 1];

A = [1, 2, 3; 4, 5, 6; 7, 8, 10];
B = [2, 0, 1; 1, -1, 3; 0, 2, 4];
C = [-1, 2, 0; 5, 3, -2; 4, -1, 1];
D = [3, 1, 0; 2, 2, 1; -1, 0, 1];

dCovMatrix_Frame2 = [A, B; C, D];
end

function [A, B, C, D] = splitBlocks_(dCovMatrix_Frame2)
A = dCovMatrix_Frame2(1:3, 1:3);
B = dCovMatrix_Frame2(1:3, 4:6);
C = dCovMatrix_Frame2(4:6, 1:3);
D = dCovMatrix_Frame2(4:6, 4:6);
end

function dRotatedCovBlock = rotateBlock_(dRotMatrix, dCovMatrix)
dRotatedCovBlock = dRotMatrix * dCovMatrix * transpose(dRotMatrix);
end
