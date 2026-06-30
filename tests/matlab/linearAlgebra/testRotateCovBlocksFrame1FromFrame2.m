function tests = testRotateCovBlocksFrame1FromFrame2
tests = functiontests(localfunctions);
end

function setupOnce(testCase)
testFolder = fileparts(mfilename('fullpath'));
rootFolder = fullfile(testFolder, '..', '..', '..');
matlabFolder = fullfile(rootFolder, 'matlab');
linearAlgebraFolder = fullfile(matlabFolder, 'linearAlgebra');
testCase.TestData.OriginalPath = path;
addpath(testFolder, matlabFolder, linearAlgebraFolder);
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

function testVectorizedFullRotationDefaultMask(testCase)
[dCovMatrix_Frame2, dRot3_Frame1FromFrame2] = makeTwoBlockInputsVectorized_(3);

dCovMatrix_Frame1 = RotateCovBlocksFrame1FromFrame2(dCovMatrix_Frame2, ...
                                                    dRot3_Frame1FromFrame2, ...
                                                    [true, true], ...
                                                    true);

expected = zeros(size(dCovMatrix_Frame2));
for page = 1:size(dCovMatrix_Frame2, 3)
    expected(:, :, page) = computeExpectedTwoBlock_(dCovMatrix_Frame2(:, :, page), ...
                                                    dRot3_Frame1FromFrame2(:, :, page), ...
                                                    [true, true]);
end

verifyEqual(testCase, dCovMatrix_Frame1, expected, "AbsTol", 1e-12);
end

function testVectorizedMixedRotationMask(testCase)
[dCovMatrix_Frame2, dRot3_Frame1FromFrame2] = makeTwoBlockInputsVectorized_(2);
bRotateBlockMask = [true, false];

dCovMatrix_Frame1 = RotateCovBlocksFrame1FromFrame2(dCovMatrix_Frame2, ...
                                                    dRot3_Frame1FromFrame2, ...
                                                    bRotateBlockMask, ...
                                                    true);

expected = zeros(size(dCovMatrix_Frame2));
for page = 1:size(dCovMatrix_Frame2, 3)
    expected(:, :, page) = computeExpectedTwoBlock_(dCovMatrix_Frame2(:, :, page), ...
                                                    dRot3_Frame1FromFrame2(:, :, page), ...
                                                    bRotateBlockMask);
end

verifyEqual(testCase, dCovMatrix_Frame1, expected, "AbsTol", 1e-12);
end

function testVectorizedNoRotationMask(testCase)
[dCovMatrix_Frame2, dRot3_Frame1FromFrame2] = makeTwoBlockInputsVectorized_(2);
bRotateBlockMask = [false, false];

dCovMatrix_Frame1 = RotateCovBlocksFrame1FromFrame2(dCovMatrix_Frame2, ...
                                                    dRot3_Frame1FromFrame2, ...
                                                    bRotateBlockMask, ...
                                                    true);

verifyEqual(testCase, dCovMatrix_Frame1, dCovMatrix_Frame2, "AbsTol", 0);
end

function testVectorizedThirdDimMismatchError(testCase)
[dCovMatrix_Frame2, dRot3_Frame1FromFrame2] = makeTwoBlockInputsVectorized_(2);
dRot3_Frame1FromFrame2 = dRot3_Frame1FromFrame2(:, :, 1);

verifyError(testCase, ...
    @() RotateCovBlocksFrame1FromFrame2(dCovMatrix_Frame2, dRot3_Frame1FromFrame2, [true, true], true), ...
    "RotateCovBlocksFrame1FromFrame2:InvalidInput");
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

function [dCovMatrix_Frame2, dRot3_Frame1FromFrame2] = makeTwoBlockInputsVectorized_(ui32NumPages)
[dBaseCovMatrix, ~] = makeTwoBlockInputs_();

dCovMatrix_Frame2 = zeros([size(dBaseCovMatrix), ui32NumPages]);
dRot3_Frame1FromFrame2 = zeros(3, 3, ui32NumPages);

for page = 1:ui32NumPages
    theta = (pi / 12) * page;
    dRot3_Frame1FromFrame2(:, :, page) = [cos(theta), -sin(theta), 0; ...
                                         sin(theta), cos(theta), 0; ...
                                         0, 0, 1];
    dCovMatrix_Frame2(:, :, page) = dBaseCovMatrix + 0.1 * page;
end
end

function [A, B, C, D] = splitBlocks_(dCovMatrix_Frame2)
A = dCovMatrix_Frame2(1:3, 1:3);
B = dCovMatrix_Frame2(1:3, 4:6);
C = dCovMatrix_Frame2(4:6, 1:3);
D = dCovMatrix_Frame2(4:6, 4:6);
end

function expected = computeExpectedTwoBlock_(dCovMatrix_Frame2, dRot3_Frame1FromFrame2, bRotateBlockMask)
[A, B, C, D] = splitBlocks_(dCovMatrix_Frame2);

if bRotateBlockMask(1) && bRotateBlockMask(2)
    expected = [rotateBlock_(dRot3_Frame1FromFrame2, A), rotateBlock_(dRot3_Frame1FromFrame2, B); ...
                rotateBlock_(dRot3_Frame1FromFrame2, C), rotateBlock_(dRot3_Frame1FromFrame2, D)];
elseif bRotateBlockMask(1) && ~bRotateBlockMask(2)
    expected = [rotateBlock_(dRot3_Frame1FromFrame2, A), dRot3_Frame1FromFrame2 * B; ...
                C * transpose(dRot3_Frame1FromFrame2), D];
elseif ~bRotateBlockMask(1) && bRotateBlockMask(2)
    expected = [A, B * transpose(dRot3_Frame1FromFrame2); ...
                dRot3_Frame1FromFrame2 * C, rotateBlock_(dRot3_Frame1FromFrame2, D)];
else
    expected = dCovMatrix_Frame2;
end
end

function dRotatedCovBlock = rotateBlock_(dRotMatrix, dCovMatrix)
dRotatedCovBlock = dRotMatrix * dCovMatrix * transpose(dRotMatrix);
end
