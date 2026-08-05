function tests = testSchurMarginalization
%% SIGNATURE
% tests = testSchurMarginalization
% -------------------------------------------------------------------------------------------------------------
%% DESCRIPTION
% Validate arbitrary-index Schur complements, full-size output placement,
% permutation invariance, and covariance/mask failure diagnostics.
% -------------------------------------------------------------------------------------------------------------
%% INPUT
% None.
% -------------------------------------------------------------------------------------------------------------
%% OUTPUT
% tests    Function-based MATLAB test suite.
% -------------------------------------------------------------------------------------------------------------
%% CHANGELOG
% 05-08-2026  Pietro Califano, Codex gpt-5.6     First implementation.
% -------------------------------------------------------------------------------------------------------------
%% DEPENDENCIES
% SchurMarginalization.
% -------------------------------------------------------------------------------------------------------------

tests = functiontests(localfunctions);
end

function setupOnce(testCase)
charTestFolder = fileparts(mfilename('fullpath'));
charRepoRoot = fullfile(charTestFolder, '..', '..', '..');
charLinearAlgebraFolder = fullfile(charRepoRoot, 'matlab', 'linearAlgebra');

testCase.TestData.charOriginalPath = path;
addpath(charLinearAlgebraFolder);
end

function teardownOnce(testCase)
path(testCase.TestData.charOriginalPath);
end

function testArbitraryMasksMatchIndependentPartitionReference(testCase)
dCovariance = BuildDenseCovariance_();
cellSelectedIdx = {[1, 2], [5, 6], [3, 4], [2, 5]};

for ui32CaseIdx = 1:numel(cellSelectedIdx)
    bSelectedStateMask = false(6, 1);
    bSelectedStateMask(cellSelectedIdx{ui32CaseIdx}) = true;

    [dRetainedConditionalCov, dSelectedConditionalCov] = ...
        SchurMarginalization(dCovariance, bSelectedStateMask);
    [dExpectedRetainedCov, dExpectedSelectedCov] = ...
        ComputePartitionReference_(dCovariance, bSelectedStateMask);

    verifyEqual(testCase, dRetainedConditionalCov, dExpectedRetainedCov, ...
                'AbsTol', 2.0e-13);
    verifyEqual(testCase, dSelectedConditionalCov, dExpectedSelectedCov, ...
                'AbsTol', 2.0e-13);
end
end

function testOutputsPreserveOriginalAllocationAndIndexPlacement(testCase)
dCovariance = BuildDenseCovariance_();
bSelectedStateMask = logical([false; true; false; false; true; false]);
bRetainedStateMask = ~bSelectedStateMask;

[dRetainedConditionalCov, dSelectedConditionalCov] = ...
    SchurMarginalization(dCovariance, bSelectedStateMask);

verifySize(testCase, dRetainedConditionalCov, size(dCovariance));
verifySize(testCase, dSelectedConditionalCov, size(dCovariance));
verifyEqual(testCase, dRetainedConditionalCov(bSelectedStateMask, :), ...
            zeros(nnz(bSelectedStateMask), size(dCovariance, 2)), 'AbsTol', 0.0);
verifyEqual(testCase, dRetainedConditionalCov(:, bSelectedStateMask), ...
            zeros(size(dCovariance, 1), nnz(bSelectedStateMask)), 'AbsTol', 0.0);
verifyEqual(testCase, dSelectedConditionalCov(bRetainedStateMask, :), ...
            zeros(nnz(bRetainedStateMask), size(dCovariance, 2)), 'AbsTol', 0.0);
verifyEqual(testCase, dSelectedConditionalCov(:, bRetainedStateMask), ...
            zeros(size(dCovariance, 1), nnz(bRetainedStateMask)), 'AbsTol', 0.0);
end

function testPermutationProducesEquivalentOriginalIndexResult(testCase)
dCovariance = BuildDenseCovariance_();
bSelectedStateMask = logical([true; false; false; true; false; true]);
ui32Permutation = [4, 1, 6, 2, 5, 3];

[dRetainedConditionalCov, dSelectedConditionalCov] = ...
    SchurMarginalization(dCovariance, bSelectedStateMask);

dPermutedCovariance = dCovariance(ui32Permutation, ui32Permutation);
bPermutedSelectedMask = bSelectedStateMask(ui32Permutation);
[dPermutedRetainedCov, dPermutedSelectedCov] = ...
    SchurMarginalization(dPermutedCovariance, bPermutedSelectedMask);

dRestoredRetainedCov = zeros(size(dCovariance));
dRestoredSelectedCov = zeros(size(dCovariance));
dRestoredRetainedCov(ui32Permutation, ui32Permutation) = dPermutedRetainedCov;
dRestoredSelectedCov(ui32Permutation, ui32Permutation) = dPermutedSelectedCov;

verifyEqual(testCase, dRestoredRetainedCov, dRetainedConditionalCov, ...
            'AbsTol', 2.0e-13);
verifyEqual(testCase, dRestoredSelectedCov, dSelectedConditionalCov, ...
            'AbsTol', 2.0e-13);
end

function testLegacyTrailingPartitionMatchesCapturedResult(testCase)
dCovarianceRoot = [2.0, 0.0, 0.0, 0.0; ...
                   0.3, 1.5, 0.0, 0.0; ...
                   0.2, -0.1, 1.2, 0.0; ...
                   0.4, 0.2, 0.1, 1.1];
dCovariance = dCovarianceRoot * transpose(dCovarianceRoot);
bSelectedStateMask = logical([false; false; true; true]);

[dRetainedConditionalCov, dSelectedConditionalCov] = ...
    SchurMarginalization(dCovariance, bSelectedStateMask);

dExpectedRetainedCov = zeros(4, 4);
dExpectedRetainedCov(1:2, 1:2) = ...
    [3.48852836709225, 0.392531438993952; ...
     0.392531438993952, 2.20179034270903];
dExpectedSelectedCov = zeros(4, 4);
dExpectedSelectedCov(3:4, 3:4) = [1.44, 0.12; 0.12, 1.22];

verifyEqual(testCase, dRetainedConditionalCov, dExpectedRetainedCov, ...
            'AbsTol', 2.0e-14);
verifyEqual(testCase, dSelectedConditionalCov, dExpectedSelectedCov, ...
            'AbsTol', 2.0e-14);
end

function testRejectsNonSquareCovariance(testCase)
verifyError(testCase, ...
    @() SchurMarginalization(ones(3, 4), logical([true; false; false])), ...
    'SchurMarginalization:InvalidCovariance');
end

function testRejectsNonFiniteCovariance(testCase)
dCovariance = eye(4);
dCovariance(2, 3) = NaN;
dCovariance(3, 2) = NaN;

verifyError(testCase, ...
    @() SchurMarginalization(dCovariance, logical([false; false; true; true])), ...
    'SchurMarginalization:InvalidCovariance');
end

function testRejectsNonSymmetricCovariance(testCase)
dCovariance = eye(4);
dCovariance(1, 2) = 0.1;

verifyError(testCase, ...
    @() SchurMarginalization(dCovariance, logical([false; false; true; true])), ...
    'SchurMarginalization:InvalidCovariance');
end

function testRejectsInvalidMasks(testCase)
dCovariance = eye(4);

verifyError(testCase, ...
    @() SchurMarginalization(dCovariance, logical([true; false; false])), ...
    'SchurMarginalization:InvalidMask');
verifyError(testCase, ...
    @() SchurMarginalization(dCovariance, false(4, 1)), ...
    'SchurMarginalization:InvalidMask');
verifyError(testCase, ...
    @() SchurMarginalization(dCovariance, true(4, 1)), ...
    'SchurMarginalization:InvalidMask');
verifyError(testCase, ...
    @() SchurMarginalization(dCovariance, [true; false; true; false] + 0.0), ...
    'SchurMarginalization:InvalidMask');
end

function testRejectsNonPositiveDefiniteSelectedBlock(testCase)
dCovariance = eye(4);
dCovariance(3:4, 3:4) = [1.0, 2.0; 2.0, 1.0];

verifyError(testCase, ...
    @() SchurMarginalization(dCovariance, logical([false; false; true; true])), ...
    'SchurMarginalization:SelectedBlockNotPositiveDefinite');
end

function testRejectsNonPositiveDefiniteRetainedBlockWhenRequested(testCase)
dCovariance = eye(4);
dCovariance(1:2, 1:2) = [1.0, 2.0; 2.0, 1.0];

verifyError(testCase, ...
    @() CallBothOutputs_(dCovariance, logical([false; false; true; true])), ...
    'SchurMarginalization:RetainedBlockNotPositiveDefinite');
end

function [dRetainedConditionalCov, dSelectedConditionalCov] = ...
        ComputePartitionReference_(dCovariance, bSelectedStateMask)
bRetainedStateMask = ~bSelectedStateMask;
dRetainedCov = dCovariance(bRetainedStateMask, bRetainedStateMask);
dSelectedCov = dCovariance(bSelectedStateMask, bSelectedStateMask);
dRetainedSelectedCov = dCovariance(bRetainedStateMask, bSelectedStateMask);

dRetainedConditionalBlock = dRetainedCov - ...
    dRetainedSelectedCov / dSelectedCov * transpose(dRetainedSelectedCov);
dSelectedConditionalBlock = dSelectedCov - ...
    transpose(dRetainedSelectedCov) / dRetainedCov * dRetainedSelectedCov;

dRetainedConditionalCov = zeros(size(dCovariance));
dSelectedConditionalCov = zeros(size(dCovariance));
dRetainedConditionalCov(bRetainedStateMask, bRetainedStateMask) = ...
    dRetainedConditionalBlock;
dSelectedConditionalCov(bSelectedStateMask, bSelectedStateMask) = ...
    dSelectedConditionalBlock;
end

function dCovariance = BuildDenseCovariance_()
dCovarianceRoot = [1.8, 0.0, 0.0, 0.0, 0.0, 0.0; ...
                   0.2, 1.5, 0.0, 0.0, 0.0, 0.0; ...
                  -0.1, 0.3, 1.4, 0.0, 0.0, 0.0; ...
                   0.4, -0.2, 0.1, 1.3, 0.0, 0.0; ...
                   0.1, 0.25, -0.15, 0.2, 1.2, 0.0; ...
                  -0.3, 0.1, 0.2, -0.1, 0.15, 1.1];
dCovariance = dCovarianceRoot * transpose(dCovarianceRoot);
end

function CallBothOutputs_(dCovariance, bSelectedStateMask)
[~, ~] = SchurMarginalization(dCovariance, bSelectedStateMask);
end
