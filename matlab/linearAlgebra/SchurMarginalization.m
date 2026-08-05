function [dRetainedConditionalCov, dSelectedConditionalCov] = ...
        SchurMarginalization(dCovMatrix, bMarginalStateMask) %#codegen
%% SIGNATURE
% [dRetainedConditionalCov, dSelectedConditionalCov] = ...
%     SchurMarginalization(dCovMatrix, bMarginalStateMask) %#codegen
% -------------------------------------------------------------------------------------------------------------
%% DESCRIPTION
% Compute both conditional covariance blocks of a joint covariance using
% Schur complements. The logical mask selects the state block being
% marginalized/conditioned upon and may describe leading, trailing, middle,
% or noncontiguous entries. Outputs preserve the complete input allocation
% and contain zeros outside their corresponding original state indices.
% -------------------------------------------------------------------------------------------------------------
%% INPUT
% dCovMatrix             Joint covariance matrix.
% bMarginalStateMask     Full-size logical mask selecting the marginalized block.
% -------------------------------------------------------------------------------------------------------------
%% OUTPUT
% dRetainedConditionalCov    Conditional covariance of the unselected states.
% dSelectedConditionalCov    Conditional covariance of the selected states.
% -------------------------------------------------------------------------------------------------------------
%% CHANGELOG
% 21-04-2024  Pietro Califano     First simple version coded.
% 03-03-2025  Pietro Califano     Update function to compute both conditional subblocks.
% 05-08-2026  Pietro Califano, Codex gpt-5.6     Generalize selection and correct the Schur API name.
% -------------------------------------------------------------------------------------------------------------
%% DEPENDENCIES
% None.
% -------------------------------------------------------------------------------------------------------------
arguments (Input)
    dCovMatrix         (:,:) double
    bMarginalStateMask (:,1)
end

arguments (Output)
    dRetainedConditionalCov (:,:) double
    dSelectedConditionalCov (:,:) double
end

ui32CovarianceSize = size(dCovMatrix, 1);
dCovarianceScale = max(1.0, norm(dCovMatrix, 'fro'));
dSymmetryTolerance = 100.0 * eps(dCovarianceScale);

% Validate the joint allocation and mask before any logical indexing. The
% selected and retained partitions must both be nonempty by contract.
if coder.target('MATLAB') || coder.target('mex')
    if size(dCovMatrix, 2) ~= ui32CovarianceSize || ...
            any(~isfinite(dCovMatrix), 'all') || ...
            norm(dCovMatrix - transpose(dCovMatrix), 'fro') > dSymmetryTolerance
        error('SchurMarginalization:InvalidCovariance', ...
            'Input covariance must be finite, square, and symmetric.');
    end

    if ~islogical(bMarginalStateMask) || ...
            numel(bMarginalStateMask) ~= ui32CovarianceSize || ...
            ~any(bMarginalStateMask) || all(bMarginalStateMask)
        error('SchurMarginalization:InvalidMask', ...
            'Selection mask must be logical, full-size, and select a proper nonempty subset.');
    end
end

% Extract the selected and retained covariance blocks using the logical mask.
bRetainedStateMask = ~bMarginalStateMask;
dSelectedCov = dCovMatrix(bMarginalStateMask, bMarginalStateMask);
dRetainedCov = dCovMatrix(bRetainedStateMask, bRetainedStateMask);
dRetainedSelectedCov = dCovMatrix(bRetainedStateMask, bMarginalStateMask);

% Factor the selected covariance and form the retained Schur complement
% without explicitly inverting the selected block.
[dSelectedCholFactor, dSelectedCholStatus] = chol(dSelectedCov, 'lower');
if dSelectedCholStatus ~= 0.0
    error('SchurMarginalization:SelectedBlockNotPositiveDefinite', ...
          'Selected covariance block must be positive definite.');
end

% Update the retained conditional block using the Cholesky factorization of the selected block.
dSelectedSolve = dSelectedCholFactor \ transpose(dRetainedSelectedCov);
dRetainedConditionalBlock = dRetainedCov - ...
    transpose(dSelectedSolve) * dSelectedSolve;

dRetainedConditionalBlock = 0.5 * ...
    (dRetainedConditionalBlock + transpose(dRetainedConditionalBlock));

dRetainedConditionalCov = zeros(size(dCovMatrix));
dSelectedConditionalCov = zeros(size(dCovMatrix));
dRetainedConditionalCov(bRetainedStateMask, bRetainedStateMask) = ...
    dRetainedConditionalBlock;

% Compute the reciprocal (selected) conditional block only when requested. This keeps
% one-output callers valid when unused fixed-allocation retained rows are zero.
if coder.const(nargout > 1)
    [dRetainedCholFactor, dRetainedCholStatus] = chol(dRetainedCov, 'lower');
    if dRetainedCholStatus ~= 0.0
        error('SchurMarginalization:RetainedBlockNotPositiveDefinite', ...
              'Retained covariance block must be positive definite when its conditional block is requested.');
    end

    dRetainedSolve = dRetainedCholFactor \ dRetainedSelectedCov;
    dSelectedConditionalBlock = dSelectedCov - ...
        transpose(dRetainedSolve) * dRetainedSolve;
        
    dSelectedConditionalBlock = 0.5 * ...
        (dSelectedConditionalBlock + transpose(dSelectedConditionalBlock));
    dSelectedConditionalCov(bMarginalStateMask, bMarginalStateMask) = ...
        dSelectedConditionalBlock;
end

end
