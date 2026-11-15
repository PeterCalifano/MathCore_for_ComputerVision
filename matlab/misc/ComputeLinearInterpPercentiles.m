function dPercentileValues = ComputeLinearInterpPercentiles(dValues, dPercentiles, ui32Dim)
%% SIGNATURE
% dPercentileValues = ComputeLinearInterpPercentiles(dValues, dPercentiles, ui32Dim)
% -------------------------------------------------------------------------------------------------------------
%% DESCRIPTION
% Compute sample percentiles by linear interpolation between order statistics, definition 7 of
% Hyndman and Fan (1996) and the default of numpy.percentile and of R quantile. For n samples sorted
% as x(1) <= ... <= x(n), percentile p sits at the 1-based position h = 1 + (n - 1)*p/100, and
%     y = x(k) + (h - k)*(x(k + 1) - x(k)),  with k = floor(h).
% Percentiles 0 and 100 return the minimum and the maximum.
% MATLAB prctile (R2024b) implements only the midpoint rule (definition 5, x(i) at 100*(i - 0.5)/n).
% The two rules differ most in sparse tails, for example in small per-run counts.
% NaN samples are omitted in each slice independently. A slice without samples, including an empty
% vector, returns NaN, as prctile does.
% -------------------------------------------------------------------------------------------------------------
%% INPUT
% dValues            [any]  Real samples, reduced along dimension ui32Dim.
% dPercentiles       [1xP]  Percentiles in [0, 100].
% ui32Dim            [1]    (Optional) Dimension of dValues to reduce. Default: first non-singleton
%                           dimension.
% -------------------------------------------------------------------------------------------------------------
%% OUTPUT
% dPercentileValues  [any]  Array of size(dValues) with size(dValues, ui32Dim) replaced by P;
%                           entry k along ui32Dim holds percentile dPercentiles(k).
% -------------------------------------------------------------------------------------------------------------
%% CHANGELOG
% 15-11-2026  Pietro Califano, Claude Opus 5.5  Add linear-interpolation percentiles.
% -------------------------------------------------------------------------------------------------------------
%% DEPENDENCIES
% None.
% -------------------------------------------------------------------------------------------------------------

arguments (Input)
    dValues {mustBeNumeric, mustBeReal}
    dPercentiles (1, :) double {mustBeReal, mustBeInRange(dPercentiles, 0, 100)}
    ui32Dim (1, 1) uint32 {mustBePositive} = FirstNonSingletonDim_(dValues)
end

% Move the reduced dimension first and flatten the remaining ones into columns.
dValues = double(dValues);
dDim = double(ui32Dim);
dOrder = [dDim, setdiff(1:max(ndims(dValues), dDim), dDim)];
dSize = size(dValues, dOrder);
dColumns = reshape(permute(dValues, dOrder), dSize(1), prod(dSize(2:end)));

% Interpolate between neighboring order statistics of the non-NaN samples of each column.
dPercentileColumns = nan(numel(dPercentiles), size(dColumns, 2));
for idCol = 1:size(dColumns, 2)

    dSorted = sort(dColumns(~isnan(dColumns(:, idCol)), idCol));
    dNumSamples = numel(dSorted);
    if dNumSamples == 0
        continue
    end

    dPosition = 1 + (dNumSamples - 1) * dPercentiles(:) / 100;
    dLowerIdx = floor(dPosition);
    dUpperIdx = min(dLowerIdx + 1, dNumSamples);
    dFraction = dPosition - dLowerIdx;

    % Interpolate only at fractional positions, so that infinite samples do not give Inf - Inf.
    dResult = dSorted(dLowerIdx);
    bInterp = dFraction > 0;
    dResult(bInterp) = dResult(bInterp) + ...
        dFraction(bInterp) .* (dSorted(dUpperIdx(bInterp)) - dResult(bInterp));
    dPercentileColumns(:, idCol) = dResult;

end

% Restore the original dimension order, with the percentiles along the reduced dimension.
dSize(1) = numel(dPercentiles);
dPercentileValues = ipermute(reshape(dPercentileColumns, dSize), dOrder);

end

%% Local functions
function ui32Dim = FirstNonSingletonDim_(dValues)
% Match the default dimension of prctile; an all-singleton input reduces along dimension 1.
ui32Dim = uint32(find(size(dValues) ~= 1, 1));
if isempty(ui32Dim)
    ui32Dim = uint32(1);
end
end
