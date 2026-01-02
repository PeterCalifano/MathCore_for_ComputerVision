function dNewIterThreshold = AdaptRansacIterThreshold(dInliersProbability, ...
                                                      dOutlierFreeProbability, ...
                                                      dSampleSize)%#codegen
arguments
    dInliersProbability     (1,1) double {mustBePositive, mustBeLessThanOrEqual(dInliersProbability, 1.0)}
    dOutlierFreeProbability (1,1) double {mustBePositive, mustBeLessThanOrEqual(dOutlierFreeProbability, 1.0)}
    dSampleSize             (1,1) double {mustBePositive}
end
%% SIGNATURE
% dNewIterThreshold = AdaptRansacIterThreshold(dInliersProbability, ...
%                                              dOutlierFreeProbability, ...
%                                              dSampleSize)%#codegen
% -------------------------------------------------------------------------------------------------------------
%% DESCRIPTION
% Computation of the iteration value based on the inliers percentage. See Hartley, Richard,
% and Andrew Zisserman. Multiple view geometry in computer vision. Cambridge university
% press, 2003, Section 4.7.1, formula 4.18.
% -------------------------------------------------------------------------------------------------------------
%% INPUTS
% dInliersProbability     (1,1) double  Inlier percentage estimate.
% dOutlierFreeProbability (1,1) double  Probability at least one sample is outlier-free.
% dSampleSize             (1,1) double  Sample size used for estimation.
% -------------------------------------------------------------------------------------------------------------
%% OUTPUTS
% dNewIterThreshold       (1,1)  New iteration count threshold
% -------------------------------------------------------------------------------------------------------------
%% CHANGELOG
% 02-01-2026    Pietro Califano     Implement from legacy code.
% -------------------------------------------------------------------------------------------------------------

% Force inlining of function
coder.inline("always");

% Apply formula
dNewIterThreshold = log10(1 - dOutlierFreeProbability) / log10( 1 - ( 1 - (1 - dInliersProbability) ) ^ dSampleSize );

end