function [dRotEntry1, dRotEntry2] = ApplyGivensRot(dEntry1, ...
                                                   dEntry2, ...
                                                   dCosTheta, ...
                                                   dSinTheta) %#codegen
arguments
    dEntry1     (1,1) double
    dEntry2     (1,1) double
    dCosTheta   (1,1) double
    dSinTheta   (1,1) double
end
%% SIGNATURE
% [dRotEntry1, dRotEntry2] = ApplyGivensRot(dEntry1, dEntry2, dCosTheta, dSinTheta) %#codegen
% -------------------------------------------------------------------------------------------------------------
%% DESCRIPTION
% Applies a 2-by-2 Givens rotation to a pair of scalar entries using the same convention adopted by the
% row-wise and column-wise Givens utilities in MathCore:
%
%   [dRotEntry1]   [ dCosTheta  -dSinTheta ] [dEntry1]
%   [dRotEntry2] = [ dSinTheta   dCosTheta ] [dEntry2]
%
% This primitive is the shared ownership point for pair-wise Givens application. Higher-level routines
% such as QR elimination and SRIF wrappers should build their algorithm-specific flow around this function.
% -------------------------------------------------------------------------------------------------------------
%% OUTPUT
% dRotEntry1    (1,1) double
% dRotEntry2    (1,1) double
% -------------------------------------------------------------------------------------------------------------
%% CHANGELOG
% 24-04-2026    Pietro Califano        First shared primitive for pair-wise Givens application.
% -------------------------------------------------------------------------------------------------------------

%% Function code
dRotEntry1 = dCosTheta * dEntry1 - dSinTheta * dEntry2;
dRotEntry2 = dSinTheta * dEntry1 + dCosTheta * dEntry2;

end
