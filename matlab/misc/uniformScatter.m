function dScatteredValues = uniformScatter(dMinVec, dMaxVec, ui32NumPoints) %#codegen
arguments
    dMinVec         (:,1) double {mustBeNumeric}
    dMaxVec         (:,1) double {mustBeNumeric}
    ui32NumPoints   (1,1) uint32 {mustBeNumeric}
end
%% PROTOTYPE
% dScatteredValues = uniformScatter(dMinVec, dMaxVec, ui32NumPoints) %#codegen
% -------------------------------------------------------------------------------------------------------------
%% DESCRIPTION
% Function randomly generating uniformly distributed 1D values given lower and upper bounds of the interval
% and the desired number of points.
% -------------------------------------------------------------------------------------------------------------
%% INPUT
% dMinVec         (:,1) double {mustBeNumeric}
% dMaxVec         (:,1) double {mustBeNumeric}
% ui32NumPoints   (1,1) uint32 {mustBeNumeric}
% -------------------------------------------------------------------------------------------------------------
%% OUTPUT
% dScatteredValues
% -------------------------------------------------------------------------------------------------------------
%% CHANGELOG
% 05-03-2024        Pietro Califano     First version coded (scalar values)
% 11-03-2025        Pietro Califano     Update of function to new coding standards  
% 15-10-2025        Pietro Califano     Improve to work with vectors, scalars and arrays of scalars
% -------------------------------------------------------------------------------------------------------------
%% DEPENDENCIES
% [-]
% -------------------------------------------------------------------------------------------------------------

%% Function code
% Get size of inputs
ui32VecSpaceDim = size(dMinVec, 1);
assert(ui32VecSpaceDim == size(dMaxVec, 1), "Max and Min vectors have different dimensions")

% Construct random points depending on inputs

% if ui32NumPoints > 1 && isscalar(dMaxVec)

dRandPoints = transpose(rand(ui32NumPoints, size(dMaxVec, 1), size(dMaxVec, 2)));
dScatteredValues = dMinVec + (dRandPoints .* (dMaxVec - dMinVec));

% end



end

