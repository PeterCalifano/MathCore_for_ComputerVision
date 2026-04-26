function [dTargetMatrix] = GivensEliminateRow(dTargetMatrix, ...
                                              ui32TargetSubscript, ...
                                              ui32AuxRowId) %#codegen
arguments
    dTargetMatrix       (:,:) double
    ui32TargetSubscript (1,2) uint32    
    ui32AuxRowId        (1,1) uint32 = ui32TargetSubscript(1) - 1; 
end
%% SIGNATURE
% [dTargetMatrix] = GivensEliminateRow(dTargetMatrix, ui32RowsIndices) %#codegen
% -------------------------------------------------------------------------------------------------------------
%% DESCRIPTION
% Function performing elimination of a matrix entry identified by ui32TargetSubscript vector (row, col),
% using Givens row-wise rotations (i.e., rotating two rows). The algorithm can be used to perform QR
% decomposition to find R in place.
% -------------------------------------------------------------------------------------------------------------
%% INPUT
% dTargetMatrix       (:,:) double
% ui32TargetSubscript (1,2) uint32 % In Golub, 2014, ith row
% ui32AuxRowId        (1,1) uint32 = ui32TargetSubscript(1) - 1; % Default: element above target row, kth row
% -------------------------------------------------------------------------------------------------------------
%% OUTPUT
% dTargetMatrix % Target matrix transformed in place
% -------------------------------------------------------------------------------------------------------------
%% CHANGELOG
% 25-02-2025    Pietro Califano     First version implemented
% -------------------------------------------------------------------------------------------------------------
%% DEPENDENCIES
% [-]
% -------------------------------------------------------------------------------------------------------------

%% Function code

% Input checks
assert(ui32AuxRowId ~= ui32TargetSubscript(1)) 
assert(ui32TargetSubscript(1) > 1, 'ERROR: Givens row-wise rotations cannot be applied to eliminate the 1st row. Target must be below the 1st pivot!')
ui32NumOfCols = size(dTargetMatrix, 2);

if coder.target("MATLAB") || coder.target("MEX") 
    % Assert sizes against indices
    % TODO
end

% Get (a,b) entries
dVal1 = dTargetMatrix(ui32AuxRowId,           ui32TargetSubscript(2));
dVal2 = dTargetMatrix(ui32TargetSubscript(1), ui32TargetSubscript(2));

% Compute Givens rotation values
[dCos, dSin] = ComputeGivensRotValues([dVal1; dVal2]);

% Transform rows
for ui32IdCol = 1:ui32NumOfCols
    
    dTmp1 = dTargetMatrix(ui32AuxRowId, ui32IdCol); 
    dTmp2 = dTargetMatrix(ui32TargetSubscript(1), ui32IdCol);

    % Rotate dTargetMatrix, equivalent to Arot = G^T * A
    [dTargetMatrix(ui32AuxRowId, ui32IdCol), ...
     dTargetMatrix(ui32TargetSubscript(1), ui32IdCol)] = ApplyGivensRot(dTmp1, dTmp2, dCos, dSin);

end

end
