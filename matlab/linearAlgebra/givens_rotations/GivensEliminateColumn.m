function [dTargetMatrix] = GivensEliminateColumn(dTargetMatrix, ...
                                                 ui32TargetSubscript, ...
                                                 ui32AuxColId) %#codegen
arguments
    dTargetMatrix       (:,:) double
    ui32TargetSubscript (1,2) uint32    
    ui32AuxColId        (1,1) uint32 = ui32TargetSubscript(2) - 1; 
end
%% SIGNATURE
% [dTargetMatrix] = GivensEliminateColumn(dTargetMatrix, ui32ColumnsIndices) %#codegen
% -------------------------------------------------------------------------------------------------------------
%% DESCRIPTION
% Function performing elimination of a matrix entry identified by ui32TargetSubscript vector (row, col),
% using Givens col-wise rotations (i.e., rotating two columns). 
% -------------------------------------------------------------------------------------------------------------
%% INPUT
% dTargetMatrix       (:,:) double
% ui32TargetSubscript (1,2) uint32 % In Golub, 2014, kth row
% ui32AuxRowId        (1,1) uint32 = ui32TargetSubscript(1) - 1; % Default: element left target col, ith col
% -------------------------------------------------------------------------------------------------------------
%% OUTPUT
% % dTargetMatrix % Target matrix transformed in place

% -------------------------------------------------------------------------------------------------------------
%% CHANGELOG
% 25-02-2025    Pietro Califano     First version implemented
% -------------------------------------------------------------------------------------------------------------
%% DEPENDENCIES
% [-]
% -------------------------------------------------------------------------------------------------------------
ui32NumOfRows = size(dTargetMatrix, 1);

if coder.target("MATLAB") || coder.target("MEX") 

    % Assert sizes against indices
    assert(ui32AuxColId ~= ui32TargetSubscript(2))
    assert(ui32TargetSubscript(2) > 1, 'ERROR: Givens column-wise rotations cannot be applied to eliminate the 1st column. Target must be below the 1st pivot!')

end

% Get (a,b) entries
dVal1 = dTargetMatrix(ui32TargetSubscript(1),           ui32AuxColId);
dVal2 = dTargetMatrix(ui32TargetSubscript(1), ui32TargetSubscript(2));

% Compute Givens rotation values
[dCos, dSin] = GivensRotVals([dVal1; dVal2]);

% Transform rows
for ui32IdRow = 1:ui32NumOfRows
    
    dTmp1 = dTargetMatrix(ui32IdRow, ui32AuxColId); 
    dTmp2 = dTargetMatrix(ui32IdRow, ui32TargetSubscript(2));

    % Rotate dTargetMatrix, equivalent to Arot = G^T * A
    dTargetMatrix(ui32IdRow, ui32AuxColId) = dCos * dTmp1 - dSin * dTmp2;
    dTargetMatrix(ui32IdRow, ui32TargetSubscript(2)) = dSin * dTmp1 + dCos * dTmp2;

end

end
