function dRotatedArrayMat3 = VectorizedRotateMat3(dArrayDCM, dArrayMat3ToRot, ui32ValidColsIdx)%#codegen
arguments (Input)
    dArrayDCM        (3,3,:) double 
    dArrayMat3ToRot  (3,:) double 
    ui32ValidColsIdx (1,:) uint32 = 1:size(dArrayDCM,3)
end
arguments (Output)
    dRotatedArrayMat3 (3,:) double 
end
%% SIGNATURE
% dRotatedArrayMat3 = VectorizedRotateMat3(dArrayDCM, dArrayMat3ToRot, ui32ValidColsIdx)%#codegen
% -------------------------------------------------------------------------------------------------------------
%% DESCRIPTION
% Function to rotate array of 3x3 matrices using vectorized operations (no for loop).
% -------------------------------------------------------------------------------------------------------------
%% INPUT
% dArrayDCM        (3,3,:) double
% dArrayMat3ToRot  (3,3,:) double
% ui32ValidColsIdx (1,:) uint32 = 1:size(dArrayVec3ToRot,2)
% -------------------------------------------------------------------------------------------------------------
%% OUTPUT
% dRotatedArrayMat3 (3,:) double 
% -------------------------------------------------------------------------------------------------------------
%% CHANGELOG
% 17-08-2025    Pietro Califano     First implementation
% -------------------------------------------------------------------------------------------------------------
%% DEPENDENCIES
% [-]
% -------------------------------------------------------------------------------------------------------------

ui32Size3rd = size(dArrayDCM,3);
if coder.target('MATLAB') || coder.target("MEX")
    assert(size(dArrayMat3ToRot,3) <= ui32Size3rd, 'ERROR: size of array of matrices to rotate must match size of array of DCMs.')
    assert(length(ui32ValidColsIdx) <= ui32Size3rd, 'ERROR: array of indices cannot be longer than the number of entries.')
end

% Initialize output array
dRotatedArrayMat3 = zeros(3, 3, length(ui32ValidColsIdx));

% Apply on left-side: R <- D*M
dRotatedArrayMat3(:,:,:) = pagemtimes( dArrayDCM(:,:,ui32ValidColsIdx), dArrayMat3ToRot(:, ui32ValidColsIdx) );
% Apply on right-side: R <- R*D^T
dRotatedArrayMat3(:,:,:) = pagemtimes( dArrayMat3ToRot(:, ui32ValidColsIdx), pagetranspose(dArrayDCM(:,:,ui32ValidColsIdx)) );
end
