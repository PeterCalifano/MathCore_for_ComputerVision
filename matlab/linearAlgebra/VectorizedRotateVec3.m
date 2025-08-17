function dRotatedArrayVec3 = VectorizedRotateVec3(dArrayDCM, dArrayVec3ToRot, ui32ValidColsIdx)%#codegen
arguments (Input)
    dArrayDCM        (3,3,:) double 
    dArrayVec3ToRot  (3,:) double 
    ui32ValidColsIdx (1,:) uint32 = 1:size(dArrayVec3ToRot,2)
end
arguments (Output)
    dRotatedArrayVec3 (3,:) double 
end
%% SIGNATURE
% dRotatedArrayVec3 = VectorizedApplyDCM(dArrayDCM, dArrayVec3ToRot, ui32ValidColsIdx)%#codegen
% -------------------------------------------------------------------------------------------------------------
%% DESCRIPTION
% Function to rotate array of vec3 using vectorized operations (no for loop).
% -------------------------------------------------------------------------------------------------------------
%% INPUT
% dArrayDCM        (3,3,:) double
% dArrayVec3ToRot  (3,:) double
% % ui32ValidColsIdx (1,:) uint32 = 1:size(dArrayVec3ToRot,2)
% -------------------------------------------------------------------------------------------------------------
%% OUTPUT
% dRotatedArrayVec3 (3,:) double 
% -------------------------------------------------------------------------------------------------------------
%% CHANGELOG
% 17-08-2025    Pietro Califano     First implementation
% -------------------------------------------------------------------------------------------------------------
%% DEPENDENCIES
% [-]
% -------------------------------------------------------------------------------------------------------------

ui32Size3rd = size(dArrayDCM,3);
if coder.target('MATLAB') || coder.target("MEX")
    assert(size(dArrayVec3ToRot,2) <= ui32Size3rd, 'ERROR: size of array of vectors to rotate must match size of array of DCMs.')
    assert(length(ui32ValidColsIdx) <= ui32Size3rd, 'ERROR: array of indices cannot be longer than the number of entries.')
end

% Initialize output array
dRotatedArrayVec3 = zeros(3, length(ui32ValidColsIdx));
% Reshape vecs for pagemtimes: [3 x 1 x N]
dTmpToRot = reshape(dArrayVec3ToRot(:, ui32ValidColsIdx), 3, 1, []);
dRotatedArrayVec3(:,:) = pagemtimes(dArrayDCM(:,:,ui32ValidColsIdx), dTmpToRot);
end
