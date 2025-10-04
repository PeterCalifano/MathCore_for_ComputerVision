function [dRotatedArrayVec3, bSizeErrorFlag] = VectorizedRotateVec3(dArrayDCM, dArrayVec3ToRot, ui32ValidColsIdx)%#codegen
arguments (Input)
    dArrayDCM        (3,3,:) double 
    dArrayVec3ToRot  (3,:) double 
    ui32ValidColsIdx (1,:) uint32 = 1:size(dArrayVec3ToRot,2)
end
arguments (Output)
    dRotatedArrayVec3 (3,:) double 
    bSizeErrorFlag    (1,1) logical
end
%% SIGNATURE
% [dRotatedArrayVec3, bSizeErrorFlag] = VectorizedApplyDCM(dArrayDCM, dArrayVec3ToRot, ui32ValidColsIdx)%#codegen
% -------------------------------------------------------------------------------------------------------------
%% DESCRIPTION
% Function to rotate array of vec3 using vectorized operations (no for loop).
% -------------------------------------------------------------------------------------------------------------
%% INPUT
% dArrayDCM             (3,3,:) double
% dArrayVec3ToRot       (3,:) double
% % ui32ValidColsIdx    (1,:) uint32 = 1:size(dArrayVec3ToRot,2)
% -------------------------------------------------------------------------------------------------------------
%% OUTPUT
% dRotatedArrayVec3 (3,:) double 
% bSizeErrorFlag    (1,1) logical
% -------------------------------------------------------------------------------------------------------------
%% CHANGELOG
% 17-08-2025    Pietro Califano     First implementation
% 04-10-2025    Pietro Califano     Fix implementation bug related to index using valid indices and assert;
%                                   extend implementation to handle different use cases (inputs)
% -------------------------------------------------------------------------------------------------------------
%% DEPENDENCIES
% [-]
% -------------------------------------------------------------------------------------------------------------

ui32Size3rd = size(dArrayDCM,3);
bSizeErrorFlag = false;

if coder.target('MATLAB') || coder.target("MEX")
    assert( size(dArrayVec3ToRot,2) <= max(ui32ValidColsIdx), 'ERROR: size of array of vectors cannot be smaller than max index.')
    assert( length(ui32ValidColsIdx) <= ui32Size3rd, 'ERROR: array of indices cannot be longer than the number of entries.')
    assert( ui32Size3rd == length(ui32ValidColsIdx) || ui32Size3rd == size(dArrayVec3ToRot,2), ...
        'ERROR: size mismatch. dArrayDCM must have 3rd size equal to number of valid indices or matching dArrayVec3ToRot size.')
end

% Initialize output array
dRotatedArrayVec3 = zeros(3, length(ui32ValidColsIdx));
% Reshape vecs for pagemtimes: [3 x 1 x N]
dTmpToRot = reshape(dArrayVec3ToRot(:, ui32ValidColsIdx), 3, 1, []);

if ui32Size3rd == length(ui32ValidColsIdx)
    % Use all DCM pages
    dRotatedArrayVec3(:,:) = pagemtimes(dArrayDCM, dTmpToRot);
    return
elseif ui32Size3rd == size(dArrayVec3ToRot,2)
    % Index DCM pages
    dRotatedArrayVec3(:,:) = pagemtimes(dArrayDCM, dTmpToRot);
    return
end

bSizeErrorFlag = true;
end
