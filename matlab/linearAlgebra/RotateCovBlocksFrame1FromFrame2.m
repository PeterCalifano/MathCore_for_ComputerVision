function [dCovMatrix_Frame1] = RotateCovBlocksFrame1FromFrame2(dCovMatrix_Frame2, ...
                                                              dRot3_Frame1FromFrame2, ...
                                                              bRotateBlockMask)%#codegen
arguments
    dCovMatrix_Frame2      (:,:) double
    dRot3_Frame1FromFrame2 (3,3) double
    bRotateBlockMask       (1,:) logical = true(1,2)
end
%% SIGNATURE
%[dCovMatrix_Frame1] = RotateCovBlocksFrame1FromFrame2(dCovMatrix_Frame2, ...
%                                                      dRot3_Frame1FromFrame2, ...
%                                                      bRotateBlockMask)%#codegen
% -------------------------------------------------------------------------------------------------------------
%% DESCRIPTION
% Function transforming a covariance matrix composed of 3x3 blocks from Frame2 to Frame1
% using the provided rotation matrix. Only the blocks marked as true in the logical mask are rotated.
% Cross-terms between rotated and non-rotated blocks are handled accordingly.
% -------------------------------------------------------------------------------------------------------------
%% INPUT
% dCovMatrix_Frame2      (:,:) double
% dRot3_Frame1FromFrame2 (3,3) double
% bRotateBlockMask       (1,:) logical = true(1,2)
% -------------------------------------------------------------------------------------------------------------
%% OUTPUT
% dCovMatrix_Frame1     (:,:) double 
% -------------------------------------------------------------------------------------------------------------
%% CHANGELOG
% 10-01-2026    Pietro Califano     First implementation of utility function
% -------------------------------------------------------------------------------------------------------------
%% DEPENDENCIES
% RotateCovBlock.m
% -------------------------------------------------------------------------------------------------------------

%% Function code
% Initialize output as input
dCovMatrix_Frame1 = dCovMatrix_Frame2;

% Assert sizes (must be 3N x 3N)
ui32NumBlocks = uint32(round(size(dCovMatrix_Frame2, 1) / 3));

assert(isequal(size(dCovMatrix_Frame2, 1), size(dCovMatrix_Frame2, 2)), ...
       "RotateCovBlocksFrame1FromFrame2:InvalidInput", ...
       "Input covariance matrix must be square.");

assert(isequal(mod(size(dCovMatrix_Frame2, 1), 3), 0), ...
       "RotateCovBlocksFrame1FromFrame2:InvalidInput", ...
       "Input covariance matrix size must be multiple of 3.");

% Assert sizes depending on index
assert(isequal(length(bRotateBlockMask), ui32NumBlocks), ...
       "RotateCovBlocksFrame1FromFrame2:InvalidInput", ...
       "Block rotation mask length must match number of 3x3 blocks in covariance matrix.");

% Rotate all sub-blocks
for k = 0:(ui32NumBlocks - 1) % Loop over rows

    % Determine rows index
    ui32RowIndex = k * 3 + 1 : (k * 3 + 3);

    for j = 0:(ui32NumBlocks - 1) % Loop over columns
        
        % Determine columns index
        ui32ColIndex = j * 3 + 1 : (j * 3 + 3);
        
        % Full rotate only if both blocks are marked for rotation
        if bRotateBlockMask(k + 1) && bRotateBlockMask(j + 1)
            dCovMatrix_Frame1( ui32RowIndex, ui32ColIndex ) = RotateCovBlock(dRot3_Frame1FromFrame2, ... ...
                                                                dCovMatrix_Frame2( ui32RowIndex, ui32ColIndex ) );
        
        elseif bRotateBlockMask(k + 1) && ~bRotateBlockMask(j + 1)
            % Only rotate rows (left-multiply)
            dCovMatrix_Frame1( ui32RowIndex, ui32ColIndex ) = dRot3_Frame1FromFrame2 * ...
                                                                dCovMatrix_Frame2( ui32RowIndex, ui32ColIndex );

        elseif ~bRotateBlockMask(k + 1) && bRotateBlockMask(j + 1)
            % Only rotate columns (right-multiply)
            dCovMatrix_Frame1( ui32RowIndex, ui32ColIndex ) = dCovMatrix_Frame2( ui32RowIndex, ui32ColIndex ) * ...
                                                                transpose(dRot3_Frame1FromFrame2);
                                                                
        end
    end
end


end