function [dTargetMatrix, dOrthogonalQ] = GivensEliminateQR(dTargetMatrix, ...
                                                            ui16ValidRowPtr,...
                                                            ui16ValidColPtr) %#codegen
arguments
    dTargetMatrix      
    ui16ValidRowPtr   (1,1) uint16 {isscalar, mustBeNumeric} = size(dTargetMatrix, 1);
    ui16ValidColPtr   (1,1) uint16 {isscalar, mustBeNumeric} = size(dTargetMatrix, 2);
end
%% SIGNATURE
% [dTargetMatrix, dOrthogonalQ] = GivensEliminateQR(dTargetMatrix, ...
%                                                   bEliminateInPlace) %#codegen
% -------------------------------------------------------------------------------------------------------------
%% DESCRIPTION
% Function performing QR decomposition using Givens Rotations.
% -------------------------------------------------------------------------------------------------------------
%% INPUT
% dTargetMatrix         (:,:) double
% bEliminateInPlace     (1,1) logical
% -------------------------------------------------------------------------------------------------------------
%% OUTPUT
% dTargetMatrix % Target matrix transformed in place
% dOrthogonalQ TODO
% -------------------------------------------------------------------------------------------------------------
%% CHANGELOG
% 26-02-2025    Pietro Califano     First version implemented and validated.
% -------------------------------------------------------------------------------------------------------------
%% DEPENDENCIES
% [-]
% -------------------------------------------------------------------------------------------------------------

if coder.target('MATLAB') || coder.target("MEX")
    % TODO add asserts

end

%% Function code

% TODO add capability of computing Q = G1 * ... * Gt
dOrthogonalQ = 0; 

% NOTE: these branches are intended as separate functions, thus self-contained. The coder should generate
% two versions according to the flag.
ui16NumOfCols = uint16(ui16ValidColPtr);

for ui16IdElimCol = uint16(1:ui16NumOfCols)

    if all(dTargetMatrix(:, ui16IdElimCol) == 0)
        % ui16LastNotZero
        % FIXME, this function does not support matrices with columns all zero. How are they handled?
        continue;
    end

    for ui16IdElimRow = uint16(ui16ValidRowPtr:-1:ui16IdElimCol+1)

        % Run GivensEliminateRow
        ui16TargetSubscript = [ui16IdElimRow, ui16IdElimCol];
        ui16AuxRowId        = ui16TargetSubscript(1) - 1;
        
        % Use "inlined" function instead of call (reduces copy of arrays)
        % dTargetMatrix = GivensEliminateRow(dTargetMatrix, ...
        %                                    ui16TargetIds, ...
        %                                    ui16AuxRowId);

        % Get (a,b) entries
        dVal1 = dTargetMatrix(ui16AuxRowId,           ui16TargetSubscript(2));
        dVal2 = dTargetMatrix(ui16TargetSubscript(1), ui16TargetSubscript(2));

        % Compute Givens rotation values
        [dCos, dSin] = ComputeGivensRotValues([dVal1; dVal2]);

        % Transform rows in place
        for ui16IdCol = 1:size(dTargetMatrix, 2)

            % if(all(dTargetMatrix(:, ui16IdCol) == 0))
            %     continue;
            % end

            dTmp1 = dTargetMatrix(ui16AuxRowId, ui16IdCol);
            dTmp2 = dTargetMatrix(ui16TargetSubscript(1), ui16IdCol);

            % Rotate dTargetMatrix, equivalent to Arot = G^T * A
            [dTargetMatrix(ui16AuxRowId, ui16IdCol), ...
                dTargetMatrix(ui16TargetSubscript(1), ui16IdCol)] = ApplyGivensRot(dTmp1, dTmp2, dCos, dSin);

        end

    end % Slide across rows
end % Slide across columns

% Remove numerical zeros
dTargetMatrix(abs(dTargetMatrix) < 1.5*eps ) = 0.0;

% TODO add capability of computing Q incrementally
return


end
