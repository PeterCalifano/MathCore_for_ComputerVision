function [dTargetMatrix, dOrthogonalQ] = GivensEliminateQR(dTargetMatrix, ...
                                                            ui16ValidRowPtr,...
                                                            ui16ValidColPtr, ...
                                                            bComputeQ) %#codegen
arguments
    dTargetMatrix      (:,:) double
    ui16ValidRowPtr   (1,1) uint16 {mustBeNumeric} = size(dTargetMatrix, 1);
    ui16ValidColPtr   (1,1) uint16 {mustBeNumeric} = size(dTargetMatrix, 2);
    bComputeQ         (1,1) logical {coder.mustBeConst} = false;
end
%% SIGNATURE
% [dTargetMatrix, dOrthogonalQ] = GivensEliminateQR(dTargetMatrix, ...
%                                                   ui16ValidRowPtr, ...
%                                                   ui16ValidColPtr, ...
%                                                   bComputeQ) %#codegen
% -------------------------------------------------------------------------------------------------------------
%% DESCRIPTION
% Function performing QR decomposition using row-wise Givens rotations over the active leading submatrix
% identified by ui16ValidRowPtr and ui16ValidColPtr. The returned matrix is the row-eliminated R factor and
% dOrthogonalQ satisfies dInputMatrix = dOrthogonalQ * dTargetMatrix within numerical precision when
% bComputeQ is true. bComputeQ must be a compile-time constant so generated code can elide Q accumulation.
% Requesting dOrthogonalQ as an output requires bComputeQ = true.
% -------------------------------------------------------------------------------------------------------------
%% INPUT
% dTargetMatrix         (:,:) double
% ui16ValidRowPtr       (1,1) uint16, last active row considered for elimination
% ui16ValidColPtr       (1,1) uint16, last active column considered for elimination
% bComputeQ             (1,1) logical, compile-time flag enabling Q accumulation
% -------------------------------------------------------------------------------------------------------------
%% OUTPUT
% dTargetMatrix         (:,:) double, R factor with eliminated entries in the active leading submatrix
% dOrthogonalQ          (:,:) double, orthogonal Q factor satisfying A = Q * R when bComputeQ is true,
%                       otherwise a zero matrix with compatible dimensions
% -------------------------------------------------------------------------------------------------------------
%% CHANGELOG
% 26-02-2025    Pietro Califano     First version implemented and validated.
% 26-04-2026    Pietro Califano     Add Q accumulation and align interface documentation.
% -------------------------------------------------------------------------------------------------------------
%% DEPENDENCIES
% ComputeGivensRotValues()
% ApplyGivensRot()
% -------------------------------------------------------------------------------------------------------------

if coder.target('MATLAB') || coder.target("MEX")
    assert(ui16ValidRowPtr <= size(dTargetMatrix, 1), ...
        'ERROR: ui16ValidRowPtr cannot exceed the number of rows in dTargetMatrix.')
    assert(ui16ValidColPtr <= size(dTargetMatrix, 2), ...
        'ERROR: ui16ValidColPtr cannot exceed the number of columns in dTargetMatrix.')
    assert(nargout <= 1 || bComputeQ, ...
        'ERROR: requesting dOrthogonalQ requires bComputeQ = true.')
end

%% Function code
if coder.const(bComputeQ)
    dOrthogonalQ = eye(size(dTargetMatrix, 1));
else
    dOrthogonalQ = zeros(size(dTargetMatrix, 1));
end

dNumOfElimCols = double(min(ui16ValidColPtr, ui16ValidRowPtr));
dValidRowPtr = double(ui16ValidRowPtr);

for idElimCol = 1:dNumOfElimCols

    if all(dTargetMatrix(:, idElimCol) == 0)
        continue;
    end

    for idElimRow = dValidRowPtr:-1:(idElimCol + 1)

        idAuxRow = idElimRow - 1;
        
        % Use "inlined" function instead of call (reduces copy of arrays)
        % dTargetMatrix = GivensEliminateRow(dTargetMatrix, ...
        %                                    ui16TargetIds, ...
        %                                    ui16AuxRowId);

        % Get (a,b) entries
        dVal1 = dTargetMatrix(idAuxRow,  idElimCol);
        dVal2 = dTargetMatrix(idElimRow, idElimCol);

        % Compute Givens rotation values
        [dCos, dSin] = ComputeGivensRotValues([dVal1; dVal2]);

        % Transform rows in place: R_next = G' * R.
        for idCol = 1:size(dTargetMatrix, 2)
            dTmp1 = dTargetMatrix(idAuxRow, idCol);
            dTmp2 = dTargetMatrix(idElimRow, idCol);

            [dTargetMatrix(idAuxRow, idCol), ...
                dTargetMatrix(idElimRow, idCol)] = ApplyGivensRot(dTmp1, dTmp2, dCos, dSin);
        end

        if coder.const(bComputeQ)
            % Accumulate the matching Q factor: Q_next = Q * G.
            for idRow = 1:size(dOrthogonalQ, 1)
                dTmp1 = dOrthogonalQ(idRow, idAuxRow);
                dTmp2 = dOrthogonalQ(idRow, idElimRow);

                [dOrthogonalQ(idRow, idAuxRow), ...
                    dOrthogonalQ(idRow, idElimRow)] = ApplyGivensRot(dTmp1, dTmp2, dCos, dSin);
            end
        end

    end % Slide across rows
end % Slide across columns

% Remove numerical zeros
for idRow = 1:size(dTargetMatrix, 1)
    for idCol = 1:size(dTargetMatrix, 2)
        if abs(dTargetMatrix(idRow, idCol)) < 0.001*eps
            dTargetMatrix(idRow, idCol) = 0.0;
        end
    end
end

end
