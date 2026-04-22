function [A_voigt] = Matrix2Voigt(A, mode)
%% PROTOTYPE
% [A_voigt] = Matrix2Voigt(A, mode)
% -------------------------------------------------------------------------------------------------------------
%% DESCRIPTION
% Converts the matrix representation of a tensor into the corresponding
% Voigt notation. Two modes of operations are possible: 1) for stress and
% deformation tensors using the common Voigt convention; 2) generic
% stacking of terms for tensors
% -------------------------------------------------------------------------------------------------------------
%% INPUT
% A: [MxN] matrix where M and N are the number of rows 
% mode: [0 or 1] select the type of operation. Default: 0, for stress and
%       deformation tensors in Structural applications
% -------------------------------------------------------------------------------------------------------------
%% OUTPUT
% A_voigt: [(M*N)x1] array representing the Voigt notation of the matrix A
% -------------------------------------------------------------------------------------------------------------
%% CONTRIBUTORS
% Pietro Califano
% -------------------------------------------------------------------------------------------------------------
%% CHANGELOG
% V1: created with 2 modes of operations, does not apply correction
%     for the deformation tensor 15/01/2022
% -------------------------------------------------------------------------------------------------------------
%% DEPENDENCIES
% None
% -------------------------------------------------------------------------------------------------------------
%% Function code
if ~exist('mode', 'var')
    mode = 0; % default mode
end

switch mode
    case 0
        A_voigt = [A(1, 1), A(2, 2), A(3, 3), A(2, 3), A(1, 3), A(1, 2)]';
    case 1
        % Stacking terms by terms moving along each row
        n_elements = length(A(1, :)) * length(A(:, 1));
        A_voigt = reshape(A', n_elements, 1);
end

  

end