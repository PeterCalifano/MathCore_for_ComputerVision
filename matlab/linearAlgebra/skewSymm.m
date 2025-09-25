function [dSkewSymmMatrix] = skewSymm(dVector) %#codegen
arguments
    dVector (3,1) double {isvector}
end
%% PROTOTYPE
% [dSkewSymmMatrix] = skewSymm(dVector)
% -------------------------------------------------------------------------------------------------------------
%% DESCRIPTION
% Function expanding a vector into a skew symmetric matrix.
% -------------------------------------------------------------------------------------------------------------
%% INPUT
% dVector
% -------------------------------------------------------------------------------------------------------------
%% OUTPUT
% dSkewSymmMatrix
% -------------------------------------------------------------------------------------------------------------
%% CHANGELOG
% 21-09-2024        Pietro Califano         Function coded for 3x3 skew symmetric matrices 
% -------------------------------------------------------------------------------------------------------------
%% DEPENDENCIES
% [-]
% -------------------------------------------------------------------------------------------------------------

%% Function code

dSkewSymmMatrix = [0,          -dVector(3), dVector(2);
                   dVector(3), 0,          -dVector(1);
                  -dVector(2), dVector(1),  0];


end

