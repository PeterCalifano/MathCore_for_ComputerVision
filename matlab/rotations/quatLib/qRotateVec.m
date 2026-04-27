function  o_vector_Ref2 = qRotateVec(i_dqRef1wrtRef2, i_vector_Ref1, i_bIS_VSRPplus) %#codegen
%% PROTOTYPE
% o_vector_Ref2 = qRotateVec(i_dqRef1wrtRef2, i_vector_Ref1, i_bIS_VSRPplus) %#codegen
% -------------------------------------------------------------------------------------------------------------
%% DESCRIPTION
% Function rotate a generic vector "i_vector_Ref1" to "o_vector_Ref2" through 
% the rotation of i_dqRef1wrtRef2 using quaternion operations. Quaternion
% convention specified by "i_bIS_VSRPplus".
% VSRPplus convention: scalar last; SVRPplus convention: scalar first.
% -------------------------------------------------------------------------------------------------------------
%% INPUT
% i_dqRef1wrtRef2
% i_vector_Ref1
% i_bIS_VSRPplus
% -------------------------------------------------------------------------------------------------------------
%% OUTPUT
% o_vector_Ref2
% -------------------------------------------------------------------------------------------------------------
%% CHANGELOG
% 13-12-2023        Pietro Califano     Function coded. Not validated. 
% -------------------------------------------------------------------------------------------------------------
%% DEPENDENCIES
% [-]
% -------------------------------------------------------------------------------------------------------------
%% Function code
% Initialize output array
o_vector_Ref2 = coder.nullcopy(zeros(3, 1));

% Compute rotation
if i_bIS_VSRPplus == true
    tmpVec = qCross( qCross( i_dqRef1wrtRef2, [i_vector_Ref1; 0] ), [-i_dqRef1wrtRef2(1:3); i_dqSCBwrtIN(4)] );
elseif i_bIS_VSRPplus == false
    tmpVec = qCross( qCross( i_dqRef1wrtRef2, [i_vector_Ref1; 0] ), [i_dqSCBwrtIN(1); -i_dqRef1wrtRef2(2:4)] );
end

% Assign output
o_vector_Ref2(1:3) = tmpVec(1:3);

%% LOCAL FUNCTION
    function q1xq2 = qCross(q1, q2) %#codegen
        %% PROTOTYPE
        % q1xq2 = qCross(q1, q2)
        % -------------------------------------------------------------------------------------------------------------
        %% DESCRIPTION
        % Computes the quaternion product between q1 and q2 using matrix Psi(q1)
        % operation: [q1 x] = [Psi(q) q]. Quaternion convention: JPL (qv, qs)
        % REFERENCE:
        % 1) Fundamentals of Spacecraft Attitude Determination and Control, Markley
        % Crassidis, 2014. Section 2, page 53.
        % -------------------------------------------------------------------------------------------------------------
        %% INPUT
        % q1: [4x1] Quaternion 1 (left side of the quat. product)
        % q2: [4x1] Quaternion 2 (right side of the quat. product)
        % -------------------------------------------------------------------------------------------------------------
        %% OUTPUT
        % q1xq2: [4x1] Quaternion cross product between q1 and q2
        % -------------------------------------------------------------------------------------------------------------
        %% CHANGELOG
        % 05-09-2023    Pietro Califano   Coded and validated for std quat. lib.
        % -------------------------------------------------------------------------------------------------------------
        %% DEPENDENCIES
        % [-]
        % -------------------------------------------------------------------------------------------------------------

        q1xq2 = [q1(4) q1(3) -q1(2) q1(1); ...
            -q1(3) q1(4) q1(1) q1(2); ...
            q1(2) -q1(1) q1(4) q1(3); ...
            -q1(1) -q1(2) -q1(3) q1(4)] * q2;

    end
end
