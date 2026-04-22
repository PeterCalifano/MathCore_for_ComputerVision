function dInvQuatSeq = qInvert(dQuatSeq, bIS_VSRPplus) %#codegen
arguments
    dQuatSeq     (4,:)
    bIS_VSRPplus (1,1) logical {islogical, isscalar} = true
end
%% PROTOTYPE
% dInvQuatSeq = qInvert(dQuatSeq, bIS_VSRPplus)
% -------------------------------------------------------------------------------------------------------------
%% DESCRIPTION
% Computes the inverse according to VSRP+ convention (right-handed). Set bIS_VSRPplus = false if SVRP+ 
% convention with scalar first is being used. 
% (SV) Scalar first, Vector last
% (P) Passive 
% (R) Successive coordinate transformations have the unmodified quaternion chain on the Right side of
%     the triple product.
% (plus) Right-Handed Rule for the imaginary numbers i, j, k. (aka Hamilton)
% REFERENCE:
% 1) Fundamentals of Spacecraft Attitude Determination and Control, Markley
% Crassidis, 2014. Section 2.
% -------------------------------------------------------------------------------------------------------------
%% INPUT
% dQuatSeq:     [4, N] Sequence of attitude quaternion with convention 
%                        determined by bConvFlag. 
% bIS_VSRPplus:   [1]  Boolean flag indicating the convention of the
%                       quaternion. 1: VSRPplus, 0: SVRPplus (as MATLAB)
% -------------------------------------------------------------------------------------------------------------
%% OUTPUT
% dInvQuatSeq: [4, N] Sequence of inverted attitude quaternions
% -------------------------------------------------------------------------------------------------------------
%% CHANGELOG
% 25-11-2023    Pietro Califano   Function coded.
% 14-12-2023    Pietro Califano   Change of convention flag for more
%                                 self-explainable code.
% -------------------------------------------------------------------------------------------------------------
%% DEPENDENCIES
% [-]
% -------------------------------------------------------------------------------------------------------------
%% Function code

if bIS_VSRPplus == true
    dInvQuatSeq = [-dQuatSeq(1:3, :); dQuatSeq(4, :)];
else
    dInvQuatSeq = [dQuatSeq(1, :); -dQuatSeq(2:4, :)];
end

end
