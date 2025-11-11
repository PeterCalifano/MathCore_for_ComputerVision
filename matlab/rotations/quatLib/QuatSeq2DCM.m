function dDCM = QuatSeq2DCM(dQuatRot, bIS_VSRPplus) %#codegen
arguments
    dQuatRot (4,:) double {ismatrix, mustBeNumeric}
    bIS_VSRPplus (1,1) logical
end
%% PROTOTYPE
% dDCM = QuatSeq2DCM(dQuatRot, bIS_VSRPplus) %#codegen
% -------------------------------------------------------------------------------------------------------------
%% DESCRIPTION
% Function converting attitude quaternion to DCM. Convention flag set to 1
% (JPL convention) by default. 
% Multiple quaternions are automatically converted to a DCM sequence.
% REFERENCE:
% 1) Section 2.9.3, page 45, Fundamentals of Spacecraft Attitude 
% Determination and Control, Markley, Crassidis, 2014
% -------------------------------------------------------------------------------------------------------------
%% INPUT
% dQuatRot:     [4, Nseq] Input quaternion to convert to DCM. Nseq is
%                              the number of quaternions in the sequence.
% bIS_VSRPplus: [1] Boolean flag indicating the convention of the
%                     quaternion. 1: VSRPplus, 0: Hamilton.
% -------------------------------------------------------------------------------------------------------------
%% OUTPUT
% dDCM: [3, 3, Nseq] Output rotation matrix sequence
% -------------------------------------------------------------------------------------------------------------
%% CHANGELOG
% 16-12-2023    Pietro Califano     Coded from reference.
% 16-07-2025    Pietro Califano     Change expected input shape
% -------------------------------------------------------------------------------------------------------------
%% DEPENDENCIES
% [-]
% -------------------------------------------------------------------------------------------------------------
%% Future upgrades
% [-]
% -------------------------------------------------------------------------------------------------------------

%% Function code
% Output array initialization
ui32Nquats = size(dQuatRot, 2);
dDCM = coder.nullcopy(zeros(3, 3, ui32Nquats));

for idV = 1:ui32Nquats
    dDCM(:, :, idV) = Quat2DCM(dQuatRot(:, idV), bIS_VSRPplus);
end


end
