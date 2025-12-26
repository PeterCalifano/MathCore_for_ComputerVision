function dQuatRot = DCM2quatSeq(dDCM, bIS_VSRPplus) %#codegen
arguments (Input)
    dDCM          (3,3,:) {mustBeNumeric} 
    bIS_VSRPplus  (1,1) {islogical, isscalar} = true
end
arguments (Output)
    dQuatRot (4,:)
end
%% PROTOTYPE
% dQuatRot = DCM2quat(dDCM, bIS_VSRPplus)%#codegen
% -------------------------------------------------------------------------------------------------------------
%% DESCRIPTION
% Function converting a DCM to Attitude quaternion. Conversion occurs
% according to JPL convention. Set bIS_JPL_CONV = false if Hamilton
% convention is being used. The conversion is made numerically optimal.
% Multiple DCM are automatically converted to a quaternion sequence.
% REFERENCE:
% 1) Section 2.9.3, page 48, Fundamentals of Spacecraft Attitude 
% Determination and Control, Markley, Crassidis, 2014
% -------------------------------------------------------------------------------------------------------------
%% INPUT
% dDCM:         [3, 3, Nseq]  Input rotation matrix to convert. Nseq is
%                              the number of quaternions in the sequence. 
% bIS_VSRPplus: [1]  Boolean flag indicating the convention of the
%                      quaternion. 1: VSRPplus, 0: SVRPplus (as MATLAB)
% -------------------------------------------------------------------------------------------------------------
%% OUTPUT
% dQuatRot:     [4, Nseq]  Output quaternion sequence
% -------------------------------------------------------------------------------------------------------------
%% CHANGELOG
% 16-12-2023    Pietro Califano     Coded from reference. Validated against
%                                   MATLAB library.
% -------------------------------------------------------------------------------------------------------------
%% DEPENDENCIES
% [-]
% -------------------------------------------------------------------------------------------------------------
%% Future upgrades
% 1) Vectorization for multiple quaternions in one call (improvement using
% vectorization over third dimension).
% -------------------------------------------------------------------------------------------------------------
%% Function code
% Quaternion output initialization
Nquats = size(dDCM, 3);
dQuatRot = coder.nullcopy(zeros(4, Nquats));

for idV = 1:Nquats
    % Compute trace of DCM
    trDCM = trace(dDCM(:, :, idV));

    [~, idMax] = max([dDCM(1,1, idV), dDCM(2,2, idV), dDCM(3,3, idV), trDCM]);

    % Adjust computation case for numerical accuracy to convert DCM to NON
    % UNITARY QUATERNION
    switch idMax
        case 4
            % Scalar part is max
            dQuatRot(:, idV) = [dDCM(2, 3, idV) - dDCM(3, 2, idV);
                dDCM(3, 1, idV) - dDCM(1, 3, idV);
                dDCM(1, 2, idV) - dDCM(2, 1, idV);
                1 + trDCM];
        case 3
            % Third vector component is max
            dQuatRot(:, idV) = [dDCM(3, 1, idV) + dDCM(1, 3, idV);
                                  dDCM(3, 2, idV) + dDCM(2, 3, idV);
                                  1 + 2*dDCM(3, 3, idV) - trDCM;
                                  dDCM(1, 2, idV) - dDCM(2, 1, idV)];
                  
        case 2
            % Second vector component is max
            dQuatRot(:, idV) = [dDCM(2, 1, idV) + dDCM(1, 2, idV);
                1 + 2*dDCM(2, 2, idV) - trDCM;
                dDCM(2, 3, idV) + dDCM(3, 2, idV);
                dDCM(3, 1, idV) - dDCM(1, 3, idV)];

        case 1
            % First vector component is max
            dQuatRot(:, idV) = [1 + 2*dDCM(1, 1, idV) - trDCM;
                dDCM(1, 2, idV) + dDCM(2, 1, idV);
                dDCM(1, 3, idV) + dDCM(3, 1, idV);
                dDCM(2, 3, idV) - dDCM(3, 2, idV)];
    end

end

% Normalize quaternion
dQuatRot = dQuatRot./vecnorm(dQuatRot, 2, 1);

% Assign quaternion components depending on quaternion convention
if bIS_VSRPplus == false
    % SVRPplus convention (false) --> SWAP
    dQuatRot = [dQuatRot(4, :);
                  dQuatRot(1, :);
                  dQuatRot(2, :);
                  dQuatRot(3, :)];

elseif bIS_VSRPplus == true
    % VSRPplus convention (false). DO NOTHING: DEFAULT
    %     dQuatRot(1) = qv1;
    %     dQuatRot(2) = qv2;
    %     dQuatRot(3) = qv3;
    %     dQuatRot(4) = qs ;
end

end


