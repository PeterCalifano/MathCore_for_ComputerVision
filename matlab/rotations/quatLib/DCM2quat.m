function dQuatRot = DCM2quat(dDCM, bIS_VSRPplus) %#codegen
arguments
    dDCM         (3,3) double
    bIS_VSRPplus (1,1) logical = false
end
%% PROTOTYPE
% dQuatRot = DCM2quat(dDCM, bIS_VSRPplus)%#codegen
% -------------------------------------------------------------------------------------------------------------
%% DESCRIPTION
% Function converting a DCM to Attitude quaternion. Conversion occurs according to VSRP+ convention
% (right-handed) Set bIS_VSRPplus = false if SVRP+ convention with scalar first is being used. 
% The conversion is made numerically optimal. Only supports one conversion per call.
% (SV) Scalar first, Vector last
% (P) Passive 
% (R) Successive coordinate transformations have the unmodified quaternion chain on the Right side of
%     the triple product.
% (plus) Right-Handed Rule for the imaginary numbers i, j, k. (aka Hamilton)
% REFERENCE:
% 1) Section 2.9.3, page 48, Fundamentals of Spacecraft Attitude
% Determination and Control, Markley, Crassidis, 2014
% -------------------------------------------------------------------------------------------------------------
%% INPUT
% dDCM:         [3, 3]  Input rotation matrix to convert
% bIS_VSRPplus: [1]  Boolean flag indicating the convention of the
%                      quaternion. 1: VSRPplus, 0: SVRPplus (as MATLAB)
% -------------------------------------------------------------------------------------------------------------
%% OUTPUT
% dQuatRot:     [4, 1]  Output quaternion
% -------------------------------------------------------------------------------------------------------------
%% CHANGELOG
% 16-12-2023    Pietro Califano     Coded from reference. Validated against
%                                   MATLAB functions.
% 27-04-2024    Pietro Califano     Modified convention of quaternion.
% 25-09-2025    Pietro Califano     Minor revision (update validations of args)
% -------------------------------------------------------------------------------------------------------------
%% DEPENDENCIES
% [-]
% -------------------------------------------------------------------------------------------------------------

%% Function code
% INPUT ASSERT CHECKS
if coder.target('MATLAB') || coder.target('MEX')
    assert( all(abs(dDCM) <= 1.0 + eps, 'all'), 'ERROR: invalid DCM. Must have a unitary norm!')
end

% Quaternion output initialization
dQuatRot = coder.nullcopy(zeros(4, 1));

% Compute trace of DCM
trDCM = trace(dDCM);

[~, idMax] = max([dDCM(1,1), dDCM(2,2), dDCM(3,3), trDCM]);

% Adjust computation case for numerical accuracy to convert DCM to NON
% UNITARY QUATERNION
switch idMax
    case 4
        % Scalar part is max
        dQuatRot = [dDCM(2, 3) - dDCM(3, 2);
            dDCM(3, 1) - dDCM(1, 3);
            dDCM(1, 2) - dDCM(2, 1);
            1 + trDCM];
    case 3
        % Third vector component is max
        dQuatRot = [dDCM(3, 1) + dDCM(1, 3);
            dDCM(3, 2) + dDCM(2, 3);
            1 + 2*dDCM(3, 3) - trDCM;
            dDCM(1, 2) - dDCM(2, 1)];

    case 2
        % Second vector component is max
        dQuatRot = [dDCM(2, 1) + dDCM(1, 2);
            1 + 2*dDCM(2, 2) - trDCM;
            dDCM(2, 3) + dDCM(3, 2);
            dDCM(3, 1) - dDCM(1, 3)];

    case 1
        % First vector component is max
        dQuatRot = [1 + 2*dDCM(1, 1) - trDCM;
            dDCM(1, 2) + dDCM(2, 1);
            dDCM(1, 3) + dDCM(3, 1);
            dDCM(2, 3) - dDCM(3, 2)];
end

% Normalize quaternion
dQuatRot = dQuatRot./norm(dQuatRot);

% Assign quaternion components depending on quaternion convention
if bIS_VSRPplus == false
    % SVRPplus convention scalar first (SWAP COMPONENTS)
    dQuatRot = [dQuatRot(4); dQuatRot(1:3)];

elseif bIS_VSRPplus == true
    % VSRPplus convention (false). DO NOTHING: DEFAULT
    %     dQuatRot(1) = qv1;
    %     dQuatRot(2) = qv2;
    %     dQuatRot(3) = qv3;
    %     dQuatRot(4) = qs ;
end

if coder.target('MATLAB') || coder.target('MEX')
    assert( all(abs(dQuatRot) <= 1.0 + eps, 'all'), 'ERROR: invalid quaternion array. Must have a unitary norm!')
end

end


