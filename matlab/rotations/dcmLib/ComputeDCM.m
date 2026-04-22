function dRotDCM = ComputeDCM(ui8AxisIDseq, dRotAngleSeq) %#codegen
%% PROTOTYPE
% dRotDCM = ComputeDCM(ui8AxisIDseq, dRotAngleSeq)
% -------------------------------------------------------------------------------------------------------------
%% DESCRIPTION
% Function computing a rotation matrix by concatenation of the selected
% single axis rotation sequence of N entries (arbitrary number).
% Angles must be in radians. Rotation order is "as written": first in the
% sequence is the last rotation to be perfomed.
% NOTE: Rotation matrix assumed Right-Hand Rule, positive clockwise in the
%       direction of the axis.
% -------------------------------------------------------------------------------------------------------------
%% INPUT
% ui8AxisIDseq: [1, N] Array sequence of axis IDs (1, 2 or 3)
% dRotAngleSeq: [1, N] Array sequence of corresponding rotation angles
% -------------------------------------------------------------------------------------------------------------
%% OUTPUT
% dRotDCM: [3, 3] Output rotation matrix
% -------------------------------------------------------------------------------------------------------------
%% CHANGELOG
% 27-11-2023    Pietro Califano     Function assembled from validated code.
% -------------------------------------------------------------------------------------------------------------
%% DEPENDENCIES
% [-]
% -------------------------------------------------------------------------------------------------------------
%% Function code
% Get number of rotations
Nr = length(dRotAngleSeq);

% Input checks
assert(Nr == length(ui8AxisIDseq), 'Input arrays length do not match.')
assert( sum(ui8AxisIDseq >= 1 & ui8AxisIDseq <= 3) == Nr , 'Invalid axis selected. Input must be: 1, 2 or 3.');

% Initialize DCM with first rotation (left-most)
dRotDCM = SingleAxisDCM(ui8AxisIDseq(end), dRotAngleSeq(end));

for idR = (Nr-1):-1:1
    % Concatenate with with next rotation matrix
    dRotDCM = SingleAxisDCM(ui8AxisIDseq(idR), dRotAngleSeq(idR)) * dRotDCM;
end

%% LOCAL FUNCTIONS
    function dSingleRotDCM = SingleAxisDCM(ui8AxisID, dRotAngle) %#codegen
        %% PROTOTYPE
        % dSingleRotDCM = SingleAxisDCM(ui8AxisID, dRotAngle)
        % -------------------------------------------------------------------------------------------------------------
        %% DESCRIPTION
        % Function computing the rotation matrix about the selected axis (1,2,3) =
        % (X,Y,Z) of the angle dRotAngle in RADIANS. Zero is returned if invalid
        % axis input is selected and a warning is issued.
        % NOTE: Rotation matrix assumed Right-Hand Rule, positive clockwise in the
        %       direction of the axis.
        % -------------------------------------------------------------------------------------------------------------
        %% INPUT
        % ui8AxisID: [1] Axis ID (1,2,3) = (X,Y,Z)
        % dRotAngle: [1] Rotation angle in RADIANS
        % -------------------------------------------------------------------------------------------------------------
        %% OUTPUT
        % dSingleRotDCM: [3, 3] DCM rotating around "AxisID" of "RotAngle" [rad]
        % -------------------------------------------------------------------------------------------------------------
        %% CHANGELOG
        % 27-11-2023    Pietro Califano     Function assembled from
        %                                   previous validated codes.
        % -------------------------------------------------------------------------------------------------------------
        %% DEPENDENCIES
        % [-]
        % -------------------------------------------------------------------------------------------------------------
        %% Function code
        switch ui8AxisID
            case 1
                % Axis 1: X
                dSingleRotDCM = [1, 0, 0; ...
                    0, cos(dRotAngle), sin(dRotAngle); ...
                    0, -sin(dRotAngle), cos(dRotAngle)];
            case 2
                % Axis 2: Y
                dSingleRotDCM = [cos(dRotAngle), 0, -sin(dRotAngle);...
                    0, 1, 0;...
                    sin(dRotAngle), 0, cos(dRotAngle)];
            case 3
                % Axis 3: Z
                dSingleRotDCM = [cos(dRotAngle), sin(dRotAngle), 0;...
                    -sin(dRotAngle), cos(dRotAngle), 0;...
                    0, 0, 1];
            otherwise
                disp('Invalid rotation axis ID. Output set to zero!')
                dSingleRotDCM = zeros(3, 3);
        end
    end

end
