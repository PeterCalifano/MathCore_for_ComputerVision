function [dDataMatrix, bIsSignSwitched, ui8howManySwitches, ...
                bSignSwitchDetectionMask] = fixQuatSignDiscontinuity(dQuat_fromAtoB) %#codegen
arguments
    dQuat_fromAtoB (:, 4) double {mustBeNumeric}
end
%% PROTOTYPE
% [dDataMatrix, bIsSignSwitched, ui8howManySwitches, bsignSwitchDetectionMask] =...
%                                                fixQuatSignDiscontinuity(dQuat_fromAtoB) %#codegen
% -------------------------------------------------------------------------------------------------------------
%% DESCRIPTION
% What the function does
% -------------------------------------------------------------------------------------------------------------
%% INPUT
% dQuat_fromAtoB
% -------------------------------------------------------------------------------------------------------------
%% OUTPUT
% dDataMatrix
% bIsSignSwitched
% ui8howManySwitches
% bsignSwitchDetectionMask
% -------------------------------------------------------------------------------------------------------------
%% CHANGELOG
% 02-05-2024    Pietro Califano     Adapted from testChebyshevInterpolation script.
% 07-12-2024    Pietro Califano     Major bug fix in discontinuity detection (last occurrence was
%                                       being fixed incorrectly in some cases)
% 18-07-2025    Pietro Califano     Improve robustness of sign switch detection
% 06-08-2025    Pietro Califano     [HOTFIX] Replace sign switch detection method with more robust one.
%                                   Previous method was failing in case any of smooth zero crossing of 
%                                   any of the components!
% -------------------------------------------------------------------------------------------------------------
%% DEPENDENCIES
% [-]
% -------------------------------------------------------------------------------------------------------------

%% Function code

% Sign discontinuity detector and fix (ATTITUDE QUATERNION ONLY)

% Initialize output variables
ui32NumOfTimes      = uint32(size(dQuat_fromAtoB, 1));
ui8howManySwitches  = uint8(0);

bSignSwitchDetectionMask = false(ui32NumOfTimes, 1);
bIsSignSwitched          = false(ui32NumOfTimes, 1);

bIsNewJump = true;
% Loop through signal timestamps to detect jumps
for idT = 1:ui32NumOfTimes

    if idT > 1
        % Check sign change of the max component of the quaternion 
        [~, ui32MaxCompIdx] = max( abs( dQuat_fromAtoB(idT,:) ) );

        % Check if product of previous timestamp and current has negative sign
        if dQuat_fromAtoB(idT, ui32MaxCompIdx) * dQuat_fromAtoB(idT-1,ui32MaxCompIdx) < 0
            % Jump detected: switch sign
            dQuat_fromAtoB(idT,:) = - dQuat_fromAtoB(idT,:);

            bIsSignSwitched(idT) = true;

            if bIsNewJump
                % Add entry in detection mask (now mainly for legacy and debug reasons)
                bSignSwitchDetectionMask(idT) = true;
                ui8howManySwitches = ui8howManySwitches + 1;
                bIsNewJump = false;
            end
        else
            bIsNewJump = true;
        end
    end

end

assert(sum(all(ischange(sign(dQuat_fromAtoB)), 2) == true) == 0, ...
    'Something may have gone wrong in fixing the discontinuity!')

% Extract three components of the quaternion
dDataMatrix = dQuat_fromAtoB';


end
