function [dDataMatrix, bIsSignSwitched, ui8howManySwitches, bsignSwitchDetectionMask] =...
    fixQuatSignDiscontinuity(dQuat_fromAtoB) %#codegen
arguments
    dQuat_fromAtoB (:, 4) double {isvector, mustBeNumeric}
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
% 02-05-2024        Pietro Califano         Adapted from testChebyshevInterpolation script.
% 07-12-2024        Pietro Califano         Major bug fix in discontinuity detection (last occurrence was
%                                           being fixed incorrectly in some cases)
% 18-07-2025        Pietro Califano         Improve robustness of sign switch detection
% -------------------------------------------------------------------------------------------------------------
%% DEPENDENCIES
% [-]
% -------------------------------------------------------------------------------------------------------------
%% Future upgrades
% [-]
% -------------------------------------------------------------------------------------------------------------
%% Function code

% Sign discontinuity detector and fix (ATTITUDE QUATERNION ONLY)
bsignSwitchDetectionMask = sign(dQuat_fromAtoB(:, 1:4)); % Use first component as sign check
% FIXME: switch detection is not robust to cases in which any of the quaternion states jumps

bsignSwitchDetectionMask = any( ischange(bsignSwitchDetectionMask), 2);
ui8howManySwitches = uint8(sum(bsignSwitchDetectionMask == true));
interpSignal = dQuat_fromAtoB;
bIsSignSwitched = false(size(dQuat_fromAtoB, 2), 1);

if ui8howManySwitches > 0
    % Get where the switches happens
    switchIdx = find(bsignSwitchDetectionMask, ui8howManySwitches);

    startIntervalsIDs = 1:2:length(switchIdx);

    for idToFix = startIntervalsIDs

        idStart = switchIdx(idToFix);

        if (idToFix == startIntervalsIDs(end) && length(startIntervalsIDs) > 1) ...
                && mod(double(ui8howManySwitches), 2) ~= 0 || isscalar(switchIdx)

            idEnd = length(interpSignal);

        elseif (idToFix == startIntervalsIDs(end) && length(startIntervalsIDs) > 1) ...
                && mod(double(ui8howManySwitches), 2) == 0 

            idEnd = switchIdx(end)-1;
            
        else
            idEnd = switchIdx(idToFix+1)-1 ;
        end

        interpSignal(idStart:idEnd, :) = -interpSignal(idStart:idEnd, :);
        bIsSignSwitched(idStart:idEnd) = true;
    end
end

assert(sum(all(ischange(sign(interpSignal)), 2) == true) == 0, ...
    'Something may have gone wrong in fixing the discontinuity!')

% Extract three components of the quaternion
% i_dDataMatrix = interpSignal(:, 1:3)';
dDataMatrix = interpSignal';


end
