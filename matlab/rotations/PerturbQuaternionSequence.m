function dQuatSequenceOut = PerturbQuaternionSequence(dQuatSequenceOrig, ...
                                                      dTimeStep, ...
                                                      dTimeConstFirstOrderGM, ...
                                                      dSigmaFirstOrderGM) %#codegen
%% SIGNATURE
% dQuatSequenceOut = PerturbQuaternionSequence(dQuatSequenceOrig, ...
%                                               dTimeStep, ...
%                                               dTimeConstFirstOrderGM, ...
%                                               dSigmaFirstOrderGM) %#codegen
% -------------------------------------------------------------------------------------------------------------
%% DESCRIPTION
% Function perturbing a sequence of quaternion using a time correlated First Order Gauss Markov process.
% This is integrated analytically with random input noise.
% -------------------------------------------------------------------------------------------------------------
%% INPUT
% dQuatSequence (4,:) double
% dTimeStep     (1,1) double {mustBePositive}
% dTau          (1,1) double {mustBePositive}
% dSigma        (1,1) double {mustBeNonnegative}
% -------------------------------------------------------------------------------------------------------------
%% OUTPUT
% dQuatSequenceOut (4,:)
% -------------------------------------------------------------------------------------------------------------
%% CHANGELOG
% 14-07-2025    Pietro Califano     First implementation.
% -------------------------------------------------------------------------------------------------------------
%% DEPENDENCIES
% [-]
% -------------------------------------------------------------------------------------------------------------
arguments (Input)
    dQuatSequenceOrig           (4,:) double
    dTimeStep                   (1,1) double {mustBePositive}
    dTimeConstFirstOrderGM      (1,1) double {mustBePositive}
    dSigmaFirstOrderGM          (1,1) double {mustBeNonnegative}
end
arguments (Output)
    dQuatSequenceOut            (4,:) double
end

% Preallocate output
ui32NumSamples      = uint32(size(dQuatSequenceOrig, 2));
dQuatSequenceOut    = zeros(4, ui32NumSamples);

% Initialize Gauss–Markov error state (small-angle approximation)
dTmpErrorAnglesPrev = zeros(3,1);

% Precompute coefficients
dAlphaIntegr     = exp(-dTimeStep / dTimeConstFirstOrderGM);
dSigmaInputConst = dSigmaFirstOrderGM * sqrt(1 - dAlphaIntegr^2);

% Loop through each sample
for idQ = 1:ui32NumSamples
    
    % Sample white noise and update error angles
    dTmpNoise       = dSigmaInputConst * randn(3,1);
    dTmpErrorAngles = dAlphaIntegr * dTmpErrorAnglesPrev + dTmpNoise;

    % Convert to quaternion error
    dQuatErr = AnglesToquatErr(dTmpErrorAngles);

    % Apply error quaternion
    dQuatSequenceOut(:,idQ) = qCross(dQuatErr, dQuatSequenceOrig(:, idQ), false);

    % Store for next iteration
    dTmpErrorAnglesPrev = dTmpErrorAngles;
end

% Normalize output quaternions to unit norm
dQuatSequenceOut(:, :) = dQuatSequenceOut(:, :) ./ vecnorm(dQuatSequenceOut, 2, 1);

end

%% Auxiliary functions
function dQuatErr = AnglesToquatErr(dAngles)
arguments
    dAngles (3,:) double {isvector, mustBeNumeric}
end

dTheta = norm(dAngles);

if dTheta > 0
    dAxis = dAngles / dTheta;
else
    dAxis = [1, 0, 0];
end

dHalfTheta = 0.5 * dTheta;
dQuatErr = [cos(dHalfTheta); dAxis * sin(dHalfTheta)];

end

