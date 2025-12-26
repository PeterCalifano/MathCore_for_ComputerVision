function dDCM = UniformlySampleSO3HaarDistr(ui32NumOfSamples)
arguments
    ui32NumOfSamples (1,1) uint32 {isscalar, mustBeNumeric} = 1
end

% Generate three uniform random numbers
dU1 = rand(1, ui32NumOfSamples);
dU2 = rand(1, ui32NumOfSamples);
dU3 = rand(1, ui32NumOfSamples);

% Compute quaternion components
dQuat1 = sqrt(1 - dU1) .* sin(2 * pi * dU2);
dQuat2 = sqrt(1 - dU1) .* cos(2 * pi * dU2);
dQuat3 = sqrt(dU1) .* sin(2 * pi * dU3);
dQuat4 = sqrt(dU1) .* cos(2 * pi * dU3);

% Convert quaternion to rotation matrix
dDCM = zeros(3,3,ui32NumOfSamples);

for idSampl = 1:ui32NumOfSamples
    dTmp1 = dQuat1(1, idSampl);
    dTmp2 = dQuat2(1, idSampl);
    dTmp3 = dQuat3(1, idSampl);
    dTmp4 = dQuat4(1, idSampl);

    % Convert quaternion to rotation matrix 
    dDCM(:,:, idSampl) = [1 - 2*(dTmp3^2 + dTmp4^2), 2*(dTmp2*dTmp3 - dTmp1*dTmp4), 2*(dTmp2*dTmp4 + dTmp1*dTmp3);
        2*(dTmp2*dTmp3 + dTmp1*dTmp4), 1 - 2*(dTmp2^2 + dTmp4^2), 2*(dTmp3*dTmp4 - dTmp1*dTmp2);
        2*(dTmp2*dTmp4 - dTmp1*dTmp3), 2*(dTmp3*dTmp4 + dTmp1*dTmp2), 1 - 2*(dTmp2^2 + dTmp3^2)];

end

end
