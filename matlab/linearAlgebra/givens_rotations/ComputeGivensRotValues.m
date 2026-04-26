function [dCosTheta, dSinTheta] = ComputeGivensRotValues(dVec2) %#codegen
arguments (Input)
    dVec2 (2,1) double
end
arguments (Output)
    dCosTheta   (1,1) double 
    dSinTheta   (1,1) double 
end
%% SIGNATURE
% [dCosTheta, dSinTheta] = ComputeGivensRotValues(dVec2) %#codegen
% -------------------------------------------------------------------------------------------------------------
%% DESCRIPTION
% Function computing the 2x2 Givens rotation matrix entries from 2x1 vector such that the second entry is nulled out when G is applied to the vector, i.e. [r; 0] = G^T * dVec2. This corresponds to theta* angle.
% -------------------------------------------------------------------------------------------------------------
%% INPUT
% dVec2 (2,1) double
% -------------------------------------------------------------------------------------------------------------
%% OUTPUT
% dCosTheta   (1,1) double
% dSinTheta   (1,1) double
% -------------------------------------------------------------------------------------------------------------
%% CHANGELOG
% 24-02-2025    Pietro Califano     First version implemented.
% -------------------------------------------------------------------------------------------------------------
%% DEPENDENCIES
% [-]
% -------------------------------------------------------------------------------------------------------------
% Allocate default values for the case b = 0
dCosTheta = 1.0;
dSinTheta = 0.0; 

if dVec2(2) == 0.0
    return;
end

% Compute G matrix entries
if abs(dVec2(2)) > abs(dVec2(1))

    % Compute ratio of dVec entries
    dTauRatio = -dVec2(1)/dVec2(2);
    
    % Compute sine of theta*
    dSinTheta = 1 / ( sqrt(1 + dTauRatio*dTauRatio) );

    % Compute cosine of theta*
    dCosTheta = dSinTheta * dTauRatio;

else

    % Compute ratio of dVec entries
    dTauRatio = -dVec2(2)/dVec2(1);
    
    % Compute sine of theta*
    dCosTheta = 1 / ( sqrt(1 + dTauRatio*dTauRatio) );

    % Compute cosine of theta*
    dSinTheta = dCosTheta * dTauRatio;

end

end
