function o_bFeasibleLoS = checkLoS2Point(i_dPos1, i_dPos2, i_dRbody) %#codegen
%% PROTOTYPE
% o_bFeasibleLoS = checkLoS2Point(i_dPos1, i_dPos2, i_dRbody) %#codegen
% -------------------------------------------------------------------------------------------------------------
%% DESCRIPTION
% Function checking for the existence of a line-of-sight joining two points at locations specified by the
% input position vectors. This may be used as preliminary check for Sun illumination or satellite to
% satellite visibility. Computationally efficient algorithm by (Alfano, 1991). See [1].
% ACHTUNG: the model assumes both the points orbiting a spherical main attractor of radius equal to i_dRbody
% The reference frame origin is assumed coincident with the centre of the latter body.
% REFERENCE
% [1] Fundamentals of Astrodynamics and Applications - D. Vallado, 4th Edition, page 308, Algorithm 35.
% -------------------------------------------------------------------------------------------------------------
%% INPUT
% i_dPos1:        [3, 1]  Position vector of Point 1 in a frame centred in target body CoM
% i_dPos2:        [3, 1]  Position vector of Point 2 in a frame centred in target body CoM
% i_dRbody:       [1]     Radius of the spherical main body
% -------------------------------------------------------------------------------------------------------------
%% OUTPUT
% o_bFeasibleLoS: [1]     Boolean checking if a LoS between P1 and P2 is feasible (true)
% -------------------------------------------------------------------------------------------------------------
%% CHANGELOG
% 14-02-2024        Pietro Califano         Function coded. Not verified.
% -------------------------------------------------------------------------------------------------------------
%% DEPENDENCIES
% [-]
% -------------------------------------------------------------------------------------------------------------
%% Future upgrades
% 1) Scaling of inputs to account for equatorial bulge (i.e. flattening)
% -------------------------------------------------------------------------------------------------------------
%% Function code

% Define default case
o_bFeasibleLoS = false;

% Compute auxiliary variables
dPos1dotPos2 = dot(i_dPos1, i_dPos2);
dPos1Norm = norm(i_dPos1);
dPos1Norm2 = dPos1Norm*dPos1Norm;

dTauMin = (dPos1Norm2 - dPos1dotPos2)/(dPos1Norm2 + norm(i_dPos2)^2 - 2*dPos1dotPos2);

% Check if TauMin determines feasible LoS
if dTauMin < 0 || dTauMin > 1
    o_bFeasibleLoS = true;
else
    cAuxVecNorm2 = (1 - dTauMin)^2 * dPos1Norm2 + dPos1dotPos2 * dTauMin;
    % Check alternative feasibility condition
    if cAuxVecNorm2 >= (i_dRbody*i_dRbody)
        o_bFeasibleLoS = true;
    end
end



end