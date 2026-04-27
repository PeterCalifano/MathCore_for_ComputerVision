function ui8ShadowType = EvaluateConeShadowing(dPosSC_W, ...
                                                dSunPos_W, ...
                                                dBodyRadius, ...
                                                dAlphaPenumbraInRad, ...
                                                dAlphaUmbraInRad) %#codegen
arguments
    dPosSC_W             (3,1) double
    dSunPos_W            (3,1) double
    dBodyRadius          (1,1) double {mustBePositive}
    dAlphaPenumbraInRad  (1,1) double {mustBePositive, mustBeLessThan(dAlphaPenumbraInRad, 1.5708)}
    dAlphaUmbraInRad     (1,1) double {mustBePositive, mustBeLessThan(dAlphaUmbraInRad, 1.5708)}
end 
%% PROTOTYPE
% ui8ShadowType = EvaluateConeShadowing(dSCpos_IN, dSunPos_IN, dRbody, dAlphaPenumbraInRad, dAlphaUmbraInRad)
% -------------------------------------------------------------------------------------------------------------
%% DESCRIPTION
% Function checking shadowing condition for a spherical planet/body orbiting the Sun. The reference frame is
% assumed centred in the body CoM; the input values (dAlphaPenumbraInRad, dAlphaUmbraInRad) depends on the
% distance from the Sun, but are relatively constant for circular orbits. See [1].
% REFERENCE
% [1] Fundamentals of Astrodynamics and Applications - D. Vallado, 4th Edition, page 301, Algorithm 34.
% -------------------------------------------------------------------------------------------------------------
%% INPUT
% dPosSC_W:             [3,1]   Spacecraft position in IN frame centred in the attractor
% dSunPos_W:            [3,1]   Position vector to the Sun in IN frame
% dBodyRadius:          [1,1]      Radius of the spherical attractor (CoM := IN origin)
% dAlphaPenumbraInRad:  [1,1]      Complementary half-angle of Penumbra shadow cone. See [1] for details
% dAlphaUmbraInRad:     [1,1]      Complementary half-angle of Umbra shadow cone. See [1] for details
% -------------------------------------------------------------------------------------------------------------
%% OUTPUT
% ui8ShadowType:   [1,1]      Integer indicating shadow type: 0: no shadow; 1: Penumbra; 2: Umbra
% -------------------------------------------------------------------------------------------------------------
%% CHANGELOG
% 14-02-2024        Pietro Califano         Function coded. Not verified.
% 02-01-2026        Pietro Califano         Update naming conventions
% -------------------------------------------------------------------------------------------------------------
%% DEPENDENCIES
% [-]
% -------------------------------------------------------------------------------------------------------------

%% Function code

% Default value assignment
ui8ShadowType = uint8(0); % No shadowing case

% Compute auxiliary variables
dSCposDotSunPos = dot(dPosSC_W, dSunPos_W);
dPosNorm = norm(dPosSC_W);

if dSCposDotSunPos < 0 % If true --> there could be shadowing
    % Check for penumbra case

    % Compute angle between directions
    dCosTheta = dSCposDotSunPos / ( dPosNorm * norm(dSunPos_W));
    dThetaAngle = acos( dCosTheta );

    % Decompose SC position vector
    dHorizSatDist = dPosNorm * dCosTheta;
    dVertSatDist = dPosNorm * sin(dThetaAngle);

    % Compute distance from Planet centre to intersection point
    dxPoint = dBodyRadius / ( sin(dAlphaPenumbraInRad) );
    dPenumbraVert = tan(dAlphaPenumbraInRad) * (dxPoint + dHorizSatDist);

    if dVertSatDist <= dPenumbraVert
        ui8ShadowType = uint8(1); % Penumbra shadow

        % Check for umbra case
        dyPoint = dBodyRadius / ( sin(dAlphaUmbraInRad) );
        dUmbraVert = tan(dAlphaUmbraInRad) * (dyPoint - dHorizSatDist);

        if dVertSatDist <= dUmbraVert
            ui8ShadowType = uint8(2); % Umbra shadow
        end
    end
end

end



