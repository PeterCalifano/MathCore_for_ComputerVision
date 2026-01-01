function [strConicData, i8ReturnId] = ComputeConicDataFromEqCoeffs(dCoeffsABCDEF, ...
                                                                   bALLOW_OPEN_CONICS)%#codegen
arguments
    dCoeffsABCDEF (1,:) double {mustBeNumeric, mustBeFinite}
    bALLOW_OPEN_CONICS (1,1) logical = true
end

%% SIGNATURE
%
% -------------------------------------------------------------------------------------------------------------
%% DESCRIPTION
% What the function does
% 
% REFERENCE:
% [1] J. A. Christian, “A Tutorial on Horizon-Based Optical Navigation and Attitude Determination With Space
% Imaging Systems,” IEEE Access, vol. 9, pp. 19819–19853, 2021, doi: 10.1109/ACCESS.2021.3051914.
% [2] Wikipedia page: "Ellipse", section: General parameters representation
% [3] J. A. Christian, Optical Navigation Using Planet s Centroid and Apparent Diameter in Image, 2015
% -------------------------------------------------------------------------------------------------------------
%% INPUT
% in1 [dim] description
% -------------------------------------------------------------------------------------------------------------
%% OUTPUT
% out1 [dim] description
% -------------------------------------------------------------------------------------------------------------
%% CHANGELOG
% 01-01-2026    Pietro Califano     Implement general purpose version from legacy ComputeHorizonConic().
% -------------------------------------------------------------------------------------------------------------
%% DEPENDENCIES
% [-]
% -------------------------------------------------------------------------------------------------------------

%% Function code
% Interesting NOTE from [2]
% The discriminant B2 – 4AC of the conic [section's quadratic equation (or equivalently the determinant 
% AC – B2/4 of the 2 × 2 matrix) and the quantity A + C (the trace of the 2 × 2 matrix) are invariant under 
% arbitrary rotations and translations of the coordinate axes, as is the determinant of the 3 × 3 matrix above.
% The constant term F and the sum D2 + E2 are invariant under rotation only.

dNUMERICAL_LIMIT = 2.5*eps;
i8ReturnId = int8(-99);

% Set inline properties
coder.inline('default');

% Initialize struct
strConicData = struct();
strConicData.dConicCentre       = zeros(1,2,'double');
strConicData.dSemiMajorAx         = 0.0;
strConicData.dSemiMinorAx         = 0.0;
strConicData.dMajorAxisAngleFromX = 0.0;
strConicData.dEccentricity        = 0.0;

coder.cstructname(strConicData, 'strConicData');

%% Auxiliary data

% Compute equation discriminant
dDiscriminant = B^2 - 4*A*C;

% Determine conic section type
if abs(dDiscriminant) < 1.5 * eps
    dDiscriminant = 0.0; % Discriminant is below machine precision, set to 0
end

dEta = 0.0;
if dDiscriminant > 0.0
    dEta = -1.0;
elseif dDiscriminant < 0.0
    dEta = 1.0;
end

% Compute common quantities
dSqrtACB = sqrt( (A-C)^2 + B^2 );

%% Compute conic data based on type
% Initialize data
dEccentricity = 0.0;
dConicCentre = zeros(1,2);

% Compute eccentricity of the Conic section
dInvDiscriminant = 1.0;
if dDiscriminant ~= 0
    dEccentricity = sqrt( (2*dSqrtACB) / (dEta * (A+C) + dSqrtACB) );
    dInvDiscriminant = 1 / dDiscriminant;
end

if dDiscriminant < -dNUMERICAL_LIMIT && ( dEccentricity >= 0 && dEccentricity < 1 )

    % Closed conic (ellipse/circle)
    dConicCentre(1) = (2*C*D - B*E) * dInvDiscriminant;
    dConicCentre(2) = (2*A*E - B*D) * dInvDiscriminant;

    % Calculate the angle of rotation
    dMajorAxisAngleFromX = 0.5 * double(atan2(-real(B), real(C-A)));

    % Calculate the semi-major and semi-minor axes
    aTmp = -sqrt( 2*(A * E^2 + C*D^2 - B*D*E + F*(dDiscriminant) ) *( A + C + dSqrtACB ) ) * dInvDiscriminant;
    bTmp = -sqrt( 2*(A * E^2 + C*D^2 - B*D*E + F*(dDiscriminant) ) *( A + C - dSqrtACB ) ) * dInvDiscriminant;

    assert(aTmp > 0 && bTmp > 0, 'Error: Conic Section considered as Closed, by either a or b turned out negative. Something may have gone wrong in the computation.')

    if aTmp > 0 && bTmp > 0
        i8ReturnId = int8(-1);
        return
    end

    if aTmp > bTmp
        dSemiMajorAx = aTmp;
        dSemiMinorAx = bTmp;

    else
        dSemiMajorAx = bTmp;
        dSemiMinorAx = aTmp;

        % Rotate angle
        dMajorAxisAngleFromX = coder.const(0.5*pi) + dMajorAxisAngleFromX;
    end


elseif dDiscriminant > dNUMERICAL_LIMIT && dEccentricity >= 1
    % Open conic (hyperbola)

    assert(bALLOW_OPEN_CONICS, 'Execution stop: Conic section is open, but bALLOW_OPEN_CONICS set to false ')
    warning('Still requires validation.')

    if not(bALLOW_OPEN_CONICS)
        i8ReturnId = int8(-98);
        return
    end

    % Build matrices Q and L
    dQ = [A, B/2; B/2, C];
    dL = [D; E];

    % Compute centre of the conic
    dConicCentre(1:2) = - dQ\dL;

    % Compute orientatin (angle of the line perpendicular to the Directrix)
    dMajorAxisAngleFromX = 0.5 * atan( B / (A-C) );

    % Find absolute values of a and b
    [~, dEigenValsMatrix] = eig(dQ);

    dEigenVals = diag(dEigenValsMatrix);
    dOneOverSqrtEig = 1 / sqrt(abs(dEigenVals(1:2)));

    dSemiMajorAx = dOneOverSqrtEig(1);
    dSemiMinorAx = dOneOverSqrtEig(2);

elseif abs(dDiscriminant) <= dNUMERICAL_LIMIT

    % Compute centre: from ∂/∂x and ∂/∂y = 0 for x^2 + y^2 + Dx + Ey + F = 0
    dConicCentre(1:2) = -real([D,E]) / (2 * real(A));

    % No rotation needed
    dMajorAxisAngleFromX = 0;

    % Radius: r = sqrt( (D^2 + E^2)/(4A^2) - F/A )
    % dTmpRadius = sqrt( real( (D^2 + E^2) / (4 * A^2) - F / A ) );
    dTmpRadius = sqrt( (D^2 + E^2 - 4 *A*F) ) / (2 * abs(A));
    
    assert( dTmpRadius > 0, ...
        'ERROR: computed circle radius is non-positive.' )

    if dTmpRadius <= 0
        i8ReturnId = int8(-2);
        return
    end

    % Semi-axes are equal
    dSemiMajorAx = dTmpRadius;
    dSemiMinorAx = dTmpRadius;
else
    % Return invalid (failed at numerical tolerance level)
    i8ReturnId = int8(-3);
    return
end

% Assign data to output
strConicData.dConicCentre(:) = dConicCentre;
strConicData.dEccentricity          = dEccentricity;
strConicData.dMajorAxisAngleFromX   = dMajorAxisAngleFromX;
strConicData.dSemiMajorAx           = dSemiMajorAx;
strConicData.dSemiMinorAx           = dSemiMinorAx;
i8ReturnId = int8(0);

end