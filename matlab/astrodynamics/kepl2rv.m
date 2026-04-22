function dxCart = kepl2rv(dxKepl, dGravParam, charUnit) %#codegen
arguments (Input)
    dxKepl      (6,1) double {mustBeFinite}
    dGravParam  (1,1) double {mustBeFinite, mustBePositive}
    charUnit    (1,1) string {mustBeMember(charUnit, ["rad", "deg"])} = "rad"
end
arguments (Output)
    dxCart      (6,1) double
end
%% PROTOTYPE
% dxCart = kepl2rv(dxKepl, dGravParam, charUnit)
% -------------------------------------------------------------------------------------------------------------
%% DESCRIPTION
% Converts Classical Keplerian orbital elements to a Cartesian position-velocity state vector in
% the inertial (ECI/ECEI) frame. Handles elliptic and hyperbolic orbits. Angles assumed in radians
% by default; pass charUnit="deg" to use degrees.
% State order: [SMA; Ecc; Incl; RAAN; ArgPeri; TrueAnomaly]
% Reference: Fundamentals of Astrodynamics and Applications, D. Vallado, Section 2.2.
% -------------------------------------------------------------------------------------------------------------
%% INPUT
% dxKepl     (6,1) double   Keplerian state [SMA; Ecc; Incl; RAAN; ArgPeri; TrueAnomaly]
%                            [km, -, rad, rad, rad, rad] (or deg if charUnit="deg")
% dGravParam (1,1) double   Gravitational parameter of central body [km^3/s^2]
% charUnit   (1,1) string   Angle unit of input elements: "rad" (default) or "deg"
% -------------------------------------------------------------------------------------------------------------
%% OUTPUT
% dxCart     (6,1) double   Cartesian state [x; y; z; vx; vy; vz] in inertial frame [km, km/s]
% -------------------------------------------------------------------------------------------------------------
%% CHANGELOG
% 23-10-2021    Pietro Califano     First version, Orbital Mechanics course @Polimi 2021/2022
% 12-07-2023    Pietro Califano     Re-worked and validated
% 01-05-2025    Pietro Califano     Refactoring and improvements
% 11-04-2026    PC, Claude Code     Fix description (was copy of PropagateKeplerianElems), add arguments Output,
%                                   fix local function variable naming to match convention.
% 22-04-2026    PC, Codex 5.4       Replace generic validator functions with shape/semantic validators,
%                                   normalize angles once, and simplify local rotation helpers to rad-only.
% -------------------------------------------------------------------------------------------------------------
%% DEPENDENCIES
% [-]
% -------------------------------------------------------------------------------------------------------------

%% Function code

%% Input pre-processing
dSMA    = dxKepl(1);
dEcc    = dxKepl(2);
dIncl   = dxKepl(3);
dRAAN   = dxKepl(4);
domega  = dxKepl(5);
dTA     = dxKepl(6);

if charUnit == "deg"
    dTA     = deg2rad(dTA);
    dIncl   = deg2rad(dIncl);
    dRAAN   = deg2rad(dRAAN);
    domega  = deg2rad(domega);
end

if dEcc > 1
    assert( dSMA < 0, 'MATLAB:assertion:failed', 'For hyperbolic orbits (dEcc > 1), semi-major axis a must be negative.');
end

%% Conversion [SMA, ECCTA] --> (R, V) in perifocal frame
% The semi-major axis a, the eccentricity e and the true anomaly defines
% the shape (so the energy) of the orbit; true anomaly identifies the
% position at the specific time t of the satellite along the orbit.
% From the first two, h can be computed. Then, r and v in perifocal coords.
% NOTE: the third components is 0 by definition.

% Compute parameter p = a*(1 - e^2) (always >0 for bound or unbound)
dSemiLatusRectum = dSMA*(1 - dEcc^2);

% Norm of the specific orbital ang. mom.
dOrbitAngMomentNorm = sqrt(dGravParam * dSemiLatusRectum);

dPositionNorm = (dOrbitAngMomentNorm.^2/dGravParam) .* (1/(1 + dEcc.*cos(dTA)));

dPositionPerifocal = dPositionNorm .* [cos(dTA); sin(dTA); 0]; % 3x1
dVelocityPerifocal = [-(dGravParam./dOrbitAngMomentNorm).*sin(dTA); (dGravParam./dOrbitAngMomentNorm ).*(dEcc + cos(dTA)); 0]; % 3x1

%% DCM Perifocal (e, p, h) --> Cartesian Inertial (I, J, K)
% Having r and v in perifocal and the Euler angles (313 seq) that defines the orbital
% plane wrt Cartesian Inertial, the transformation is done by means of the
% corresponding DCM (ZXZ). Inclination rotation matches the XY plane.

% Order of rotation: 1st RAAN, 2nd inclination, 3rd omega
dECEI2perifocal = Rot3(domega) * Rot1(dIncl) * Rot3(dRAAN);

dxCart = zeros(6, 1);
dxCart(1:3) = dECEI2perifocal' * dPositionPerifocal;
dxCart(4:6) = dECEI2perifocal' * dVelocityPerifocal;

%% LOCAL functions
%%% Rotation about 1st axis
    function [dRotMat1] = Rot1(dRotAngle)
        dSinAngle = sin(dRotAngle);
        dCosAngle = cos(dRotAngle);

        dRotMat1 = [1, 0, 0; ...
            0, dCosAngle, dSinAngle; ...
            0, -dSinAngle, dCosAngle];
    end

%%% Rotation about 3rd axis
    function [dRotMat3] = Rot3(dRotAngle)
        dSinAngle = sin(dRotAngle);
        dCosAngle = cos(dRotAngle);

        dRotMat3 = [dCosAngle, dSinAngle, 0; ...
            -dSinAngle, dCosAngle, 0; ...
            0, 0, 1];
    end

end
