function dxCart = kepl2rv(dxKepl, dGravParam, charUnit) %#codegen
arguments
    dxKepl      (6,1) double {isvector, mustBeNumeric}
    dGravParam  (1,1) double {isscalar, mustBeNumeric, mustBeGreaterThan(dGravParam, 0)}
    charUnit    (1,:) string {mustBeMember(charUnit, ["rad", "deg"])} = "rad"
end
%% PROTOTYPE
% dxCart = kepl2rv(dxKepl, dGravParam, charUnit) %#codegen
% -------------------------------------------------------------------------------------------------------------
%% DESCRIPTION
% General purpose Keplerian propagator in Classical Keplerian elements. It propagates the initial state over
% the specified time grid of time intervals. The function handles all types of keplerian orbits.
% 1) Fundamentals of Astrodynamics and Applications - D. Vallado. Section
%    2.2, Kepler's Problem. Algorithms 3 and 4.
% -------------------------------------------------------------------------------------------------------------
%% INPUT
% arguments
%     dxKepl      (6,1) double {isvector, mustBeNumeric}
%     dGravParam  (1,1) double {isscalar, mustBeNumeric, mustBeGreaterThan(dGravParam, 0)}
%     charUnit    (1,:) string {mustBeMember(charUnit, ["rad", "deg"])} = "rad"
% end
% -------------------------------------------------------------------------------------------------------------
%% OUTPUT
% dxCart: [6, 1] State in cartesian coordinates 
% -------------------------------------------------------------------------------------------------------------
%% CHANGELOG
% 23-10-2021    Pietro Califano     First version, Orbital Mechanics course @Polimi 2021/2022
% 12-07-2023    Pietro Califano     Re-worked and validated
% 01-05-2025    Pietro Califano     Refactoring and improvements
% -------------------------------------------------------------------------------------------------------------
%% DEPENDENCIES
% [-]
% -------------------------------------------------------------------------------------------------------------
%% Future upgrades
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
dECEI2perifocal = Rot3(domega, "rad") * Rot1(dIncl, "rad") * Rot3(dRAAN, "rad");

dxCart = zeros(6, 1);
dxCart(1:3) = dECEI2perifocal' * dPositionPerifocal; 
dxCart(4:6) = dECEI2perifocal' * dVelocityPerifocal;

%% LOCAL functions
%%% Rotation about 1st axis
    function [Rot_Mat1] = Rot1(dRotAngle, charUnit)
        if charUnit == "rad"

            Rot_Mat1 = [1, 0, 0; ...
                0, cos(dRotAngle), sin(dRotAngle); ...
                0, -sin(dRotAngle), cos(dRotAngle)];

        elseif charUnit == "deg"

            Rot_Mat1 = [1, 0, 0; ...
                0, cosd(dRotAngle), sind(dRotAngle); ...
                0, -sind(dRotAngle), cosd(dRotAngle)];
        end

    end

%%% Rotation about 3rd axis
    function [Rot_Mat3] = Rot3(rot_angle, unit)

        if unit == "rad"

            Rot_Mat3 = [cos(rot_angle), sin(rot_angle), 0;...
                -sin(rot_angle), cos(rot_angle), 0;...
                0, 0, 1];

        elseif unit == "deg"

            Rot_Mat3 = [cosd(rot_angle), sind(rot_angle), 0;...
                -sind(rot_angle), cosd(rot_angle), 0;...
                0, 0, 1];

        end
    end

end
