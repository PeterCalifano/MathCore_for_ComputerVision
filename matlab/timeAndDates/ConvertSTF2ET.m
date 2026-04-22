function dEphemerisTime = ConvertSTF2ET(dSTF) %#codegen
arguments
    dSTF (:,1) double {mustBeNumeric}
end
%% PROTOTYPE
% dEphemerisTime = ConvertSTF2ET(dSTF)
% -------------------------------------------------------------------------------------------------------------
%% DESCRIPTION
% Converts Standard Tyvak Format (STF, GPS days) to SPICE Ephemeris Time (ET).
% Inverse of convert_et2stf.
% -------------------------------------------------------------------------------------------------------------
%% INPUT
% dSTF             (:,1) GPS days in Standard Tyvak Format [days]
% -------------------------------------------------------------------------------------------------------------
%% OUTPUT
% dEphemerisTime   (:,1) SPICE Ephemeris Time [s past J2000.0]
% -------------------------------------------------------------------------------------------------------------
%% CHANGELOG
% 12-04-2026    Pietro Califano     Modernized interface, imported to MathCore from RCS-1.
% -------------------------------------------------------------------------------------------------------------
%% DEPENDENCIES
% [-]
% -------------------------------------------------------------------------------------------------------------

%% Function code
dEphemerisTime = 32.184 + 19 + 86400 * (dSTF - 7300.50015046139);

end
