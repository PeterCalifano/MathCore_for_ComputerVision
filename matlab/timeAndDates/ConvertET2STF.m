function dSTF = ConvertET2STF(dEphemerisTime) %#codegen
arguments
    dEphemerisTime (:,1) double {mustBeNumeric}
end
%% PROTOTYPE
% dSTF = ConvertET2STF(dEphemerisTime)
% -------------------------------------------------------------------------------------------------------------
%% DESCRIPTION
% Converts SPICE Ephemeris Time (ET) to Standard Tyvak Format (STF, GPS days).
% Conversion chain:
%   ET  = TAI + 32.184       [s]
%   GPS = TAI - 19           [s]
%   GPS = ET  - 32.184 - 19  [s]
% An additional constant accounts for the epoch difference:
%   ET epoch:  J2000.0 (12:00:00 TDB, 1 Jan 2000)
%   GPS epoch: 00:00:00 UTC, 6 Jan 1980
% -------------------------------------------------------------------------------------------------------------
%% INPUT
% dEphemerisTime   (:,1) SPICE Ephemeris Time [s past J2000.0]
% -------------------------------------------------------------------------------------------------------------
%% OUTPUT
% dSTF             (:,1) GPS days in Standard Tyvak Format [days]
% -------------------------------------------------------------------------------------------------------------
%% CHANGELOG
% 12-04-2026    Pietro Califano     Modernized interface, imported to MathCore from RCS-1.
% -------------------------------------------------------------------------------------------------------------
%% DEPENDENCIES
% [-]
% -------------------------------------------------------------------------------------------------------------

%% Function code

dDenom = coder.const(1 / 86400); % seconds per day
dSTF = (dEphemerisTime - 32.184 - 19) * dDenom + 7300.50015046139;

end
