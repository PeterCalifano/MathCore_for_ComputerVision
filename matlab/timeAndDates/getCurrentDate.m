function [o_dDayOfYear, o_dUTtimeSecNow, o_dYear] = getCurrentDate(i_dDate0, i_dElapsedSeconds) %#codegen
%% PROTOTYPE
% [o_dDayOfYear, o_dUTtimeSecNow, o_dYear] = getCurrentDate(i_dDate0, i_dElapsedSeconds)
% -------------------------------------------------------------------------------------------------------------
%% DESCRIPTION
% Function updating the date tailored for MATLAB NRLMSISE Atmosphere model
% requested inputs, given the simulation elapsed time in seconds and the
% input date as [YYYY, MM, DD, hh, mm, ss].
% -------------------------------------------------------------------------------------------------------------
%% INPUT
% i_dDate0:          [6x1] Array of double specifying the initial date at t0: [YYYY, MM, DD, hh, mm, ss]
% i_dElapsedSeconds: [1] Number elapsed seconds from simulation start (from t0). Fractional increments allowed. 
% -------------------------------------------------------------------------------------------------------------
%% OUTPUT
% o_dDayOfYear:    [1] Current day of the year counting from January 1,
%                      assuming 365 days in a year.
% o_dUTtimeSecNow: [1] Number of seconds from Midnight of current day (UT)
% o_dYear:         [1] Current year of the simulation
% -------------------------------------------------------------------------------------------------------------
%% CHANGELOG
% 20-01-2024     Pietro Califano     Function coded and verified at
%                                    day and year change.
% -------------------------------------------------------------------------------------------------------------
%% DEPENDENCIES
% [-]
% -------------------------------------------------------------------------------------------------------------
%% Future upgrades
% 1) Increase in accuracy (e.g. accounting for actual year duration)
% -------------------------------------------------------------------------------------------------------------
%% Function code
% Assumed Entries in i_dDate0
% year    = i_dDate0(1);
% month   = i_dDate0(2);
% day     = i_dDate0(3);
% hour    = i_dDate0(4);
% minutes = i_dDate0(5);
% second  = i_dDate0(6);

% Constant variables
nDaysInMonths = [31 28 31 30 31 30 31 31 30 31 30 31];
secondsInDay = 3600*24;

% Compute UT time in current DOY
o_dUTtimeSecNow = (i_dDate0(6) + 60*i_dDate0(5) + 3600* i_dDate0(4)) + i_dElapsedSeconds;

% Compute "day overflow" and adjust UT accordingly
DeltaDays = floor(o_dUTtimeSecNow/secondsInDay);

if DeltaDays > 0
    o_dUTtimeSecNow = o_dUTtimeSecNow - secondsInDay * DeltaDays;
end

% Compute the current day of the year
o_dDayOfYear = sum( nDaysInMonths(1:(i_dDate0(2))-1) ) + i_dDate0(3) + DeltaDays;

% Check if "year overflow" occurred
if o_dDayOfYear >= 365
    o_dDayOfYear = o_dDayOfYear - 365;
    o_dYear = i_dDate0(1) + 1;
else
    % Get Year of Date0
    o_dYear = i_dDate0(1);
end


end
