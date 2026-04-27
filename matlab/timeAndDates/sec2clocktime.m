function [hours, minutes, seconds] = sec2clocktime(time_in_seconds, display_flag, var_name)
% PROTOTYPE
% [hours, minutes, seconds] = sec2clocktime(time_in_seconds, display_flag, var_name)
% -------------------------------------------------------------------
% DESCRIPTION
% Converts time expressed in second into hours-minutes-seconds format
% -------------------------------------------------------------------
% INPUT
%    time_in_seconds    [1]         time in [s] to convert 
%    var_name           [string]    name of the converted time variable
%    display_flag       [0 or 1]    flag to enable displaying of the output
% -------------------------------------------------------------------------
% OUTPUT
% [-]
% -------------------------------------------------------------------------
% CONTRIBUTORS
%    01-12-2021    Pietro Califano    First Version (tested)
%    15-03-2023    Pietro Califano    Added display flag
% -------------------------------------------------------------------------
% DEPENDENCIES
% [-]
% -------------------------------------------------------------------------
% Future upgrades
% [-]
% -------------------------------------------------------------------------

if ~exist('var_name', 'var')
    string = 'Time = ';
elseif exist('var_name', 'var')
    string = var_name;
end

if time_in_seconds < 60
    hours = 0;
    minutes = 0;
    seconds = time_in_seconds;

    if display == 1
        disp([string, num2str(seconds), ' [s]']);
    end

elseif time_in_seconds > 60
    % Convert to hours
    time_in_seconds = time_in_seconds/3600;
    % Get hours
    hours = floor(time_in_seconds);
    % Get minutes
    minutes = floor((time_in_seconds - hours)*60);
    % Get seconds
    tempsec = (time_in_seconds - hours)*60 - minutes;
    seconds = round((tempsec) * 60);

    if display_flag == 1
        fprintf('\n%s\t%2.0d h %2.0d min %2.0d s\n', string, hours, minutes, seconds);
    end
end

end