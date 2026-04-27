function [fDot, fDotDot] = FDMkernel(samplesPoints, samplesTimegrid, type) %#codegen
%% PROTOTYPE
% [fDot, fDotDot] = FDMkernel(samplesPoints, samplesTimegrid, type)
% -------------------------------------------------------------------------------------------------------------
%% DESCRIPTION
% Function computing the approximate 1st and 2nd order derivative from
% sample points and scalar domain (intended for time derivatives).
% -------------------------------------------------------------------------------------------------------------
%% INPUT
% samplesPoints
% samplesTimegrid
% type: [scalar] 1: Forward, 2: Central
% -------------------------------------------------------------------------------------------------------------
%% OUTPUT
% fDot
% fDotDot
% -------------------------------------------------------------------------------------------------------------
%% CHANGELOG
% 12-06-2023  Pietro Califano  FDM Forward (1st order) and Central (2nd 
%                              order method coded and validated with  
%                              analytical function 
% 03-07-2023  Pietro Califano  FDM Forward up to 6th order accuracy;
%                              Central of 8th order accuracy
% -------------------------------------------------------------------------------------------------------------
%% DEPENDENCIES
%  [-]
% -------------------------------------------------------------------------------------------------------------
[stateLen, Npoints] = size(samplesPoints);

y = zeros(Npoints, 1);
y(:, 1) = samplesPoints;

dt = samplesTimegrid(2) - samplesTimegrid(1);

% Static allocation
fDot = zeros(stateLen, 1);

if nargout > 1
    fDotDot = zeros(stateLen, 1);
end

% Input checks
assert(Npoints == length(samplesTimegrid), 'Number of sample points and timetags are not matched.');

if type == 1 % Apply Forward difference

    % Storage cell
    FDcoeff = {0;
        [-1, 1];
        [-3/2, 2, -1/2];
        [-11/6, 3 - 3/2, 1/3];
        [-25/12, 4, -3, 4/3, -1/4];
        [-137/60, 5, -5, 10/3, -5/4, 1/5];
        [-49/20, 6, -15/2, 20/3, -15/4, 6/5, -1/6]};

    % Get required coefficients
    coeff = FDcoeff{Npoints};

    % Apply Forward Difference
    if Npoints <= 7
        fDot = (1/dt) * (coeff * y)';
    else
        error('Forward Difference with more than 7 points not implemented');
    end

elseif type == 2 % Apply central difference

    CDcoeff = {0;
        [-1/2, 0, 1/2];
        [1/280, -4/105, 1/5, -4/5, 0, ]};

    if Npoints == 3 || Npoints == 9

        if Npoints == 3
            coeff = CDcoeff{2};
        elseif Npoints == 9
            coeff = CDcoeff{3};
        end

        fDot = 1/dt .* (coeff * y)';

        if nargout > 1 && Npoints == 3
            fDotDot = ( y(3) - 2*y(2) + y(1) )/(dt^2);
        end

    else
        error('Central Difference with n° points not equal to 3 or 9 not implemented')
    end
end

end