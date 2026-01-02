function [xKepl, figKeplElems] = plotKeplElem(X_SC, mu, timeGrid, angleUnit, figKeplElems)
%% PROTOTYPE
% [xKepl, figKeplElems] = plotKeplElem(X_SC, mu, timeGrid, angleUnit)
% -------------------------------------------------------------------------------------------------------------
%% DESCRIPTION
% Function plotting the history of the Osculating Keplerian Elements
% computed from the input timeseries of the SC state vector. A timegrid can
% be added as optional input for the plotting routine. The unit of measure
% for the angle in output can be selected as either "deg" or "rad".
% -------------------------------------------------------------------------------------------------------------
%% INPUT
% X_SC
% mu
% timeGrid
% angleUnit
% -------------------------------------------------------------------------------------------------------------
%% OUTPUT
% xKepl,
% figKeplElems
% -------------------------------------------------------------------------------------------------------------
%% CHANGELOG
% 20-08-2023    Pietro Califano     Coded from CRE_UnitTest script. First
%                                   version tested.
% -------------------------------------------------------------------------------------------------------------
%% DEPENDENCIES
% rv2kepl()
% -------------------------------------------------------------------------------------------------------------

%% Function code
keplName = {'SMA', 'ECC', 'INCL', 'RAAN', 'OMEGA', 'TA'};
unitsPlot = {' [m]', ' [-]', ' [deg]', ' [deg]', ' [deg]', ' [deg]'};


% Define mu, timegrid, angle unit if not given
    if nargin < 4
        angleUnit = 'rad'; % Default angle units for output
        if nargin < 3
            timeGrid = 0:size(X_SC, 2); % Default timegrid
            if nargin < 2
                mu = 398600; % Earth gravitational parameter [km^3/s^2]
            end
        end
    end


% Convert from Cartesian to Osculating Keplerian elements
xKepl = CartTraj2KeplTraj(X_SC, mu, angleUnit);

% Plot

if nargin < 5
    figKeplElems = figure('WindowState', 'maximized', 'Name', 'Keplerian Elements');
else
    figure(figKeplElems);
    hold on;
end

cmap = turbo(6);

for i = 1:6
    subplot(2, 3, i);
    plot(timeGrid, xKepl(i, :), '-', 'Color', cmap(i, :), 'LineWidth', 1.05, ...
        'DisplayName', keplName{i});
    % Plot options
    grid minor
    ax_gca = gca;
    ax_gca.XAxisLocation = 'bottom';
    ax_gca.YAxisLocation = 'left';
    ax_gca.XMinorTick = 'on';
    ax_gca.YMinorTick = 'on';
    ax_gca.LineWidth = 1.04;
    axis padded
    % Labels 
    ylabel(strcat(keplName{i}, unitsPlot{i}));
    xlabel('Time [TU]')
    xlim([timeGrid(1), timeGrid(end)]);
end

sgtitle('Time Trajectory of Osculating Keplerian elements')
hold off;

%% LOCAL FUNCTION
    function xKepl = CartTraj2KeplTraj(xCart, mu, angleUnit)

        % Static allocation
        Nstates = size(xCart, 2);
        xKepl = zeros(6, Nstates);

        isParActive = not(isempty(gcp('nocreate')));

        if not(isParActive)
            % Cycle through trajectory (serial)
            for id = 1:Nstates
                xKepltemp = rv2kepl(xCart(:, id), mu);
                if ~isreal(xKepl)
                    error('One or more Keplerian elements not real!');
                end
                xKepl(:, id) = xKepltemp;

            end
        else
            % Cycle through trajectory (parallel)
            parfor id = 1:Nstates
                xKepltemp  = rv2kepl(xCart(:, id), mu);
                if ~isreal(xKepltemp)                 
                    error('One or more Keplerian elements not real!');
                end
                xKepl(:, id) = xKepltemp;
            end
        end

        if nargin > 3 && strcmpi(angleUnit, 'deg')
            xKepl(3:6, :) = rad2deg(xKepl(3:6, :));
        end

    end
end