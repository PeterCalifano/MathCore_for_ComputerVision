function [xState, evalcounter, starter_fk] = PropagateRK2Heun(RHS, tspan, x0)
%% PROTOTYPE
% [xState, evalcounter, starter_fk] = RK2Heun(RHS, tspan, state0)
% -------------------------------------------------------------------------------------------------------------
%% DESCRIPTION
% Implementation of Runge-Kutta 2 fixed-step integrator. Butcher's table of
% the method taken from MSAS 2022/2023 course lectures at PoliMi. (TO
% REVIEW, NOT WORKING)
% -------------------------------------------------------------------------------------------------------------
%% INPUT
% RHS
% tspan
% x0
% -------------------------------------------------------------------------------------------------------------
%% OUTPUT
% xState
% evalcounter
% starter_fk
% -------------------------------------------------------------------------------------------------------------
%% CHANGELOG
% 14-08-2023    Pietro Califano     RK2 Heun code retrieved and adapted 
%                                   from function coded during MSAS course
%                                   2022/2023 at PoliMi
% -------------------------------------------------------------------------------------------------------------
%% DEPENDENCIES
% [-]
% -------------------------------------------------------------------------------------------------------------
%% Function code
ifnotstarter = 1;

h_fixed = tspan(2) - tspan(1);
h_vec = h_fixed*ones(length(tspan), 1);

% Modify last stepsize to solve up to the correct interval
if tspan(end) - tspan(end-1) ~= h_fixed
    h_vec(end) = tspan(end) - tspan(end-1);
end

if nargout > 2
    ifnotstarter = 0; % If the method is used as starter for multi-step,
    % fk is evaluated at all point of tspan. Otherwise last cycle is skipped.
    if timestep == 1
        starter_fk = zeros(length(tspan), length(fk));
    end
    starter_fk(timestep, :) = fk;
end

% ADD CHECK on tspan and h

% ADD CHECK on x0 and RHS size and adjust input to RHS appropriately to
% handle multivariable functions

% Static memory allocation
xk = nan(length(tspan)+1, length(x0));
evalcounter = 0;
% Initial condition
xk(1, :) = x0;

% h = tspan(2) - tspan(1);

for idt = 1:length(tspan)

    fk = RHS(xk(idt, :), tspan(idt));
    evalcounter = evalcounter + 1;

    if nargout > 2
        starter_fk(idt, :) = fk;
    end

    % Predictor step
    x_p = xk(idt, :) + h_vec(idt)*fk;
    % Corrector step to get x at k+1 step
    xk(idt + 1, :) = xk(idt, :) + h_vec(idt)* (0.5*fk + 0.5* RHS(x_p, tspan(idt)+h_vec(idt)));
    evalcounter = evalcounter + 1;
end
xState = xk(1:end-1, :);

end
