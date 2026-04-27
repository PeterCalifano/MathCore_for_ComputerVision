function [xState, evalcounter, starter_fk] = PropagateRK4(RHS, tspan, state0) %#codegen
%% PROTOTYPE
% [state, evalcounter, starter_fk] = RK4(RHS, tspan, state0)
% -------------------------------------------------------------------------------------------------------------
%% DESCRIPTION
% Implementation of Runge-Kutta 4 fixed-step integrator. Butcher's table of
% the method taken from MSAS 2022/2023 course lectures at PoliMi.
% RHS must be a function handle declared as @(t, x) RHS(t, x, params).
% -------------------------------------------------------------------------------------------------------------
%% INPUT
% RHS
% tspan
% state0
% -------------------------------------------------------------------------------------------------------------
%% OUTPUT
% xState
% evalcounter
% starter_fk
% -------------------------------------------------------------------------------------------------------------
%% CHANGELOG
% 14-08-2023    Pietro Califano     RK4 code retrieved and adapted from
%                                   function coded during MSAS course 
%                                   2022/2023 at PoliMi. Verified.
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
end

%% RK4 Integrator scheme
%         RK4timer = tic;
evalcounter = 0;

time_vec = tspan; %et0:h:etf;
%         h = tspan(2) - tspan(1); % Fixed step derived from tspan

alfa = [1/2, 1/2, 1, 1]';

beta = [1/2, 0, 0, 0;
    0, 1/2, 0, 0;
    0, 0, 1, 0;
    1/6, 1/3, 1/3, 1/6];

xState = zeros(length(state0), length(time_vec));
% Initial condition assignment
xState(:, 1) = state0; %x0;

for timestep = 1:length(time_vec) - ifnotstarter

    x_k = xState(:, 1);
    t_k = time_vec(timestep);

    % 1st Predictor stage
    fk = RHS(t_k, x_k);
    evalcounter = evalcounter + 1;
    x_P1 = x_k + alfa(1)* h_vec(timestep) * fk;  %f(x_k, t_k)

    if nargout > 2
        if timestep == 1
            starter_fk = zeros(length(tspan), length(fk));
        end

        starter_fk(timestep, :) = fk;
    end

    % 2nd Predictor stage
    fP1 = RHS(t_k + alfa(1)*h_vec(timestep), x_P1);
    evalcounter = evalcounter + 1;
    x_P2 = x_k + alfa(2)*h_vec(timestep) *fP1;  %f(x_P1, t_k + alfa(1)*h)

    % 3rd Predictor stage
    fP2 = RHS(t_k + alfa(2)*h_vec(timestep), x_P2);
    evalcounter = evalcounter + 1;
    x_P3 = x_k + alfa(3)*h_vec(timestep) *fP2; %f(x_P2, t_k + alfa(2)*h)

    % 4th Corrector stage
    fP3 = RHS(t_k + alfa(3)*h_vec(timestep), x_P3);
    evalcounter = evalcounter + 1;

    x_next = x_k + alfa(4)*h_vec(timestep) * (beta(4, 1)*fk  + beta(4, 2)* fP1 + beta(4, 3)* fP2 + beta(4, 4)* fP3);

    xState(:, timestep + 1) = x_next;

end

if ifnotstarter == 0
    xState = xState(:, 1:end-1);
end

end

