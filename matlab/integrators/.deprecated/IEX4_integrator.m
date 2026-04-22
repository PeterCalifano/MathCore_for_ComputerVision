function [x, evalcounter] = IEX4_integrator(RHS, tspan, x0) %#codegen
%% PROTOTYPE
% -------------------------------------------------------------------------------------------------------------
%% DESCRIPTION
% -------------------------------------------------------------------------------------------------------------
%% INPUT
% -------------------------------------------------------------------------------------------------------------
%% OUTPUT
% -------------------------------------------------------------------------------------------------------------
%% CHANGELOG
% Date, User, brief summary of the modification
% -------------------------------------------------------------------------------------------------------------
%% DEPENDENCIES
% -------------------------------------------------------------------------------------------------------------
%% Function code

IEX4_timer = tic;
% Integrator parameters and counters
alfa = [-1/6, 4, -27/2, 32/3]';
h_fixed = tspan(2) - tspan(1);
h = [1, 1/2, 1/3, 1/4] * h_fixed;

options = optimoptions('fsolve', 'Display', 'off', 'FunctionTolerance', 1e-12, 'UseParallel', false);

% Static memory allocation
x = nan(length(tspan), length(x0));
x(1, :) = x0;

evalcounter = 0;

for timestep = 1:length(tspan)-1

    x_k = x(timestep, :)';
    t_k = tspan(timestep);

    % 1st Predictor stage
    % Note: stages evaluations are independent: can be parallelized

    % Implicit non-linear equation to solve
%     f = @(xP1) xP1 - x_k - h(1) * RHS(t_k + h(1), xP1);
    [xP1, ~, ~, output] = fsolve(@(xP1) xP1 - x_k - h(1) * RHS(t_k + h(1), xP1), x_k, options);  
    evalcounter = evalcounter + output.funcCount;

    % 2nd Predictor stage
%     f = @(xP2_mid) xP2_mid - x_k - h(2) * RHS(t_k + h(2), xP2_mid);
    [xP2_a, ~, ~, output] = fsolve(@(xP2_a) xP2_a - x_k - h(2) * RHS(t_k + h(2), xP2_a), x_k, options);
    evalcounter = evalcounter + output.funcCount;

%     f = @(xP2) xP2 - xP2_a - h(2) * RHS(t_k + h(1), xP2); % Nbp_sys(t_k + h(1), xP2, bodies, Rframe, Oframe);
    [xP2, ~, ~, output] = fsolve(@(xP2) xP2 - xP2_a - h(2) * RHS(t_k + h(1), xP2), xP2_a, options);  
    evalcounter = evalcounter + output.funcCount;

    % 3th Predictor stage

%     f = @(xP3_a) xP3_a - x_k - h(3) * RHS(t_k + h(3), xP3_a);
    [xP3_a, ~, ~, output] = fsolve(@(xP3_a) xP3_a - x_k - h(3) * RHS(t_k + h(3), xP3_a), x_k, options); 
    evalcounter = evalcounter + output.funcCount;

%     f = @(xP3_b) xP3_b - xP3_a - h(3) * RHS(t_k + h(3), xP3_b);
    [xP3_b, ~, ~, output] = fsolve(@(xP3_b) xP3_b - xP3_a - h(3) * RHS(t_k + h(3), xP3_b), xP3_a, options);
    evalcounter = evalcounter + output.funcCount;

%     f = @(xP3) xP3 - xP3_b - h(3) * RHS(t_k + h(1), xP3); % Nbp_sys(t_k + h(1), xP3, bodies, Rframe, Oframe);
    [xP3, ~, ~, output] = fsolve(@(xP3) xP3 - xP3_b - h(3) * RHS(t_k + h(1), xP3), xP3_b, options);
    evalcounter = evalcounter + output.funcCount;

    % 4th Predictor stage

%     f = @(xP4_a) xP4_a - x_k - h(4) * RHS(t_k + h(4), xP4_a);
    [xP4_a, ~, ~, output]  = fsolve(@(xP4_a) xP4_a - x_k - h(4) * RHS(t_k + h(4), xP4_a), x_k, options); 
    evalcounter = evalcounter + output.funcCount;

%     f = @(xP4_b) xP4_b - xP3_a - h(4) * RHS(t_k + h(4), xP4_b);
    [xP4_b, ~, ~, output] = fsolve(@(xP4_b) xP4_b - xP4_a - h(4) * RHS(t_k + h(4), xP4_b), xP4_a, options);
    evalcounter = evalcounter + output.funcCount;

%     f = @(xP4_c) xP4_c - xP4_b - h(4) * RHS(t_k + h(4), xP4_c);
    [xP4_c, ~, ~, output] = fsolve(@(xP4_c) xP4_c - xP4_b - h(4) * RHS(t_k + h(4), xP4_c), xP4_b, options);
    evalcounter = evalcounter + output.funcCount;

%     f = @(xP4) xP4 - xP4_c - h(4) * RHS(t_k + h(1), xP4); % Nbp_sys(t_k + h(1), xP4, bodies, Rframe, Oframe);
    [xP4, ~, ~, output] = fsolve(@(xP4) xP4 - xP4_c - h(4) * RHS(t_k + h(1), xP4), xP4_c, options);
    evalcounter = evalcounter + output.funcCount;

    % 5th Corrector stage
    x_next = alfa(1)*xP1 + alfa(2)*xP2 + alfa(3)*xP3 + alfa(4)*xP4;

    x(timestep + 1, :) = x_next;
end

toc(IEX4_timer)

end