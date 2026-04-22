function [time_vec, x, evalcounter] = BIn_theta_integrator(RHS, tspan, x0, order, h, theta) %#codegen
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
% Input checks
 
% Static memory allocation
evalcounter = 0;
xk = nan(length(tspan), length(x0)); % Allocation memory for solution
xk(1, :) = x0;

switch order

    case 1 % BI1_theta
        
        for idt = 1:length(tspan)-1

            tk = tspan(idt);
            fk = RHS(xk(idt, :), tk);
            evalcounter = evalcounter + 1;
            
            % Compute reference point with Forward Euler to be used as guess solution
            xleft_ktheta = xk(idt, :) + theta*h *fk;

            opts = optimoptions('fsolve', 'Display', 'off', 'FunctionTolerance', fsolve_tol);
            % Function handle definition at each step may be avoided by
            % providing directly the function to fsolve

%             f = @(x_next) x_next - (1-theta)*h*RHS(x_next, tk+h) - xleft_ktheta;
            [xk(idt+1), ~, exitflag, output] = fsolve(@(x_next) x_next - (1-theta)*h*RHS(x_next, tk+h) - xleft_ktheta, xleft_ktheta, opts);
            evalcounter = evalcounter + output.funcCount;

            %             evalcounter = evalcounter + 1;

        end
        

    case 4
        err('BRK4 not yet implemented')

end

time_vec = tspan;
x = xk;

end