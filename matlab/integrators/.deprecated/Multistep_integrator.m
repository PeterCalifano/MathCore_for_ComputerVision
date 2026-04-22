function [time_vec, x, evalcounter] = Multistep_integrator(RHS, tspan, x0, class, order)
%% PROTOTYPE
% [time_vec, x, evalcounter] = Multistep_integrator(RHS, tspan, x0, class, order)
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
%% DEVNOTE: DEBUG STILL ON GOING. VALIDATION TO DO

%% Function code
% ACHTUNG: RHS MUST BE SPECIFIED WITH VARS AS (x, t)
% Check inputs

% Scheme selection
if strcmpi('AB', class)
    class = 0;
elseif strcmpi('AM', class)
    class = 1;
elseif strcmpi('ABM', class)
    class = 2;
elseif strcmpi('BDF', class)
    class = 3;
else
    err('Specified method name not found or available')
end

% Schemes parameters and counters
h = tspan(2) - tspan(1);
starting_idt = 1;
starter_counter = 0;
evalcounter = 0;

% Static memory allocation

xk = nan(length(tspan), length(x0)); % Allocation memory for solution

fk = nan(length(tspan), length(RHS(x0, tspan(1))));
evalcounter = evalcounter + 1;

xk(1, :) = x0;

if class == 0 || class == 1 || class == 2 % AM/AB/ABM

    switch order
        case 3  
            [x_starter, starter_counter, starter_fk] = RKN_integrator(RHS, tspan(1:3), x0, 4);
            xk(1:3, :) = x_starter;
            fk(1:3, :) = starter_fk;

            if class == 1
                starting_idt = 2;
            else
                starting_idt = 3;
            end
            
        case 4  
            [x_starter, starter_counter, starter_fk] = RKN_integrator(RHS, tspan(1:4), x0, 4);
            xk(1:3, :) = x_starter;
            fk(1:3, :) = starter_fk;
            starting_idt = 4;
    end

elseif class == 3 % BDF

    switch order
        case 3 % BDF3
            [x_starter, starter_counter, starter_fk] = RKN_integrator(RHS, tspan(1:3), x0, 4);
            xk(1:3, :) = x_starter;
            fk(1:3, :) = starter_fk;
            
            starting_idt = 3;
    end
end

evalcounter = evalcounter + starter_counter;

switch class

    case 0 % 'AB'

        for idt = starting_idt:length(tspan)-1

            tk = tspan(idt);          
            fk(idt, :) = RHS(xk(idt, :), tk);
            evalcounter = evalcounter + 1;

            if order == 3
                % AB3
                xk(idt+1, :) = xk(idt, :) + 1/12 * h*(23*fk(idt, :) - 16*fk(idt-1, :) + 5*fk(idt-2, :));

            elseif order == 4
                % AB4
                xk(idt+1, :) = xk(idt, :) + 1/24 * h*(55*fk(idt, :) - 59*fk(idt-1, :) + 37*fk(idt-2) - 9*fk(idt-3, :));

            end
        end

    case 1 %'AM'

        opts = optimoptions('fsolve', 'Display', 'off', 'FunctionTolerance', fsolve_tol);

        for idt = starting_idt:length(tspan)-1

            tk = tspan(idt);

            fk(idt, :) = RHS(xk(idt, :), tk);
            evalcounter = evalcounter + 1;

            if order == 3
                % AM3
                [xk(idt+1, :), ~, exitflag, output] = fsolve(@(x_next) x_next - xk(idt, :) - (1/12 * h* ( 5*(RHS(x_next, tk+h))' + 8*fk(idt, :) - fk(idt-1, :) ) ) , xk(idt, :), opts);
                evalcounter = evalcounter + output.funcCount;

            elseif order == 4
                % AM4
                [xk(idt+1, :), ~, exitflag, output] = fsolve(@(x_next) x_next - xk(idt, :) - (1/24 * h*(9*(RHS(x_next, tk+h))' + 19*fk(idt, :) - 5*fk(idt-1, :) + fk(idt-2, :))) , xk(idt, :), opts);
                evalcounter = evalcounter + output.funcCount;
            end
        end

    case 2 %'ABM'
        for idt = starting_idt:length(tspan)-1

            tk = tspan(idt);
            % COMPUTE idt-1 and idt-2 with RK3 or 4
            fk(idt, :) = RHS(xk(idt, :), tk);
            evalcounter = evalcounter + 1;

            if order == 3
                % ABM3
                % Predictor stage - AB3
                xp = xk(idt, :) + 1/12 * h*(23*fk(idt, :) - 16*fk(idt-1, :) + 5*fk(idt-2, :));
                % Corrector stage - AM3
                xk(idt+1, :) = xk(idt, :) + 1/12 * h*(5*(RHS(xp, tk+h))' + 8*fk(idt, :) - fk(idt-1, :));
                evalcounter = evalcounter + 1;
            elseif order == 4
                % ABM4
                % Predictor stage - AB4
                xp = xk(idt, :) + 1/24 * h*(55*fk(idt, :) - 59*fk(idt-1, :) + 37*fk(idt-2) - 9*fk(idt-3, :));
                % Corrector stage - AM4
                xk(idt+1, :) = xk(idt, :) + 1/24 * h*(9*(RHS(xp, tk+h))' + 19*fk(idt, :) - 5*fk(idt-1, :) + fk(idt-2, :));
                evalcounter = evalcounter + 1;
            end
        end

    case 3 %'BDF'

        % BDF3: xk(idt+1, :) = 18/11 * xk(idt, :) - 9/11 * xk(idt-1, :) + 2/11 * xk(idt-2, :) + 6/11* h * RHS(x_next, tk+h);

        opts = optimoptions('fsolve', 'Display', 'off', 'FunctionTolerance', fsolve_tol);

        for idt = starting_idt:length(tspan)-1

            tk = tspan(idt);

            [xk(idt+1, :), ~, exitflag, output] = fsolve(@(x_next) x_next - (18/11 * xk(idt, :) - 9/11 * xk(idt-1, :) + 2/11 * xk(idt-2, :) + 6/11* h * (RHS(x_next, tk+h))') , xk(idt, :), opts);
            evalcounter = evalcounter + output.funcCount;
        end

end

time_vec = tspan;
x = xk;

