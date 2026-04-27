function [xSol, feval, info] = zerofinderNewtonSym(F, x0, tol, maxIter)
%% PROTOTYPE
% [xSol, feval, info] = zerofinderNewtonSym(F, x0, tol, maxIter)
% -------------------------------------------------------------------------------------------------------------
%% DESCRIPTION
% Function implementing Newton's method for zero finding problems of the
% form: F = 0, where F is a multivariable vectorial function of dimension M 
% with domain X in R^N (variable x: [Nx1]). Symbolic toolbox is used to
% compute the jacobian of the function analytically. Function handles 
% definining F can be given as input and are automatically converted.
% -------------------------------------------------------------------------------------------------------------
%% INPUT
% fcn: [MATLAB object] Function handle or MATLAB symbolic function
% x0: [Nx1] Initial guess vector
% -------------------------------------------------------------------------------------------------------------
%% OUTPUT
% xSol: [Nx1] Solution vector at feval <= tol
% feval: [Mx1] Function value at solution vector
% -------------------------------------------------------------------------------------------------------------
%% CHANGELOG
% 13-08-2023    Pietro Califano    Retrieved and adapted from exercise 2
%                                  MSAS Assignement 1 (2022/2023). Correct
%                                  functioning verified with generic input
%                                  size function and variables.
% -------------------------------------------------------------------------------------------------------------
%% DEPENDENCIES
% MATLAB Symbolic Toolbox
% -------------------------------------------------------------------------------------------------------------
%% Function code
% Convert function handle to symbolic if needed
nDims = size(x0, 1);

if isa(F, "function_handle")
    x = sym('x', [nDims, 1]);
    Fsym = F(x);
else
    Fsym = F;
    % Get symbolic variables in expression
    x = transpose(symvar(F)); 
end

if nargin < 4
    if nargin < 3
        disp('Default error tolerance setting: 1e-8')
        tol = 1e-8;
    end
    disp('Default max iteration number: 100')
    maxIter = 100;
end

% Compute Symbolic jacobian
Jf = jacobian(Fsym, x);

% Initialize cycle counter
cycle_iter = 0;
fcnEval = 0;
% Initialize error with large value
err = 1e2 * ones(nDims, 1);

% Evaluate function at initial guess x0 vector to start iter
xk = x0;
Fk = eval(subs(Fsym, x, x0));

fcnEval = fcnEval + 1;

while sum(err > tol) ~= 0
    if cycle_iter == maxIter
        warning('Max number of iteration reached without convergence.')
        break;
    end

    % Execute iteration to compute xNew vector 
    % Equation: xNext = xk - [Jf(xk)]^-1 f(xk)
    xNew = eval(xk - subs(Jf, x, xk)\Fk);
    % Evaluate function at new xNew vector
    FxNew = eval(subs(Fsym, x, xNew));
    fcnEval = fcnEval + 1;

    % Update cycle counter
    cycle_iter = cycle_iter + 1;
    % Evaluate current error
    err = abs(FxNew);

    % Replace xk, fk with current solution
    xk = xNew;
    Fk = FxNew;
    
end

% Store solution
xSol = xk;
feval = Fk;

if nargout > 2
    info.fcnEval = fcnEval;
    info.cycle_iter = cycle_iter;
    info.errorEnd = err;
end

end