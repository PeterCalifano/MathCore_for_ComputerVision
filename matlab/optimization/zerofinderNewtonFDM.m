function [xSol, feval, info] = zerofinderNewtonFDM(F, x0, tol, maxIter) %#codegen
%% PROTOTYPE
% [xSol, feval, info] = zerofinderNewtonFDM(F, x0, tol, maxIter)
% -------------------------------------------------------------------------------------------------------------
%% DESCRIPTION
% Function implementing Newton's method for zero finding problems of the
% form: F = 0, where F is a multivariable vectorial function of dimension M 
% with domain X in R^N (variable x: [Nx1]). Finite differencing (Forward)
% is used to estimate the jacobian of the function.
% -------------------------------------------------------------------------------------------------------------
%% INPUT
% fcn: [MATLAB object] Function handle or MATLAB function
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
% [-]
% -------------------------------------------------------------------------------------------------------------
%% Function code
nDims = size(x0, 1);

if nargin < 4
    if nargin < 3
        disp('Default error tolerance setting: 1e-8')
        tol = 1e-8;
    end
    disp('Default max iteration number: 100')
    maxIter = 100;
end

% Initialize cycle counter
cycle_iter = 0;
fcnEval = 0;
% Initialize error with large value
err = 1e2 * ones(nDims, 1);

% Evaluate function at initial guess x0 vector to start iter
xk = x0;
Fk = F(x0);

Mdims = size(Fk, 1);

while sum(err > tol) ~= 0
    if cycle_iter == maxIter
        warning('Max number of iteration reached without convergence.')
        break;
    end

    % Initialize Jacobian matrix and perturbation vector
    Jf = zeros(Mdims, nDims);
    dX = zeros(nDims, 1);

    for id = 1:nDims
        % Compute perturbation vector dX
        dX(id) = sqrt(eps)*max(1, abs(xk(id)));
        % Evaluate function at perturbed state
        FxPert = F(xk + dX);
        fcnEval = fcnEval + 1;
        % Evaluate Jacobian entry
        Jf(:, id) = (FxPert - Fk)./dX(id);        
    end

    % Execute iteration to compute xNew vector 
    % Equation: xNext = xk - [Jf(xk)]^-1 f(xk)
    xNew = xk - Jf\Fk;
    % Evaluate function at new xNew vector
    FxNew = F(xNew);
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