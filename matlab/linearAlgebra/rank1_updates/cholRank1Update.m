function dLnew = cholRank1Update(dL, dU, dBeta, dAlpha) %#codegen
arguments
dL     (:, :) double {ismatrix, mustBeNumeric}
dU     (:, :) double {ismatrix, mustBeNumeric}
dBeta  (1, 1) double {mustBeNumeric, isscalar} = 1.0
dAlpha (1, 1) double {mustBeNumeric, isscalar} = 1.0
end
%% PROTOTYPE
% dLnew = cholRank1Update(dL, dU, dBeta, dAlpha)
% -------------------------------------------------------------------------------------------------------------
%% DESCRIPTION
% Function computing the Cholesky Rank 1 Update on the columns of matrix
% dU using "Triangular Rank-One update". Given a matrix Sigma with
% Cholesky (lower) factor L, a vector v and two real scalar alpha and beta,
% the algorithm computes the updated Cholesky factor Lnew such that:
% Lnew = chol(SigmaNew), where SigmaNew = alpha*Sigma + beta * v*v';
% REFERENCE:
% 1) A More Efficient Rank-one Covariance Matrix Update for Evolution 
% Strategies, Krause, Igel, 2015
% -------------------------------------------------------------------------------------------------------------
%% INPUT
% dL     (:, :)               Original Cholesky factor L (can be either lower or upper, automatically handled)
% dU     (:, :)               Update Matrix with update vectors v as columns
% dBeta  (1, 1) double = 1.0  Scaling coefficient of Update matrix
% dAlpha (1, 1) double = 1.0  Scaling coefficient of cholesky factor
% -------------------------------------------------------------------------------------------------------------
%% OUTPUT
% dLnew
% -------------------------------------------------------------------------------------------------------------
%% CHANGELOG
% 17-09-2023    Pietro Califano     First prototype coded. Debug needed:
%                                   slightly output with respect to chol().
%                                   Error is in the order of 10^-2 on the
%                                   diagonal whereas other elements are
%                                   below machine precision in double.
% 04-02-2024    Pietro Califano     Code validated against MATLAB cholUpdate
% -------------------------------------------------------------------------------------------------------------
%% DEPENDENCIES
% [-]
% -------------------------------------------------------------------------------------------------------------
%% Function code

% Check if input Cholesky factor is lower triangular, else transpose
UPPER = logical(istriu(dL));

if UPPER
    dL = transpose(dL);
end

% Check if scaling is present
if nargin < 4
    dAlpha = 1.0;
    if nargin < 3
        dBeta = 1.0;
    end
end

% Determine number of columns of U
NrowColU = uint16(size(dU));
NrowU = NrowColU(1);
NcolU = NrowColU(2);

ui16N = uint16(size(dL, 1));
assert(NrowU == ui16N, 'Dimension of rows of L and U does not match!')

% Allocate output matrix
dLnew = zeros(size(dL));

% Input: Cholesky factor dL of size N
alphaSqrt = sqrt(dAlpha);

for idCol = 1:NcolU

    % Triangular Rank1 Update
    v = dU(:, idCol);
    
    if all(v == 0)
        continue;
    end

    w = v; % Initialize rank-1 vector
    b = 1.0; % Initialize b coefficient

    for idj = 1:ui16N
        
        dSqrtEntry = dAlpha * dL(idj, idj)^2 + (dBeta/b) * w(idj)^2;
        assert(dSqrtEntry >= 0, 'ERROR: detected negative value under square root while updating Cholesky factor.')

        % Compute entry of new Cholesky factor
        dLnew(idj, idj) = sqrt( dSqrtEntry );
        
        gamma = dAlpha * dL(idj, idj)^2 * b + dBeta * w(idj)^2;

        for idk = idj+1:ui16N     
            % Update kth value in "update" vector
            w(idk) = w(idk) - (w(idj)/dL(idj, idj)) * alphaSqrt * dL(idk, idj);
            % Compute entry of new Cholesky factor
            dLnew(idk, idj) = alphaSqrt * (dLnew(idj, idj)/dL(idj, idj)) * dL(idk, idj) +...
                (dLnew(idj, idj) * dBeta * w(idj)) * w(idk)/gamma;         
        end

        % Update auxiliary variable b
        b = b + dBeta * w(idj)^2/(dAlpha * dL(idj, idj)^2);

    end
end

% If input was Upper, provide output as Upper Cholesky factor
if UPPER
    dLnew = transpose(dLnew);
end
