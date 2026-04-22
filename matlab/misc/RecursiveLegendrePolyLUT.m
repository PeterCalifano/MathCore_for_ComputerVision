function Plm_LUT = RecursiveLegendrePolyLUT(gamma, dgamma, lMax) %#codegen
%% PROTOTYPE
% Plm_LUT = RecursiveLegendrePolyLUT(gamma, dgamma, lMax)
% -------------------------------------------------------------------------------------------------------------
%% DESCRIPTION
% ACHTUNG: IMPLEMENTATION NOT DEBUGGE. ISSUE IS PRESENT.
% Function computing a Look-Up table (lower diagonal array) containing all
% the Legendre Polynomials up to order lMax, for input gamma. dgamma must
% be the derivative of gamma for the generation to be consistent. Recursive
% formulas are used. Coding and verification specifically for SHE model. 
% The generation goes up to lMax + 2 because of the need of SHE model.
%
% Reference:
% 1) Fundamentals Of Astrodynamics And Applications, Vallado (chapters 8.6, 8.7)
%
% Note 1: the function is tailored for use with Exterior SHE model by
% providing gamma and dgamma (derivative of gamma) equal to sin(Lat) and 
% cos(Lat) as inputs. However, it is completely general provided that gamma
% and dgamma are consistent.
%
% Note 2: The memory storage required by the funciton may be optimized by
% replacing zeros() with sparse(). The gain is in the order of hundreds of 
% KB at most, though.
% -------------------------------------------------------------------------------------------------------------
%% INPUT
% gamma:  [1] Evaluation point of the Polynomial. 
%               For SHE, gamma = sin(Phi_latitude)
% dgamma: [1] Derivative of gamma. For SHE, dgamma = cos(Phi_latitude)
% lMax:   [1] Maximum degree of the Polynomial to compute
% -------------------------------------------------------------------------------------------------------------
%% OUTPUT
% Plm_LUT: [lMax+2, lMax+2] Lower triangular matrix of the Derived Legendre
%           Polynomials of gamma. Degree l along row, order m along column.
% -------------------------------------------------------------------------------------------------------------
%% CHANGELOG
% 30-11-2023    Pietro Califano    Documented from previous code 
%                                  (July 2023). Known issue: Plm not
%                                  correctly generated.
% -------------------------------------------------------------------------------------------------------------
%% DEPENDENCIES
% [-]
% -------------------------------------------------------------------------------------------------------------
%% Function code

% Output array is lower diagonal
% Plm_LUT: [(l+2)x(m+2)]

% Initialize LUT array
Plm_LUT = zeros(lMax+2, lMax+2);

% Recursion Start-up values
Plm_LUT(1, 1) = 1; % l=0, m=0
Plm_LUT(2, 1) = gamma;
Plm_LUT(2, 2) = dgamma;

% Recursion loop start
ldeg = 2;
for idl = 3:lMax + 2

    % Compute Legendre Polynomial of Zonal Harmonics (m=0)
    Plm_LUT(idl, 1) = 1/ldeg * ( (2*ldeg-1)*Plm_LUT(idl-1, 1) - (ldeg-1)*Plm_LUT(idl-2, 1));

    for idM = 2:idl-1 % for m that goes from 1 to l-1
        Plm_LUT(idl, idM) = Plm_LUT(idl-2, idM) + (2*ldeg-1)* dgamma * Plm_LUT(idl-1, idM-1);
    end

    % Compute Legendre Polynomial for m=l (diagonal of array
    Plm_LUT(idl, idl) = (2*ldeg-1)* dgamma * Plm_LUT(idl-1, idM-1);

    % Update value of l degree
    ldeg = ldeg + 1;
end

flagNoFailure = istril(Plm_LUT);
assert(flagNoFailure, 'Plm_LUT generation failed.')
end