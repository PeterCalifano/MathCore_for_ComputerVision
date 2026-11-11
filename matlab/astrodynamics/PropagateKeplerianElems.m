function [dxStates] = PropagateKeplerianElems(dxStateKeplInit, ...
                                                dGravParam, ...
                                                dTimeGrid, ...
                                                kwargs) %#codegen
arguments (Input)
    dxStateKeplInit (6,1) {mustBeNumeric}
    dGravParam      (1,1) {mustBeGreaterThan(dGravParam, 0)}
    dTimeGrid       (1,:) {mustBeNumeric}
end
arguments (Input)
    kwargs.ui32MaxIter   (1,1) uint32 {mustBeNumeric} = 50
    kwargs.dRadAngleTol  (1,1) double {mustBeNumeric} = 1e-9; % [rad]
    kwargs.bConvert2Cart (1,1) logical = false
end
%% SIGNATURE
% [dxStates] = PropagateKeplerianElems(dxStateKeplInit, ...
%                                      dGravParam, ...
%                                      dTimeGrid, ...
%                                      kwargs) %#codegen
% -------------------------------------------------------------------------------------------------------------
%% DESCRIPTION
% General purpose Keplerian propagator in Classical Keplerian elements. It propagates the initial state over
% the specified time grid of time intervals. The function handles all types of keplerian orbits. The state
% carries the true anomaly; each step converts it to the eccentric (elliptic) or hyperbolic anomaly, advances
% the mean anomaly, solves Kepler's equation by Newton-Raphson and converts back to the true anomaly.
% 1) Fundamentals of Astrodynamics and Applications - D. Vallado. Section
%    2.2, Kepler's Problem. Algorithms 3 and 4.
% -------------------------------------------------------------------------------------------------------------
%% INPUT
% arguments
%     dxStateKeplInit (6,1) {mustBeNumeric}
%     dGravParam      (1,1) {mustBeGreaterThan(dGravParam, 0)}
%     dTimeGrid       (1,:) {mustBeNumeric}
% end
% arguments
%     kwargs.ui32MaxIter   (1,1) uint32 {mustBeNumeric} = 50
%     kwargs.dRadAngleTol  (1,1) double {mustBeNumeric} = 1e-9; % [rad]
%     kwargs.bConvert2Cart (1,1) logical = false
% end
% -------------------------------------------------------------------------------------------------------------
%% OUTPUT
% dxStates:
%   - If dTimeGrid is a scalar horizon dt > 0: [6,2] with states at t=[0,dt].
%   - If dTimeGrid is scalar dt <= 0: [6,1] initial state only.
%   - If dTimeGrid is a strictly increasing vector of absolute times: [6,Ntimes].
%   States are in keplerian or cartesian coordinates based on kwargs.bConvert2Cart.
% -------------------------------------------------------------------------------------------------------------
%% CHANGELOG
% 19-10-2023    Pietro Califano     First prototype coded. VALIDATED.
% 01-05-2025    Pietro Califano     Major code refactoring for usage in Sequences dataset generation. Unit
%                                   test case (class): testPropagateKeplerianElems
% 14-02-2026   Pietro Califano      Update implementation and change of input interface, fix bugs
% 28-10-2026   Pietro Califano, Claude Fable 5.1   Fix: Kepler's equation was evaluated with the true anomaly
%                                   in place of the eccentric/hyperbolic anomaly (wrong for every non-circular
%                                   orbit); Newton solvers warn whenever the tolerance is not reached.
% 09-11-2026   Pietro Califano, Claude Opus 5.5   Fix: the hyperbolic solver starts from an upper bound of the
%                                   root, min(asinh(|M|/(e-1)), (6|M|/e)^(1/3)) with the sign of M (the former
%                                   start failed from |M| of about 51); a non-finite iterate raises an error.
% 11-11-2026   Pietro Califano, Claude Opus 5.5   Readability pass.
% -------------------------------------------------------------------------------------------------------------
%% DEPENDENCIES
% [-]
% -------------------------------------------------------------------------------------------------------------

%% Function code

dTimeGrid = double(dTimeGrid(:)');
ui32NumInputTimes = uint32(numel(dTimeGrid));
ui32MaxIter = kwargs.ui32MaxIter;
dEpsTime = 1e-12;
dEcc = dxStateKeplInit(2);

assert(dEcc >= 0, 'MATLAB:assertion:failed', ...
    sprintf('ERROR: eccentricity value must be non-negative. Found value: %.6g', dEcc))

if ui32NumInputTimes == 1
    dPropagationHorizon = dTimeGrid(1);
    assert(dPropagationHorizon >= 0.0, "Scalar propagation horizon must be non-negative.");

    if dPropagationHorizon <= dEpsTime
        dxStates = dxStateKeplInit;
        if kwargs.bConvert2Cart == true
            dxStates = kepl2rv(dxStates, dGravParam, "rad");
        end
        return;
    end

    dEvalTimeGrid = [0.0, dPropagationHorizon];
else
    assert(issorted(dTimeGrid) && all(diff(dTimeGrid) > 0), ...
        'Time grid must be strictly increasing when passed as a vector of absolute times.');
    dEvalTimeGrid = dTimeGrid;
end

ui32TimegridSize = uint32(numel(dEvalTimeGrid));

% Allocate output equal to input
dxStates = zeros(6, double(ui32TimegridSize));
dxStates(1:6, 1) = dxStateKeplInit; % Copy initial value

% Compute common quantities
dMeanAngVel = sqrt(dGravParam / abs(dxStateKeplInit(1))^3); % Mean angular velocity of the orbit

% Propagation loop
for idT = 2:double(ui32TimegridSize)

    % Initialize propagator
    dDeltaTime_ = dEvalTimeGrid(idT) - dEvalTimeGrid(idT - 1);
    dxStates(1:6, idT) = dxStates(1:6, idT - 1);

    if abs(dDeltaTime_) <= dEpsTime
        continue;
    end

    % Call single step propagator dispatcher
    dxStates(1:6, idT) = PropagateSingleStep_(dDeltaTime_, ...
                                               dxStates(1:6, idT), ...
                                               dMeanAngVel, ...
                                               kwargs.dRadAngleTol, ...
                                               ui32MaxIter);
end

if kwargs.bConvert2Cart == true
    % Convert all states to cartesian coordinates
    for idS = 1:size(dxStates, 2)
        dxStates(:, idS) = kepl2rv(dxStates(:, idS), dGravParam, "rad");
    end
end

end
%% LOCAL COMPONENTS
%%% Single dispatch function
function dKeplStateNext = PropagateSingleStep_(dDeltaTime, dKeplStatePrev, dMeanAngVel, dRadAngleTol, ui32MaxIter)
arguments (Input)
    dDeltaTime
    dKeplStatePrev
    dMeanAngVel
    dRadAngleTol
    ui32MaxIter
end

% Cache useful values
persistent dSqrt1PlusEcc dEccPrev

if isempty(dEccPrev)
    % Initialize cached values (for faster computation)
    dEccPrev = dKeplStatePrev(2);
    dSqrt1PlusEcc = sqrt(1 + dKeplStatePrev(2));

elseif abs(dEccPrev - dKeplStatePrev(2)) > 1E-14

    % Eccentricity has changed, recompute values
    dEccPrev = dKeplStatePrev(2);
    dSqrt1PlusEcc = sqrt(1 + dKeplStatePrev(2));
end

dEccThreshold = 1E-10;
dKeplStateNext = dKeplStatePrev;

% Get state values
dEcc = dKeplStatePrev(2);
dTrueAnom = dKeplStatePrev(6);

% Compute initial mean Anomaly based on the orbit type. Kepler's equation takes the eccentric anomaly E
% (elliptic, M = E - e sin E) or the hyperbolic anomaly H (hyperbolic, M = e sinh H - H), never the true
% anomaly nu carried by the state: E and H come from the half-angle relations
% tan(E/2) = sqrt((1-e)/(1+e)) tan(nu/2) and sinh(H) = sqrt(e^2-1) sin(nu) / (1 + e cos(nu)) (Vallado, 2.2).
if dEcc >= 1 + dEccThreshold
    dHyperbolicAnom = asinh(sqrt(dEcc^2 - 1) * sin(dTrueAnom) / (1 + dEcc * cos(dTrueAnom)));
    dMeanAnom = dEcc * sinh(dHyperbolicAnom) - dHyperbolicAnom;
    i8OrbitType = int8(0); % Hyperbolic

elseif dEcc >= dEccThreshold && dEcc <= 1 - dEccThreshold
    dEccentricAnom = 2 * atan2(sqrt(1 - dEcc) * sin(dTrueAnom / 2), dSqrt1PlusEcc * cos(dTrueAnom / 2));
    dMeanAnom = dEccentricAnom - dEcc * sin(dEccentricAnom);
    i8OrbitType = int8(1); % Elliptical

elseif dEcc <= 0 + dEccThreshold
    dMeanAnom = dTrueAnom;
    i8OrbitType = int8(2); % Circular
else
    error('Eccentricity very close to 0 or 1. Case not handled.')
end

% Update mean anomaly for propagation
dMeanAnomNext = dMeanAnom + dMeanAngVel * dDeltaTime;

% Dispatch call to orbit-specific solver
switch i8OrbitType
    case 0
        %%% HYPERBOLIC ORBIT
        dHyperbolicAnomNext = SolveHyperbolicAnomaly_(dMeanAnomNext, dEcc, dRadAngleTol, ui32MaxIter);

        % Compute corresponding True Anomaly and assign output state
        dKeplStateNext(6) = 2 * atan2(dSqrt1PlusEcc * sinh(dHyperbolicAnomNext / 2), ...
                                      sqrt(dEcc - 1) * cosh(dHyperbolicAnomNext / 2));

    case 1
        %%% ELLIPTICAL ORBIT
        dEccentricAnomNext = SolveEllipticalAnomaly_(dMeanAnomNext, dEcc, dRadAngleTol, ui32MaxIter);

        % Compute corresponding True Anomaly and assign output state
        dKeplStateNext(6) = 2 * atan2(dSqrt1PlusEcc * sin(dEccentricAnomNext / 2), ...
                                      sqrt(1 - dEcc) * cos(dEccentricAnomNext / 2));
        dKeplStateNext(6) = wrapTo2Pi(dKeplStateNext(6));

    case 2
        %%% CIRCULAR ORBIT
        % Circular: true anomaly equals mean anomaly
        dKeplStateNext(6) = dMeanAnomNext;
        dKeplStateNext(6) = wrapTo2Pi(dKeplStateNext(6));

    otherwise
        error('Orbit type with eccentricity value %.6g not handled by current implementation!', dEcc)
end

end

%%% Hyperbolic case
function dHyperbolicAnomNext = SolveHyperbolicAnomaly_(dMeanAnom, dEcc, dRadAngleTol, ui32MaxIter)
arguments (Input)
    dMeanAnom
    dEcc
    dRadAngleTol
    ui32MaxIter
end

% Newton-Raphson loop for H
dErrorVal = 1.0; % Error initialization for zero-finding cycle
ui32IterCount = uint32(0);

% Initial guess: an upper bound of the root (in magnitude), so that Newton-Raphson on the convex, increasing
% f(H) = e sinh(H) - H - M (H >= 0) converges monotonically and never overshoots. For H >= 0,
% e sinh(H) - H >= (e - 1) sinh(H) (since sinh(H) >= H) and e sinh(H) - H >= e H^3 / 6, so both
% asinh(|M| / (e - 1)) and (6 |M| / e)^(1/3) bound the root from above. The negative branch follows from
% the odd symmetry H(-M) = -H(M). The former start (M -/+ e or M / (e - 1)) stopped unconverged from
% |M| of about 51 and overflowed sinh above |M| of about 710; asinh(M / e) alone overflows sinh for
% near-parabolic orbits at small |M| (first step about M / (e - 1)).
dAbsMeanAnom = abs(dMeanAnom);
dHtmp = sign(dMeanAnom) * min(asinh(dAbsMeanAnom / (dEcc - 1)), nthroot(6 * dAbsMeanAnom / dEcc, 3));

% Newton-Raphson cycle
dHyperbolicAnomNext = dHtmp;
while dErrorVal > dRadAngleTol && ui32IterCount < ui32MaxIter

    % Execute iteration of Newton-Raphson to solve for H anomaly
    dHyperbolicAnomNext = dHtmp + (dMeanAnom - dEcc * sinh(dHtmp) + dHtmp) ./ (dEcc * cosh(dHtmp) - 1);

    % Stop on a non-finite iterate (non-finite mean anomaly or sinh overflow) instead of returning NaN
    if ~isfinite(dHyperbolicAnomNext)
        error('PropagateKeplerianElems:nonFiniteAnomaly', ...
            'Hyperbolic anomaly iterate is not finite (M = %.6g, e = %.6g).', dMeanAnom, dEcc);
    end

    % Evaluate error for convergence check
    dErrorVal = abs(dHyperbolicAnomNext - dHtmp);
    dHtmp = dHyperbolicAnomNext;

    ui32IterCount = ui32IterCount + uint32(1);

end

if dErrorVal > dRadAngleTol
    warning("Max number %d of iterations reached with error %6g.", ui32MaxIter, dErrorVal);
end

end

%%% Elliptical case
function dEccentricAnomNext = SolveEllipticalAnomaly_(dMeanAnom, dEcc, dRadAngleTol, ui32MaxIter)
arguments (Input)
    dMeanAnom
    dEcc
    dRadAngleTol
    ui32MaxIter
end

ui32IterCount = uint32(0);
dErrorValue = 1; % Error for Newton-Raphson cycle

% Initial guess generation
% dEtmp = dMeanAnom + sign(sin(dMeanAnom))*0.8*dEcc;

if (dMeanAnom > -pi && dMeanAnom < 0) || dMeanAnom > pi
    dEtmp = dMeanAnom - dEcc;
else
    dEtmp = dMeanAnom + dEcc;
end

% Newton-Raphson cycle
dEccentricAnomNext = dEtmp;
while dErrorValue > dRadAngleTol && ui32IterCount < ui32MaxIter

    % Execute iteration of Newton-Raphson to solve for E anomaly
    dEccentricAnomNext = dEtmp + (dMeanAnom - dEtmp + dEcc * sin(dEtmp)) ./ (1 - dEcc * cos(dEtmp));

    % Evaluate error for convergence check
    dErrorValue = abs(dEccentricAnomNext - dEtmp);
    dEtmp = dEccentricAnomNext;

    ui32IterCount = ui32IterCount + uint32(1);

end

if dErrorValue > dRadAngleTol
    warning("Max number %d of iterations reached with error %6g.", ui32MaxIter, dErrorValue);
end

end

