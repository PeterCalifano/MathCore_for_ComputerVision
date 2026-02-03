function [dxStates] = PropagateKeplerianElems(dxStateKeplInit, ...
                                                dGravParam, ...
                                                dTimeGrid, ...
                                                kwargs) %#codegen
arguments
    dxStateKeplInit (6,1) {mustBeNumeric}
    dGravParam      (1,1) {mustBeGreaterThan(dGravParam, 0)}
    dTimeGrid       (1,:) {mustBeNumeric}
end
arguments
    kwargs.ui32MaxIter   (1,1) uint32 {mustBeNumeric} = 50
    kwargs.dRadAngleTol  (1,1) double {mustBeNumeric} = 1e-9; % [rad]
    kwargs.bConvert2Cart (1,1) logical {islogical} = false
end
%% PROTOTYPE
% [dxStates] = PropagateKeplerianElems(dxStateKeplInit, ...
%                                      dGravParam, ...
%                                      dDeltaTime, ...
%                                      kwargs) %#codegen
% -------------------------------------------------------------------------------------------------------------
%% DESCRIPTION
% General purpose Keplerian propagator in Classical Keplerian elements. It propagates the initial state over
% the specified time grid of time intervals. The function handles all types of keplerian orbits.
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
%     kwargs.bConvert2Cart (1,1) logical {islogical} = false
% end
% -------------------------------------------------------------------------------------------------------------
%% OUTPUT
% dxStates: [6, Ntimes] States in keplerian or cartesian coordinates on the specified timegrid
% -------------------------------------------------------------------------------------------------------------
%% CHANGELOG
% 19-10-2023    Pietro Califano     First prototype coded. VALIDATED.
% 01-05-2025    Pietro Califano     Major code refactoring for usage in Sequences dataset generation. Unit
%                                   test case (class): testPropagateKeplerianElems
% -------------------------------------------------------------------------------------------------------------
%% DEPENDENCIES
% [-]
% -------------------------------------------------------------------------------------------------------------

%% Function code

ui32TimegridSize = length(dTimeGrid);
ui32MaxIter = kwargs.ui32MaxIter;
dEPS_TIME = 1e-12;
dEcc = dxStateKeplInit(2);

assert(dEcc >= 0, 'MATLAB:assertion:failed', ...
    sprintf('ERROR: eccentricity value must be non-negative. Found value: %.6g', dEcc))

% Allocate output equal to input
dxStates = zeros(6, ui32TimegridSize);
dxStates(1:6, 1) = dxStateKeplInit; % Copy initial value

% Check timegrid validity (t0 < t1 < ... < tN)
if ui32TimegridSize > 2
    assert(issorted(dTimeGrid) && all(diff(dTimeGrid) > 0), 'Time grid must be strictly increasing');

elseif dTimeGrid(1) <= dEPS_TIME
    dxStates(1:6, 1) = dxStateKeplInit; % Copy initial value
    dxStates = dxStates(1:6, 1);
    return
end

% Compute common quantities
dMeanAngVel = sqrt(dGravParam / abs(dxStateKeplInit(1))^3 ); % Mean angular velocity of the orbit

% Propagation loop
for idT = 2:ui32TimegridSize
    
    % Initialize propagator
    dDeltaTime_ = dTimeGrid(idT) - dTimeGrid(idT-1);
    dxStates(1:6, idT) = dxStates(1:6, idT-1);

    if abs(dDeltaTime_) <= dEPS_TIME
        continue;
    end

    % Call single step propagator dispatcher
    dxStates(1:6, idT) = PropagateSingleStep(dDeltaTime_, ...
                                              dxStates(1:6, idT), ...
                                              dMeanAngVel, ...
                                              kwargs.dRadAngleTol, ...
                                              ui32MaxIter);
end

if kwargs.bConvert2Cart == true
    % Convert all states to cartesian coordinates
    for idS = 1:size(dxStates,2)
        dxStates(:,idS) = kepl2rv(dxStates(:,idS), dGravParam, "rad");
    end
end

end
%% LOCAL COMPONENTS
%%% Single dispatch function
function dKeplStateNext = PropagateSingleStep(dDeltaTime, dKeplStatePrev, dMeanAngVel, dRadAngleTol, ui32MaxIter)
arguments
    dDeltaTime
    dKeplStatePrev
    dMeanAngVel
    dRadAngleTol
    ui32MaxIter
end

% Cache useful values
persistent dSqrt1PlusEcc dEccPrev 

if isempty(dEccPrev) 
    % Initialize values
    dEccPrev = dKeplStatePrev(2);
    dSqrt1PlusEcc  = sqrt(1 + dKeplStatePrev(2));
elseif abs(dEccPrev - dKeplStatePrev(2)) > 1E-14
    % Eccentricity has changes, recompute values
    dSqrt1PlusEcc  = sqrt(1 + dKeplStatePrev(2));
end

dEccThreshold = 1E-10;
dKeplStateNext = dKeplStatePrev;

% Get state values
dSMA        = dKeplStatePrev(1);
dEcc        = dKeplStatePrev(2);
dTrueAnom   = dKeplStatePrev(6);

% Compute initial mean Anomaly based on the orbit type
if dEcc >= 1 + dEccThreshold
    dMeanAnom = dEcc * sinh(dTrueAnom) - dTrueAnom;
    i8OrbitType = int8(0); % Hyperbolic

elseif dEcc >= dEccThreshold && dEcc <= 1 - dEccThreshold
    dMeanAnom = dTrueAnom - dEcc * sin(dTrueAnom);
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
        dHyperbolicAnomNext = SolveHyperbolicAnomaly(dMeanAnomNext, dEcc, dRadAngleTol, ui32MaxIter);
        
        % Compute corresponding True Anomaly and assign output state
        dKeplStateNext(6) = 2 * atan2( dSqrt1PlusEcc * sinh(dHyperbolicAnomNext/2), sqrt(dEcc-1) * cosh(dHyperbolicAnomNext/2) );

    case 1
        %%% ELLIPTICAL ORBIT
        dEccentricAnomNext = SolveEllipticalAnomaly(dMeanAnomNext, dEcc, dRadAngleTol, ui32MaxIter);

        % Compute corresponding True Anomaly and assign output state
        dKeplStateNext(6) = 2 * atan2( dSqrt1PlusEcc * sin(dEccentricAnomNext/2), sqrt(1-dEcc) * cos(dEccentricAnomNext/2));
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
function dHyperbolicAnomNext = SolveHyperbolicAnomaly(dMeanAnom, dEcc, dRadAngleTol, ui32MaxIter)
arguments
    dMeanAnom
    dEcc
    dRadAngleTol
    ui32MaxIter
end

% Newton-Raphson loop for H
dErrorVal = 1.0; % Error initialization for zero-finding cycle
ui32IterCount = uint32(0);

% Initial guess generation
% dHtmp = asinh(dMeanAnom/dEcc);

if dEcc < 1.6
    if (dMeanAnom > -pi && Mnext < 0) || dMeanAnom > pi
        dHtmp = dMeanAnom - dEcc;
    else
        dHtmp = dMeanAnom + dEcc;
    end
else
    if dEcc < 3.6 && abs(dMeanAnom) > pi
        dHtmp = dMeanAnom - sign(dMeanAnom) * dEcc;
    else
        dHtmp = dMeanAnom/(dEcc-1);
    end
end

% Newton-Raphson cycle
while dErrorVal > dRadAngleTol && ui32IterCount <= ui32MaxIter

    % Execute iteration of Newton-Raphson to solve for H anomaly
    dHyperbolicAnomNext = dHtmp + ( dMeanAnom - dEcc*sinh(dHtmp) + dHtmp )./(dEcc * cosh(dHtmp) - 1);

    % Evaluate error for convergence check
    dErrorVal = abs(dHyperbolicAnomNext - dHtmp);
    dHtmp = dHyperbolicAnomNext;

    ui32IterCount = ui32IterCount + uint32(1);

end

if ui32IterCount == ui32MaxIter
    warning("Max number %d of iterations reached with error %6g.", ui32MaxIter, dErrorVal);
    return
end

end

%%% Elliptical case
function dEccentricAnomNext = SolveEllipticalAnomaly(dMeanAnom, dEcc, dRadAngleTol, ui32MaxIter)
arguments
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
while dErrorValue > dRadAngleTol && ui32IterCount <= ui32MaxIter

    % Execute iteration of Newton-Raphson to solve for H anomaly
    dEccentricAnomNext = dEtmp + ( dMeanAnom - dEtmp + dEcc * sin(dEtmp) )./(1 - dEcc*cos(dEtmp));

    % Evaluate error for convergence check
    dErrorValue = abs(dEccentricAnomNext - dEtmp);
    dEtmp = dEccentricAnomNext;

    ui32IterCount = ui32IterCount + uint32(1);

    if ui32IterCount == ui32MaxIter
        warning("Max number %d of iterations reached with error %6g.", ui32MaxIter, dErrorVal);
        return
    end
end

end



