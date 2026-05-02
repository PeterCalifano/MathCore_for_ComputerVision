function [strFitStats, dChbvInterpVector] = checkFitChbvPoly(ui32PolyDeg, ...
                                                            dInterpDomain, ...
                                                            dChbvCoeffs, ...
                                                            dDataMatrix, ...
                                                            dDomainLB, ...
                                                            dDomainUB, ...
                                                            bIS_ATT_QUAT, ...
                                                            dSwitchIntervals, ...
                                                            bEnableErrorThrow, ...
                                                            dPercRelErrorTol) %#codegen
arguments
    ui32PolyDeg         (1,1) uint32
    dInterpDomain       (:,1) double
    dChbvCoeffs         (:,1) double
    dDataMatrix         (:,:) double
    dDomainLB           (1,1) double
    dDomainUB           (1,1) double
    bIS_ATT_QUAT        (1,1) logical
    dSwitchIntervals    (:,:) double = []
    bEnableErrorThrow   (1,1) logical = true
    dPercRelErrorTol    (1,1) double {mustBeNumeric} = 0.1
end
%% PROTOTYPE
% [strfitStats, dChbvInterpVector] = checkFitChbvPoly(ui32PolyDeg, ...
%                                                     dInterpDomain, ...
%                                                     dChbvCoeffs, ...
%                                                     dDataMatrix, ...
%                                                     dDomainLB, ...
%                                                     dDomainUB, ...
%                                                     bIS_ATT_QUAT, ...
%                                                     dSwitchIntervals, ...
%                                                     bEnableErrorThrow, ...
%                                                     dPercRelErrorTol) %#codegen
% -------------------------------------------------------------------------------------------------------------
%% DESCRIPTION
% Function performing fitting check for Chebyshev interpolation functions. It uses randomly picked input
% sample points to perform verification that the interpolation has been fitted correctly. The function uses
% the dot product to evaluate quaternions instead of subtraction.
% -------------------------------------------------------------------------------------------------------------
%% INPUT
% ui32PolyDeg         (1, 1) uint32
% dInterpDomain       (:, 1) double
% dChbvCoeffs         (:, 1) double
% dDataMatrix         (:, :) double
% dDomainLB           (1, 1) double
% dDomainUB           (1, 1) double
% bIS_ATT_QUAT        (1, 1) logical
% dSwitchIntervals    (:, :) double = []
% bEnableErrorThrow   (1,1) logical = true
% dPercRelErrorTol    (1,1) double {mustBeNumeric} = 0.1
% -------------------------------------------------------------------------------------------------------------
%% OUTPUT
% strfitStats
% dChbvInterpVector
% -------------------------------------------------------------------------------------------------------------
%% CHANGELOG
% 08-05-2024    Pietro Califano     Function adapted from testing script.
% 06-08-2025    Pietro Califano     Add error tolerance check and error throwing
% -------------------------------------------------------------------------------------------------------------
%% DEPENDENCIES
% evalAttQuatChbvPolyWithCoeffs()
% evalChbvPolyWithCoeffs()
% -------------------------------------------------------------------------------------------------------------

%% Function code
if bIS_ATT_QUAT == true
    ui8OutputSize = 4; % HARDCODED for specialization
    assert( size(dDataMatrix, 1) == ui8OutputSize );
else
    ui8OutputSize = size(dDataMatrix, 1);
    dSwitchIntervals = [];
end

assert( size(dDataMatrix, 2) == length(dInterpDomain) );

% Evaluation at test points
ui32Npoints = 5000;
ui32Npoints = ui32Npoints - 2;

testpointsIDs = sort( randi( length(dInterpDomain), ui32Npoints, 1 ), 'ascend' );

TestPoints_Time = [dInterpDomain(1); dInterpDomain(testpointsIDs); dInterpDomain(end)];
TestPoints_Labels = [dDataMatrix(:, 1),...
                dDataMatrix(:, testpointsIDs), ...
                dDataMatrix(:, end)];

dChbvInterpVector = zeros(ui8OutputSize, length(TestPoints_Time));

evalRunTime = zeros(length(TestPoints_Time), 1);

for idP = 1:length(TestPoints_Time)
    dEvalPoint = TestPoints_Time(idP);

    if bIS_ATT_QUAT == true 
        tic
        dChbvInterpVector(:, idP) = evalAttQuatChbvPolyWithCoeffs(ui32PolyDeg, ui8OutputSize, ...
            dEvalPoint, dChbvCoeffs, dSwitchIntervals, dDomainLB, dDomainUB);

    else
        tic
        dChbvInterpVector(:, idP) = evalChbvPolyWithCoeffs(ui32PolyDeg, ui8OutputSize, ...
            dEvalPoint, dChbvCoeffs, dDomainLB, dDomainUB);
    end

    evalRunTime(idP) = toc;
end
fprintf("\nAverage interpolant evaluation time: %4.4g [s]\n", mean(evalRunTime))

% Error evaluation
strFitStats = struct();


if bIS_ATT_QUAT
    strFitStats.dAbsErrVec = abs(abs(dot(dChbvInterpVector, TestPoints_Labels, 1)) - 1);

    strFitStats.dMaxAbsErr = max(strFitStats.dAbsErrVec, [], 'all');
    strFitStats.dAvgAbsErr = mean(strFitStats.dAbsErrVec, 2);

    fprintf('Max absolute difference of (q1-dot-q2 - 1): %4.4g [-]\n', strFitStats.dMaxAbsErr);
    fprintf('Average absolute difference of (q1-dot-q2 - 1): %4.4g [-]\n', strFitStats.dAvgAbsErr);

    if bEnableErrorThrow
        assert( all([strFitStats.dAvgAbsErr, strFitStats.dMaxAbsErr] < 1E-3), ...
            ['ERROR: fitting validation failed to meet tolerances for attitude quaternion. ' ...
            'Found distance greater than 1E-3 at sampling nodes.']);
    end
else
    strFitStats.dAbsErrVec = abs(dChbvInterpVector - TestPoints_Labels);
    strFitStats.dRelErrVec = strFitStats.dAbsErrVec./vecnorm(TestPoints_Labels, 2, 1);

    strFitStats.dMaxAbsErr = max(strFitStats.dAbsErrVec, [], 'all');
    strFitStats.dAvgAbsErr = mean(strFitStats.dAbsErrVec, 2);

    strFitStats.dMaxRelErr = 100*max(strFitStats.dRelErrVec, [], 'all');
    strFitStats.dAvgRelErr = 100*mean(strFitStats.dRelErrVec, 2);

    % Printing
    fprintf('Max absolute error: %4.4g [-]\n', strFitStats.dMaxAbsErr);

    % Print average absolute errors per component on one line (wraps for large N)
    fprintf('Average absolute errors (per component): ');
    fprintf('%4.4g ', strFitStats.dAvgAbsErr);
    fprintf('[-]\n');

    fprintf('\nMax relative error: %4.4g [%%]\n', strFitStats.dMaxRelErr);
    fprintf('Average relative errors (per component): ');
    fprintf('%4.4g ', strFitStats.dAvgRelErr);
    fprintf('[%%]\n');

    if bEnableErrorThrow
        assert( all([strFitStats.dMaxRelErr, max(strFitStats.dAvgRelErr)] <= dPercRelErrorTol), ...
            ['ERROR: fitting validation failed to meet relative tolerances. ' ...
            'Found relative error greater than 0.1% at sampling nodes.']);
    end
end

end
