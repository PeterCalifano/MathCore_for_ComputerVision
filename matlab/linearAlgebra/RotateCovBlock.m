function [dRotatedCovBlock] = RotateCovBlock(dRotMatrix, dCovMatrix)
arguments
    dRotMatrix  (3,3) double
    dCovMatrix  (3,3) double
end
%% SIGNATURE
% [objPerfMonitor] = ComputeSmoothedFullSolution(objGraphEstimator, objModeManager, objPerfMonitor)
% -------------------------------------------------------------------------------------------------------------
%% DESCRIPTION
% Function computing the entire smoothed solution of pose, velocity and target map from estimator object.
% -------------------------------------------------------------------------------------------------------------
%% INPUT
% objGraphEstimator       (1,1) {mustBeA(objGraphEstimator, ["gtsam.ISAM2", "CGraphEstimator"])} % GTSAM/nav-backend
% objModeManager          (1,1) {mustBeA(objModeManager, "CModeManager")} % In nav-backend/SimulationGears_for_SpaceNav
% objPerfMonitor          (1,1) {mustBeA(objPerfMonitor, "SPerfMonitor")} % In nav-backend
% -------------------------------------------------------------------------------------------------------------
%% OUTPUT
% objPerfMonitor          (1,1) {mustBeA(objPerfMonitor, "SPerfMonitor")} % In nav-backend
% -------------------------------------------------------------------------------------------------------------
%% CHANGELOG
% 12-09-2025    Pietro Califano     First implementation.
% 24-11-2025    Pietro Califano     Extend and split between incrementally smoothed and post-processed batch
%                                   solution (by further optimizing using DogLeg); add keyframes time index
% -------------------------------------------------------------------------------------------------------------
%% DEPENDENCIES
% gtsam library implementation
% -------------------------------------------------------------------------------------------------------------


dRotatedCovBlock = dRotMatrix * dCovMatrix * transpose(dRotMatrix);

end
