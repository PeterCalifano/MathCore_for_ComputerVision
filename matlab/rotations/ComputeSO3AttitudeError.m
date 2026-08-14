function [dErrorRotationVector, dTheta, dErrorDCM] = ComputeSO3AttitudeError(dEstimatedDCM, dTrueDCM) %#codegen
%% SIGNATURE
% [dErrorRotationVector, dTheta, dErrorDCM] = ComputeSO3AttitudeError(dEstimatedDCM, dTrueDCM)
% -------------------------------------------------------------------------------------------------------------
%% DESCRIPTION
% Compute the principal truth-relative attitude error using R_error = R_estimated * R_true'. Both inputs must
% represent the same directed frame relation at the same epoch.
% -------------------------------------------------------------------------------------------------------------
%% INPUT
% dEstimatedDCM          Estimated direction cosine matrix.
% dTrueDCM               True direction cosine matrix for the same relation and epoch.
% -------------------------------------------------------------------------------------------------------------
%% OUTPUT
% dErrorRotationVector   Principal truth-relative logarithmic error vector [rad].
% dTheta                 Principal attitude-error angle [rad].
% dErrorDCM              Truth-relative error direction cosine matrix.
% -------------------------------------------------------------------------------------------------------------
%% CHANGELOG
% 13-08-2026  Pietro Califano, Codex gpt-5.6     First implementation.
% -------------------------------------------------------------------------------------------------------------
%% DEPENDENCIES
% LogMap_SO3toR3.
% -------------------------------------------------------------------------------------------------------------

arguments (Input)
    dEstimatedDCM (3,3) double {mustBeFinite}
    dTrueDCM (3,3) double {mustBeFinite}
end
arguments (Output)
    dErrorRotationVector (3,1) double
    dTheta (1,1) double
    dErrorDCM (3,3) double
end

dErrorDCM = dEstimatedDCM * transpose(dTrueDCM);
dErrorRotationVector = LogMap_SO3toR3(dErrorDCM);
dTheta = norm(dErrorRotationVector);

end
