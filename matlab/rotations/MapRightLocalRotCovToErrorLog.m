function [dErrorCovariance, dErrorJacobian] = MapRightLocalRotCovToErrorLog(dEstimatedDCM, ...
    dTrueDCM, dRightLocalCovariance) %#codegen
%% SIGNATURE
% [dErrorCovariance, dErrorJacobian] = MapRightLocalRotCovToErrorLog( ...
%     dEstimatedDCM, dTrueDCM, dRightLocalCovariance)
% -------------------------------------------------------------------------------------------------------------
%% DESCRIPTION
% First-order map a right-local rotation covariance, R_estimated * Exp(delta), into the coordinates of
% Log(R_estimated * R_true'). The transform includes both the truth-frame adjoint and inverse right Jacobian at the
% realized error; it is therefore not generally a mean-frame covariance rotation.
% -------------------------------------------------------------------------------------------------------------
%% INPUT
% dEstimatedDCM          Estimated direction cosine matrix.
% dTrueDCM               True direction cosine matrix for the same relation and epoch.
% dRightLocalCovariance  Covariance of the right-local perturbation delta [rad^2].
% -------------------------------------------------------------------------------------------------------------
%% OUTPUT
% dErrorCovariance       First-order covariance of the truth-relative logarithmic error [rad^2].
% dErrorJacobian         Jacobian of the error logarithm with respect to the right-local perturbation.
% -------------------------------------------------------------------------------------------------------------
%% CHANGELOG
% 13-08-2026  Pietro Califano, Codex gpt-5.6     First implementation.
% -------------------------------------------------------------------------------------------------------------
%% DEPENDENCIES
% ComputeSO3AttitudeError, skewSymm.
% -------------------------------------------------------------------------------------------------------------

arguments (Input)
    dEstimatedDCM (3,3) double {mustBeFinite}
    dTrueDCM (3,3) double {mustBeFinite}
    dRightLocalCovariance (3,3) double {mustBeFinite}
end
arguments (Output)
    dErrorCovariance (3,3) double
    dErrorJacobian (3,3) double
end

dSymmetryTolerance = 1.0e-12 * max(1.0, norm(dRightLocalCovariance, 'fro'));
if norm(dRightLocalCovariance - transpose(dRightLocalCovariance), 'fro') > dSymmetryTolerance
    error('MapRightLocalRotCovToErrorLog:InvalidCovariance', ...
        'The right-local covariance must be symmetric.');
end

[dErrorRotationVector, dErrorTheta] = ComputeSO3AttitudeError(dEstimatedDCM, dTrueDCM);
dErrorSkew = skewSymm(dErrorRotationVector);

if dErrorTheta < 1.0e-6
    % Use a second-order Taylor expansion of the inverse right Jacobian for small angles
    dThetaSquared = dErrorTheta * dErrorTheta;
    dSecondOrderCoefficient = 1.0 / 12.0 + dThetaSquared / 720.0 + ...
        dThetaSquared * dThetaSquared / 30240.0;
else
    % Use the exact inverse right Jacobian for larger angles
    dSecondOrderCoefficient = 1.0 / (dErrorTheta * dErrorTheta) - ...
        (1.0 + cos(dErrorTheta)) / (2.0 * dErrorTheta * sin(dErrorTheta));
end

% Compute the inverse right Jacobian and the error covariance
dInverseRightJacobian = eye(3) + 0.5 .* dErrorSkew + ...
    dSecondOrderCoefficient .* (dErrorSkew * dErrorSkew);

dErrorJacobian = dInverseRightJacobian * dTrueDCM;
dErrorCovariance = dErrorJacobian * dRightLocalCovariance * transpose(dErrorJacobian);
dErrorCovariance = 0.5 .* (dErrorCovariance + transpose(dErrorCovariance));

end
