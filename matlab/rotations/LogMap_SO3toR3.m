function dPhiVec = LogMap_SO3toR3(dRot3) %#codegen
%% SIGNATURE
% dPhiVec = LogMap_SO3toR3(dRot3)
% -------------------------------------------------------------------------------------------------------------
%% DESCRIPTION
% Compute the principal SO(3) logarithm as a rotation vector. Small angles use the skew part directly; rotations near
% pi recover the axis from the symmetric part to avoid division by a vanishing sine.
% -------------------------------------------------------------------------------------------------------------
%% INPUT
% dRot3       Proper orthogonal 3-by-3 direction cosine matrix.
% -------------------------------------------------------------------------------------------------------------
%% OUTPUT
% dPhiVec     Principal rotation vector with norm in [0, pi] [rad].
% -------------------------------------------------------------------------------------------------------------
%% CHANGELOG
% 01-05-2026  Pietro Califano     First prototype.
% 13-08-2026  Pietro Califano, Codex gpt-5.6     Stabilize the principal logarithm near pi.
% -------------------------------------------------------------------------------------------------------------
%% DEPENDENCIES
% [-]
% -------------------------------------------------------------------------------------------------------------

arguments (Input)
    dRot3 (3,3) double {mustBeFinite}
end
arguments (Output)
    dPhiVec (3,1) double
end

% Convert rotation matrix to rotation vector
dCosTheta = min(1.0, max(-1.0, 0.5 * (trace(dRot3) - 1.0)));
dTheta = acos(dCosTheta);
dSkewDifference = dRot3 - transpose(dRot3);
dSinThetaAxis = 0.5 .* [dSkewDifference(3,2); dSkewDifference(1,3); dSkewDifference(2,1)];

if dTheta < 1.0e-10
    dPhiVec = dSinThetaAxis;
    return
end

if pi - dTheta < 1.0e-6
    % Near pi the skew part loses the axis. Recover its largest component
    % from (R + I)/2, then use symmetric off-diagonal terms for the rest.
    dAxisOuterProduct = 0.5 .* (dRot3 + eye(3));
    dAxisDiagonal = max(diag(dAxisOuterProduct), 0.0);
    [~, ui32LargestAxisIndex] = max(dAxisDiagonal);
    dRotationAxis = zeros(3,1);
    dLargestAxisComponent = sqrt(dAxisDiagonal(ui32LargestAxisIndex));

    if dLargestAxisComponent > 1.0e-12
        dRotationAxis(ui32LargestAxisIndex) = dLargestAxisComponent;
        switch ui32LargestAxisIndex
            case 1
                dRotationAxis(2) = (dRot3(1,2) + dRot3(2,1)) / (4.0 * dLargestAxisComponent);
                dRotationAxis(3) = (dRot3(1,3) + dRot3(3,1)) / (4.0 * dLargestAxisComponent);
            case 2
                dRotationAxis(1) = (dRot3(1,2) + dRot3(2,1)) / (4.0 * dLargestAxisComponent);
                dRotationAxis(3) = (dRot3(2,3) + dRot3(3,2)) / (4.0 * dLargestAxisComponent);
            otherwise
                dRotationAxis(1) = (dRot3(1,3) + dRot3(3,1)) / (4.0 * dLargestAxisComponent);
                dRotationAxis(2) = (dRot3(2,3) + dRot3(3,2)) / (4.0 * dLargestAxisComponent);
        end
    else
        dRotationAxis = [1.0; 0.0; 0.0];
    end

    dRotationAxis = dRotationAxis ./ norm(dRotationAxis);
    if norm(dSinThetaAxis) > 1.0e-12 && dot(dRotationAxis, dSinThetaAxis) < 0.0
        dRotationAxis = -dRotationAxis;
    end

    dPhiVec = dTheta .* dRotationAxis;
    return
end

% For angles away from zero and pi, the axis is recovered from the skew part.
dPhiVec = (dTheta / sin(dTheta)) .* dSinThetaAxis;

end
