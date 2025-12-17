function dPhiVec = LogMap_SO3toR3(dRot3)%#codegen
arguments
    dRot3 (3,3) double {mustBeFinite}
end
% LogMap_SO3toR3  Log map SO(3)->R^3, returns a 3x1 rotation error vector dPhi.
% Robust for small angles. Uses: log(R) = (theta/(2*sin(theta))) * (R - R').

dTmpTrace = trace(dRot3);
dTmpCosTheta    = 0.5 * (dTmpTrace - 1.0);
dTmpCosTheta    = min(1.0, max(-1.0, dTmpCosTheta)); % Clamp cos(theta)
dTmpTheta       = acos(dTmpCosTheta);

dTmpSkew = dRot3 - transpose(dRot3);
dTmpVee  = 0.5 * [dTmpSkew(3,2); dTmpSkew(1,3); dTmpSkew(2,1)]; % vee( (R-R')/2 )

if dTmpTheta < 1e-10
    % Small-angle: dPhi ≈ vee(R - I) at first order
    dPhiVec = dTmpVee;
else
    % Full Log map
    dTmpScale = dTmpTheta / sin(dTmpTheta);
    dPhiVec = dTmpScale * dTmpVee;
end
end
