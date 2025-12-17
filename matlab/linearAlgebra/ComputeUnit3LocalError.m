function dLocalVec2 = ComputeUnit3LocalError(dUnit3Base, dUnit3Other)
% ComputeUnit3LocalError  Implements GTSAM Unit3::localCoordinates.
% Returns 2x1 tangent vector error at dUnit3Base.
arguments (Input)
    dUnit3Base  (3,1) double {mustBeFinite}
    dUnit3Other (3,1) double {mustBeFinite}
end
arguments (Output)
    dLocalVec2 (2,1) double
end

dTmpA = dUnit3Base  ./ norm(dUnit3Base);
dTmpB = dUnit3Other ./ norm(dUnit3Other);

% Build a deterministic tangent basis B = [b1 b2]
[dOrthogAx1, dOrthogAx2] = Build3dOrthonormalBasis(dTmpA); % 3x2
dTmpOrthogB = [dOrthogAx1, dOrthogAx2];

dTmpX = dot(dTmpA, dTmpB);
dTmpX = min(1.0, max(-1.0, dTmpX));
dTmpZ = 1.0 - dTmpX*dTmpX;

if dTmpZ < eps
    % Very small error case
    if dTmpX > 0.0
        % First-order expansion at x=1
        dTmpY = 1.0 - (dTmpX - 1.0) / 3.0;
    else
        dLocalVec2 = [pi; 0.0];
        return
    end
else
    % Full mapping
    dTmpY = acos(dTmpX) / sqrt(dTmpZ); % y = theta/sin(theta)
end

% Compute full error (Sphere Log map)
dLocalVec2 = transpose(dTmpOrthogB) * (dTmpY * (dTmpB - dTmpX * dTmpA));

end
