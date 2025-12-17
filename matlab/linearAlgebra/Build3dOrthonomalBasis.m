function [dOrthogAx1, dOrthogAx2, dUnitAx3] = Build3dOrthonomalBasis(dVec3)%#codegen
arguments
    dVec3 (3,1) double
end

dUnitAx3 = dVec3 / norm(dVec3);

% Get auxiliary axes
if abs(dUnitAx3(3)) < 0.9
    dAuxAxis = [0; 0; 1];
else
    dAuxAxis = [1; 0; 0];
end

% Compute axes 1 (X)
dOrthogAx1 = cross(dUnitAx3, dAuxAxis);
dOrthogAx1 = dOrthogAx1 ./ norm(dOrthogAx1);

% Compute axes 2 (Y)
dOrthogAx2 = cross(dUnitAx3, dOrthogAx1);
dOrthogAx2 = dOrthogAx2 ./ norm(dOrthogAx2);

end
