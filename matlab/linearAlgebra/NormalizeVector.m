function dUnitVector = NormalizeVector(dVector, dMinNorm)%#codegen
arguments
    dVector  (3,1) double
    dMinNorm (1,1) double {mustBePositive} = eps
end
% Normalizes a vector and with zero norm safe division.

coder.inline('always');
if coder.const(nargin < 2 || isempty(dMinNorm))
    dMinNorm = eps;
end

dNorm = norm(dVector);
if dNorm < dMinNorm
    dUnitVector = zeros(3,1);
    if coder.target('MATLAB')
        warning('NormalizeVector:ZeroNorm', 'Vector norm is below minimum threshold. Returning zero vector.');
    end
else
    dUnitVector = dVector ./ dNorm;
end

end
