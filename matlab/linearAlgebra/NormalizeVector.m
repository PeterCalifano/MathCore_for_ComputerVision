function dUnitVector = NormalizeVector(dVector, dMinNorm)%#codegen
arguments
    dVector  (3,1) 
    dMinNorm (1,1) {mustBePositive, coder.mustBeConst} = eps
end
% Normalizes a vector and with zero norm safe division.

coder.inline('always');

if coder.const(nargin < 2 || isempty(dMinNorm))
    dMinNorm = eps;
end

dNorm = max(norm(dVector), dMinNorm);
assert(isfinite(dNorm) && dNorm > dMinNorm);
dUnitVector = dVector ./ dNorm;

end
