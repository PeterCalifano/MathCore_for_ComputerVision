% TEST SETUP
dVector = [3.0; 4.0; 0.0];

%% test nominal normalization
dUnitVector = NormalizeVector(dVector, 1.0e-12);
assert(norm(dUnitVector - [0.6; 0.8; 0.0]) < 1.0e-14);
assert(abs(norm(dUnitVector) - 1.0) < 1.0e-14);

%% test default tolerance
dUnitVectorDefault = NormalizeVector(dVector);
assert(norm(dUnitVectorDefault - dUnitVector) < 1.0e-14);

%% test near-zero vector fallback
lastwarn('', '');
dZeroUnitVector = NormalizeVector([0.0; 0.0; 0.0], 1.0e-12);
[~, charWarningId] = lastwarn;
assert(norm(dZeroUnitVector) == 0.0);
assert(strcmp(charWarningId, 'NormalizeVector:ZeroNorm'), ...
       'NormalizeVector must warn when returning zero for near-zero vectors.');
