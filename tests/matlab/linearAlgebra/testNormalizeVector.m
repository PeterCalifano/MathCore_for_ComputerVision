% TEST SETUP
dVector = [3.0; 4.0; 0.0];

%% test nominal normalization
dUnitVector = NormalizeVector(dVector, 1.0e-12);
assert(norm(dUnitVector - [0.6; 0.8; 0.0]) < 1.0e-14);
assert(abs(norm(dUnitVector) - 1.0) < 1.0e-14);

%% test default tolerance
dUnitVectorDefault = NormalizeVector(dVector);
assert(norm(dUnitVectorDefault - dUnitVector) < 1.0e-14);

%% test near-zero vector failure
bThrew = false;
try
    NormalizeVector([0.0; 0.0; 0.0], 1.0e-12);
catch
    bThrew = true;
end
assert(bThrew, 'NormalizeVector must fail for near-zero vectors.');
