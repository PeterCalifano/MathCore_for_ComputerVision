function tests = testAttitudePropagation
tests = functiontests(localfunctions);
end

function setupOnce(testCase)
testFolder = fileparts(mfilename('fullpath'));
rootFolder = fullfile(testFolder, '..', '..', '..');
matlabFolder = fullfile(rootFolder, 'matlab');
rotationsFolder = fullfile(matlabFolder, 'rotations');
quatFolder = fullfile(rotationsFolder, 'quatLib');
testCase.TestData.OriginalPath = path;
addpath(testFolder, matlabFolder, rotationsFolder, quatFolder);
end

function teardownOnce(testCase)
path(testCase.TestData.OriginalPath);
end

function testConstantAngularVelocityPropagation(testCase)
objQuatIntegr = CQuatKinematicsIntegrator();
dQuat0 = [1.0; 0.0; 0.0; 0.0];
dOmega = [0.0; 0.0; pi];
dTimegrid = 0.0:0.1:1.0;
dTimestep = 0.1;

[dQuatEnd, dTimegridOut, dQuatSeq] = objQuatIntegr.integrate(dQuat0, dOmega, dTimegrid, 'rkmk4', dTimestep);

verifyEqual(testCase, dTimegridOut, dTimegrid, "AbsTol", 1.0e-12);
verifyEqual(testCase, dQuatEnd, [0.0; 0.0; 0.0; 1.0], "AbsTol", 1.0e-10);
verifyEqual(testCase, vecnorm(dQuatSeq, 2, 1), ones(1, numel(dTimegrid)), "AbsTol", 1.0e-12);
end

function testTorqueFreeSymmetricTopPropagation(testCase)
dInertia = diag([2.0, 2.0, 1.0]);
dQuat0 = [1.0; 0.0; 0.0; 0.0];
dOmega0 = [0.1; 0.0; 1.0];
dTorque = zeros(3, 1);
dTimegrid = 0.0:0.1:2.0;
objIntegrator = CRigidBodyDynamicsIntegrator(dInertia);

[dQuatSeq, dOmegaSeq] = objIntegrator.integrate(dTimegrid, dQuat0, dOmega0, dTorque, 0.1, 'rk4_rkmk4', true, 1.0);

dOmegaPrecession = (dInertia(3, 3) - dInertia(1, 1)) / dInertia(1, 1) * dOmega0(3);
dTransverseAmplitude = dOmega0(1);

verifyEqual(testCase, vecnorm(dQuatSeq, 2, 1), ones(1, numel(dTimegrid)), "AbsTol", 1.0e-10);
for idTime = 1:numel(dTimegrid)
    dTime = dTimegrid(idTime);
    verifyEqual(testCase, dOmegaSeq(1, idTime), dTransverseAmplitude * cos(dOmegaPrecession * dTime), "AbsTol", 1.0e-4);
    verifyEqual(testCase, dOmegaSeq(2, idTime), dTransverseAmplitude * sin(dOmegaPrecession * dTime), "AbsTol", 1.0e-4);
    verifyEqual(testCase, dOmegaSeq(3, idTime), dOmega0(3), "AbsTol", 1.0e-12);
end
end
