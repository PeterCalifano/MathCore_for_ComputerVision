function tests = testPlotAttitude
tests = functiontests(localfunctions);
end

function setupOnce(testCase)
testFolder = fileparts(mfilename('fullpath'));
rootFolder = fullfile(testFolder, '..', '..', '..');
matlabFolder = fullfile(rootFolder, 'matlab');
rotationsFolder = fullfile(matlabFolder, 'rotations');
quatFolder = fullfile(rotationsFolder, 'quatLib');
testCase.TestData.OriginalPath = path;
testCase.TestData.DefaultFigureVisible = get(groot, 'DefaultFigureVisible');
addpath(testFolder, matlabFolder, rotationsFolder, quatFolder);
set(groot, 'DefaultFigureVisible', 'off');
end

function teardownOnce(testCase)
close all force;
set(groot, 'DefaultFigureVisible', testCase.TestData.DefaultFigureVisible);
path(testCase.TestData.OriginalPath);
end

function testPlotAttitudeQuatSyntheticSequence(testCase)
dQuatSeq = [1.0, cos(pi/4), cos(pi/4);
            0.0, 0.0,       0.0;
            0.0, 0.0,       0.0;
            0.0, sin(pi/4), -sin(pi/4)];
dOriginPos = zeros(3, size(dQuatSeq, 2));

[fig, dDCM_Target2Fixed] = PlotAttitudeQuat(dQuatSeq, dOriginPos, false, false, false, 0.0);
cleanupObj = onCleanup(@() close(fig)); %#ok<NASGU>

verifyTrue(testCase, isgraphics(fig, 'figure'));
verifySize(testCase, dDCM_Target2Fixed, [3, 3, size(dQuatSeq, 2)]);
for idQuat = 1:size(dQuatSeq, 2)
    verifyEqual(testCase, dDCM_Target2Fixed(:, :, idQuat), Quat2DCM(dQuatSeq(:, idQuat), false), "AbsTol", 1.0e-12);
end
end
