close all
clear
clc

%% testVectorSampling
% TODO this is currently identical to scalar multisampling!
ui32NumPoints = 10;
dMinVec = 1:10;
dMaxVec = dMinVec + 1;

dScatteredValues = uniformScatter(dMinVec, dMaxVec, ui32NumPoints);

return


%% testScalarMultisampling
ui32NumPoints = 10;
dMinVec = 0;
dMaxVec = 10;

dScatteredValues = uniformScatter(dMinVec, dMaxVec, ui32NumPoints);

return

%% testScalarMultisamplingAsArray
ui32NumPoints = 100;
dMinVec = rand(1,10);
dMaxVec = dMinVec + 10;

dScatteredValues = uniformScatter(dMinVec, dMaxVec, ui32NumPoints);

assert(all(dScatteredValues <= dMaxVec', 'all'));
assert(all(dScatteredValues >= dMinVec', 'all'));

return