close all
clear
clc


%% TEST: AVAR function by GPT 
%% SCRIPT NAME
% -------------------------------------------------------------------------------------------------------------
%% DESCRIPTION
% What the script does
% -------------------------------------------------------------------------------------------------------------
%% NEEDED FROM BASE WORKSPACE
% in1 [dim] description
% -------------------------------------------------------------------------------------------------------------
%% OUT TO BASE WORKSPACE
% out1 [dim] description
% -------------------------------------------------------------------------------------------------------------
%% CHANGELOG
% 19-06-2024      Pietro Califano     Testing of Allan Variance function 
% -------------------------------------------------------------------------------------------------------------
%% DEPENDENCIES
% -------------------------------------------------------------------------------------------------------------
%% Future upgrades

% Load data
addpath("Gyro_Samples/");

gyroDataStruct = load('Gyro_Samples/gyro1_1kHz.mat');
gyroDataFields = fieldnames(gyroDataStruct);
% gyroDataMat = table2array(gyroDataStruct.(gyroDataFields{1}));


%% Plot signals
timeStamps = gyroDataMat(:, 1); % [ms]
deltaTimes = gyroDataMat(:, 2); % [ms]
accData  = gyroDataMat(:, 3:5); % [m/s^2]
rateData = gyroDataMat(:, 6:8); % [deg/s]

cellName = {'Acc. X', 'Acc. Y', 'Acc. Z'};

% Accelerations
figure;
for id = 1:3

    subplot(3,1,id)
    plot(timeStamps, accData(:,id));

    title(cellName{id})
    xlabel('Time [s]')
    ylabel('Acceleration [m/s^2]')
    DefaultPlotOpts();

end

% Angular rates
cellName = {'Rate X', 'Rate Y', 'Rate Z'};

figure;
for id = 1:3

    subplot(3,1,id)
    plot(timeStamps, rateData(:,id));

    title(cellName{id})
    xlabel('Time [s]')
    ylabel('Angular rate [deg/s]')
    DefaultPlotOpts();

end


%% Perform Allan Variance
i_strGyrosSetupData.dSamplingFrequency = 1000 * 1./round(mean(deltaTimes(2:end))); % [Hz]
i_dGyroSampleData = [accData, rateData];

% Perform AVAR 
[o_dAllanVarData, o_uiMaxClusterSize, o_dTauAveragingTime] = ApplyAVAR(i_dGyroSampleData, i_strGyrosSetupData);

fprintf('\n');

smoothedAllanVar = smoothdata(o_dAllanVarData, 1, "sgolay", 50);

% Compute Allan Standard Deviation at specified Averaging Window times
i_dTauAveragingTime = [0.01, 0.1];
[o_dAllanStdDev_Target, o_dAllanVarData_Target] = GetAVARatAvgPeriod(i_dTauAveragingTime, i_dGyroSampleData, i_strGyrosSetupData);

% Plot the results
NstartSmoothed = 500;
figure;
loglog(o_dTauAveragingTime, sqrt(o_dAllanVarData(:, 1:3)));
hold on;
loglog(o_dTauAveragingTime(NstartSmoothed:end), sqrt(smoothedAllanVar(NstartSmoothed:end, 1:3)));

grid on;
xlabel('\tau (s)');
ylabel('\sigma(\tau) (m/s^2)');
title('Allan Deviation: Accelerations');
DefaultPlotOpts()
% legend(arrayfun(@(x) sprintf('Gyro %d', x), 1:M, 'UniformOutput', false));

figure;
loglog(o_dTauAveragingTime, sqrt(o_dAllanVarData(:, 4:6)));
hold on;
loglog(o_dTauAveragingTime(NstartSmoothed:end), sqrt(smoothedAllanVar(NstartSmoothed:end, 4:6)));
grid on;
xlabel('\tau (s)');
ylabel('\sigma(\tau) (deg/s)');
title('Allan Deviation: Angular rates');
DefaultPlotOpts()
