function [o_dAllanVarData, o_ui32ClusterSizes, o_dTauAveragingTime] = ApplyAVAR(i_dGyroSampleData, i_strGyrosSetupData) %#codegen
%% PROTOTYPE
% [o_dAllanVarData, o_uiMaxClusterSize, o_dTauAveragingTime] = ApplyAVAR(i_dGyroSampleData, i_strGyrosSetupData)
% -------------------------------------------------------------------------------------------------------------
%% DESCRIPTION
% What the function does
% -------------------------------------------------------------------------------------------------------------
%% INPUT
% in1 [dim] description
% Name1                     []
% Name2                     []
% Name3                     []
% Name4                     []
% Name5                     []
% Name6                     []
% -------------------------------------------------------------------------------------------------------------
%% OUTPUT
% out1 [dim] description
% Name1                     []
% Name2                     []
% Name3                     []
% Name4                     []
% Name5                     []
% Name6                     []
% -------------------------------------------------------------------------------------------------------------
%% CHANGELOG
% 19-06-2024          Pietro Califano         Prototype function to perform Allan Variance
% -------------------------------------------------------------------------------------------------------------
%% DEPENDENCIES
% [-]
% -------------------------------------------------------------------------------------------------------------
%% Future upgrades
% [-]
% -------------------------------------------------------------------------------------------------------------
%% Function code

dSamplingFrequency = i_strGyrosSetupData.dSamplingFrequency; % [Hz]
% i_strGyrosSetupData

% Allan Variance Calculation for Gyroscope Signals
% Ensure you have i_dGyroSampleData (NxM matrix, where N is number of samples and M is number of signals)
% and dSamplingFrequency (sampling frequency in Hz) variables in your workspace.

% Parameters
[Nsamples, Msignals] = size(i_dGyroSampleData); % N: number of samples, Msignals: number of signals
dSamplingPeriod = 1 / dSamplingFrequency;           % Basic sample interval (sampling period)

% Maximum number of cluster sizes
MaxClusterSize = floor(Nsamples / 2);

deltaClusterSize = 2;

o_ui32ClusterSizes = 1:deltaClusterSize:MaxClusterSize;

% Prepare to store Allan Variance results
o_dAllanVarData = zeros(length(o_ui32ClusterSizes), Msignals);

% Loop over each signal (column)
for j = 1:Msignals

    fprintf('Analyzing signal ID %d\n', j);

    % Extract the j-th gyroscope signal
    y = i_dGyroSampleData(:, j);
    
    % Calculate Allan Variance for different cluster sizes
    for idCluster = 1:length(o_ui32ClusterSizes)
        
        mClusterSize = o_ui32ClusterSizes(idCluster);
        % MaxClusterSize = mClusterSize * dSamplingPeriod;
        numClusters = floor(Nsamples / mClusterSize);
        
        % Compute cluster means
        % clusterMeans = arrayfun(@(i) mean(y((i-1)*Msignals+1:i*Msignals)), 1:numClusters);
        
        % Initialize array to store cluster means
        clusterMeans = zeros(numClusters, 1);

        % Compute cluster means using a loop
        for i = 1:numClusters
            clusterMeans(i) = mean(y((i-1)*mClusterSize+1:i*mClusterSize));
        end

        % Compute differences between consecutive cluster means
        diffClusterMeans = diff(clusterMeans);
        
        % Compute Allan Variance
        o_dAllanVarData(idCluster, j) = 0.5 * mean(diffClusterMeans.^2);
    end
end

% Compute corresponding tau values
o_dTauAveragingTime = o_ui32ClusterSizes * dSamplingPeriod;





end
