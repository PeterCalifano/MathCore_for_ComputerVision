function [o_dAllanStdDev, dAllanVarData] = GetAVARatAvgPeriod(i_dTauAveragingTime, i_dGyroSampleData, i_strGyrosSetupData) %#codegen
%% PROTOTYPE
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
% 19-06-2024          Pietro Califano         Prototype function to compute Allan Variance at specified times
% -------------------------------------------------------------------------------------------------------------
%% DEPENDENCIES
% [-]
% -------------------------------------------------------------------------------------------------------------
%% Future upgrades
% [-]
% -------------------------------------------------------------------------------------------------------------
%% Function code
[Nsamples, Msignals] = size(i_dGyroSampleData); % Nsamples: number of samples, Msignals: number of signals

dSamplingFrequency = i_strGyrosSetupData.dSamplingFrequency; % [Hz]
dSamplingPeriod = 1 / dSamplingFrequency;   

% Calculate the corresponding cluster sizes (mValues values)
for idTAU = 1:length(i_dTauAveragingTime)
    mValues = max(1, round(i_dTauAveragingTime(idTAU)) / dSamplingPeriod); % mValues should always be > 1
end

% Prepare to store Allan Variance results
dAllanVarData = zeros(length(mValues), Msignals);

% Loop over each signal (column)
for j = 1:Msignals
    fprintf('Analyzing signal ID %d\n', j);
    % Extract the j-th gyroscope signal
    y = i_dGyroSampleData(:, j);
    
    % Calculate Allan Variance for specified tau values
    for k = 1:length(Msignals)

        mSamples = mValues(k);

        if mSamples > Nsamples / 2
            error('Tau value too large for the given dataset');
        end
        
        numClusters = floor(Nsamples / mSamples);
        
        % Compute cluster means
        % clusterMeans = arrayfun(@(i) mean(y((i-1)*Msignals+1:i*Msignals)), 1:numClusters);
        
        % Initialize array to store cluster means
        clusterMeans = zeros(numClusters, 1);

        % Compute cluster means using a loop
        for i = 1:numClusters
            clusterMeans(i) = mean(y((i-1)*mSamples+1:i*mSamples));
        end

        % Compute differences between consecutive cluster means
        diffClusterMeans = diff(clusterMeans);
        
        % Compute Allan Variance
        dAllanVarData(k, j) = 0.5 * mean(diffClusterMeans.^2);
    end
end

% Compute Allan Deviation (square root of Allan Variance)
o_dAllanStdDev = sqrt(dAllanVarData);


end

