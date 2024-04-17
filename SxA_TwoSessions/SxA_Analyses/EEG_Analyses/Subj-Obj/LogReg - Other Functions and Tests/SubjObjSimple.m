% Simple Subj-Obj Power Analysis
clear
clc
subj=17:22;
timep=[250 750]; %(data aligned to warning signal, target at 800)
electrodes=[25:30 62:64]; % occipital
alpharange=[8:12];
%% Calculate TF on a trial basis
% Parameters
triggercodes={71;72;73}; % Warning Signals per condition
%triggercodes=[71 72 73]; % All Warning Signals (rhythm, interval, irregular)
timerange=[-200 1500];
paddingLength=500; % in ms
TF_trial_wavFreqs=1:40; % Range taken from Elmira, Assaf uses 1:30 in plos bio

for s=1:length(subj)
% Load Data
cd 'C:\Users\cbruckmann\Documents\PhD Projects\Proj1 - StructurexAwareness\Two Session Version - Code and Data\Data Analysis\EEG\EEG Data'
loadfilename=sprintf('EEG_SxA_Subj%i_Session2_pp.mat',subj(s));
savefilename=sprintf('EEG_SxA_Subj%i_Results_New.mat',subj(s));
load(loadfilename)

% Extract Info and Data from File
data=SDATA.data;
artifacts=SDATA.metadata.artifacts;
triggers=SDATA.events.triggerChannel;
srate=SDATA.info.sampling_rate;

% Segmentation
for c=1:size(triggercodes,1) % For conditions
    [segmentedData, TF_isNotArtifact, timeVec]=segmentContEEGdata(triggercodes{c}, timerange+[-paddingLength paddingLength], data, triggers, artifacts, srate);

    % Check how many trials are artifact-free
    sprintf('Proportion of artifact-free trials: %.2f', mean(TF_isNotArtifact))

    % Remove trials with artifacts
    artrej=1;
    if artrej==1
        segmentedData=segmentedData(:,:,TF_isNotArtifact==1); 
    end

    % TF Analysis
    tic
    parfor el=1:width(data)
        [wvlt_amp, ~] = morletwave(TF_trial_wavFreqs, 12, squeeze(segmentedData(:,el,:))', srate, 0, 'waitbar', 'on'); % time points x frequnencies x trials
        %inducedMat=squeeze(mean(wvlt_amp,3)); % average over trials
        condResults(:,:,:,el)=wvlt_amp; % frequencies x time points x trials x electrodes
    end
    toc

    % Remove Padding
    condResults=condResults(:,timeVec>=timerange(1) & timeVec<=timerange(2) ,:,1:64);
    timeVec=timeVec(timeVec>=timerange(1) & timeVec<=timerange(2));
    TF_Results_Trial{c}=condResults; % Time Points, Frequencies, Electrodes
    TF_trial_timeVec{c}=timeVec;
    TF_NotArtifact{c}=TF_isNotArtifact;
    clear condResults segmentedData timeVec
end

cd 'C:\Users\cbruckmann\Documents\PhD Projects\Proj1 - StructurexAwareness\Two Session Version - Code and Data\Data Analysis\SxA_EEG_Analyses_Current\Results'
save(savefilename, 'TF_Results_Trial', 'TF_trial_timeVec', 'TF_trial_wavFreqs','TF_NotArtifact',"-v7.3",'-append');
disp('TF Results Saved')
end
%% Calculate average power of 500ms pre target window (smearing?)
for s=1:length(subj)
% Load TF Data 
cd 'C:\Users\cbruckmann\Documents\PhD Projects\Proj1 - StructurexAwareness\Two Session Version - Code and Data\Data Analysis\SxA_EEG_Analyses_Current\Results'

loadfilename=sprintf('EEG_SxA_Subj%i_Results_New.mat',subj(s));
load(loadfilename,'TF_Results_Trial','TF_trial_timeVec','TF_NotArtifact')

TF_res{s}=TF_Results_Trial; % Size: frequencies timepoints trials electrodes
TF_time{s}=TF_trial_timeVec; 
TF_ArtVectors{s}=TF_isNotArtifact;
end

% Convert to power and average across 500ms pre target window
for s=1:length(subj)
    for c=1:3
        % Select Data
        data=TF_res{1, s}{1, c};
        timevec=TF_time{1, s}{1, c};
        idx=timevec>=timep(1)&timevec<=timep(2);
        data=mean(data(alpharange,idx,:,electrodes),4); % average across electrodes (freq,time, trials)

        % Convert to Power
        data_power=(abs(data)).^2;

        % Average Across all Time points and frequencies
        data_power=squeeze(mean(data_power,2));
        power_avg{c}=mean(squeeze(data_power),1); % average alpha power pe trial across electrodes and alpha frequencies
    end
    power_group{s}=power_avg; %(subject, condition, power)
    clear power_avg data_power
end

%% Save group level average power
cd 'C:\Users\cbruckmann\Documents\PhD Projects\Proj1 - StructurexAwareness\Two Session Version - Code and Data\Data Analysis\SxA_EEG_Analyses_Current\Results'
save("GroupLevelAlphaPower", 'power_group','TF_ArtVectors');

%% Load Behavioural Data
cd 'C:\Users\cbruckmann\Documents\PhD Projects\Proj1 - StructurexAwareness\Two Session Version - Code and Data\Data Analysis\Behaviour'
for s=1:length(subj)
    loadfilename=sprintf('SxA_ResultsSubject%i_Total',subj(s));
    load(loadfilename,'alldataclean')
    behaviour_group{s}=alldataclean;
    clear alldataclean
end
cd 'C:\Users\cbruckmann\Documents\PhD Projects\Proj1 - StructurexAwareness\Two Session Version - Code and Data\Data Analysis\SxA_EEG_Analyses_Current\Results'
save('behaviour_group');

%% Correlate
clearvars -except subj
clc

% Load
cd 'C:\Users\cbruckmann\Documents\PhD Projects\Proj1 - StructurexAwareness\SxA_TwoSessions\Data Analysis\SxA_EEG_Analyses_Current\Results'
load GroupLevelAlphaPower
load behaviour_group

% For each subject, create a table with power (averaged across time points, electrodes and time points),subj,obj for each trial
for s=1:length(subj)
    for c=1:3
    logidx=logical(behaviour_group{s}{:,'Condition'}==c);
    obj_resp=behaviour_group{s}{logidx,'Correct/Incorrect'};
    subj_resp=behaviour_group{s}{logidx,'Binary Visibility'};
    al_power=power_group{s}{c};
    all_data{s}{c}=[obj_resp subj_resp al_power];
    end
end