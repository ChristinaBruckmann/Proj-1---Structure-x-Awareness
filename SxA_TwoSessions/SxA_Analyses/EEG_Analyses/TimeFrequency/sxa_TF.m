% TF Analysis
function []=sxa_TF(subj)
%% TF Analysis
% Input: Subject Number, generate plots (1/0)
disp('Starting Time Frequency Analysis')
singletrials=1; %(save single trials?)
% Parameters
triggercodes={71;72;73}; % Warning Signals per condition
%triggercodes=[71 72 73]; % All Warning Signals (rhythm, interval, irregular)
timerange=[-200 1500];
paddingLength=500; % in ms
TF_wavFreqs=1:40; % Range taken from Elmira, Assaf uses 1:30 in plos bio
%alpharange=8:12;

% Load Data
cd 'C:\Users\cbruckmann\Documents\PhD Projects\Proj1 - StructurexAwareness\SxA_TwoSessions\Data\EEG Preprocessed'
loadfilename=sprintf('EEG_SxA_Subj%i_Session2_pp.mat',subj);
savefilename=sprintf('EEG_SxA_Subj%i_Results_New.mat',subj);
load(loadfilename)

% Extract Info and Data from File
data=SDATA.data;
artifacts=SDATA.metadata.artifacts;
triggers=SDATA.events.triggerChannel;
srate=SDATA.info.sampling_rate;

% Segmentation
for c=1:size(triggercodes,1) % For conditions
    [segmentedData, isNotArtifact, timeVec]=segmentContEEGdata(triggercodes{c}, timerange+[-paddingLength paddingLength], data, triggers, artifacts, srate);

    % Check how many trials are artifact-free
    sprintf('Proportion of artifact-free trials: %.2f', mean(isNotArtifact))

    % Remove trials with artifacts
    artrej=1;
    if artrej==1
        segmentedData=segmentedData(:,:,isNotArtifact==1); 
    end

    % TF Analysis
    tic
    if ~singletrials
        parfor el=1:width(data)
            [wvlt_amp, ~] = morletwave(TF_wavFreqs, 12, squeeze(segmentedData(:,el,:))', srate, 0, 'waitbar', 'on'); % time points x frequnencies x trials
            inducedMat=squeeze(mean(wvlt_amp,3)); % average over trials
            condResults(:,:,el)=inducedMat'; % reformat to: time points x frequencies x electrodes
        end
    else
        parfor el=1:width(data)
            [wvlt_amp, ~] = morletwave(TF_wavFreqs, 12, squeeze(segmentedData(:,el,:))', srate, 0, 'waitbar', 'on'); % time points x frequnencies x trials
            inducedMat=wvlt_amp; % do not average over trials
            condResults(:,:,:,el)=permute(inducedMat,[2,1,3]); % reformat to: time points x frequencies x trials x electrodes
        end
    end
    toc

    % Remove Padding
    if ~singletrials
        condResults=condResults(timeVec>=timerange(1) & timeVec<=timerange(2) ,:,:);
    else
        condResults=condResults(timeVec>=timerange(1) & timeVec<=timerange(2) ,:,:,:);
    end

    timeVec=timeVec(timeVec>=timerange(1) & timeVec<=timerange(2));
    TF_Results{c}=condResults; % Time Points, Frequencies, Trials, Electrodes
    TF_timeVecTotal{c}=timeVec;
    TF_NotArtifact{c}=isNotArtifact;
    clear condResults segmentedData timeVec subj
end
cd 'C:\Users\cbruckmann\Documents\PhD Projects\Proj1 - StructurexAwareness\SxA_TwoSessions\Data Analysis\SxA_EEG_Analyses_Current\Results'
if ~singletrials
    save(savefilename, 'TF_Results', 'TF_timeVecTotal', 'TF_wavFreqs');
else
    TF_Results_Trial=TF_Results;
    TF_trial_timeVec=TF_timeVecTotal;
    save(savefilename, 'TF_Results_Trial','TF_trial_timeVec','TF_NotArtifact','-v7.3');
end
disp('TF Results Saved')
end